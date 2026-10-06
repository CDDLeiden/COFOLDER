"""Bias training dataset assembly and ordering."""

from __future__ import annotations

from rdkit import DataStructs
from pathlib import Path
import pandas as pd

from ._common import (
    BIAS_PAIRING_COMPONENT_1_ONLY,
    BIAS_PAIRING_COMPONENT_2_ONLY,
    BIAS_PAIRING_LIGAND_ONLY,
    BIAS_PAIRING_PAIRED,
    BIAS_PAIRING_PROTEIN_ONLY,
    BIAS_PLOT_LIGAND_THRESHOLD,
    BIAS_PLOT_SEQUENCE_THRESHOLD,
    BIAS_TRAINING_DATA_COLUMNS,
    REFERENCE_PATH_COLUMNS,
    SAME_TYPE_PAIR_DATA_COLUMNS,
    _canonical_reference_source,
    _norm_id,
    _normalize_smiles,
)

from ._similarity import (
    _morgan_fp_from_smiles,
)


def _reference_key_from_row(row: pd.Series, *, is_ligand: bool) -> str:
    source = _canonical_reference_source(row.get("source"), "custom")
    pdb_id = row.get("pdb_id")
    if pd.notna(pdb_id) and str(pdb_id).strip() and source == "public":
        return f"{source}:pdb:{_norm_id(pdb_id)}"

    complex_id = row.get("complex_id")
    if pd.notna(complex_id) and str(complex_id).strip():
        dataset = str(row.get("dataset_name") or "custom").strip()
        return f"{source}:complex:{dataset}:{str(complex_id).strip()}"

    for column in REFERENCE_PATH_COLUMNS:
        value = row.get(column)
        if pd.notna(value) and str(value).strip():
            return f"{source}:path:{Path(str(value)).expanduser().resolve()}"

    if pd.notna(pdb_id) and str(pdb_id).strip():
        return f"{source}:pdb:{_norm_id(pdb_id)}"

    dataset_name = row.get("dataset_name")
    if pd.notna(dataset_name) and str(dataset_name).strip():
        return f"{source}:dataset:{str(dataset_name).strip()}:{'ligand' if is_ligand else 'protein'}:{row.name}"

    return f"{source}:{'ligand' if is_ligand else 'protein'}:{row.name}"


def _collapse_pair_value(
    values: list[object], *, mixed_label: str | None = None
) -> str | None:
    normalized = [
        str(value).strip() for value in values if pd.notna(value) and str(value).strip()
    ]
    if not normalized:
        return None
    unique_values = list(dict.fromkeys(normalized))
    if len(unique_values) == 1:
        return unique_values[0]
    if mixed_label is not None:
        return mixed_label
    return " | ".join(sorted(unique_values))


def _reference_label_from_merged_row(row: pd.Series) -> str | None:
    for value in (
        row.get("pdb_id_protein"),
        row.get("pdb_id_ligand"),
        row.get("source_structure_path_protein"),
        row.get("source_reference_path_protein"),
        row.get("source_structure_path_ligand"),
        row.get("source_reference_path_ligand"),
        row.get("dataset_name_protein"),
        row.get("dataset_name_ligand"),
    ):
        if pd.isna(value) or not str(value).strip():
            continue
        text = str(value).strip()
        candidate = (
            Path(text).stem
            if "/" in text or text.endswith((".pdb", ".cif", ".sdf", ".csv"))
            else text
        )
        if candidate:
            return candidate
    return None


def _prepare_reference_side(
    view: pd.DataFrame, *, is_ligand: bool, suffix: str
) -> pd.DataFrame:
    prepared = view.copy()
    if "reference_key" not in prepared.columns:
        prepared["reference_key"] = prepared.apply(
            lambda row: _reference_key_from_row(row, is_ligand=is_ligand),
            axis=1,
        )
    return prepared.rename(
        columns=lambda column: (
            column if column == "reference_key" else f"{column}_{suffix}"
        )
    )


def _combine_prefixed_rows(
    protein_row: pd.Series | None,
    ligand_row: pd.Series | None,
) -> pd.Series:
    combined: dict[str, object] = {}
    if protein_row is not None:
        combined.update(protein_row.to_dict())
    if ligand_row is not None:
        combined.update(ligand_row.to_dict())
    return pd.Series(combined)


def _matching_reference_rows(
    reference_view: pd.DataFrame, reference_key: object
) -> pd.DataFrame:
    if reference_view.empty or pd.isna(reference_key) or not str(reference_key).strip():
        return reference_view.iloc[0:0].copy()
    return reference_view[reference_view["reference_key"] == reference_key].copy()


def _direct_overlap_reference_keys(
    left: pd.DataFrame,
    right: pd.DataFrame,
) -> set[str]:
    if left.empty or right.empty:
        return set()
    left_keys = {
        str(value).strip()
        for value in left.get("reference_key", pd.Series(dtype=object)).dropna()
        if str(value).strip()
    }
    right_keys = {
        str(value).strip()
        for value in right.get("reference_key", pd.Series(dtype=object)).dropna()
        if str(value).strip()
    }
    return left_keys & right_keys


def _finalize_bias_training_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=BIAS_TRAINING_DATA_COLUMNS)
    ordered = df.copy()
    if "sequence_similarity_pairwise" not in ordered.columns:
        ordered["sequence_similarity_pairwise"] = pd.NA
    if "plot_sequence_similarity_pairwise" not in ordered.columns:
        ordered["plot_sequence_similarity_pairwise"] = pd.NA
    ligand_only = ordered["pairing_status"].astype(str).eq(BIAS_PAIRING_LIGAND_ONLY)
    ordered.loc[ligand_only, "sequence_similarity"] = pd.NA
    ordered.loc[ligand_only, "sequence_similarity_pairwise"] = pd.NA
    ordered.loc[ligand_only, "plot_sequence_similarity"] = pd.NA
    ordered.loc[ligand_only, "plot_sequence_similarity_pairwise"] = pd.NA
    mmseqs = pd.to_numeric(ordered["sequence_similarity"], errors="coerce")
    pairwise = pd.to_numeric(ordered["sequence_similarity_pairwise"], errors="coerce")
    if "sequence_similarity_method" not in ordered.columns:
        ordered["sequence_similarity_method"] = pd.NA
    methods = ordered["sequence_similarity_method"].astype("string")
    missing_method = methods.isna() | methods.str.strip().eq("")
    methods = methods.mask(missing_method & mmseqs.notna(), "mmseqs_pident")
    methods = methods.mask(
        missing_method & mmseqs.isna() & pairwise.notna(), "pairwise_aligner"
    )
    methods = methods.mask(
        missing_method & mmseqs.isna() & pairwise.isna(), "unavailable"
    )
    methods = methods.mask(ligand_only, "unavailable")
    ordered["sequence_similarity_method"] = methods
    methods = methods.fillna("")
    invalid = (
        (mmseqs.notna() & pairwise.notna())
        | (mmseqs.notna() & methods.ne("mmseqs_pident"))
        | (pairwise.notna() & methods.ne("pairwise_aligner"))
        | (mmseqs.isna() & pairwise.isna() & methods.ne("unavailable"))
    )
    if invalid.any():
        raise ValueError(
            "Bias training rows must keep MMseqs pident and PairwiseAligner "
            "similarities in their method-specific columns"
        )
    ordered["plot_sequence_similarity"] = pd.to_numeric(
        ordered["plot_sequence_similarity"],
        errors="coerce",
    ).clip(lower=0.0, upper=1.0)
    ordered["plot_sequence_similarity_pairwise"] = pd.to_numeric(
        ordered.get(
            "plot_sequence_similarity_pairwise",
            pd.Series(pd.NA, index=ordered.index, dtype="Float64"),
        ),
        errors="coerce",
    ).clip(lower=0.0, upper=1.0)
    ordered["plot_ecfp_similarity"] = pd.to_numeric(
        ordered["plot_ecfp_similarity"],
        errors="coerce",
    ).clip(lower=0.0, upper=1.0)
    ordered = _sort_bias_training_dataset(ordered)
    for column in BIAS_TRAINING_DATA_COLUMNS:
        if column not in ordered.columns:
            ordered[column] = pd.NA
    return ordered.loc[:, BIAS_TRAINING_DATA_COLUMNS].reset_index(drop=True)


def _reference_pdb_id_from_bias_row(row: pd.Series) -> str | None:
    for value in (
        row.get("pdb_id"),
        row.get("protein_pdb_id"),
        row.get("ligand_pdb_id"),
        row.get("reference_label"),
    ):
        if pd.isna(value):
            continue
        token = str(value).strip().upper()
        if len(token) == 4 and token.isalnum():
            return token
    return None


def _reference_pdb_id_from_same_type_row(row: pd.Series) -> str | None:
    for value in (
        row.get("reference_pdb_id"),
        row.get("component_1_pdb_id"),
        row.get("component_2_pdb_id"),
        row.get("reference_label"),
    ):
        if pd.isna(value):
            continue
        token = str(value).strip().upper()
        if len(token) == 4 and token.isalnum():
            return token
    return None


def _reference_ligand_similarity_rows(
    *,
    pdb_id: str,
    query_smiles: str,
    ligand_reference_df: pd.DataFrame | None,
    ligand_reference_index: dict[str, list[dict[str, str]]] | None = None,
    exclude_ligand_ids: set[str] | None = None,
) -> list[dict[str, object]]:
    query_fp = _morgan_fp_from_smiles(query_smiles)
    if query_fp is None:
        return []

    excluded_ids = {
        _norm_id(ligand_id)
        for ligand_id in (exclude_ligand_ids or set())
        if str(ligand_id).strip()
    }
    candidate_rows: list[dict[str, str]]
    normalized_pdb_id = str(pdb_id).strip().upper()
    if ligand_reference_index is not None:
        candidate_rows = ligand_reference_index.get(normalized_pdb_id, [])
    else:
        if (
            ligand_reference_df is None
            or ligand_reference_df.empty
            or "pdb_id" not in ligand_reference_df.columns
        ):
            return []
        subset = ligand_reference_df[
            ligand_reference_df["pdb_id"].astype(str).str.upper() == normalized_pdb_id
        ].copy()
        if subset.empty or "smiles" not in subset.columns:
            return []
        subset["smiles"] = subset["smiles"].apply(_normalize_smiles)
        subset = subset[subset["smiles"].astype(str).str.strip().astype(bool)].copy()
        if subset.empty:
            return []
        subset = subset.drop_duplicates(
            subset=["ligand_id", "smiles"], keep="first"
        ).reset_index(drop=True)
        candidate_rows = [
            {
                "ligand_id": str(ligand_row.get("ligand_id", "")).strip(),
                "smiles": str(ligand_row.get("smiles", "")).strip(),
            }
            for _, ligand_row in subset.iterrows()
        ]

    rows: list[dict[str, object]] = []
    for ligand_row in candidate_rows:
        ligand_id = str(ligand_row.get("ligand_id", "")).strip()
        if ligand_id and _norm_id(ligand_id) in excluded_ids:
            continue
        smiles = str(ligand_row.get("smiles", "")).strip()
        if not smiles:
            continue
        ligand_fp = _morgan_fp_from_smiles(smiles)
        if ligand_fp is None:
            continue
        rows.append(
            {
                "ligand_id": ligand_id,
                "smiles": smiles,
                "ecfp_similarity": float(
                    DataStructs.TanimotoSimilarity(query_fp, ligand_fp)
                ),
            }
        )
    return rows


def _build_protein_lookup_index(
    protein_lookup_view: pd.DataFrame | None,
) -> dict[str, list[dict[str, object]]]:
    if (
        protein_lookup_view is None
        or protein_lookup_view.empty
        or "pdb_id" not in protein_lookup_view.columns
    ):
        return {}

    ordered = protein_lookup_view.copy()
    if "sequence_similarity" in ordered.columns:
        ordered["sequence_similarity"] = pd.to_numeric(
            ordered["sequence_similarity"],
            errors="coerce",
        )
    if "sequence_similarity_pairwise" in ordered.columns:
        ordered["sequence_similarity_pairwise"] = pd.to_numeric(
            ordered["sequence_similarity_pairwise"],
            errors="coerce",
        )
    else:
        ordered["sequence_similarity_pairwise"] = pd.NA
    ordered = ordered.sort_values(
        by=[
            "pdb_id",
            "sequence_similarity",
            "sequence_similarity_pairwise",
            "source",
            "dataset_name",
            "sequence",
        ],
        ascending=[True, False, False, True, True, True],
        na_position="last",
    )

    index: dict[str, list[dict[str, object]]] = {}
    for _, protein_row in ordered.iterrows():
        pdb_id = protein_row.get("pdb_id")
        if pd.isna(pdb_id) or not str(pdb_id).strip():
            continue
        index.setdefault(str(pdb_id).strip().upper(), []).append(protein_row.to_dict())
    return index


def _reference_protein_similarity_rows(
    *,
    pdb_id: str,
    protein_lookup_index: dict[str, list[dict[str, object]]] | None,
) -> list[dict[str, object]]:
    if not protein_lookup_index:
        return []
    return [
        dict(row) for row in protein_lookup_index.get(str(pdb_id).strip().upper(), [])
    ]


def _build_ligand_reference_index(
    ligand_reference_df: pd.DataFrame | None,
) -> dict[str, list[dict[str, str]]]:
    if ligand_reference_df is None or ligand_reference_df.empty:
        return {}
    if (
        "pdb_id" not in ligand_reference_df.columns
        or "smiles" not in ligand_reference_df.columns
    ):
        return {}

    subset = ligand_reference_df.copy()
    subset["pdb_id"] = subset["pdb_id"].astype(str).str.strip().str.upper()
    subset["ligand_id"] = (
        subset.get("ligand_id", pd.Series(dtype=object)).astype(str).str.strip()
    )
    subset["smiles"] = subset["smiles"].apply(_normalize_smiles)
    subset = subset[
        subset["pdb_id"].astype(bool)
        & subset["smiles"].astype(str).str.strip().astype(bool)
    ].copy()
    if subset.empty:
        return {}
    subset = subset.drop_duplicates(
        subset=["pdb_id", "ligand_id", "smiles"], keep="first"
    )

    index: dict[str, list[dict[str, str]]] = {}
    for pdb_id, group in subset.groupby("pdb_id", sort=False):
        index[str(pdb_id)] = [
            {
                "ligand_id": str(row.get("ligand_id", "")).strip(),
                "smiles": str(row.get("smiles", "")).strip(),
            }
            for _, row in group.iterrows()
        ]
    return index


def _component_similarity_column(component_type: str) -> str:
    return "sequence_similarity" if component_type == "protein" else "ecfp_similarity"


def _component_threshold(component_type: str) -> float:
    return (
        BIAS_PLOT_SEQUENCE_THRESHOLD
        if component_type == "protein"
        else BIAS_PLOT_LIGAND_THRESHOLD
    )


def _normalize_component_similarity(value: float | int | None, component_type: str):
    if pd.isna(value):
        return pd.NA
    numeric_value = float(value)
    if component_type == "protein":
        return numeric_value / 100.0
    return numeric_value


def _reference_label_from_component_rows(
    component_1_row: pd.Series | None,
    component_2_row: pd.Series | None,
) -> str | None:
    for row, suffix in (
        (component_1_row, "component_1"),
        (component_2_row, "component_2"),
    ):
        if row is None:
            continue
        for value in (
            row.get(f"pdb_id_{suffix}"),
            row.get(f"source_structure_path_{suffix}"),
            row.get(f"source_reference_path_{suffix}"),
            row.get(f"dataset_name_{suffix}"),
        ):
            if pd.isna(value) or not str(value).strip():
                continue
            text = str(value).strip()
            candidate = (
                Path(text).stem
                if "/" in text or text.endswith((".pdb", ".cif", ".sdf", ".csv"))
                else text
            )
            if candidate:
                return candidate
    return None


def _same_type_pair_row_from_sources(
    *,
    query_pair_id: str,
    component_1_id: str,
    component_2_id: str,
    component_type: str,
    reference_key: object,
    reference_label: str | None,
    pairing_status: str,
    component_1_row: pd.Series | None,
    component_2_row: pd.Series | None,
    forced_component_1_similarity: float | None = None,
    forced_component_2_similarity: float | None = None,
) -> dict[str, object]:
    similarity_column = _component_similarity_column(component_type)
    component_1_reference_key = (
        component_1_row.get("reference_key") if component_1_row is not None else None
    )
    component_2_reference_key = (
        component_2_row.get("reference_key") if component_2_row is not None else None
    )
    component_1_similarity = pd.to_numeric(
        pd.Series(
            [
                (
                    forced_component_1_similarity
                    if forced_component_1_similarity is not None
                    else (
                        component_1_row.get(f"{similarity_column}_component_1")
                        if component_1_row is not None
                        else pd.NA
                    )
                )
            ]
        ),
        errors="coerce",
    ).iloc[0]
    component_2_similarity = pd.to_numeric(
        pd.Series(
            [
                (
                    forced_component_2_similarity
                    if forced_component_2_similarity is not None
                    else (
                        component_2_row.get(f"{similarity_column}_component_2")
                        if component_2_row is not None
                        else pd.NA
                    )
                )
            ]
        ),
        errors="coerce",
    ).iloc[0]
    return {
        "query_pair_id": query_pair_id,
        "component_1_id": component_1_id,
        "component_1_type": component_type,
        "component_2_id": component_2_id,
        "component_2_type": component_type,
        "reference_key": reference_key,
        "reference_label": reference_label,
        "reference_pdb_id": _collapse_pair_value(
            [
                (
                    component_1_row.get("pdb_id_component_1")
                    if component_1_row is not None
                    else None
                ),
                (
                    component_2_row.get("pdb_id_component_2")
                    if component_2_row is not None
                    else None
                ),
            ]
        ),
        "pairing_status": pairing_status,
        "source": _collapse_pair_value(
            [
                (
                    component_1_row.get("source_component_1")
                    if component_1_row is not None
                    else None
                ),
                (
                    component_2_row.get("source_component_2")
                    if component_2_row is not None
                    else None
                ),
            ],
            mixed_label="mixed",
        ),
        "dataset_name": _collapse_pair_value(
            [
                (
                    component_1_row.get("dataset_name_component_1")
                    if component_1_row is not None
                    else None
                ),
                (
                    component_2_row.get("dataset_name_component_2")
                    if component_2_row is not None
                    else None
                ),
            ]
        ),
        "component_1_reference_key": component_1_reference_key,
        "component_1_pdb_id": (
            component_1_row.get("pdb_id_component_1")
            if component_1_row is not None
            else pd.NA
        ),
        "component_1_release_date": (
            component_1_row.get("release_date_component_1")
            if component_1_row is not None
            else pd.NA
        ),
        "component_1_source": (
            component_1_row.get("source_component_1")
            if component_1_row is not None
            else pd.NA
        ),
        "component_1_dataset_name": (
            component_1_row.get("dataset_name_component_1")
            if component_1_row is not None
            else pd.NA
        ),
        "component_1_source_structure_path": (
            component_1_row.get("source_structure_path_component_1")
            if component_1_row is not None
            else pd.NA
        ),
        "component_1_source_reference_path": (
            component_1_row.get("source_reference_path_component_1")
            if component_1_row is not None
            else pd.NA
        ),
        "component_2_reference_key": component_2_reference_key,
        "component_2_pdb_id": (
            component_2_row.get("pdb_id_component_2")
            if component_2_row is not None
            else pd.NA
        ),
        "component_2_release_date": (
            component_2_row.get("release_date_component_2")
            if component_2_row is not None
            else pd.NA
        ),
        "component_2_source": (
            component_2_row.get("source_component_2")
            if component_2_row is not None
            else pd.NA
        ),
        "component_2_dataset_name": (
            component_2_row.get("dataset_name_component_2")
            if component_2_row is not None
            else pd.NA
        ),
        "component_2_source_structure_path": (
            component_2_row.get("source_structure_path_component_2")
            if component_2_row is not None
            else pd.NA
        ),
        "component_2_source_reference_path": (
            component_2_row.get("source_reference_path_component_2")
            if component_2_row is not None
            else pd.NA
        ),
        "component_1_similarity": (
            float(component_1_similarity) if pd.notna(component_1_similarity) else pd.NA
        ),
        "component_2_similarity": (
            float(component_2_similarity) if pd.notna(component_2_similarity) else pd.NA
        ),
        "plot_component_1_similarity": _normalize_component_similarity(
            component_1_similarity, component_type
        ),
        "plot_component_2_similarity": _normalize_component_similarity(
            component_2_similarity, component_type
        ),
        "_component_1_ligand_id": (
            component_1_row.get("ligand_id_component_1")
            if component_1_row is not None
            else pd.NA
        ),
        "_component_1_smiles": (
            component_1_row.get("smiles_component_1")
            if component_1_row is not None
            else pd.NA
        ),
        "_component_2_ligand_id": (
            component_2_row.get("ligand_id_component_2")
            if component_2_row is not None
            else pd.NA
        ),
        "_component_2_smiles": (
            component_2_row.get("smiles_component_2")
            if component_2_row is not None
            else pd.NA
        ),
    }


def _same_type_pair_sort_order() -> dict[str, int]:
    return {
        BIAS_PAIRING_PAIRED: 0,
        BIAS_PAIRING_COMPONENT_1_ONLY: 1,
        BIAS_PAIRING_COMPONENT_2_ONLY: 2,
    }


def _sort_same_type_pair_dataset(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df.copy()
    ordered = df.copy()
    ordered["_pairing_order"] = (
        ordered["pairing_status"].map(_same_type_pair_sort_order()).fillna(99)
    )
    return ordered.sort_values(
        by=[
            "query_pair_id",
            "_pairing_order",
            "source",
            "dataset_name",
            "plot_component_1_similarity",
            "plot_component_2_similarity",
            "reference_label",
        ],
        ascending=[True, True, True, True, False, False, True],
        na_position="last",
    ).drop(columns="_pairing_order")


def _finalize_same_type_pair_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=SAME_TYPE_PAIR_DATA_COLUMNS)
    ordered = df.copy()
    ordered["plot_component_1_similarity"] = pd.to_numeric(
        ordered["plot_component_1_similarity"],
        errors="coerce",
    ).clip(lower=0.0, upper=1.0)
    ordered["plot_component_2_similarity"] = pd.to_numeric(
        ordered["plot_component_2_similarity"],
        errors="coerce",
    ).clip(lower=0.0, upper=1.0)
    ordered = ordered.drop_duplicates().reset_index(drop=True)
    ordered = _sort_same_type_pair_dataset(ordered)
    for column in SAME_TYPE_PAIR_DATA_COLUMNS:
        if column not in ordered.columns:
            ordered[column] = pd.NA
    return ordered.loc[:, SAME_TYPE_PAIR_DATA_COLUMNS].reset_index(drop=True)


def _bias_training_row_from_sources(
    *,
    query_protein_chain_id: str | None,
    query_ligand_chain_id: str | None,
    reference_key: object,
    reference_label: str | None,
    pairing_status: str,
    protein_row: pd.Series | None,
    ligand_row: pd.Series | None,
    protein_unique_counts: dict[object, int],
    ligand_id_unique_counts: dict[object, int],
    ligand_smiles_unique_counts: dict[object, int],
    forced_ecfp_similarity: float | None = None,
) -> dict[str, object]:
    protein_chain = query_protein_chain_id or None
    ligand_chain = query_ligand_chain_id or None
    protein_reference_key = (
        protein_row.get("reference_key") if protein_row is not None else None
    )
    ligand_reference_key = (
        ligand_row.get("reference_key") if ligand_row is not None else None
    )
    protein_source = (
        protein_row.get("source_protein") if protein_row is not None else None
    )
    ligand_source = ligand_row.get("source_ligand") if ligand_row is not None else None

    sequence_similarity = pd.to_numeric(
        pd.Series(
            [
                (
                    protein_row.get("sequence_similarity_protein")
                    if protein_row is not None
                    else pd.NA
                )
            ]
        ),
        errors="coerce",
    ).iloc[0]
    sequence_similarity_pairwise = pd.to_numeric(
        pd.Series(
            [
                (
                    protein_row.get("sequence_similarity_pairwise_protein")
                    if protein_row is not None
                    else pd.NA
                )
            ]
        ),
        errors="coerce",
    ).iloc[0]
    if pd.notna(sequence_similarity) and pd.notna(sequence_similarity_pairwise):
        raise ValueError(
            "Protein reference row contains both MMseqs pident and PairwiseAligner similarity"
        )
    sequence_similarity_method = (
        protein_row.get("sequence_similarity_method_protein")
        if protein_row is not None
        else pd.NA
    )
    if (
        pd.isna(sequence_similarity_method)
        or not str(sequence_similarity_method).strip()
    ):
        if pd.notna(sequence_similarity):
            sequence_similarity_method = "mmseqs_pident"
        elif pd.notna(sequence_similarity_pairwise):
            sequence_similarity_method = "pairwise_aligner"
        else:
            sequence_similarity_method = "unavailable"
    if (
        (
            pd.notna(sequence_similarity)
            and sequence_similarity_method != "mmseqs_pident"
        )
        or (
            pd.notna(sequence_similarity_pairwise)
            and sequence_similarity_method != "pairwise_aligner"
        )
        or (
            pd.isna(sequence_similarity)
            and pd.isna(sequence_similarity_pairwise)
            and sequence_similarity_method != "unavailable"
        )
    ):
        raise ValueError(
            "Protein similarity value does not match sequence_similarity_method"
        )
    ecfp_similarity = pd.to_numeric(
        pd.Series(
            [
                (
                    forced_ecfp_similarity
                    if forced_ecfp_similarity is not None
                    else (
                        ligand_row.get("ecfp_similarity_ligand")
                        if ligand_row is not None
                        else pd.NA
                    )
                )
            ]
        ),
        errors="coerce",
    ).iloc[0]
    plot_sequence_similarity = (
        float(sequence_similarity) / 100.0 if pd.notna(sequence_similarity) else pd.NA
    )
    plot_sequence_similarity_pairwise = (
        float(sequence_similarity_pairwise) / 100.0
        if pd.notna(sequence_similarity_pairwise)
        else pd.NA
    )
    plot_ecfp_similarity = (
        float(ecfp_similarity) if pd.notna(ecfp_similarity) else pd.NA
    )

    protein_is_public_pdb_reference = pd.notna(protein_reference_key) and str(
        protein_reference_key
    ).startswith("public:pdb:")
    ligand_is_public_pdb_reference = pd.notna(ligand_reference_key) and str(
        ligand_reference_key
    ).startswith("public:pdb:")
    if (
        protein_row is not None
        and protein_is_public_pdb_reference
        and protein_unique_counts.get(protein_reference_key, 0) > 1
    ):
        sequence_value = pd.NA
    else:
        sequence_value = (
            protein_row.get("sequence_protein") if protein_row is not None else pd.NA
        )
    if (
        ligand_row is not None
        and ligand_is_public_pdb_reference
        and ligand_id_unique_counts.get(ligand_reference_key, 0) > 1
    ):
        ligand_id_value = pd.NA
    else:
        ligand_id_value = (
            ligand_row.get("ligand_id_ligand") if ligand_row is not None else pd.NA
        )
    if (
        ligand_row is not None
        and ligand_is_public_pdb_reference
        and ligand_smiles_unique_counts.get(ligand_reference_key, 0) > 1
    ):
        smiles_value = pd.NA
    else:
        smiles_value = (
            ligand_row.get("smiles_ligand") if ligand_row is not None else pd.NA
        )

    return {
        "query_pair_id": f"{protein_chain or 'none'}__{ligand_chain or 'none'}",
        "query_protein_chain_id": protein_chain,
        "query_ligand_chain_id": ligand_chain,
        "reference_key": reference_key,
        "reference_label": reference_label,
        "pairing_status": pairing_status,
        "source": _collapse_pair_value(
            [protein_source, ligand_source], mixed_label="mixed"
        ),
        "dataset_name": _collapse_pair_value(
            [
                (
                    protein_row.get("dataset_name_protein")
                    if protein_row is not None
                    else None
                ),
                (
                    ligand_row.get("dataset_name_ligand")
                    if ligand_row is not None
                    else None
                ),
            ]
        ),
        "complex_id": _collapse_pair_value(
            [
                (
                    protein_row.get("complex_id_protein")
                    if protein_row is not None
                    else None
                ),
                ligand_row.get("complex_id_ligand") if ligand_row is not None else None,
            ]
        ),
        "pdb_id": _collapse_pair_value(
            [
                protein_row.get("pdb_id_protein") if protein_row is not None else None,
                ligand_row.get("pdb_id_ligand") if ligand_row is not None else None,
            ]
        ),
        "protein_pdb_id": (
            protein_row.get("pdb_id_protein") if protein_row is not None else pd.NA
        ),
        "ligand_pdb_id": (
            ligand_row.get("pdb_id_ligand") if ligand_row is not None else pd.NA
        ),
        "protein_release_date": (
            protein_row.get("release_date_protein")
            if protein_row is not None
            else pd.NA
        ),
        "ligand_release_date": (
            ligand_row.get("release_date_ligand") if ligand_row is not None else pd.NA
        ),
        "protein_source": protein_source,
        "protein_dataset_name": (
            protein_row.get("dataset_name_protein")
            if protein_row is not None
            else pd.NA
        ),
        "protein_source_structure_path": (
            protein_row.get("source_structure_path_protein")
            if protein_row is not None
            else pd.NA
        ),
        "protein_source_reference_path": (
            protein_row.get("source_reference_path_protein")
            if protein_row is not None
            else pd.NA
        ),
        "ligand_source": ligand_source,
        "ligand_dataset_name": (
            ligand_row.get("dataset_name_ligand") if ligand_row is not None else pd.NA
        ),
        "ligand_source_structure_path": (
            ligand_row.get("source_structure_path_ligand")
            if ligand_row is not None
            else pd.NA
        ),
        "ligand_source_reference_path": (
            ligand_row.get("source_reference_path_ligand")
            if ligand_row is not None
            else pd.NA
        ),
        "sequence_similarity": (
            float(sequence_similarity) if pd.notna(sequence_similarity) else pd.NA
        ),
        "sequence_similarity_pairwise": (
            float(sequence_similarity_pairwise)
            if pd.notna(sequence_similarity_pairwise)
            else pd.NA
        ),
        "sequence_similarity_method": sequence_similarity_method,
        "ecfp_similarity": (
            float(ecfp_similarity) if pd.notna(ecfp_similarity) else pd.NA
        ),
        "plot_sequence_similarity": plot_sequence_similarity,
        "plot_sequence_similarity_pairwise": plot_sequence_similarity_pairwise,
        "plot_ecfp_similarity": plot_ecfp_similarity,
        "sequence": sequence_value,
        "ligand_id": ligand_id_value,
        "smiles": smiles_value,
    }


def _build_bias_training_rows(
    *,
    query_protein_chain_id: str | None,
    query_ligand_chain_id: str | None,
    protein_view: pd.DataFrame,
    ligand_view: pd.DataFrame,
    protein_lookup_view: pd.DataFrame | None = None,
    ligand_lookup_view: pd.DataFrame | None = None,
) -> list[dict[str, object]]:
    left = _prepare_reference_side(protein_view, is_ligand=False, suffix="protein")
    right = _prepare_reference_side(ligand_view, is_ligand=True, suffix="ligand")
    left_lookup = _prepare_reference_side(
        protein_lookup_view if protein_lookup_view is not None else protein_view,
        is_ligand=False,
        suffix="protein",
    )
    right_lookup = _prepare_reference_side(
        ligand_lookup_view if ligand_lookup_view is not None else ligand_view,
        is_ligand=True,
        suffix="ligand",
    )

    # Public references are PDB-level aggregates in the legacy bias builder.
    # When multiple protein or ligand rows exist for the same PDB, we keep the
    # paired numeric bias rows but avoid presenting sequence/ligand identifiers
    # as if each cross-product row were a validated one-to-one reference pair.
    protein_unique_counts = (
        left_lookup.groupby("reference_key")["sequence_protein"]
        .nunique(dropna=True)
        .to_dict()
        if "sequence_protein" in left_lookup.columns
        else {}
    )
    ligand_id_unique_counts = (
        right_lookup.groupby("reference_key")["ligand_id_ligand"]
        .nunique(dropna=True)
        .to_dict()
        if "ligand_id_ligand" in right_lookup.columns
        else {}
    )
    ligand_smiles_unique_counts = (
        right_lookup.groupby("reference_key")["smiles_ligand"]
        .nunique(dropna=True)
        .to_dict()
        if "smiles_ligand" in right_lookup.columns
        else {}
    )

    rows: list[dict[str, object]] = []
    overlap_keys = _direct_overlap_reference_keys(left, right)

    for _, row in left.iterrows():
        reference_key = row.get("reference_key")
        matched_ligands = (
            _matching_reference_rows(right, reference_key)
            if reference_key in overlap_keys
            else right.iloc[0:0].copy()
        )
        if not matched_ligands.empty:
            for _, ligand_match in matched_ligands.iterrows():
                combined_row = _combine_prefixed_rows(row, ligand_match)
                rows.append(
                    _bias_training_row_from_sources(
                        query_protein_chain_id=query_protein_chain_id,
                        query_ligand_chain_id=query_ligand_chain_id,
                        reference_key=reference_key,
                        reference_label=_reference_label_from_merged_row(combined_row),
                        pairing_status=BIAS_PAIRING_PAIRED,
                        protein_row=row,
                        ligand_row=ligand_match,
                        protein_unique_counts=protein_unique_counts,
                        ligand_id_unique_counts=ligand_id_unique_counts,
                        ligand_smiles_unique_counts=ligand_smiles_unique_counts,
                    )
                )
            continue
        rows.append(
            _bias_training_row_from_sources(
                query_protein_chain_id=query_protein_chain_id,
                query_ligand_chain_id=query_ligand_chain_id,
                reference_key=reference_key,
                reference_label=_reference_label_from_merged_row(row),
                pairing_status=BIAS_PAIRING_PROTEIN_ONLY,
                protein_row=row,
                ligand_row=None,
                protein_unique_counts=protein_unique_counts,
                ligand_id_unique_counts=ligand_id_unique_counts,
                ligand_smiles_unique_counts=ligand_smiles_unique_counts,
                forced_ecfp_similarity=(
                    0.0 if query_ligand_chain_id is not None else None
                ),
            )
        )

    for _, row in right.iterrows():
        reference_key = row.get("reference_key")
        if reference_key in overlap_keys:
            continue
        rows.append(
            _bias_training_row_from_sources(
                query_protein_chain_id=query_protein_chain_id,
                query_ligand_chain_id=query_ligand_chain_id,
                reference_key=reference_key,
                reference_label=_reference_label_from_merged_row(row),
                pairing_status=BIAS_PAIRING_LIGAND_ONLY,
                protein_row=None,
                ligand_row=row,
                protein_unique_counts=protein_unique_counts,
                ligand_id_unique_counts=ligand_id_unique_counts,
                ligand_smiles_unique_counts=ligand_smiles_unique_counts,
            )
        )

    return rows


def _build_same_type_pair_rows(
    *,
    query_pair_id: str,
    component_1_id: str,
    component_2_id: str,
    component_type: str,
    component_1_view: pd.DataFrame,
    component_2_view: pd.DataFrame,
    component_1_lookup_view: pd.DataFrame | None = None,
    component_2_lookup_view: pd.DataFrame | None = None,
) -> list[dict[str, object]]:
    is_ligand = component_type == "ligand"
    left = _prepare_reference_side(
        component_1_view, is_ligand=is_ligand, suffix="component_1"
    )
    right = _prepare_reference_side(
        component_2_view, is_ligand=is_ligand, suffix="component_2"
    )
    left_lookup = _prepare_reference_side(
        (
            component_1_lookup_view
            if component_1_lookup_view is not None
            else component_1_view
        ),
        is_ligand=is_ligand,
        suffix="component_1",
    )
    right_lookup = _prepare_reference_side(
        (
            component_2_lookup_view
            if component_2_lookup_view is not None
            else component_2_view
        ),
        is_ligand=is_ligand,
        suffix="component_2",
    )

    rows: list[dict[str, object]] = []
    overlap_keys = _direct_overlap_reference_keys(left, right)

    for _, row in left.iterrows():
        reference_key = row.get("reference_key")
        matches = (
            _matching_reference_rows(right, reference_key)
            if reference_key in overlap_keys
            else (
                right.iloc[0:0].copy()
                if component_type == "ligand"
                else _matching_reference_rows(right_lookup, reference_key)
            )
        )
        if not matches.empty:
            for _, match in matches.iterrows():
                combined_row = _combine_prefixed_rows(row, match)
                rows.append(
                    _same_type_pair_row_from_sources(
                        query_pair_id=query_pair_id,
                        component_1_id=component_1_id,
                        component_2_id=component_2_id,
                        component_type=component_type,
                        reference_key=reference_key,
                        reference_label=_reference_label_from_component_rows(
                            combined_row, combined_row
                        ),
                        pairing_status=BIAS_PAIRING_PAIRED,
                        component_1_row=row,
                        component_2_row=match,
                    )
                )
            continue
        rows.append(
            _same_type_pair_row_from_sources(
                query_pair_id=query_pair_id,
                component_1_id=component_1_id,
                component_2_id=component_2_id,
                component_type=component_type,
                reference_key=reference_key,
                reference_label=_reference_label_from_component_rows(row, None),
                pairing_status=BIAS_PAIRING_COMPONENT_1_ONLY,
                component_1_row=row,
                component_2_row=None,
                forced_component_2_similarity=0.0,
            )
        )

    for _, row in right.iterrows():
        reference_key = row.get("reference_key")
        if reference_key in overlap_keys:
            continue
        matches = (
            left.iloc[0:0].copy()
            if component_type == "ligand"
            else _matching_reference_rows(left_lookup, reference_key)
        )
        if not matches.empty:
            for _, match in matches.iterrows():
                combined_row = _combine_prefixed_rows(match, row)
                rows.append(
                    _same_type_pair_row_from_sources(
                        query_pair_id=query_pair_id,
                        component_1_id=component_1_id,
                        component_2_id=component_2_id,
                        component_type=component_type,
                        reference_key=reference_key,
                        reference_label=_reference_label_from_component_rows(
                            combined_row, combined_row
                        ),
                        pairing_status=BIAS_PAIRING_PAIRED,
                        component_1_row=match,
                        component_2_row=row,
                    )
                )
            continue
        rows.append(
            _same_type_pair_row_from_sources(
                query_pair_id=query_pair_id,
                component_1_id=component_1_id,
                component_2_id=component_2_id,
                component_type=component_type,
                reference_key=reference_key,
                reference_label=_reference_label_from_component_rows(None, row),
                pairing_status=BIAS_PAIRING_COMPONENT_2_ONLY,
                component_1_row=None,
                component_2_row=row,
                forced_component_1_similarity=0.0,
            )
        )

    return rows


def _build_bias_training_dataset(
    *,
    protein_view: pd.DataFrame,
    ligand_views: dict[str, pd.DataFrame],
    protein_lookup_view: pd.DataFrame | None = None,
    ligand_lookup_views: dict[str, pd.DataFrame] | None = None,
) -> pd.DataFrame:
    protein_chain_ids = sorted(
        {
            str(value).strip()
            for value in protein_view.get(
                "query_chain_id", pd.Series(dtype=object)
            ).dropna()
            if str(value).strip()
        }
    )
    ligand_chain_ids = sorted(
        str(chain_id).strip() for chain_id in ligand_views if str(chain_id).strip()
    )

    if not protein_chain_ids and not ligand_chain_ids:
        return pd.DataFrame(columns=BIAS_TRAINING_DATA_COLUMNS)

    rows: list[dict[str, object]] = []
    protein_iter = protein_chain_ids or [None]
    ligand_iter = ligand_chain_ids or [None]

    for protein_chain_id in protein_iter:
        protein_subset = (
            protein_view[
                protein_view["query_chain_id"].astype(str) == protein_chain_id
            ].copy()
            if protein_chain_id is not None and "query_chain_id" in protein_view.columns
            else protein_view.iloc[0:0].copy()
        )
        protein_lookup_subset = (
            protein_lookup_view[
                protein_lookup_view["query_chain_id"].astype(str) == protein_chain_id
            ].copy()
            if (
                protein_lookup_view is not None
                and protein_chain_id is not None
                and "query_chain_id" in protein_lookup_view.columns
            )
            else protein_subset.copy()
        )
        for ligand_chain_id in ligand_iter:
            ligand_subset = (
                ligand_views.get(ligand_chain_id, pd.DataFrame()).copy()
                if ligand_chain_id is not None
                else pd.DataFrame()
            )
            ligand_lookup_subset = (
                ligand_lookup_views.get(ligand_chain_id, pd.DataFrame()).copy()
                if ligand_lookup_views is not None and ligand_chain_id is not None
                else ligand_subset.copy()
            )
            rows.extend(
                _build_bias_training_rows(
                    query_protein_chain_id=protein_chain_id,
                    query_ligand_chain_id=ligand_chain_id,
                    protein_view=protein_subset,
                    ligand_view=ligand_subset,
                    protein_lookup_view=protein_lookup_subset,
                    ligand_lookup_view=ligand_lookup_subset,
                )
            )

    bias_training_df = pd.DataFrame(rows)
    return _finalize_bias_training_dataframe(bias_training_df)


def _pairing_sort_order() -> dict[str, int]:
    return {
        BIAS_PAIRING_PAIRED: 0,
        BIAS_PAIRING_PROTEIN_ONLY: 1,
        BIAS_PAIRING_LIGAND_ONLY: 2,
    }


def _sort_bias_training_dataset(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df.copy()
    ordered = df.copy()
    ordered["_pairing_order"] = (
        ordered["pairing_status"].map(_pairing_sort_order()).fillna(99)
    )
    return ordered.sort_values(
        by=[
            "query_protein_chain_id",
            "query_ligand_chain_id",
            "_pairing_order",
            "source",
            "dataset_name",
            "plot_sequence_similarity",
            "plot_ecfp_similarity",
            "reference_label",
        ],
        ascending=[True, True, True, True, True, False, False, True],
        na_position="last",
    ).drop(columns="_pairing_order")


def _safe_pair_file_token(value: object) -> str:
    text = str(value).strip()
    if not text:
        return "unknown"
    safe = "".join(character if character.isalnum() else "_" for character in text)
    safe = safe.strip("_")
    while "__" in safe:
        safe = safe.replace("__", "_")
    return safe or "unknown"


def _protein_query_groups(
    chain_df: pd.DataFrame,
    protein_queries: dict[str, str | None],
) -> list[dict[str, object]]:
    protein_chain_ids = {
        str(value).strip().upper()
        for value in chain_df.loc[
            chain_df["ENTITY_TYPE"].astype(str) == "protein", "CHAIN_ID"
        ].astype(str)
        if str(value).strip()
    }
    grouped: dict[str, dict[str, object]] = {}
    for raw_chain_id, query_sequence in protein_queries.items():
        chain_id = str(raw_chain_id).strip().upper()
        if not chain_id or chain_id not in protein_chain_ids:
            continue
        identity_key = (
            f"sequence:{str(query_sequence).strip()}"
            if query_sequence
            else f"chain:{chain_id}"
        )
        bucket = grouped.setdefault(
            identity_key,
            {
                "chain_ids": [],
                "label": None,
                "query_value": str(query_sequence) if query_sequence else None,
            },
        )
        bucket["chain_ids"].append(chain_id)

    groups: list[dict[str, object]] = []
    for bucket in grouped.values():
        chain_ids = sorted(set(bucket["chain_ids"]))
        if not chain_ids:
            continue
        groups.append(
            {
                "chain_ids": chain_ids,
                "label": chain_ids[0],
                "query_value": bucket.get("query_value"),
            }
        )
    return sorted(groups, key=lambda group: str(group["label"]))


def _ligand_query_groups(
    chain_df: pd.DataFrame,
    ligand_queries: dict[str, str | None],
) -> list[dict[str, object]]:
    ligand_rows = chain_df[chain_df["ENTITY_TYPE"].astype(str) == "ligand"].copy()
    if ligand_rows.empty:
        return []

    grouped: dict[str, dict[str, object]] = {}
    for _, row in ligand_rows.iterrows():
        chain_id = str(row.get("CHAIN_ID", "")).strip().upper()
        if not chain_id:
            continue
        molecule_id = str(row.get("ligand_molecule_id", "")).strip()
        query_smiles = ligand_queries.get(chain_id)
        if molecule_id and molecule_id != "UNKNOWN_LIGAND":
            identity_key = f"molecule:{molecule_id.upper()}"
            identity_label = molecule_id.upper()
        elif query_smiles and str(query_smiles).strip():
            identity_key = f"smiles:{str(query_smiles).strip()}"
            identity_label = str(query_smiles).strip()
        else:
            identity_key = f"chain:{chain_id}"
            identity_label = chain_id
        bucket = grouped.setdefault(
            identity_key,
            {
                "chain_ids": [],
                "identity_label": identity_label,
                "query_value": str(query_smiles).strip() if query_smiles else None,
            },
        )
        bucket["chain_ids"].append(chain_id)

    groups: list[dict[str, object]] = []
    for bucket in grouped.values():
        chain_ids = sorted(set(bucket["chain_ids"]))
        if not chain_ids:
            continue
        label = chain_ids[0]
        if len(chain_ids) > 1:
            identity_label = str(bucket.get("identity_label", "")).strip()
            if identity_label and all(
                character.isalnum() or character in {"_", "-"}
                for character in identity_label
            ):
                label = identity_label
        groups.append(
            {
                "chain_ids": chain_ids,
                "label": label,
                "query_value": bucket.get("query_value"),
            }
        )
    return sorted(groups, key=lambda group: str(group["label"]))


def _all_query_groups(
    chain_df: pd.DataFrame,
    protein_queries: dict[str, str | None],
    ligand_queries: dict[str, str | None],
) -> list[dict[str, object]]:
    groups: list[dict[str, object]] = []
    for group in _protein_query_groups(chain_df, protein_queries):
        groups.append(
            {
                "component_type": "protein",
                "chain_ids": group["chain_ids"],
                "label": group["label"],
                "query_value": group.get("query_value"),
            }
        )
    for group in _ligand_query_groups(chain_df, ligand_queries):
        groups.append(
            {
                "component_type": "ligand",
                "chain_ids": group["chain_ids"],
                "label": group["label"],
                "query_value": group.get("query_value"),
            }
        )
    return groups


def _dedupe_query_group_view(view: pd.DataFrame) -> pd.DataFrame:
    if view.empty:
        return view.copy()
    subset = [column for column in view.columns if column != "query_chain_id"]
    if not subset:
        return view.copy()
    return view.drop_duplicates(subset=subset).reset_index(drop=True)


def _subset_group_protein_view(
    *,
    view: pd.DataFrame,
    chain_ids: set[str],
) -> pd.DataFrame:
    if view.empty or "query_chain_id" not in view.columns:
        return view.iloc[0:0].copy()
    subset = view[view["query_chain_id"].astype(str).str.upper().isin(chain_ids)].copy()
    return _dedupe_query_group_view(subset)


def _subset_group_ligand_view(
    *,
    views: dict[str, pd.DataFrame],
    chain_ids: set[str],
) -> pd.DataFrame:
    frames = [
        views.get(chain_id, pd.DataFrame()).copy() for chain_id in sorted(chain_ids)
    ]
    frames = [frame for frame in frames if not frame.empty]
    if not frames:
        return pd.DataFrame()
    subset = pd.concat(frames, ignore_index=True)
    return _dedupe_query_group_view(subset)


def _combined_mixed_pair_bias_training_dataframe(
    pair_artifacts: list[dict[str, object]],
) -> pd.DataFrame:
    mixed_frames = [
        artifact["dataframe"]
        for artifact in pair_artifacts
        if artifact.get("pair_type") == "mixed"
        and isinstance(artifact.get("dataframe"), pd.DataFrame)
        and not artifact["dataframe"].empty
    ]
    if not mixed_frames:
        return pd.DataFrame(columns=BIAS_TRAINING_DATA_COLUMNS)
    return _finalize_bias_training_dataframe(pd.concat(mixed_frames, ignore_index=True))


__all__ = [
    "_all_query_groups",
    "_bias_training_row_from_sources",
    "_build_bias_training_dataset",
    "_build_bias_training_rows",
    "_build_ligand_reference_index",
    "_build_protein_lookup_index",
    "_build_same_type_pair_rows",
    "_collapse_pair_value",
    "_combine_prefixed_rows",
    "_combined_mixed_pair_bias_training_dataframe",
    "_component_similarity_column",
    "_component_threshold",
    "_dedupe_query_group_view",
    "_direct_overlap_reference_keys",
    "_finalize_bias_training_dataframe",
    "_finalize_same_type_pair_dataframe",
    "_ligand_query_groups",
    "_matching_reference_rows",
    "_normalize_component_similarity",
    "_pairing_sort_order",
    "_prepare_reference_side",
    "_protein_query_groups",
    "_reference_key_from_row",
    "_reference_label_from_component_rows",
    "_reference_label_from_merged_row",
    "_reference_ligand_similarity_rows",
    "_reference_pdb_id_from_bias_row",
    "_reference_pdb_id_from_same_type_row",
    "_reference_protein_similarity_rows",
    "_safe_pair_file_token",
    "_same_type_pair_row_from_sources",
    "_same_type_pair_sort_order",
    "_sort_bias_training_dataset",
    "_sort_same_type_pair_dataset",
    "_subset_group_ligand_view",
    "_subset_group_protein_view",
]
