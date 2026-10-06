"""Bias dataset and plot artifact materialization."""

from __future__ import annotations

from pathlib import Path
from itertools import combinations
import pandas as pd
from cofolder.modules.analytics.plots import plot_bias_reference_overlap
from cofolder.modules.analytics.plots import plot_reference_overlap_scatter

from ._common import (
    BIAS_LIGAND_PAIR_FILE_STEM,
    BIAS_PLOT_FILE_STEM,
    BIAS_PLOT_LIGAND_THRESHOLD,
    BIAS_PLOT_SEQUENCE_THRESHOLD,
    BIAS_PROTEIN_PAIR_FILE_STEM,
    _order_ligand_training_columns,
)

from ._datasets import (
    _all_query_groups,
    _build_bias_training_rows,
    _build_same_type_pair_rows,
    _component_threshold,
    _finalize_same_type_pair_dataframe,
    _safe_pair_file_token,
    _subset_group_ligand_view,
    _subset_group_protein_view,
)

from ._similarity import (
    _best_reference_row,
    _build_ligand_training_views,
    _build_protein_training_view,
)

from ._enrichment import (
    _enrich_mixed_bias_dataset_with_pdb_backfill,
    _enrich_same_type_ligand_pair_dataset_with_pdb_backfill,
)


def _materialize_reference_landscape_summary(
    *,
    chain_df: pd.DataFrame,
    proteins_df: pd.DataFrame,
    ligands_df: pd.DataFrame,
    protein_queries: dict[str, str | None],
    ligand_queries: dict[str, str | None],
    output_dir: Path,
) -> None:
    summary_rows: list[dict[str, object]] = []

    protein_chain_ids = set(
        chain_df.loc[
            chain_df["ENTITY_TYPE"].astype(str) == "protein", "CHAIN_ID"
        ].astype(str)
    )
    for chain_id, query_seq in protein_queries.items():
        if str(chain_id) not in protein_chain_ids or not query_seq:
            continue
        public_df = proteins_df[proteins_df["source"].astype(str) == "public"].copy()
        custom_df = proteins_df[proteins_df["source"].astype(str) != "public"].copy()
        public_best = _best_reference_row(public_df, query_seq, is_ligand=False)
        custom_best = _best_reference_row(custom_df, query_seq, is_ligand=False)
        overall_best = _best_reference_row(proteins_df, query_seq, is_ligand=False)
        summary_rows.append(
            _build_reference_landscape_row(
                chain_id=str(chain_id),
                entity_type="protein",
                public_best=public_best,
                custom_best=custom_best,
                overall_best=overall_best,
            )
        )

    ligand_chain_ids = set(
        chain_df.loc[
            chain_df["ENTITY_TYPE"].astype(str) == "ligand", "CHAIN_ID"
        ].astype(str)
    )
    for chain_id, query_smiles in ligand_queries.items():
        if str(chain_id) not in ligand_chain_ids or not query_smiles:
            continue
        ligand_rows = chain_df[chain_df["CHAIN_ID"].astype(str) == str(chain_id)].copy()
        mol_id = None
        if "ligand_molecule_id" in ligand_rows.columns and not ligand_rows.empty:
            raw_mol_id = ligand_rows["ligand_molecule_id"].iloc[0]
            if pd.notna(raw_mol_id) and str(raw_mol_id).strip():
                mol_id = str(raw_mol_id)
        public_df = ligands_df[ligands_df["source"].astype(str) == "public"].copy()
        custom_df = ligands_df[ligands_df["source"].astype(str) != "public"].copy()
        public_best = _best_reference_row(
            public_df, query_smiles, is_ligand=True, mol_id=mol_id
        )
        custom_best = _best_reference_row(
            custom_df, query_smiles, is_ligand=True, mol_id=mol_id
        )
        overall_best = _best_reference_row(
            ligands_df, query_smiles, is_ligand=True, mol_id=mol_id
        )
        summary_rows.append(
            _build_reference_landscape_row(
                chain_id=str(chain_id),
                entity_type="ligand",
                public_best=public_best,
                custom_best=custom_best,
                overall_best=overall_best,
            )
        )

    pd.DataFrame(summary_rows).to_csv(
        output_dir / "reference_landscape_summary.csv", index=False
    )


def _build_reference_landscape_row(
    *,
    chain_id: str,
    entity_type: str,
    public_best: pd.Series | None,
    custom_best: pd.Series | None,
    overall_best: pd.Series | None,
) -> dict[str, object]:
    nearest_overall_source = (
        str(overall_best.get("source")).strip()
        if overall_best is not None and pd.notna(overall_best.get("source"))
        else None
    )
    public_similarity = (
        float(public_best["best_similarity"])
        if public_best is not None and pd.notna(public_best.get("best_similarity"))
        else None
    )
    custom_similarity = (
        float(custom_best["best_similarity"])
        if custom_best is not None and pd.notna(custom_best.get("best_similarity"))
        else None
    )
    return {
        "query_chain_id": chain_id,
        "entity_type": entity_type,
        "nearest_public_pdb_id": (
            public_best.get("pdb_id") if public_best is not None else None
        ),
        "nearest_public_dataset_name": (
            public_best.get("dataset_name") if public_best is not None else None
        ),
        "nearest_public_similarity": public_similarity,
        "nearest_public_similarity_method": (
            public_best.get("best_similarity_method")
            if public_best is not None
            else None
        ),
        "nearest_custom_pdb_id": (
            custom_best.get("pdb_id") if custom_best is not None else None
        ),
        "nearest_custom_dataset_name": (
            custom_best.get("dataset_name") if custom_best is not None else None
        ),
        "nearest_custom_similarity": custom_similarity,
        "nearest_custom_similarity_method": (
            custom_best.get("best_similarity_method")
            if custom_best is not None
            else None
        ),
        "nearest_overall_source": nearest_overall_source,
        "nearest_overall_pdb_id": (
            overall_best.get("pdb_id") if overall_best is not None else None
        ),
        "nearest_overall_dataset_name": (
            overall_best.get("dataset_name") if overall_best is not None else None
        ),
        "nearest_overall_similarity": (
            float(overall_best["best_similarity"])
            if overall_best is not None
            and pd.notna(overall_best.get("best_similarity"))
            else None
        ),
        "nearest_overall_similarity_method": (
            overall_best.get("best_similarity_method")
            if overall_best is not None
            else None
        ),
        "custom_changes_nearest_reference": bool(
            nearest_overall_source == "custom"
            and (
                public_similarity is None
                or (
                    custom_similarity is not None
                    and custom_similarity > public_similarity
                )
            )
        ),
    }


def _materialize_bias_training_views(
    chain_df: pd.DataFrame,
    ligands_df: pd.DataFrame,
    proteins_df: pd.DataFrame,
    protein_queries: dict[str, str | None],
    ligand_queries: dict[str, str | None],
    boltz_cache_path: Path,
    output_dir: Path,
    components_cif_path: Path | None = None,
) -> tuple[
    pd.DataFrame, dict[str, pd.DataFrame], pd.DataFrame, dict[str, pd.DataFrame]
]:
    """Write bias training files for debug/inspection under results/bias_train."""
    output_dir.mkdir(parents=True, exist_ok=True)
    for stale_path in output_dir.glob("ligand_training_data_*.csv"):
        stale_path.unlink()

    proteins_out = _build_protein_training_view(proteins_df, protein_queries)
    protein_lookup_out = _build_protein_training_view(
        proteins_df,
        protein_queries,
        minimum_similarity=None,
        top_n=None,
    )
    proteins_out.to_csv(output_dir / "protein_training_data.csv", index=False)
    _order_ligand_training_columns(ligands_df).to_csv(
        output_dir / "ligand_training_data.csv",
        index=False,
    )

    ligand_views = _build_ligand_training_views(
        chain_df=chain_df,
        ligands_df=ligands_df,
        ligand_queries=ligand_queries,
        boltz_cache_path=boltz_cache_path,
        components_cif_path=components_cif_path,
    )
    ligand_lookup_views = _build_ligand_training_views(
        chain_df=chain_df,
        ligands_df=ligands_df,
        ligand_queries=ligand_queries,
        boltz_cache_path=boltz_cache_path,
        components_cif_path=components_cif_path,
        minimum_similarity=None,
        fallback_top_n=None,
    )

    lig_rows = chain_df[chain_df["ENTITY_TYPE"].astype(str) == "ligand"].copy()
    if not lig_rows.empty and "ligand_molecule_id" not in lig_rows.columns:
        return proteins_out, ligand_views, protein_lookup_out, ligand_lookup_views

    for chain_id, ligand_view in ligand_views.items():
        out_path = output_dir / f"ligand_training_data_{chain_id}.csv"
        ligand_view.drop(columns=["query_chain_id"], errors="ignore").to_csv(
            out_path, index=False
        )

    return proteins_out, ligand_views, protein_lookup_out, ligand_lookup_views


def _skip_plot_message() -> str:
    return (
        "Bias reference-overlap scatter plot skipped because the merged plotting dataset "
        "did not contain rows with both protein and ligand similarity axes. "
        "Bias artifacts remain reference-overlap diagnostics only and do not predict "
        "affinity, structural confidence, or cofolding success.\n"
    )


def _write_bias_plot_artifacts(
    *,
    bias_training_df: pd.DataFrame,
    output_dir: Path,
    file_stem: str,
    sequence_threshold: float = BIAS_PLOT_SEQUENCE_THRESHOLD,
    ligand_threshold: float = BIAS_PLOT_LIGAND_THRESHOLD,
) -> list[Path]:
    skipped_plot_path = output_dir / f"{file_stem}.skipped.txt"
    plot_paths = plot_bias_reference_overlap(
        bias_training_df,
        output_dir=output_dir,
        file_stem=file_stem,
        sequence_threshold=sequence_threshold,
        ligand_threshold=ligand_threshold,
    )
    if plot_paths:
        if skipped_plot_path.exists():
            skipped_plot_path.unlink()
        return plot_paths

    for suffix in ("png", "pdf"):
        (output_dir / f"{file_stem}.{suffix}").unlink(missing_ok=True)
    skipped_plot_path.write_text(_skip_plot_message(), encoding="utf-8")
    return []


def _write_same_type_plot_artifacts(
    *,
    pair_df: pd.DataFrame,
    output_dir: Path,
    file_stem: str,
    component_type: str,
    threshold: float | None = None,
) -> list[Path]:
    skipped_plot_path = output_dir / f"{file_stem}.skipped.txt"
    threshold = (
        _component_threshold(component_type) if threshold is None else float(threshold)
    )
    component_label = "Protein" if component_type == "protein" else "Ligand"
    plot_paths = plot_reference_overlap_scatter(
        pair_df,
        output_dir=output_dir,
        file_stem=file_stem,
        x_col="plot_component_1_similarity",
        y_col="plot_component_2_similarity",
        x_threshold=threshold,
        y_threshold=threshold,
        x_label=f"{component_label} reference overlap (component 1)",
        y_label=f"{component_label} reference overlap (component 2)",
        title=f"Bias {component_type}-{component_type} reference-overlap diagnostic",
        query_1_col="component_1_id",
        query_2_col="component_2_id",
    )
    if plot_paths:
        if skipped_plot_path.exists():
            skipped_plot_path.unlink()
        return plot_paths

    for suffix in ("png", "pdf"):
        (output_dir / f"{file_stem}.{suffix}").unlink(missing_ok=True)
    skipped_plot_path.write_text(_skip_plot_message(), encoding="utf-8")
    return []


def _clear_pair_specific_artifacts(output_dir: Path) -> None:
    patterns = (
        "bias_training_data_*.csv",
        f"{BIAS_PLOT_FILE_STEM}_*.png",
        f"{BIAS_PLOT_FILE_STEM}_*.pdf",
        f"{BIAS_PLOT_FILE_STEM}_*.skipped.txt",
        "bias_protein_pair_data_*.csv",
        f"{BIAS_PROTEIN_PAIR_FILE_STEM}_*.png",
        f"{BIAS_PROTEIN_PAIR_FILE_STEM}_*.pdf",
        f"{BIAS_PROTEIN_PAIR_FILE_STEM}_*.skipped.txt",
        "bias_ligand_pair_data_*.csv",
        f"{BIAS_LIGAND_PAIR_FILE_STEM}_*.png",
        f"{BIAS_LIGAND_PAIR_FILE_STEM}_*.pdf",
        f"{BIAS_LIGAND_PAIR_FILE_STEM}_*.skipped.txt",
    )
    for pattern in patterns:
        for stale_path in output_dir.glob(pattern):
            stale_path.unlink()


def _mixed_pair_dataset(
    *,
    protein_label: str,
    ligand_label: str,
    protein_query_sequence: str | None,
    ligand_query_smiles: str | None,
    protein_view: pd.DataFrame,
    ligand_view: pd.DataFrame,
    protein_lookup_view: pd.DataFrame,
    ligand_lookup_view: pd.DataFrame,
    boltz_cache_path: Path,
    components_cif_path: Path | None,
    ligand_reference_df: pd.DataFrame | None = None,
    ligand_reference_index: dict[str, list[dict[str, str]]] | None = None,
    progress_output_path: Path | None = None,
    progress_label: str | None = None,
) -> pd.DataFrame:
    rows = _build_bias_training_rows(
        query_protein_chain_id=protein_label,
        query_ligand_chain_id=ligand_label,
        protein_view=protein_view,
        ligand_view=ligand_view,
        protein_lookup_view=protein_lookup_view,
        ligand_lookup_view=ligand_lookup_view,
    )
    return _enrich_mixed_bias_dataset_with_pdb_backfill(
        pd.DataFrame(rows),
        protein_queries={protein_label: protein_query_sequence},
        ligand_queries={ligand_label: ligand_query_smiles},
        protein_lookup_view=protein_lookup_view,
        boltz_cache_path=boltz_cache_path,
        components_cif_path=components_cif_path,
        ligand_reference_df=ligand_reference_df,
        ligand_reference_index=ligand_reference_index,
        progress_output_path=progress_output_path,
        progress_label=progress_label,
    )


def _same_type_pair_dataset(
    *,
    component_1_label: str,
    component_2_label: str,
    component_type: str,
    component_1_view: pd.DataFrame,
    component_2_view: pd.DataFrame,
    component_1_lookup_view: pd.DataFrame,
    component_2_lookup_view: pd.DataFrame,
    ligand_queries: dict[str, str | None] | None = None,
    boltz_cache_path: Path | None = None,
    components_cif_path: Path | None = None,
    ligand_reference_df: pd.DataFrame | None = None,
    ligand_reference_index: dict[str, list[dict[str, str]]] | None = None,
    progress_output_path: Path | None = None,
    progress_label: str | None = None,
) -> pd.DataFrame:
    rows = _build_same_type_pair_rows(
        query_pair_id=f"{component_1_label}__{component_2_label}",
        component_1_id=component_1_label,
        component_2_id=component_2_label,
        component_type=component_type,
        component_1_view=component_1_view,
        component_2_view=component_2_view,
        component_1_lookup_view=component_1_lookup_view,
        component_2_lookup_view=component_2_lookup_view,
    )
    pair_df = pd.DataFrame(rows)
    if component_type == "ligand":
        pair_df = _enrich_same_type_ligand_pair_dataset_with_pdb_backfill(
            pair_df,
            ligand_queries=ligand_queries or {},
            boltz_cache_path=(
                boltz_cache_path
                if boltz_cache_path is not None
                else Path("~/.boltz").expanduser()
            ),
            components_cif_path=components_cif_path,
            ligand_reference_df=ligand_reference_df,
            ligand_reference_index=ligand_reference_index,
            progress_output_path=progress_output_path,
            progress_label=progress_label,
        )
    return _finalize_same_type_pair_dataframe(pair_df)


def _pair_specific_bias_training_datasets(
    *,
    chain_df: pd.DataFrame,
    protein_queries: dict[str, str | None],
    ligand_queries: dict[str, str | None],
    ligand_reference_df: pd.DataFrame,
    ligand_reference_index: dict[str, list[dict[str, str]]] | None,
    protein_view: pd.DataFrame,
    ligand_views: dict[str, pd.DataFrame],
    protein_lookup_view: pd.DataFrame,
    ligand_lookup_views: dict[str, pd.DataFrame],
    boltz_cache_path: Path,
    components_cif_path: Path | None,
    output_dir: Path,
) -> list[dict[str, object]]:
    all_groups = _all_query_groups(chain_df, protein_queries, ligand_queries)
    datasets: list[dict[str, object]] = []

    for component_1, component_2 in combinations(all_groups, 2):
        component_1_type = str(component_1["component_type"])
        component_2_type = str(component_2["component_type"])
        component_1_label = str(component_1["label"])
        component_2_label = str(component_2["label"])
        component_1_chain_ids = set(component_1["chain_ids"])
        component_2_chain_ids = set(component_2["chain_ids"])

        if component_1_type != component_2_type:
            protein_group = (
                component_1 if component_1_type == "protein" else component_2
            )
            ligand_group = component_2 if component_2_type == "ligand" else component_1
            protein_label = str(protein_group["label"])
            ligand_label = str(ligand_group["label"])
            pair_suffix = f"{_safe_pair_file_token(protein_label)}__{_safe_pair_file_token(ligand_label)}"
            data_file = f"bias_training_data_{pair_suffix}.csv"
            pair_df = _mixed_pair_dataset(
                protein_label=protein_label,
                ligand_label=ligand_label,
                protein_query_sequence=(
                    str(protein_group.get("query_value"))
                    if protein_group.get("query_value") is not None
                    else None
                ),
                ligand_query_smiles=(
                    str(ligand_group.get("query_value"))
                    if ligand_group.get("query_value") is not None
                    else None
                ),
                protein_view=_subset_group_protein_view(
                    view=protein_view,
                    chain_ids=set(protein_group["chain_ids"]),
                ),
                ligand_view=_subset_group_ligand_view(
                    views=ligand_views,
                    chain_ids=set(ligand_group["chain_ids"]),
                ),
                protein_lookup_view=_subset_group_protein_view(
                    view=protein_lookup_view,
                    chain_ids=set(protein_group["chain_ids"]),
                ),
                ligand_lookup_view=_subset_group_ligand_view(
                    views=ligand_lookup_views,
                    chain_ids=set(ligand_group["chain_ids"]),
                ),
                boltz_cache_path=boltz_cache_path,
                components_cif_path=components_cif_path,
                ligand_reference_df=ligand_reference_df,
                ligand_reference_index=ligand_reference_index,
                progress_output_path=output_dir / data_file,
                progress_label=f"{protein_label}__{ligand_label}",
            )
            datasets.append(
                {
                    "pair_type": "mixed",
                    "pair_label": f"{protein_label}__{ligand_label}",
                    "pair_suffix": pair_suffix,
                    "dataframe": pair_df,
                    "data_file": data_file,
                    "plot_file_stem": f"{BIAS_PLOT_FILE_STEM}_{pair_suffix}",
                }
            )
            continue

        pair_suffix = f"{_safe_pair_file_token(component_1_label)}__{_safe_pair_file_token(component_2_label)}"
        component_1_view = (
            _subset_group_protein_view(
                view=protein_view, chain_ids=component_1_chain_ids
            )
            if component_1_type == "protein"
            else _subset_group_ligand_view(
                views=ligand_views, chain_ids=component_1_chain_ids
            )
        )
        component_2_view = (
            _subset_group_protein_view(
                view=protein_view, chain_ids=component_2_chain_ids
            )
            if component_2_type == "protein"
            else _subset_group_ligand_view(
                views=ligand_views, chain_ids=component_2_chain_ids
            )
        )
        component_1_lookup = (
            _subset_group_protein_view(
                view=protein_lookup_view, chain_ids=component_1_chain_ids
            )
            if component_1_type == "protein"
            else _subset_group_ligand_view(
                views=ligand_lookup_views, chain_ids=component_1_chain_ids
            )
        )
        component_2_lookup = (
            _subset_group_protein_view(
                view=protein_lookup_view, chain_ids=component_2_chain_ids
            )
            if component_2_type == "protein"
            else _subset_group_ligand_view(
                views=ligand_lookup_views, chain_ids=component_2_chain_ids
            )
        )

        if component_1_type == "protein":
            data_file = f"bias_protein_pair_data_{pair_suffix}.csv"
            plot_file_stem = f"{BIAS_PROTEIN_PAIR_FILE_STEM}_{pair_suffix}"
            pair_type = "protein_pair"
            progress_output_path = None
        else:
            data_file = f"bias_ligand_pair_data_{pair_suffix}.csv"
            plot_file_stem = f"{BIAS_LIGAND_PAIR_FILE_STEM}_{pair_suffix}"
            pair_type = "ligand_pair"
            progress_output_path = output_dir / data_file

        pair_df = _same_type_pair_dataset(
            component_1_label=component_1_label,
            component_2_label=component_2_label,
            component_type=component_1_type,
            component_1_view=component_1_view,
            component_2_view=component_2_view,
            component_1_lookup_view=component_1_lookup,
            component_2_lookup_view=component_2_lookup,
            ligand_queries=ligand_queries,
            boltz_cache_path=boltz_cache_path,
            components_cif_path=components_cif_path,
            ligand_reference_df=ligand_reference_df,
            ligand_reference_index=ligand_reference_index,
            progress_output_path=progress_output_path,
            progress_label=f"{component_1_label}__{component_2_label}",
        )
        datasets.append(
            {
                "pair_type": pair_type,
                "pair_label": f"{component_1_label}__{component_2_label}",
                "pair_suffix": pair_suffix,
                "dataframe": pair_df,
                "data_file": data_file,
                "plot_file_stem": plot_file_stem,
            }
        )

    return datasets


__all__ = [
    "_build_reference_landscape_row",
    "_clear_pair_specific_artifacts",
    "_materialize_bias_training_views",
    "_materialize_reference_landscape_summary",
    "_mixed_pair_dataset",
    "_pair_specific_bias_training_datasets",
    "_same_type_pair_dataset",
    "_skip_plot_message",
    "_write_bias_plot_artifacts",
    "_write_same_type_plot_artifacts",
]
