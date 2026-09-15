"""Pure protein and ligand similarity calculations for bias analytics."""

from __future__ import annotations

from rdkit import Chem
from rdkit import DataStructs
from pathlib import Path
from functools import lru_cache
import pandas as pd

from ._common import (
    ALIGNER,
    BIAS_LIGAND_VIEW_MIN_SIMILARITY,
    BIAS_PROTEIN_VIEW_MIN_SIMILARITY,
    MORGAN_FP,
    _norm_id,
    _normalize_smiles,
    _order_ligand_training_columns,
    _order_protein_training_columns,
)

from ._references import (
    _smiles_from_ccd_cache,
    _smiles_from_components_cif,
)


def _pairwise_sequence_identity_percent(query: str, target: str) -> float:
    if not query or not target:
        return 0.0
    score = ALIGNER.score(query, target)
    denom = max(len(query), len(target))
    if denom == 0:
        return 0.0
    return float(100.0 * float(score) / float(denom))


def _best_protein_hit(query_seq: str, proteins_df: pd.DataFrame) -> float | None:
    if "sequence_similarity" in proteins_df.columns:
        sims = pd.to_numeric(
            proteins_df["sequence_similarity"], errors="coerce"
        ).dropna()
        if len(sims):
            return float(sims.max())
    return None


def _best_pairwise_protein_hit(
    query_seq: str, proteins_df: pd.DataFrame
) -> float | None:
    if "sequence_similarity_pairwise" in proteins_df.columns:
        sims = pd.to_numeric(
            proteins_df["sequence_similarity_pairwise"], errors="coerce"
        ).dropna()
        if len(sims):
            return float(sims.max())
    best_score = None
    for _, row in proteins_df.iterrows():
        if pd.notna(row.get("sequence_similarity")):
            continue
        score = _pairwise_sequence_identity_percent(query_seq, str(row["sequence"]))
        if best_score is None or score > best_score:
            best_score = score
    return best_score


def _morgan_fp_from_smiles(smiles: str):
    normalized_smiles = _normalize_smiles(smiles)
    if not normalized_smiles:
        return None
    return _morgan_fp_from_normalized_smiles(normalized_smiles)


@lru_cache(maxsize=65536)
def _morgan_fp_from_normalized_smiles(normalized_smiles: str):
    mol = Chem.MolFromSmiles(normalized_smiles)
    if mol is None:
        return None
    return MORGAN_FP.GetFingerprint(mol)


def _best_ligand_hit(query_smiles: str, ligands_df: pd.DataFrame) -> float | None:
    query_fp = _morgan_fp_from_smiles(query_smiles)
    if query_fp is None:
        return None

    best_sim = None
    for _, row in ligands_df.iterrows():
        fp = _morgan_fp_from_smiles(str(row["smiles"]))
        if fp is None:
            continue
        sim = float(DataStructs.TanimotoSimilarity(query_fp, fp))
        if best_sim is None or sim > best_sim:
            best_sim = sim
    return best_sim


def _reference_rows_for_query(frame: pd.DataFrame, chain_id: str) -> pd.DataFrame:
    """Select chain-specific rows plus reusable rows without a query identity."""
    if "query_chain_id" not in frame.columns:
        return frame.copy()
    identities = frame["query_chain_id"].fillna("").astype(str).str.strip().str.upper()
    selected = identities.eq(str(chain_id).strip().upper()) | identities.eq("")
    return frame.loc[selected].copy()


def _protein_similarity_series(query_seq: str, proteins_df: pd.DataFrame) -> pd.Series:
    if proteins_df.empty:
        return pd.Series(dtype=float)
    if "sequence_similarity" in proteins_df.columns:
        return pd.to_numeric(proteins_df["sequence_similarity"], errors="coerce")
    return pd.Series(pd.NA, index=proteins_df.index, dtype="Float64")


def _protein_pairwise_similarity_series(
    query_seq: str,
    proteins_df: pd.DataFrame,
) -> pd.Series:
    if proteins_df.empty:
        return pd.Series(dtype=float)
    if "sequence_similarity_pairwise" in proteins_df.columns:
        pairwise = pd.to_numeric(
            proteins_df["sequence_similarity_pairwise"], errors="coerce"
        )
    else:
        pairwise = pd.Series(pd.NA, index=proteins_df.index, dtype="Float64")
    mmseqs = _protein_similarity_series(query_seq, proteins_df)
    missing = pairwise.isna() & mmseqs.isna()
    if missing.any() and "sequence" in proteins_df.columns:
        pairwise.loc[missing] = proteins_df.loc[missing, "sequence"].apply(
            lambda sequence: _pairwise_sequence_identity_percent(
                query_seq, str(sequence)
            )
        )
    return pd.to_numeric(pairwise, errors="coerce")


def _ligand_similarity_series(
    query_smiles: str,
    ligands_df: pd.DataFrame,
    *,
    mol_id: str | None = None,
) -> pd.Series:
    if ligands_df.empty:
        return pd.Series(dtype=float)
    query_fp = _morgan_fp_from_smiles(query_smiles)
    if query_fp is None:
        return pd.Series(
            [pd.NA] * len(ligands_df), index=ligands_df.index, dtype="Float64"
        )

    similarities = []
    for _, row in ligands_df.iterrows():
        row_smiles = str(row["smiles"])
        row_ligand_id = (
            _norm_id(row["ligand_id"])
            if pd.notna(row.get("ligand_id")) and str(row.get("ligand_id")).strip()
            else None
        )
        if row_smiles == str(query_smiles):
            similarities.append(1.0)
            continue
        if mol_id and row_ligand_id == _norm_id(mol_id):
            precomputed = pd.to_numeric(
                pd.Series([row.get("ecfp_similarity")]), errors="coerce"
            ).iloc[0]
            if pd.notna(precomputed):
                similarities.append(float(precomputed))
                continue
        fp = _morgan_fp_from_smiles(str(row["smiles"]))
        if fp is None:
            similarities.append(pd.NA)
        else:
            similarities.append(float(DataStructs.TanimotoSimilarity(query_fp, fp)))
    return pd.to_numeric(
        pd.Series(similarities, index=ligands_df.index), errors="coerce"
    )


def _best_reference_row(
    df: pd.DataFrame,
    query_value: str | None,
    *,
    is_ligand: bool,
    mol_id: str | None = None,
) -> pd.Series | None:
    if df.empty or not query_value:
        return None
    if is_ligand:
        similarities = _ligand_similarity_series(query_value, df, mol_id=mol_id)
        similarity_methods = pd.Series("ecfp4_tanimoto", index=df.index, dtype="string")
    else:
        mmseqs_similarities = _protein_similarity_series(query_value, df)
        pairwise_similarities = _protein_pairwise_similarity_series(query_value, df)
        similarities = mmseqs_similarities.combine_first(pairwise_similarities)
        similarity_methods = pd.Series(pd.NA, index=df.index, dtype="string")
        similarity_methods.loc[mmseqs_similarities.notna()] = "mmseqs_pident"
        similarity_methods.loc[
            mmseqs_similarities.isna() & pairwise_similarities.notna()
        ] = "pairwise_aligner"
    if similarities.dropna().empty:
        return None
    best_idx = similarities.fillna(float("-inf")).idxmax()
    best_row = df.loc[best_idx].copy()
    best_row["best_similarity"] = float(similarities.loc[best_idx])
    best_row["best_similarity_method"] = similarity_methods.loc[best_idx]
    return best_row


def _build_protein_training_view(
    proteins_df: pd.DataFrame,
    protein_queries: dict[str, str | None],
    *,
    minimum_similarity: float | None = BIAS_PROTEIN_VIEW_MIN_SIMILARITY,
    top_n: int | None = 100,
) -> pd.DataFrame:
    prot_rows: list[pd.DataFrame] = []
    for chain_id, query_seq in protein_queries.items():
        if not query_seq:
            continue
        sub = _reference_rows_for_query(proteins_df, chain_id)
        sub["sequence_similarity"] = _protein_similarity_series(str(query_seq), sub)
        sub["sequence_similarity_pairwise"] = _protein_pairwise_similarity_series(
            str(query_seq), sub
        )
        sub["sequence_similarity"] = pd.to_numeric(
            sub["sequence_similarity"], errors="coerce"
        ).clip(
            lower=0.0,
            upper=100.0,
        )
        sub["sequence_similarity_pairwise"] = pd.to_numeric(
            sub["sequence_similarity_pairwise"], errors="coerce"
        ).clip(lower=0.0, upper=100.0)
        sub["sequence_similarity_method"] = pd.Series(
            pd.NA, index=sub.index, dtype="string"
        )
        sub.loc[sub["sequence_similarity"].notna(), "sequence_similarity_method"] = (
            "mmseqs_pident"
        )
        sub.loc[
            sub["sequence_similarity"].isna()
            & sub["sequence_similarity_pairwise"].notna(),
            "sequence_similarity_method",
        ] = "pairwise_aligner"
        sub["_effective_sequence_similarity"] = sub["sequence_similarity"].copy()
        missing_effective_similarity = sub["_effective_sequence_similarity"].isna()
        sub.loc[missing_effective_similarity, "_effective_sequence_similarity"] = (
            sub.loc[missing_effective_similarity, "sequence_similarity_pairwise"]
        )
        if minimum_similarity is not None:
            sub = sub[sub["_effective_sequence_similarity"] > minimum_similarity].copy()
        sub = sub.sort_values("_effective_sequence_similarity", ascending=False)
        if top_n is not None:
            sub = sub.head(top_n).copy()
        sub = sub.drop(columns="_effective_sequence_similarity")
        sub["query_chain_id"] = str(chain_id).strip()
        prot_rows.append(sub)

    if prot_rows:
        proteins_out = pd.concat(prot_rows, ignore_index=True)
        proteins_out = proteins_out.drop_duplicates(
            subset=[
                "pdb_id",
                "release_date",
                "sequence",
                "query_chain_id",
                "source",
                "dataset_name",
                "source_structure_path",
                "source_reference_path",
            ]
        )
        proteins_out = proteins_out.sort_values(
            by=[
                "query_chain_id",
                "sequence_similarity",
                "sequence_similarity_pairwise",
                "source",
                "pdb_id",
            ],
            ascending=[True, False, False, True, True],
            na_position="last",
        )
    else:
        proteins_out = proteins_df.iloc[0:0].copy()
        proteins_out["sequence_similarity"] = pd.Series(dtype=float)
        proteins_out["sequence_similarity_pairwise"] = pd.Series(dtype=float)
        proteins_out["sequence_similarity_method"] = pd.Series(dtype="string")
        proteins_out["query_chain_id"] = pd.Series(dtype=str)

    return _order_protein_training_columns(proteins_out)


def _build_ligand_training_views(
    chain_df: pd.DataFrame,
    ligands_df: pd.DataFrame,
    ligand_queries: dict[str, str | None],
    boltz_cache_path: Path,
    components_cif_path: Path | None = None,
    *,
    minimum_similarity: float | None = BIAS_LIGAND_VIEW_MIN_SIMILARITY,
    fallback_top_n: int | None = 100,
) -> dict[str, pd.DataFrame]:
    lig_rows = chain_df[chain_df["ENTITY_TYPE"].astype(str) == "ligand"].copy()
    if lig_rows.empty or "ligand_molecule_id" not in lig_rows.columns:
        return {}

    views: dict[str, pd.DataFrame] = {}
    ordered = (
        lig_rows[["CHAIN_ID", "ligand_molecule_id"]]
        .dropna()
        .drop_duplicates()
        .sort_values(["CHAIN_ID", "ligand_molecule_id"])
    )
    for _, row in ordered.iterrows():
        chain_id = str(row["CHAIN_ID"]).strip()
        mol_id = str(row["ligand_molecule_id"]).strip()
        if not chain_id or not mol_id:
            continue

        query_smiles = ligand_queries.get(chain_id)
        if not query_smiles:
            query_smiles = _smiles_from_components_cif(mol_id, components_cif_path)
        if not query_smiles:
            query_smiles = _smiles_from_ccd_cache(mol_id, boltz_cache_path)

        if not query_smiles:
            ref = ligands_df[
                ligands_df.get("ligand_id").astype(str) == _norm_id(mol_id)
            ]
            if not ref.empty:
                query_smiles = str(ref.iloc[0]["smiles"])

        if not query_smiles:
            empty_view = ligands_df.iloc[0:0].copy()
            empty_view["query_chain_id"] = pd.Series(dtype=str)
            views[chain_id] = _order_ligand_training_columns(empty_view)
            continue

        sub = _reference_rows_for_query(ligands_df, chain_id)
        sub["ecfp_similarity"] = _ligand_similarity_series(
            query_smiles, sub, mol_id=mol_id
        )
        sub = sub.dropna(subset=["ecfp_similarity"]).copy()
        sub = sub.sort_values("ecfp_similarity", ascending=False, na_position="last")
        if minimum_similarity is not None:
            above = sub[sub["ecfp_similarity"] > minimum_similarity].copy()
            if not above.empty:
                sub = above
            elif fallback_top_n is not None:
                sub = sub.head(fallback_top_n).copy()
        elif fallback_top_n is not None:
            sub = sub.head(fallback_top_n).copy()
        sub["query_chain_id"] = chain_id
        views[chain_id] = _order_ligand_training_columns(sub)

    return views


__all__ = [
    "_best_ligand_hit",
    "_best_pairwise_protein_hit",
    "_best_protein_hit",
    "_best_reference_row",
    "_build_ligand_training_views",
    "_build_protein_training_view",
    "_ligand_similarity_series",
    "_morgan_fp_from_normalized_smiles",
    "_morgan_fp_from_smiles",
    "_pairwise_sequence_identity_percent",
    "_protein_pairwise_similarity_series",
    "_protein_similarity_series",
    "_reference_rows_for_query",
]
