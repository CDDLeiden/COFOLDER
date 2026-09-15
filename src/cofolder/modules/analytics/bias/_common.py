"""Shared constants and normalization helpers for bias analytics."""

from __future__ import annotations

from cofolder.modules.analytics.bias_database import PROTEIN_SIMILARITY_CUTOFF
from cofolder.modules.analytics.bias_database import PROTEIN_SIMILARITY_PERCENT_CUTOFF
from Bio.Align import PairwiseAligner
from pathlib import Path
from datetime import date
import pandas as pd
from rdkit.Chem import rdFingerprintGenerator

ALIGNER = PairwiseAligner(mode="global")

ALIGNER.match_score = 1.0

ALIGNER.mismatch_score = 0.0

ALIGNER.open_gap_score = 0.0

ALIGNER.extend_gap_score = 0.0

MORGAN_FP = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

CORE_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"

CORE_NONPOLY_URL = (
    "https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{pdb_id}/{entity_id}"
)

FASTA_URL = "https://www.rcsb.org/fasta/entry/{pdb_id}/download"

PROVENANCE_COLUMNS = [
    "source",
    "dataset_name",
    "complex_id",
    "source_component_id",
    "source_structure_path",
    "source_reference_path",
]

PROTEIN_BASE_COLUMNS = [
    "pdb_id",
    "release_date",
    "sequence",
    "sequence_similarity",
    "sequence_similarity_pairwise",
    "sequence_similarity_method",
]

LIGAND_BASE_COLUMNS = [
    "pdb_id",
    "release_date",
    "ligand_id",
    "smiles",
    "ecfp_similarity",
]

CUSTOM_PROTEIN_REQUIRED_COLUMNS = {"sequence"}

CUSTOM_LIGAND_REQUIRED_COLUMNS = {"smiles"}

REFERENCE_PATH_COLUMNS = ("source_structure_path", "source_reference_path")

BIAS_PROTEIN_VIEW_MIN_SIMILARITY = PROTEIN_SIMILARITY_PERCENT_CUTOFF

BIAS_LIGAND_VIEW_MIN_SIMILARITY = 0.35

BIAS_PLOT_SEQUENCE_THRESHOLD = PROTEIN_SIMILARITY_CUTOFF

BIAS_PLOT_LIGAND_THRESHOLD = 0.35

BIAS_PLOT_FILE_STEM = "bias_reference_overlap_scatter"

BIAS_PROTEIN_PAIR_FILE_STEM = "bias_protein_pair_scatter"

BIAS_LIGAND_PAIR_FILE_STEM = "bias_ligand_pair_scatter"

BIAS_PLOT_SKIPPED_FILE = f"{BIAS_PLOT_FILE_STEM}.skipped.txt"

BIAS_PROGRESS_LOG_EVERY = 10

BIAS_PROGRESS_WRITE_EVERY = 50

BIAS_PAIRING_PAIRED = "paired"

BIAS_PAIRING_PROTEIN_ONLY = "protein_only"

BIAS_PAIRING_LIGAND_ONLY = "ligand_only"

BIAS_PAIRING_COMPONENT_1_ONLY = "component_1_only"

BIAS_PAIRING_COMPONENT_2_ONLY = "component_2_only"

BIAS_TRAINING_DATA_COLUMNS = [
    "query_pair_id",
    "query_protein_chain_id",
    "query_ligand_chain_id",
    "reference_key",
    "reference_label",
    "pairing_status",
    "source",
    "dataset_name",
    "complex_id",
    "pdb_id",
    "protein_pdb_id",
    "ligand_pdb_id",
    "protein_release_date",
    "ligand_release_date",
    "protein_source",
    "protein_dataset_name",
    "protein_source_structure_path",
    "protein_source_reference_path",
    "ligand_source",
    "ligand_dataset_name",
    "ligand_source_structure_path",
    "ligand_source_reference_path",
    "sequence_similarity",
    "sequence_similarity_pairwise",
    "sequence_similarity_method",
    "ecfp_similarity",
    "plot_sequence_similarity",
    "plot_sequence_similarity_pairwise",
    "plot_ecfp_similarity",
    "sequence",
    "ligand_id",
    "smiles",
]

SAME_TYPE_PAIR_DATA_COLUMNS = [
    "query_pair_id",
    "component_1_id",
    "component_1_type",
    "component_2_id",
    "component_2_type",
    "reference_key",
    "reference_label",
    "reference_pdb_id",
    "pairing_status",
    "source",
    "dataset_name",
    "component_1_reference_key",
    "component_1_pdb_id",
    "component_1_release_date",
    "component_1_source",
    "component_1_dataset_name",
    "component_1_source_structure_path",
    "component_1_source_reference_path",
    "component_2_reference_key",
    "component_2_pdb_id",
    "component_2_release_date",
    "component_2_source",
    "component_2_dataset_name",
    "component_2_source_structure_path",
    "component_2_source_reference_path",
    "component_1_similarity",
    "component_2_similarity",
    "plot_component_1_similarity",
    "plot_component_2_similarity",
]


def _empty_protein_reference_frame() -> pd.DataFrame:
    return pd.DataFrame(columns=PROTEIN_BASE_COLUMNS + PROVENANCE_COLUMNS)


def _empty_ligand_reference_frame() -> pd.DataFrame:
    return pd.DataFrame(columns=LIGAND_BASE_COLUMNS + PROVENANCE_COLUMNS)


def _order_ligand_training_columns(df: pd.DataFrame) -> pd.DataFrame:
    preferred = [
        "query_chain_id",
        "pdb_id",
        "release_date",
        "source",
        "dataset_name",
        "source_structure_path",
        "source_reference_path",
        "ligand_id",
        "ecfp_similarity",
        "smiles",
    ]
    cols = [c for c in preferred if c in df.columns] + [
        c for c in df.columns if c not in preferred
    ]
    return df.loc[:, cols].copy()


def _order_protein_training_columns(df: pd.DataFrame) -> pd.DataFrame:
    preferred = [
        "query_chain_id",
        "pdb_id",
        "release_date",
        "source",
        "dataset_name",
        "source_structure_path",
        "source_reference_path",
        "sequence_similarity",
        "sequence_similarity_pairwise",
        "sequence_similarity_method",
        "sequence",
    ]
    cols = [c for c in preferred if c in df.columns] + [
        c for c in df.columns if c not in preferred
    ]
    return df.loc[:, cols].copy()


def _parse_iso_date(value: str) -> date | None:
    try:
        return date.fromisoformat(str(value)[:10])
    except Exception:
        return None


def _is_before_cutoff(value: str, cutoff: date) -> bool:
    d = _parse_iso_date(value)
    return bool(d and d < cutoff)


def _norm_id(value) -> str:
    return str(value).strip().upper()


def _canonical_reference_source(value, default: str) -> str:
    text = str(value).strip().lower() if value is not None else ""
    if text in {"", "nan", "none"}:
        return default
    if text in {"public", "pdb"}:
        return "public"
    if text in {"custom", "private", "non-public", "non_public"}:
        return "custom"
    return text


def _normalize_reference_path_fields(
    df: pd.DataFrame, input_path: Path
) -> pd.DataFrame:
    for column in REFERENCE_PATH_COLUMNS:
        if column not in df.columns:
            df[column] = pd.NA
            continue
        normalized_values = []
        for value in df[column]:
            if pd.isna(value) or not str(value).strip():
                normalized_values.append(pd.NA)
                continue
            raw_path = Path(str(value).strip()).expanduser()
            resolved_path = (
                raw_path if raw_path.is_absolute() else (input_path.parent / raw_path)
            )
            if not resolved_path.exists():
                raise ValueError(
                    f"{input_path.name} column '{column}' references a missing file: {resolved_path}"
                )
            if not resolved_path.is_file():
                raise ValueError(
                    f"{input_path.name} column '{column}' is not a file: {resolved_path}"
                )
            normalized_values.append(str(resolved_path))
        df[column] = normalized_values
    return df


def _normalize_smiles(smiles: str | None) -> str:
    if smiles is None:
        return ""
    normalized = str(smiles).strip()
    while len(normalized) >= 2 and (
        (normalized[0] == '"' and normalized[-1] == '"')
        or (normalized[0] == "'" and normalized[-1] == "'")
    ):
        normalized = normalized[1:-1].strip()
    return normalized.strip(";").strip()


def _path_cache_token(path: Path | None) -> str:
    if path is None:
        return ""
    return str(Path(path).expanduser().resolve())


__all__ = [
    "ALIGNER",
    "BIAS_LIGAND_PAIR_FILE_STEM",
    "BIAS_LIGAND_VIEW_MIN_SIMILARITY",
    "BIAS_PAIRING_COMPONENT_1_ONLY",
    "BIAS_PAIRING_COMPONENT_2_ONLY",
    "BIAS_PAIRING_LIGAND_ONLY",
    "BIAS_PAIRING_PAIRED",
    "BIAS_PAIRING_PROTEIN_ONLY",
    "BIAS_PLOT_FILE_STEM",
    "BIAS_PLOT_LIGAND_THRESHOLD",
    "BIAS_PLOT_SEQUENCE_THRESHOLD",
    "BIAS_PLOT_SKIPPED_FILE",
    "BIAS_PROGRESS_LOG_EVERY",
    "BIAS_PROGRESS_WRITE_EVERY",
    "BIAS_PROTEIN_PAIR_FILE_STEM",
    "BIAS_PROTEIN_VIEW_MIN_SIMILARITY",
    "BIAS_TRAINING_DATA_COLUMNS",
    "CORE_ENTRY_URL",
    "CORE_NONPOLY_URL",
    "CUSTOM_LIGAND_REQUIRED_COLUMNS",
    "CUSTOM_PROTEIN_REQUIRED_COLUMNS",
    "FASTA_URL",
    "LIGAND_BASE_COLUMNS",
    "MORGAN_FP",
    "PROTEIN_BASE_COLUMNS",
    "PROVENANCE_COLUMNS",
    "REFERENCE_PATH_COLUMNS",
    "SAME_TYPE_PAIR_DATA_COLUMNS",
    "_canonical_reference_source",
    "_empty_ligand_reference_frame",
    "_empty_protein_reference_frame",
    "_is_before_cutoff",
    "_norm_id",
    "_normalize_reference_path_fields",
    "_normalize_smiles",
    "_order_ligand_training_columns",
    "_order_protein_training_columns",
    "_parse_iso_date",
    "_path_cache_token",
]
