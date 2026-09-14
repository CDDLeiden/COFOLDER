"""Lightweight bias metrics against protein/ligand training references."""

from __future__ import annotations

from collections.abc import Callable
from contextlib import nullcontext
from datetime import date
from functools import lru_cache
from itertools import combinations
import json
import logging
from pathlib import Path
import pickle
import subprocess
import tempfile
from urllib.request import urlopen

from Bio.Align import PairwiseAligner
import gemmi
import pandas as pd
from pandas.errors import EmptyDataError
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

from cofolder.modules.analytics.bias_training import _resolve_mmseqs_bin
from cofolder.modules.analytics.bias_database import (
    PROTEIN_SIMILARITY_CUTOFF,
    PROTEIN_SIMILARITY_PERCENT_CUTOFF,
    parse_bias_release_policy,
)
from cofolder.modules.analytics.plots import plot_bias_reference_overlap, plot_reference_overlap_scatter
from cofolder.modules.utils.timing import DebugTimingCollector

ALIGNER = PairwiseAligner(mode="global")
ALIGNER.match_score = 1.0
ALIGNER.mismatch_score = 0.0
ALIGNER.open_gap_score = 0.0
ALIGNER.extend_gap_score = 0.0
MORGAN_FP = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
CORE_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
CORE_NONPOLY_URL = "https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{pdb_id}/{entity_id}"
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
LIGAND_BASE_COLUMNS = ["pdb_id", "release_date", "ligand_id", "smiles", "ecfp_similarity"]
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
    cols = [c for c in preferred if c in df.columns] + [c for c in df.columns if c not in preferred]
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
    cols = [c for c in preferred if c in df.columns] + [c for c in df.columns if c not in preferred]
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


def _normalize_reference_path_fields(df: pd.DataFrame, input_path: Path) -> pd.DataFrame:
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
            resolved_path = raw_path if raw_path.is_absolute() else (input_path.parent / raw_path)
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


def _finalize_reference_frame(
    df: pd.DataFrame,
    *,
    input_path: Path,
    default_source: str,
    default_dataset_name: str,
    is_ligand: bool,
) -> pd.DataFrame:
    if df.empty and len(df.columns) == 0:
        df = _empty_ligand_reference_frame() if is_ligand else _empty_protein_reference_frame()

    base_columns = LIGAND_BASE_COLUMNS if is_ligand else PROTEIN_BASE_COLUMNS
    for column in base_columns:
        if column not in df.columns:
            df[column] = pd.NA

    if "source" not in df.columns:
        df["source"] = default_source
    else:
        df["source"] = df["source"].apply(lambda value: _canonical_reference_source(value, default_source))

    if "dataset_name" not in df.columns:
        df["dataset_name"] = default_dataset_name
    else:
        df["dataset_name"] = df["dataset_name"].apply(
            lambda value: str(value).strip() if pd.notna(value) and str(value).strip() else default_dataset_name
        )

    df = _normalize_reference_path_fields(df, input_path)

    if is_ligand:
        if "ligand_id" not in df.columns:
            df["ligand_id"] = None
        else:
            df["ligand_id"] = df["ligand_id"].apply(
                lambda x: _norm_id(x) if pd.notna(x) and str(x).strip() else None
            )
        if "smiles" in df.columns:
            df["smiles"] = df["smiles"].apply(_normalize_smiles)
        if "ecfp_similarity" in df.columns:
            df["ecfp_similarity"] = pd.to_numeric(df["ecfp_similarity"], errors="coerce")
    else:
        if "sequence" in df.columns:
            df["sequence"] = df["sequence"].astype(str)
        df["sequence_similarity"] = pd.to_numeric(df["sequence_similarity"], errors="coerce")
        df["sequence_similarity_pairwise"] = pd.to_numeric(
            df["sequence_similarity_pairwise"], errors="coerce"
        )
        both_present = (
            df["sequence_similarity"].notna()
            & df["sequence_similarity_pairwise"].notna()
        )
        if both_present.any():
            raise ValueError(
                "Protein reference rows cannot contain both sequence_similarity "
                "(MMseqs pident) and sequence_similarity_pairwise"
            )
        method = df["sequence_similarity_method"].astype("string")
        missing_method = method.isna() | method.str.strip().eq("")
        method = method.mask(
            missing_method & df["sequence_similarity"].notna(),
            "mmseqs_pident",
        )
        method = method.mask(
            missing_method
            & df["sequence_similarity"].isna()
            & df["sequence_similarity_pairwise"].notna(),
            "pairwise_aligner",
        )
        method = method.mask(
            missing_method
            & df["sequence_similarity"].isna()
            & df["sequence_similarity_pairwise"].isna(),
            "unavailable",
        )
        invalid_mmseqs = df["sequence_similarity"].notna() & method.ne("mmseqs_pident")
        invalid_pairwise = (
            df["sequence_similarity_pairwise"].notna()
            & method.ne("pairwise_aligner")
        )
        invalid_unavailable = (
            df["sequence_similarity"].isna()
            & df["sequence_similarity_pairwise"].isna()
            & method.ne("unavailable")
        )
        if (
            invalid_mmseqs.any()
            or invalid_pairwise.any()
            or invalid_unavailable.any()
        ):
            raise ValueError(
                "Protein similarity values do not match sequence_similarity_method"
            )
        df["sequence_similarity_method"] = method

    return df


def _smiles_from_ccd_cache(ccd_id: str, boltz_cache_path: Path) -> str | None:
    """Resolve CCD ID to SMILES from Boltz mol cache (.boltz/mols/<CCD>.pkl)."""
    ccd = _norm_id(ccd_id)
    pkl = boltz_cache_path / "mols" / f"{ccd}.pkl"
    if not pkl.exists():
        return None
    try:
        with pkl.open("rb") as fh:
            mol = pickle.load(fh)
        if mol is None:
            return None
        if isinstance(mol, Chem.Mol):
            return Chem.MolToSmiles(mol)
    except Exception:
        return None
    return None


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


@lru_cache(maxsize=8)
def _components_cif_smiles_index(components_cif_token: str) -> dict[str, str]:
    if not components_cif_token:
        return {}

    components_cif_path = Path(components_cif_token)
    if not components_cif_path.exists():
        return {}

    try:
        doc = gemmi.cif.read_file(str(components_cif_path))
    except Exception:
        return {}

    descriptor_priority = {
        "SMILES_CANONICAL": 0,
        "SMILES": 1,
    }
    best_by_ccd: dict[str, tuple[int, str]] = {}

    for block in doc:
        block_id = _norm_id(block.find_value("_chem_comp.id") or block.name)
        if not block_id:
            continue

        table = block.find(
            "_pdbx_chem_comp_descriptor.",
            ["comp_id", "type", "program", "descriptor"],
        )
        for row in table:
            if len(row) != 4:
                continue
            row_ccd_id = _norm_id(row[0])
            if row_ccd_id != block_id:
                continue
            descriptor_type = str(row[1]).strip().upper()
            descriptor = _normalize_smiles(str(row[3]))
            if not descriptor:
                continue
            rank = descriptor_priority.get(descriptor_type, 99)
            current = best_by_ccd.get(row_ccd_id)
            if current is None or rank < current[0]:
                best_by_ccd[row_ccd_id] = (rank, descriptor)

    return {
        ccd_id: descriptor
        for ccd_id, (_, descriptor) in best_by_ccd.items()
    }


def _smiles_from_components_cif(ccd_id: str, components_cif_path: Path | None) -> str | None:
    if components_cif_path is None or not components_cif_path.exists():
        return None

    target_id = _norm_id(ccd_id)
    return _components_cif_smiles_index(_path_cache_token(components_cif_path)).get(target_id)


def _load_public_protein_training(path: Path, cutoff: date | None, top_n: int = 100) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {"pdb_id", "release_date", "sequence"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Protein training data missing required columns: {sorted(missing)}")

    release_mask = (
        pd.Series(True, index=df.index)
        if cutoff is None
        else df["release_date"].map(lambda x: _is_before_cutoff(x, cutoff)).fillna(False).astype(bool)
    )
    out = df.loc[release_mask, :].copy()
    if "sequence" not in out.columns:
        out = out.reindex(columns=df.columns)
    out = _finalize_reference_frame(
        out,
        input_path=path,
        default_source="public",
        default_dataset_name="public",
        is_ligand=False,
    )

    if top_n > 0 and len(out) > top_n:
        if "sequence_similarity" in out.columns:
            out = out.sort_values("sequence_similarity", ascending=False).head(top_n).copy()
        else:
            out = out.head(top_n).copy()

    return out


def _load_custom_protein_training(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    missing = CUSTOM_PROTEIN_REQUIRED_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(
            "Custom protein reference data missing required columns: "
            f"{sorted(missing)}. Supported optional provenance columns: {list(PROVENANCE_COLUMNS)}"
        )
    return _finalize_reference_frame(
        df.copy(),
        input_path=path,
        default_source="custom",
        default_dataset_name=path.stem,
        is_ligand=False,
    )


def _read_ligand_training_from_sdf(path: Path) -> pd.DataFrame:
    rows = []
    supplier = Chem.SDMolSupplier(str(path), removeHs=False)
    for mol in supplier:
        if mol is None:
            continue
        smiles = Chem.MolToSmiles(mol)
        pdb_id = mol.GetProp("pdb_id") if mol.HasProp("pdb_id") else None
        release_date = mol.GetProp("release_date") if mol.HasProp("release_date") else None
        ligand_id = mol.GetProp("ligand_id") if mol.HasProp("ligand_id") else None
        source = mol.GetProp("source") if mol.HasProp("source") else None
        dataset_name = mol.GetProp("dataset_name") if mol.HasProp("dataset_name") else None
        source_structure_path = (
            mol.GetProp("source_structure_path") if mol.HasProp("source_structure_path") else None
        )
        source_reference_path = (
            mol.GetProp("source_reference_path") if mol.HasProp("source_reference_path") else None
        )
        ecfp_similarity = mol.GetProp("ecfp_similarity") if mol.HasProp("ecfp_similarity") else None
        if smiles:
            rows.append(
                {
                    "pdb_id": str(pdb_id) if pdb_id is not None else None,
                    "release_date": str(release_date) if release_date is not None else None,
                    "smiles": str(smiles),
                    "ligand_id": str(ligand_id) if ligand_id is not None else None,
                    "source": str(source) if source is not None else None,
                    "dataset_name": str(dataset_name) if dataset_name is not None else None,
                    "source_structure_path": (
                        str(source_structure_path) if source_structure_path is not None else None
                    ),
                    "source_reference_path": (
                        str(source_reference_path) if source_reference_path is not None else None
                    ),
                    "ecfp_similarity": str(ecfp_similarity) if ecfp_similarity is not None else None,
                }
            )
    return pd.DataFrame(rows)


def _load_public_ligand_training(path: Path, cutoff: date | None) -> pd.DataFrame:
    if path.suffix.lower() == ".sdf":
        df = _read_ligand_training_from_sdf(path)
    else:
        try:
            df = pd.read_csv(path)
        except EmptyDataError:
            df = _empty_ligand_reference_frame()

    if df.empty and len(df.columns) == 0:
        df = _empty_ligand_reference_frame()

    required = {"pdb_id", "release_date", "smiles"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Ligand training data missing required columns: {sorted(missing)}")

    release_mask = (
        pd.Series(True, index=df.index)
        if cutoff is None
        else df["release_date"].map(lambda x: _is_before_cutoff(x, cutoff)).fillna(False).astype(bool)
    )
    out = df.loc[release_mask, :].copy()
    return _finalize_reference_frame(
        out,
        input_path=path,
        default_source="public",
        default_dataset_name="public",
        is_ligand=True,
    )


def _load_custom_ligand_training(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".sdf":
        df = _read_ligand_training_from_sdf(path)
    else:
        try:
            df = pd.read_csv(path)
        except EmptyDataError:
            df = _empty_ligand_reference_frame()

    if df.empty and len(df.columns) == 0:
        df = _empty_ligand_reference_frame()

    missing = CUSTOM_LIGAND_REQUIRED_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(
            "Custom ligand reference data missing required columns: "
            f"{sorted(missing)}. Supported optional provenance columns/properties: {list(PROVENANCE_COLUMNS)}"
        )

    return _finalize_reference_frame(
        df.copy(),
        input_path=path,
        default_source="custom",
        default_dataset_name=path.stem,
        is_ligand=True,
    )


def _concat_reference_frames(frames: list[pd.DataFrame], *, is_ligand: bool) -> pd.DataFrame:
    non_empty_frames = [frame for frame in frames if frame is not None]
    if not non_empty_frames:
        return _empty_ligand_reference_frame() if is_ligand else _empty_protein_reference_frame()
    return pd.concat(non_empty_frames, ignore_index=True, sort=False)


def _load_ligand_training(
    path: Path,
    cutoff: date | None,
    *,
    custom_reference_path: Path | None = None,
) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    if path is not None:
        frames.append(_load_public_ligand_training(path, cutoff))
    if custom_reference_path is not None:
        frames.append(_load_custom_ligand_training(custom_reference_path))
    return _concat_reference_frames(frames, is_ligand=True)


def _load_protein_training(
    path: Path,
    cutoff: date | None,
    top_n: int = 100,
    *,
    custom_reference_path: Path | None = None,
) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    if path is not None:
        frames.append(_load_public_protein_training(path, cutoff, top_n=top_n))
    if custom_reference_path is not None:
        frames.append(_load_custom_protein_training(custom_reference_path))
    return _concat_reference_frames(frames, is_ligand=False)


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
        sims = pd.to_numeric(proteins_df["sequence_similarity"], errors="coerce").dropna()
        if len(sims):
            return float(sims.max())
    return None


def _best_pairwise_protein_hit(query_seq: str, proteins_df: pd.DataFrame) -> float | None:
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


def _fetch_json(url: str, timeout: int = 30) -> dict | list:
    with urlopen(url, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8"))


def _fetch_text(url: str, timeout: int = 30) -> str:
    with urlopen(url, timeout=timeout) as response:
        return response.read().decode("utf-8", errors="replace")


def _extract_nonpoly_comp_id(nonpoly_payload: dict) -> str | None:
    if not isinstance(nonpoly_payload, dict):
        return None
    container = nonpoly_payload.get("rcsb_nonpolymer_entity_container_identifiers", {})
    if isinstance(container, dict):
        comp = container.get("nonpolymer_comp_id")
        if comp:
            return str(comp).upper()
    ent = nonpoly_payload.get("pdbx_entity_nonpoly", {})
    if isinstance(ent, dict):
        comp = ent.get("comp_id")
        if comp:
            return str(comp).upper()
    return None


@lru_cache(maxsize=4096)
def _entry_ligand_ids(pdb_id: str, timeout: int = 30) -> tuple[str, ...]:
    try:
        payload = _fetch_json(CORE_ENTRY_URL.format(pdb_id=pdb_id), timeout=timeout)
    except Exception:
        return tuple()
    ids_obj = payload.get("rcsb_entry_container_identifiers", {}) if isinstance(payload, dict) else {}
    entity_ids = ids_obj.get("non_polymer_entity_ids", []) if isinstance(ids_obj, dict) else []
    out: set[str] = set()
    for ent_id in entity_ids:
        try:
            nonpoly = _fetch_json(
                CORE_NONPOLY_URL.format(pdb_id=pdb_id, entity_id=ent_id),
                timeout=timeout,
            )
        except Exception:
            continue
        comp = _extract_nonpoly_comp_id(nonpoly)
        if comp:
            out.add(comp)
    return tuple(sorted(out))


@lru_cache(maxsize=4096)
def _entry_fasta_sequences(pdb_id: str, timeout: int = 30) -> tuple[str, ...]:
    try:
        txt = _fetch_text(FASTA_URL.format(pdb_id=pdb_id), timeout=timeout)
    except Exception:
        return tuple()
    sequences: list[str] = []
    buffer: list[str] = []
    for raw_line in txt.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if buffer:
                sequences.append("".join(buffer))
                buffer = []
            continue
        buffer.append(line)
    if buffer:
        sequences.append("".join(buffer))
    return tuple(sequence for sequence in sequences if sequence)


@lru_cache(maxsize=65536)
def _resolve_pdb_ligand_smiles(
    ccd_id: str,
    boltz_cache_token: str,
    components_cif_token: str,
) -> str:
    normalized_ccd_id = _norm_id(ccd_id)
    smiles = _components_cif_smiles_index(components_cif_token).get(normalized_ccd_id)
    if smiles:
        return smiles
    if boltz_cache_token:
        cached = _smiles_from_ccd_cache(normalized_ccd_id, Path(boltz_cache_token))
        if cached:
            return _normalize_smiles(cached)
    return ""


@lru_cache(maxsize=4096)
def _cached_pdb_ligand_entries(
    pdb_id: str,
    boltz_cache_token: str,
    components_cif_token: str,
    timeout: int = 30,
) -> tuple[tuple[str, str], ...]:
    entries: list[tuple[str, str]] = []
    for ligand_id in _entry_ligand_ids(pdb_id, timeout=timeout):
        smiles = _resolve_pdb_ligand_smiles(
            ligand_id,
            boltz_cache_token,
            components_cif_token,
        )
        if not smiles:
            continue
        entries.append((str(ligand_id).strip(), smiles))
    return tuple(entries)


@lru_cache(maxsize=8192)
def _cached_pdb_protein_similarity_rows(
    pdb_id: str,
    query_sequence: str,
    timeout: int = 30,
) -> tuple[tuple[str, float], ...]:
    if not query_sequence:
        return tuple()
    target_sequences = tuple(
        sequence
        for sequence in _entry_fasta_sequences(pdb_id, timeout=timeout)
        if sequence
    )
    if not target_sequences:
        return tuple()

    resolved_mmseqs = _resolved_mmseqs_bin_token()
    if not resolved_mmseqs:
        return tuple()

    with tempfile.TemporaryDirectory(prefix="cofolder-bias-mmseqs-") as tmpdir:
        tmp_root = Path(tmpdir)
        query_fasta = tmp_root / "query.fasta"
        target_fasta = tmp_root / "target.fasta"
        _write_fasta_records(query_fasta, [("query", str(query_sequence))])
        _write_fasta_records(
            target_fasta,
            [
                (f"{str(pdb_id).strip().upper()}_{index}", sequence)
                for index, sequence in enumerate(target_sequences, start=1)
            ],
        )

        qdb = tmp_root / "query_db"
        tdb = tmp_root / "target_db"
        rdb = tmp_root / "result_db"
        out_tsv = tmp_root / "result.tsv"
        search_tmp = tmp_root / "search_tmp"
        search_tmp.mkdir(parents=True, exist_ok=True)

        commands = [
            [resolved_mmseqs, "createdb", str(query_fasta), str(qdb)],
            [resolved_mmseqs, "createdb", str(target_fasta), str(tdb)],
            [
                resolved_mmseqs,
                "search",
                str(qdb),
                str(tdb),
                str(rdb),
                str(search_tmp),
                "--threads",
                "1",
                "--max-seqs",
                str(max(1, len(target_sequences))),
            ],
            [
                resolved_mmseqs,
                "convertalis",
                str(qdb),
                str(tdb),
                str(rdb),
                str(out_tsv),
                "--format-output",
                "target,pident,tseq",
            ],
        ]

        try:
            for command in commands:
                subprocess.run(
                    command,
                    check=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    timeout=max(30, int(timeout)),
                )
        except (OSError, subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
            logging.getLogger(__name__).warning(
                "MMseqs protein fallback failed for pdb_id=%s: %s",
                str(pdb_id).strip().upper(),
                exc,
            )
            return tuple()

        similarity_by_sequence: dict[str, float] = {}
        if out_tsv.exists():
            try:
                hits = pd.read_csv(
                    out_tsv,
                    sep="\t",
                    header=None,
                    names=["target", "pident", "tseq"],
                )
            except EmptyDataError:
                hits = pd.DataFrame(columns=["target", "pident", "tseq"])
            if not hits.empty:
                hits["pident"] = pd.to_numeric(hits["pident"], errors="coerce").fillna(0.0)
                hits["tseq"] = hits["tseq"].astype(str)
                for _, hit in hits.iterrows():
                    sequence = str(hit["tseq"])
                    similarity = float(hit["pident"])
                    current = similarity_by_sequence.get(sequence)
                    if current is None or similarity > current:
                        similarity_by_sequence[sequence] = similarity

    return tuple(
        (sequence, float(similarity_by_sequence.get(sequence, 0.0)))
        for sequence in target_sequences
    )


@lru_cache(maxsize=1)
def _resolved_mmseqs_bin_token() -> str:
    return _resolve_mmseqs_bin() or ""


def _write_fasta_records(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for record_id, sequence in records:
            handle.write(f">{record_id}\n{sequence}\n")


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
        return pd.Series([pd.NA] * len(ligands_df), index=ligands_df.index, dtype="Float64")

    similarities = []
    for _, row in ligands_df.iterrows():
        row_smiles = str(row["smiles"])
        row_ligand_id = _norm_id(row["ligand_id"]) if pd.notna(row.get("ligand_id")) and str(row.get("ligand_id")).strip() else None
        if row_smiles == str(query_smiles):
            similarities.append(1.0)
            continue
        if mol_id and row_ligand_id == _norm_id(mol_id):
            precomputed = pd.to_numeric(pd.Series([row.get("ecfp_similarity")]), errors="coerce").iloc[0]
            if pd.notna(precomputed):
                similarities.append(float(precomputed))
                continue
        fp = _morgan_fp_from_smiles(str(row["smiles"]))
        if fp is None:
            similarities.append(pd.NA)
        else:
            similarities.append(float(DataStructs.TanimotoSimilarity(query_fp, fp)))
    return pd.to_numeric(pd.Series(similarities, index=ligands_df.index), errors="coerce")


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
        similarity_methods = pd.Series(
            "ecfp4_tanimoto", index=df.index, dtype="string"
        )
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

    protein_chain_ids = set(chain_df.loc[chain_df["ENTITY_TYPE"].astype(str) == "protein", "CHAIN_ID"].astype(str))
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

    ligand_chain_ids = set(chain_df.loc[chain_df["ENTITY_TYPE"].astype(str) == "ligand", "CHAIN_ID"].astype(str))
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
        public_best = _best_reference_row(public_df, query_smiles, is_ligand=True, mol_id=mol_id)
        custom_best = _best_reference_row(custom_df, query_smiles, is_ligand=True, mol_id=mol_id)
        overall_best = _best_reference_row(ligands_df, query_smiles, is_ligand=True, mol_id=mol_id)
        summary_rows.append(
            _build_reference_landscape_row(
                chain_id=str(chain_id),
                entity_type="ligand",
                public_best=public_best,
                custom_best=custom_best,
                overall_best=overall_best,
            )
        )

    pd.DataFrame(summary_rows).to_csv(output_dir / "reference_landscape_summary.csv", index=False)


def _build_reference_landscape_row(
    *,
    chain_id: str,
    entity_type: str,
    public_best: pd.Series | None,
    custom_best: pd.Series | None,
    overall_best: pd.Series | None,
) -> dict[str, object]:
    nearest_overall_source = (
        str(overall_best.get("source")).strip() if overall_best is not None and pd.notna(overall_best.get("source")) else None
    )
    public_similarity = (
        float(public_best["best_similarity"]) if public_best is not None and pd.notna(public_best.get("best_similarity")) else None
    )
    custom_similarity = (
        float(custom_best["best_similarity"]) if custom_best is not None and pd.notna(custom_best.get("best_similarity")) else None
    )
    return {
        "query_chain_id": chain_id,
        "entity_type": entity_type,
        "nearest_public_pdb_id": public_best.get("pdb_id") if public_best is not None else None,
        "nearest_public_dataset_name": public_best.get("dataset_name") if public_best is not None else None,
        "nearest_public_similarity": public_similarity,
        "nearest_public_similarity_method": (
            public_best.get("best_similarity_method")
            if public_best is not None
            else None
        ),
        "nearest_custom_pdb_id": custom_best.get("pdb_id") if custom_best is not None else None,
        "nearest_custom_dataset_name": custom_best.get("dataset_name") if custom_best is not None else None,
        "nearest_custom_similarity": custom_similarity,
        "nearest_custom_similarity_method": (
            custom_best.get("best_similarity_method")
            if custom_best is not None
            else None
        ),
        "nearest_overall_source": nearest_overall_source,
        "nearest_overall_pdb_id": overall_best.get("pdb_id") if overall_best is not None else None,
        "nearest_overall_dataset_name": overall_best.get("dataset_name") if overall_best is not None else None,
        "nearest_overall_similarity": (
            float(overall_best["best_similarity"])
            if overall_best is not None and pd.notna(overall_best.get("best_similarity"))
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
                or (custom_similarity is not None and custom_similarity > public_similarity)
            )
        ),
    }


def _extract_system_queries(
    sys_obj,
    ligands_df: pd.DataFrame,
    boltz_cache_path: Path,
    components_cif_path: Path | None = None,
):
    sequences = sys_obj.find_value(key="sequences") or []
    protein_by_chain = {}
    ligand_by_chain = {}
    unresolved_ccd_by_chain: dict[str, tuple[str, ...]] = {}
    ligand_fallback_by_id = {}

    for _, row in ligands_df.dropna(subset=["ligand_id", "smiles"]).iterrows():
        ligand_fallback_by_id[_norm_id(row["ligand_id"])] = str(row["smiles"])

    for seq in sequences:
        if not isinstance(seq, dict):
            continue
        if "protein" in seq:
            p = seq["protein"] or {}
            ids = p.get("id")
            if ids is not None:
                if not isinstance(ids, list):
                    ids = [ids]
                sequence = p.get("sequence") or p.get("fasta")
                for cid in ids:
                    protein_by_chain[str(cid).strip()] = str(sequence) if sequence else None

        if "ligand" in seq:
            ligand_data = seq["ligand"] or {}
            ids = ligand_data.get("id")
            if ids is None:
                continue
            if not isinstance(ids, list):
                ids = [ids]

            smiles = ligand_data.get("smiles")
            ccd = ligand_data.get("ccd")
            ccd_ids = ccd if isinstance(ccd, list) else ([ccd] if ccd else [])
            ccd_ids = [_norm_id(x) for x in ccd_ids if str(x).strip()]
            for cid in ids:
                chain_id = str(cid).strip()
                if smiles:
                    ligand_by_chain[chain_id] = str(smiles)
                    continue
                if not ccd_ids:
                    ligand_by_chain[chain_id] = None
                    continue

                mapped_smiles = None
                for ccd_id in ccd_ids:
                    mapped_smiles = ligand_fallback_by_id.get(ccd_id)
                    if not mapped_smiles:
                        mapped_smiles = _smiles_from_components_cif(ccd_id, components_cif_path)
                    if not mapped_smiles:
                        mapped_smiles = _smiles_from_ccd_cache(ccd_id, boltz_cache_path)
                    if mapped_smiles:
                        break

                ligand_by_chain[chain_id] = mapped_smiles
                if mapped_smiles is None:
                    unresolved_ccd_by_chain[chain_id] = tuple(ccd_ids)

    return protein_by_chain, ligand_by_chain, unresolved_ccd_by_chain


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
        sub["sequence_similarity"] = pd.to_numeric(sub["sequence_similarity"], errors="coerce").clip(
            lower=0.0,
            upper=100.0,
        )
        sub["sequence_similarity_pairwise"] = pd.to_numeric(
            sub["sequence_similarity_pairwise"], errors="coerce"
        ).clip(lower=0.0, upper=100.0)
        sub["sequence_similarity_method"] = pd.Series(pd.NA, index=sub.index, dtype="string")
        sub.loc[sub["sequence_similarity"].notna(), "sequence_similarity_method"] = "mmseqs_pident"
        sub.loc[
            sub["sequence_similarity"].isna() & sub["sequence_similarity_pairwise"].notna(),
            "sequence_similarity_method",
        ] = "pairwise_aligner"
        sub["_effective_sequence_similarity"] = sub["sequence_similarity"].copy()
        missing_effective_similarity = sub["_effective_sequence_similarity"].isna()
        sub.loc[missing_effective_similarity, "_effective_sequence_similarity"] = sub.loc[
            missing_effective_similarity, "sequence_similarity_pairwise"
        ]
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
            ref = ligands_df[ligands_df.get("ligand_id").astype(str) == _norm_id(mol_id)]
            if not ref.empty:
                query_smiles = str(ref.iloc[0]["smiles"])

        if not query_smiles:
            empty_view = ligands_df.iloc[0:0].copy()
            empty_view["query_chain_id"] = pd.Series(dtype=str)
            views[chain_id] = _order_ligand_training_columns(empty_view)
            continue

        sub = _reference_rows_for_query(ligands_df, chain_id)
        sub["ecfp_similarity"] = _ligand_similarity_series(query_smiles, sub, mol_id=mol_id)
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


def _collapse_pair_value(values: list[object], *, mixed_label: str | None = None) -> str | None:
    normalized = [str(value).strip() for value in values if pd.notna(value) and str(value).strip()]
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
        candidate = Path(text).stem if "/" in text or text.endswith((".pdb", ".cif", ".sdf", ".csv")) else text
        if candidate:
            return candidate
    return None


def _prepare_reference_side(view: pd.DataFrame, *, is_ligand: bool, suffix: str) -> pd.DataFrame:
    prepared = view.copy()
    if "reference_key" not in prepared.columns:
        prepared["reference_key"] = prepared.apply(
            lambda row: _reference_key_from_row(row, is_ligand=is_ligand),
            axis=1,
        )
    return prepared.rename(
        columns=lambda column: column if column == "reference_key" else f"{column}_{suffix}"
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


def _matching_reference_rows(reference_view: pd.DataFrame, reference_key: object) -> pd.DataFrame:
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
    pairwise = pd.to_numeric(
        ordered["sequence_similarity_pairwise"], errors="coerce"
    )
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


def _pdb_ligand_similarity_rows(
    *,
    pdb_id: str,
    query_smiles: str,
    boltz_cache_path: Path,
    components_cif_path: Path | None,
    exclude_ligand_ids: set[str] | None = None,
    timeout: int = 30,
) -> list[dict[str, object]]:
    query_fp = _morgan_fp_from_smiles(query_smiles)
    if query_fp is None:
        return []
    excluded_ids = {
        _norm_id(ligand_id)
        for ligand_id in (exclude_ligand_ids or set())
        if str(ligand_id).strip()
    }
    boltz_cache_token = _path_cache_token(boltz_cache_path)
    components_cif_token = _path_cache_token(components_cif_path)
    rows: list[dict[str, object]] = []
    for ligand_id, smiles in _cached_pdb_ligand_entries(
        str(pdb_id).strip().upper(),
        boltz_cache_token,
        components_cif_token,
        timeout=timeout,
    ):
        if _norm_id(ligand_id) in excluded_ids:
            continue
        ligand_fp = _morgan_fp_from_smiles(smiles)
        if ligand_fp is None:
            continue
        rows.append(
            {
                "ligand_id": ligand_id,
                "smiles": smiles,
                "ecfp_similarity": float(DataStructs.TanimotoSimilarity(query_fp, ligand_fp)),
            }
        )
    return rows


def _pdb_protein_similarity_rows(
    *,
    pdb_id: str,
    query_sequence: str,
    timeout: int = 30,
) -> list[dict[str, object]]:
    if not query_sequence:
        return []
    return [
        {
            "sequence": sequence,
            "sequence_similarity": similarity,
            "sequence_similarity_pairwise": pd.NA,
            "sequence_similarity_method": "mmseqs_pident",
        }
        for sequence, similarity in _cached_pdb_protein_similarity_rows(
            str(pdb_id).strip().upper(),
            str(query_sequence),
            timeout=timeout,
        )
    ]


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
        if ligand_reference_df is None or ligand_reference_df.empty or "pdb_id" not in ligand_reference_df.columns:
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
        subset = subset.drop_duplicates(subset=["ligand_id", "smiles"], keep="first").reset_index(drop=True)
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
                "ecfp_similarity": float(DataStructs.TanimotoSimilarity(query_fp, ligand_fp)),
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
        dict(row)
        for row in protein_lookup_index.get(str(pdb_id).strip().upper(), [])
    ]


def _build_ligand_reference_index(
    ligand_reference_df: pd.DataFrame | None,
) -> dict[str, list[dict[str, str]]]:
    if ligand_reference_df is None or ligand_reference_df.empty:
        return {}
    if "pdb_id" not in ligand_reference_df.columns or "smiles" not in ligand_reference_df.columns:
        return {}

    subset = ligand_reference_df.copy()
    subset["pdb_id"] = subset["pdb_id"].astype(str).str.strip().str.upper()
    subset["ligand_id"] = subset.get("ligand_id", pd.Series(dtype=object)).astype(str).str.strip()
    subset["smiles"] = subset["smiles"].apply(_normalize_smiles)
    subset = subset[
        subset["pdb_id"].astype(bool)
        & subset["smiles"].astype(str).str.strip().astype(bool)
    ].copy()
    if subset.empty:
        return {}
    subset = subset.drop_duplicates(subset=["pdb_id", "ligand_id", "smiles"], keep="first")

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


def _write_progress_checkpoint(
    *,
    rows: list[dict[str, object]],
    pending_rows: list[dict[str, object]],
    output_path: Path | None,
    finalize: Callable[[pd.DataFrame], pd.DataFrame],
    progress_label: str,
    completed_lookups: int,
    total_lookups: int,
) -> None:
    logger = logging.getLogger(__name__)
    if output_path is None:
        logger.info(
            "Bias lookup progress for %s: %d/%d lookups completed",
            progress_label,
            completed_lookups,
            total_lookups,
        )
        return

    checkpoint_df = finalize(pd.DataFrame([*rows, *pending_rows]))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    checkpoint_df.to_csv(output_path, index=False)
    logger.info(
        "Bias lookup progress for %s: %d/%d lookups completed; wrote checkpoint to %s",
        progress_label,
        completed_lookups,
        total_lookups,
        output_path,
    )


def _log_bias_lookup_progress(
    *,
    progress_label: str,
    completed_lookups: int,
    total_lookups: int,
) -> None:
    logging.getLogger(__name__).info(
        "Bias lookup progress for %s: %d/%d lookups completed",
        progress_label,
        completed_lookups,
        total_lookups,
    )


def _flatten_resolved_bias_rows(
    resolved_rows: dict[int, list[dict[str, object]]],
) -> list[dict[str, object]]:
    flattened: list[dict[str, object]] = []
    for index in sorted(resolved_rows):
        flattened.extend(resolved_rows[index])
    return flattened


def _pending_bias_rows(
    *,
    resolved_rows: dict[int, list[dict[str, object]]],
    input_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    return [
        row_dict
        for index, row_dict in enumerate(input_rows)
        if index not in resolved_rows
    ]


def _enrich_mixed_bias_row_with_pdb_backfill(
    *,
    row_dict: dict[str, object],
    protein_queries: dict[str, str | None],
    ligand_queries: dict[str, str | None],
    protein_lookup_index: dict[str, list[dict[str, object]]] | None = None,
    boltz_cache_path: Path,
    components_cif_path: Path | None,
    ligand_reference_df: pd.DataFrame | None = None,
    ligand_reference_index: dict[str, list[dict[str, str]]] | None = None,
    timeout: int = 30,
) -> list[dict[str, object]]:
    row = pd.Series(row_dict)
    pdb_id = _reference_pdb_id_from_bias_row(row)
    if not pdb_id:
        return [row_dict]

    generated_rows: list[dict[str, object]] = []
    if (
        row.get("pairing_status") == BIAS_PAIRING_PROTEIN_ONLY
        and pd.isna(row.get("ligand_pdb_id"))
    ):
        ligand_query_id = str(row.get("query_ligand_chain_id")).strip()
        ligand_query_smiles = ligand_queries.get(ligand_query_id)
        if ligand_query_smiles:
            ligand_matches = _reference_ligand_similarity_rows(
                pdb_id=pdb_id,
                query_smiles=ligand_query_smiles,
                ligand_reference_df=ligand_reference_df,
                ligand_reference_index=ligand_reference_index,
            )
            if not ligand_matches:
                ligand_matches = _pdb_ligand_similarity_rows(
                    pdb_id=pdb_id,
                    query_smiles=ligand_query_smiles,
                    boltz_cache_path=boltz_cache_path,
                    components_cif_path=components_cif_path,
                    timeout=timeout,
                )
            for ligand_match in ligand_matches:
                updated_row = dict(row_dict)
                updated_row["pairing_status"] = BIAS_PAIRING_PAIRED
                updated_row["pdb_id"] = pdb_id
                updated_row["ligand_pdb_id"] = pdb_id
                updated_row["ligand_release_date"] = updated_row.get("protein_release_date", pd.NA)
                updated_row["ligand_source"] = updated_row.get("ligand_source", pd.NA)
                if pd.isna(updated_row["ligand_source"]) or not str(updated_row["ligand_source"]).strip():
                    updated_row["ligand_source"] = "public"
                updated_row["ligand_dataset_name"] = updated_row.get("ligand_dataset_name", pd.NA)
                if pd.isna(updated_row["ligand_dataset_name"]) or not str(updated_row["ligand_dataset_name"]).strip():
                    updated_row["ligand_dataset_name"] = "public"
                updated_row["ligand_id"] = ligand_match["ligand_id"]
                updated_row["smiles"] = ligand_match["smiles"]
                updated_row["ecfp_similarity"] = ligand_match["ecfp_similarity"]
                updated_row["plot_ecfp_similarity"] = float(ligand_match["ecfp_similarity"])
                updated_row["source"] = _collapse_pair_value(
                    [updated_row.get("protein_source"), updated_row.get("ligand_source")],
                    mixed_label="mixed",
                )
                updated_row["dataset_name"] = _collapse_pair_value(
                    [updated_row.get("protein_dataset_name"), updated_row.get("ligand_dataset_name")]
                )
                generated_rows.append(updated_row)

    if (
        row.get("pairing_status") == BIAS_PAIRING_LIGAND_ONLY
        and pd.isna(row.get("protein_pdb_id"))
    ):
        protein_query_id = str(row.get("query_protein_chain_id")).strip()
        protein_query_sequence = protein_queries.get(protein_query_id)
        if protein_query_sequence:
            protein_matches = _reference_protein_similarity_rows(
                pdb_id=pdb_id,
                protein_lookup_index=protein_lookup_index,
            )
            used_protein_lookup = bool(protein_matches)
            if not protein_matches:
                if not _resolved_mmseqs_bin_token():
                    logging.getLogger(__name__).warning(
                        "MMseqs protein fallback unavailable for mixed bias row %s "
                        "(query_protein_chain_id=%s, pdb_id=%s); leaving row ligand_only",
                        str(row.get("query_pair_id")).strip(),
                        protein_query_id,
                        pdb_id,
                    )
                    return [row_dict]
                protein_matches = _pdb_protein_similarity_rows(
                    pdb_id=pdb_id,
                    query_sequence=protein_query_sequence,
                    timeout=timeout,
                )

            for protein_match in protein_matches:
                sequence_similarity = pd.to_numeric(
                    pd.Series([protein_match.get("sequence_similarity")]),
                    errors="coerce",
                ).iloc[0]
                if pd.isna(sequence_similarity):
                    continue
                updated_row = dict(row_dict)
                updated_row["pairing_status"] = BIAS_PAIRING_PAIRED
                updated_row["pdb_id"] = pdb_id
                updated_row["protein_pdb_id"] = pdb_id
                updated_row["protein_release_date"] = protein_match.get(
                    "release_date",
                    updated_row.get("ligand_release_date", pd.NA),
                )
                updated_row["protein_source"] = protein_match.get(
                    "source",
                    updated_row.get("protein_source", pd.NA),
                )
                if pd.isna(updated_row["protein_source"]) or not str(updated_row["protein_source"]).strip():
                    updated_row["protein_source"] = "public"
                updated_row["protein_dataset_name"] = protein_match.get(
                    "dataset_name",
                    updated_row.get("protein_dataset_name", pd.NA),
                )
                if pd.isna(updated_row["protein_dataset_name"]) or not str(updated_row["protein_dataset_name"]).strip():
                    updated_row["protein_dataset_name"] = "public"
                updated_row["protein_source_structure_path"] = protein_match.get(
                    "source_structure_path",
                    updated_row.get("protein_source_structure_path", pd.NA),
                )
                updated_row["protein_source_reference_path"] = protein_match.get(
                    "source_reference_path",
                    updated_row.get("protein_source_reference_path", pd.NA),
                )
                updated_row["sequence"] = protein_match["sequence"]
                updated_row["sequence_similarity"] = float(sequence_similarity)
                updated_row["sequence_similarity_pairwise"] = pd.NA
                updated_row["sequence_similarity_method"] = "mmseqs_pident"
                updated_row["plot_sequence_similarity"] = float(sequence_similarity) / 100.0
                updated_row["plot_sequence_similarity_pairwise"] = pd.NA
                updated_row["source"] = _collapse_pair_value(
                    [updated_row.get("protein_source"), updated_row.get("ligand_source")],
                    mixed_label="mixed",
                )
                updated_row["dataset_name"] = _collapse_pair_value(
                    [updated_row.get("protein_dataset_name"), updated_row.get("ligand_dataset_name")]
                )
                if not used_protein_lookup and float(sequence_similarity) > BIAS_PROTEIN_VIEW_MIN_SIMILARITY:
                    logging.getLogger(__name__).warning(
                        "Protein MMseqs fallback exceeded reference threshold for mixed bias row %s "
                        "(query_protein_chain_id=%s, pdb_id=%s, sequence_similarity=%.3f); "
                        "PDB absent from protein_training_data reference view",
                        str(row.get("query_pair_id")).strip(),
                        protein_query_id,
                        pdb_id,
                        float(sequence_similarity),
                    )
                generated_rows.append(updated_row)

    return generated_rows or [row_dict]


def _enrich_same_type_ligand_pair_row_with_pdb_backfill(
    *,
    row_dict: dict[str, object],
    ligand_queries: dict[str, str | None],
    boltz_cache_path: Path,
    components_cif_path: Path | None,
    ligand_reference_df: pd.DataFrame | None = None,
    ligand_reference_index: dict[str, list[dict[str, str]]] | None = None,
    timeout: int = 30,
) -> list[dict[str, object]]:
    row = pd.Series(row_dict)
    pdb_id = _reference_pdb_id_from_same_type_row(row)
    if not pdb_id:
        return [row_dict]

    generated_rows: list[dict[str, object]] = []
    if row.get("pairing_status") == BIAS_PAIRING_COMPONENT_1_ONLY:
        component_1_similarity = pd.to_numeric(
            pd.Series([row.get("component_1_similarity")]),
            errors="coerce",
        ).iloc[0]
        component_2_query_id = str(row.get("component_2_id")).strip()
        component_2_query_smiles = ligand_queries.get(component_2_query_id)
        if (
            pd.notna(component_1_similarity)
            and float(component_1_similarity) > BIAS_LIGAND_VIEW_MIN_SIMILARITY
            and component_2_query_smiles
        ):
            exclude_ligand_ids = (
                {str(row.get("_component_1_ligand_id")).strip()}
                if pd.notna(row.get("_component_1_ligand_id"))
                and str(row.get("_component_1_ligand_id")).strip()
                else set()
            )
            ligand_matches = _reference_ligand_similarity_rows(
                pdb_id=pdb_id,
                query_smiles=component_2_query_smiles,
                ligand_reference_df=ligand_reference_df,
                ligand_reference_index=ligand_reference_index,
                exclude_ligand_ids=exclude_ligand_ids,
            )
            if not ligand_matches:
                ligand_matches = _pdb_ligand_similarity_rows(
                    pdb_id=str(pdb_id),
                    query_smiles=str(component_2_query_smiles),
                    boltz_cache_path=boltz_cache_path,
                    components_cif_path=components_cif_path,
                    exclude_ligand_ids=set(exclude_ligand_ids),
                    timeout=timeout,
                )
            for ligand_match in ligand_matches:
                updated_row = dict(row_dict)
                updated_row["pairing_status"] = BIAS_PAIRING_PAIRED
                updated_row["reference_pdb_id"] = pdb_id
                if pd.isna(updated_row.get("component_2_reference_key")) or not str(updated_row.get("component_2_reference_key")).strip():
                    updated_row["component_2_reference_key"] = updated_row.get("reference_key")
                updated_row["component_2_pdb_id"] = pdb_id
                updated_row["component_2_release_date"] = updated_row.get("component_1_release_date", pd.NA)
                if pd.isna(updated_row.get("component_2_source")) or not str(updated_row.get("component_2_source")).strip():
                    updated_row["component_2_source"] = "public"
                if pd.isna(updated_row.get("component_2_dataset_name")) or not str(updated_row.get("component_2_dataset_name")).strip():
                    updated_row["component_2_dataset_name"] = "public"
                updated_row["component_2_similarity"] = ligand_match["ecfp_similarity"]
                updated_row["plot_component_2_similarity"] = float(ligand_match["ecfp_similarity"])
                updated_row["_component_2_ligand_id"] = ligand_match["ligand_id"]
                updated_row["_component_2_smiles"] = ligand_match["smiles"]
                updated_row["source"] = _collapse_pair_value(
                    [updated_row.get("component_1_source"), updated_row.get("component_2_source")],
                    mixed_label="mixed",
                )
                updated_row["dataset_name"] = _collapse_pair_value(
                    [updated_row.get("component_1_dataset_name"), updated_row.get("component_2_dataset_name")]
                )
                generated_rows.append(updated_row)

    if row.get("pairing_status") == BIAS_PAIRING_COMPONENT_2_ONLY:
        component_2_similarity = pd.to_numeric(
            pd.Series([row.get("component_2_similarity")]),
            errors="coerce",
        ).iloc[0]
        component_1_query_id = str(row.get("component_1_id")).strip()
        component_1_query_smiles = ligand_queries.get(component_1_query_id)
        if (
            pd.notna(component_2_similarity)
            and float(component_2_similarity) > BIAS_LIGAND_VIEW_MIN_SIMILARITY
            and component_1_query_smiles
        ):
            exclude_ligand_ids = (
                {str(row.get("_component_2_ligand_id")).strip()}
                if pd.notna(row.get("_component_2_ligand_id"))
                and str(row.get("_component_2_ligand_id")).strip()
                else set()
            )
            ligand_matches = _reference_ligand_similarity_rows(
                pdb_id=pdb_id,
                query_smiles=component_1_query_smiles,
                ligand_reference_df=ligand_reference_df,
                ligand_reference_index=ligand_reference_index,
                exclude_ligand_ids=exclude_ligand_ids,
            )
            if not ligand_matches:
                ligand_matches = _pdb_ligand_similarity_rows(
                    pdb_id=str(pdb_id),
                    query_smiles=str(component_1_query_smiles),
                    boltz_cache_path=boltz_cache_path,
                    components_cif_path=components_cif_path,
                    exclude_ligand_ids=set(exclude_ligand_ids),
                    timeout=timeout,
                )
            for ligand_match in ligand_matches:
                updated_row = dict(row_dict)
                updated_row["pairing_status"] = BIAS_PAIRING_PAIRED
                updated_row["reference_pdb_id"] = pdb_id
                if pd.isna(updated_row.get("component_1_reference_key")) or not str(updated_row.get("component_1_reference_key")).strip():
                    updated_row["component_1_reference_key"] = updated_row.get("reference_key")
                updated_row["component_1_pdb_id"] = pdb_id
                updated_row["component_1_release_date"] = updated_row.get("component_2_release_date", pd.NA)
                if pd.isna(updated_row.get("component_1_source")) or not str(updated_row.get("component_1_source")).strip():
                    updated_row["component_1_source"] = "public"
                if pd.isna(updated_row.get("component_1_dataset_name")) or not str(updated_row.get("component_1_dataset_name")).strip():
                    updated_row["component_1_dataset_name"] = "public"
                updated_row["component_1_similarity"] = ligand_match["ecfp_similarity"]
                updated_row["plot_component_1_similarity"] = float(ligand_match["ecfp_similarity"])
                updated_row["_component_1_ligand_id"] = ligand_match["ligand_id"]
                updated_row["_component_1_smiles"] = ligand_match["smiles"]
                updated_row["source"] = _collapse_pair_value(
                    [updated_row.get("component_1_source"), updated_row.get("component_2_source")],
                    mixed_label="mixed",
                )
                updated_row["dataset_name"] = _collapse_pair_value(
                    [updated_row.get("component_1_dataset_name"), updated_row.get("component_2_dataset_name")]
                )
                generated_rows.append(updated_row)

    return generated_rows or [row_dict]


def _enrich_mixed_bias_dataset_with_pdb_backfill(
    df: pd.DataFrame,
    *,
    protein_queries: dict[str, str | None],
    ligand_queries: dict[str, str | None],
    protein_lookup_view: pd.DataFrame | None = None,
    boltz_cache_path: Path,
    components_cif_path: Path | None,
    ligand_reference_df: pd.DataFrame | None = None,
    ligand_reference_index: dict[str, list[dict[str, str]]] | None = None,
    progress_output_path: Path | None = None,
    progress_label: str | None = None,
    timeout: int = 30,
) -> pd.DataFrame:
    if df.empty:
        return _finalize_bias_training_dataframe(df)

    input_rows = [row.to_dict() for _, row in df.iterrows()]
    total_lookups = 0
    for row_dict in input_rows:
        row = pd.Series(row_dict)
        pdb_id = _reference_pdb_id_from_bias_row(row)
        if (
            row.get("pairing_status") == BIAS_PAIRING_PROTEIN_ONLY
            and pd.isna(row.get("ligand_pdb_id"))
            and pdb_id
            and ligand_queries.get(str(row.get("query_ligand_chain_id")).strip())
        ):
            total_lookups += 1
        elif (
            row.get("pairing_status") == BIAS_PAIRING_LIGAND_ONLY
            and pd.isna(row.get("protein_pdb_id"))
            and pdb_id
            and protein_queries.get(str(row.get("query_protein_chain_id")).strip())
        ):
            total_lookups += 1

    progress_name = progress_label or "mixed_bias_backfill"
    protein_lookup_index = _build_protein_lookup_index(protein_lookup_view)
    _log_bias_lookup_progress(
        progress_label=progress_name,
        completed_lookups=0,
        total_lookups=total_lookups,
    )

    enriched_rows: list[dict[str, object]] = []
    completed_lookups = 0
    for index, row_dict in enumerate(input_rows):
        row = pd.Series(row_dict)
        pdb_id = _reference_pdb_id_from_bias_row(row)
        lookup_attempted = False

        if (
            row.get("pairing_status") == BIAS_PAIRING_PROTEIN_ONLY
            and pd.isna(row.get("ligand_pdb_id"))
            and pdb_id
            and ligand_queries.get(str(row.get("query_ligand_chain_id")).strip())
        ):
            lookup_attempted = True
        elif (
            row.get("pairing_status") == BIAS_PAIRING_LIGAND_ONLY
            and pd.isna(row.get("protein_pdb_id"))
            and pdb_id
            and protein_queries.get(str(row.get("query_protein_chain_id")).strip())
        ):
            lookup_attempted = True

        enriched_rows.extend(
            _enrich_mixed_bias_row_with_pdb_backfill(
                row_dict=row_dict,
                protein_queries=protein_queries,
                ligand_queries=ligand_queries,
                protein_lookup_index=protein_lookup_index,
                boltz_cache_path=boltz_cache_path,
                components_cif_path=components_cif_path,
                ligand_reference_df=ligand_reference_df,
                ligand_reference_index=ligand_reference_index,
                timeout=timeout,
            )
        )

        if lookup_attempted:
            completed_lookups += 1
            if (
                completed_lookups % BIAS_PROGRESS_LOG_EVERY == 0
                or completed_lookups == total_lookups
            ):
                _log_bias_lookup_progress(
                    progress_label=progress_name,
                    completed_lookups=completed_lookups,
                    total_lookups=total_lookups,
                )
            if (
                completed_lookups % BIAS_PROGRESS_WRITE_EVERY == 0
                or completed_lookups == total_lookups
            ):
                _write_progress_checkpoint(
                    rows=enriched_rows,
                    pending_rows=input_rows[index + 1 :],
                    output_path=progress_output_path,
                    finalize=_finalize_bias_training_dataframe,
                    progress_label=progress_name,
                    completed_lookups=completed_lookups,
                    total_lookups=total_lookups,
                )

    return _finalize_bias_training_dataframe(pd.DataFrame(enriched_rows))


def _component_similarity_column(component_type: str) -> str:
    return "sequence_similarity" if component_type == "protein" else "ecfp_similarity"


def _component_threshold(component_type: str) -> float:
    return BIAS_PLOT_SEQUENCE_THRESHOLD if component_type == "protein" else BIAS_PLOT_LIGAND_THRESHOLD


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
    for row, suffix in ((component_1_row, "component_1"), (component_2_row, "component_2")):
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
            candidate = Path(text).stem if "/" in text or text.endswith((".pdb", ".cif", ".sdf", ".csv")) else text
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
                forced_component_1_similarity
                if forced_component_1_similarity is not None
                else (
                    component_1_row.get(f"{similarity_column}_component_1")
                    if component_1_row is not None
                    else pd.NA
                )
            ]
        ),
        errors="coerce",
    ).iloc[0]
    component_2_similarity = pd.to_numeric(
        pd.Series(
            [
                forced_component_2_similarity
                if forced_component_2_similarity is not None
                else (
                    component_2_row.get(f"{similarity_column}_component_2")
                    if component_2_row is not None
                    else pd.NA
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
                component_1_row.get("pdb_id_component_1") if component_1_row is not None else None,
                component_2_row.get("pdb_id_component_2") if component_2_row is not None else None,
            ]
        ),
        "pairing_status": pairing_status,
        "source": _collapse_pair_value(
            [
                component_1_row.get("source_component_1") if component_1_row is not None else None,
                component_2_row.get("source_component_2") if component_2_row is not None else None,
            ],
            mixed_label="mixed",
        ),
        "dataset_name": _collapse_pair_value(
            [
                component_1_row.get("dataset_name_component_1") if component_1_row is not None else None,
                component_2_row.get("dataset_name_component_2") if component_2_row is not None else None,
            ]
        ),
        "component_1_reference_key": component_1_reference_key,
        "component_1_pdb_id": component_1_row.get("pdb_id_component_1") if component_1_row is not None else pd.NA,
        "component_1_release_date": (
            component_1_row.get("release_date_component_1") if component_1_row is not None else pd.NA
        ),
        "component_1_source": component_1_row.get("source_component_1") if component_1_row is not None else pd.NA,
        "component_1_dataset_name": (
            component_1_row.get("dataset_name_component_1") if component_1_row is not None else pd.NA
        ),
        "component_1_source_structure_path": (
            component_1_row.get("source_structure_path_component_1") if component_1_row is not None else pd.NA
        ),
        "component_1_source_reference_path": (
            component_1_row.get("source_reference_path_component_1") if component_1_row is not None else pd.NA
        ),
        "component_2_reference_key": component_2_reference_key,
        "component_2_pdb_id": component_2_row.get("pdb_id_component_2") if component_2_row is not None else pd.NA,
        "component_2_release_date": (
            component_2_row.get("release_date_component_2") if component_2_row is not None else pd.NA
        ),
        "component_2_source": component_2_row.get("source_component_2") if component_2_row is not None else pd.NA,
        "component_2_dataset_name": (
            component_2_row.get("dataset_name_component_2") if component_2_row is not None else pd.NA
        ),
        "component_2_source_structure_path": (
            component_2_row.get("source_structure_path_component_2") if component_2_row is not None else pd.NA
        ),
        "component_2_source_reference_path": (
            component_2_row.get("source_reference_path_component_2") if component_2_row is not None else pd.NA
        ),
        "component_1_similarity": float(component_1_similarity) if pd.notna(component_1_similarity) else pd.NA,
        "component_2_similarity": float(component_2_similarity) if pd.notna(component_2_similarity) else pd.NA,
        "plot_component_1_similarity": _normalize_component_similarity(component_1_similarity, component_type),
        "plot_component_2_similarity": _normalize_component_similarity(component_2_similarity, component_type),
        "_component_1_ligand_id": (
            component_1_row.get("ligand_id_component_1") if component_1_row is not None else pd.NA
        ),
        "_component_1_smiles": (
            component_1_row.get("smiles_component_1") if component_1_row is not None else pd.NA
        ),
        "_component_2_ligand_id": (
            component_2_row.get("ligand_id_component_2") if component_2_row is not None else pd.NA
        ),
        "_component_2_smiles": (
            component_2_row.get("smiles_component_2") if component_2_row is not None else pd.NA
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
    ordered["_pairing_order"] = ordered["pairing_status"].map(_same_type_pair_sort_order()).fillna(99)
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
    protein_reference_key = protein_row.get("reference_key") if protein_row is not None else None
    ligand_reference_key = ligand_row.get("reference_key") if ligand_row is not None else None
    protein_source = protein_row.get("source_protein") if protein_row is not None else None
    ligand_source = ligand_row.get("source_ligand") if ligand_row is not None else None

    sequence_similarity = pd.to_numeric(
        pd.Series(
            [protein_row.get("sequence_similarity_protein") if protein_row is not None else pd.NA]
        ),
        errors="coerce",
    ).iloc[0]
    sequence_similarity_pairwise = pd.to_numeric(
        pd.Series(
            [
                protein_row.get("sequence_similarity_pairwise_protein")
                if protein_row is not None
                else pd.NA
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
    if pd.isna(sequence_similarity_method) or not str(sequence_similarity_method).strip():
        if pd.notna(sequence_similarity):
            sequence_similarity_method = "mmseqs_pident"
        elif pd.notna(sequence_similarity_pairwise):
            sequence_similarity_method = "pairwise_aligner"
        else:
            sequence_similarity_method = "unavailable"
    if (
        pd.notna(sequence_similarity)
        and sequence_similarity_method != "mmseqs_pident"
    ) or (
        pd.notna(sequence_similarity_pairwise)
        and sequence_similarity_method != "pairwise_aligner"
    ) or (
        pd.isna(sequence_similarity)
        and pd.isna(sequence_similarity_pairwise)
        and sequence_similarity_method != "unavailable"
    ):
        raise ValueError(
            "Protein similarity value does not match sequence_similarity_method"
        )
    ecfp_similarity = pd.to_numeric(
        pd.Series(
            [
                forced_ecfp_similarity
                if forced_ecfp_similarity is not None
                else (ligand_row.get("ecfp_similarity_ligand") if ligand_row is not None else pd.NA)
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
    plot_ecfp_similarity = float(ecfp_similarity) if pd.notna(ecfp_similarity) else pd.NA

    protein_is_public_pdb_reference = (
        pd.notna(protein_reference_key) and str(protein_reference_key).startswith("public:pdb:")
    )
    ligand_is_public_pdb_reference = (
        pd.notna(ligand_reference_key) and str(ligand_reference_key).startswith("public:pdb:")
    )
    if (
        protein_row is not None
        and protein_is_public_pdb_reference
        and protein_unique_counts.get(protein_reference_key, 0) > 1
    ):
        sequence_value = pd.NA
    else:
        sequence_value = protein_row.get("sequence_protein") if protein_row is not None else pd.NA
    if (
        ligand_row is not None
        and ligand_is_public_pdb_reference
        and ligand_id_unique_counts.get(ligand_reference_key, 0) > 1
    ):
        ligand_id_value = pd.NA
    else:
        ligand_id_value = ligand_row.get("ligand_id_ligand") if ligand_row is not None else pd.NA
    if (
        ligand_row is not None
        and ligand_is_public_pdb_reference
        and ligand_smiles_unique_counts.get(ligand_reference_key, 0) > 1
    ):
        smiles_value = pd.NA
    else:
        smiles_value = ligand_row.get("smiles_ligand") if ligand_row is not None else pd.NA

    return {
        "query_pair_id": f"{protein_chain or 'none'}__{ligand_chain or 'none'}",
        "query_protein_chain_id": protein_chain,
        "query_ligand_chain_id": ligand_chain,
        "reference_key": reference_key,
        "reference_label": reference_label,
        "pairing_status": pairing_status,
        "source": _collapse_pair_value([protein_source, ligand_source], mixed_label="mixed"),
        "dataset_name": _collapse_pair_value(
            [
                protein_row.get("dataset_name_protein") if protein_row is not None else None,
                ligand_row.get("dataset_name_ligand") if ligand_row is not None else None,
            ]
        ),
        "complex_id": _collapse_pair_value(
            [
                protein_row.get("complex_id_protein") if protein_row is not None else None,
                ligand_row.get("complex_id_ligand") if ligand_row is not None else None,
            ]
        ),
        "pdb_id": _collapse_pair_value(
            [
                protein_row.get("pdb_id_protein") if protein_row is not None else None,
                ligand_row.get("pdb_id_ligand") if ligand_row is not None else None,
            ]
        ),
        "protein_pdb_id": protein_row.get("pdb_id_protein") if protein_row is not None else pd.NA,
        "ligand_pdb_id": ligand_row.get("pdb_id_ligand") if ligand_row is not None else pd.NA,
        "protein_release_date": protein_row.get("release_date_protein") if protein_row is not None else pd.NA,
        "ligand_release_date": ligand_row.get("release_date_ligand") if ligand_row is not None else pd.NA,
        "protein_source": protein_source,
        "protein_dataset_name": protein_row.get("dataset_name_protein") if protein_row is not None else pd.NA,
        "protein_source_structure_path": protein_row.get("source_structure_path_protein") if protein_row is not None else pd.NA,
        "protein_source_reference_path": protein_row.get("source_reference_path_protein") if protein_row is not None else pd.NA,
        "ligand_source": ligand_source,
        "ligand_dataset_name": ligand_row.get("dataset_name_ligand") if ligand_row is not None else pd.NA,
        "ligand_source_structure_path": ligand_row.get("source_structure_path_ligand") if ligand_row is not None else pd.NA,
        "ligand_source_reference_path": ligand_row.get("source_reference_path_ligand") if ligand_row is not None else pd.NA,
        "sequence_similarity": float(sequence_similarity) if pd.notna(sequence_similarity) else pd.NA,
        "sequence_similarity_pairwise": (
            float(sequence_similarity_pairwise)
            if pd.notna(sequence_similarity_pairwise)
            else pd.NA
        ),
        "sequence_similarity_method": sequence_similarity_method,
        "ecfp_similarity": float(ecfp_similarity) if pd.notna(ecfp_similarity) else pd.NA,
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
        left_lookup.groupby("reference_key")["sequence_protein"].nunique(dropna=True).to_dict()
        if "sequence_protein" in left_lookup.columns
        else {}
    )
    ligand_id_unique_counts = (
        right_lookup.groupby("reference_key")["ligand_id_ligand"].nunique(dropna=True).to_dict()
        if "ligand_id_ligand" in right_lookup.columns
        else {}
    )
    ligand_smiles_unique_counts = (
        right_lookup.groupby("reference_key")["smiles_ligand"].nunique(dropna=True).to_dict()
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
    left = _prepare_reference_side(component_1_view, is_ligand=is_ligand, suffix="component_1")
    right = _prepare_reference_side(component_2_view, is_ligand=is_ligand, suffix="component_2")
    left_lookup = _prepare_reference_side(
        component_1_lookup_view if component_1_lookup_view is not None else component_1_view,
        is_ligand=is_ligand,
        suffix="component_1",
    )
    right_lookup = _prepare_reference_side(
        component_2_lookup_view if component_2_lookup_view is not None else component_2_view,
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
                        reference_label=_reference_label_from_component_rows(combined_row, combined_row),
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
                        reference_label=_reference_label_from_component_rows(combined_row, combined_row),
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


def _enrich_same_type_ligand_pair_dataset_with_pdb_backfill(
    df: pd.DataFrame,
    *,
    ligand_queries: dict[str, str | None],
    boltz_cache_path: Path,
    components_cif_path: Path | None,
    ligand_reference_df: pd.DataFrame | None = None,
    ligand_reference_index: dict[str, list[dict[str, str]]] | None = None,
    progress_output_path: Path | None = None,
    progress_label: str | None = None,
    timeout: int = 30,
) -> pd.DataFrame:
    if df.empty:
        return df.copy()

    input_rows = [row.to_dict() for _, row in df.iterrows()]
    total_lookups = 0
    for row_dict in input_rows:
        row = pd.Series(row_dict)
        pdb_id = _reference_pdb_id_from_same_type_row(row)
        if row.get("pairing_status") == BIAS_PAIRING_COMPONENT_1_ONLY and pdb_id:
            component_1_similarity = pd.to_numeric(
                pd.Series([row.get("component_1_similarity")]),
                errors="coerce",
            ).iloc[0]
            component_2_query_id = str(row.get("component_2_id")).strip()
            if (
                pd.notna(component_1_similarity)
                and float(component_1_similarity) > BIAS_LIGAND_VIEW_MIN_SIMILARITY
                and ligand_queries.get(component_2_query_id)
            ):
                total_lookups += 1
        elif row.get("pairing_status") == BIAS_PAIRING_COMPONENT_2_ONLY and pdb_id:
            component_2_similarity = pd.to_numeric(
                pd.Series([row.get("component_2_similarity")]),
                errors="coerce",
            ).iloc[0]
            component_1_query_id = str(row.get("component_1_id")).strip()
            if (
                pd.notna(component_2_similarity)
                and float(component_2_similarity) > BIAS_LIGAND_VIEW_MIN_SIMILARITY
                and ligand_queries.get(component_1_query_id)
            ):
                total_lookups += 1

    progress_name = progress_label or "ligand_pair_backfill"
    _log_bias_lookup_progress(
        progress_label=progress_name,
        completed_lookups=0,
        total_lookups=total_lookups,
    )

    enriched_rows: list[dict[str, object]] = []
    completed_lookups = 0
    for index, row_dict in enumerate(input_rows):
        row = pd.Series(row_dict)
        pdb_id = _reference_pdb_id_from_same_type_row(row)
        lookup_attempted = False

        if row.get("pairing_status") == BIAS_PAIRING_COMPONENT_1_ONLY and pdb_id:
            component_1_similarity = pd.to_numeric(
                pd.Series([row.get("component_1_similarity")]),
                errors="coerce",
            ).iloc[0]
            component_2_query_id = str(row.get("component_2_id")).strip()
            if (
                pd.notna(component_1_similarity)
                and float(component_1_similarity) > BIAS_LIGAND_VIEW_MIN_SIMILARITY
                and ligand_queries.get(component_2_query_id)
            ):
                lookup_attempted = True
        elif row.get("pairing_status") == BIAS_PAIRING_COMPONENT_2_ONLY and pdb_id:
            component_2_similarity = pd.to_numeric(
                pd.Series([row.get("component_2_similarity")]),
                errors="coerce",
            ).iloc[0]
            component_1_query_id = str(row.get("component_1_id")).strip()
            if (
                pd.notna(component_2_similarity)
                and float(component_2_similarity) > BIAS_LIGAND_VIEW_MIN_SIMILARITY
                and ligand_queries.get(component_1_query_id)
            ):
                lookup_attempted = True

        enriched_rows.extend(
            _enrich_same_type_ligand_pair_row_with_pdb_backfill(
                row_dict=row_dict,
                ligand_queries=ligand_queries,
                boltz_cache_path=boltz_cache_path,
                components_cif_path=components_cif_path,
                ligand_reference_df=ligand_reference_df,
                ligand_reference_index=ligand_reference_index,
                timeout=timeout,
            )
        )

        if lookup_attempted:
            completed_lookups += 1
            if (
                completed_lookups % BIAS_PROGRESS_LOG_EVERY == 0
                or completed_lookups == total_lookups
            ):
                _log_bias_lookup_progress(
                    progress_label=progress_name,
                    completed_lookups=completed_lookups,
                    total_lookups=total_lookups,
                )
            if (
                completed_lookups % BIAS_PROGRESS_WRITE_EVERY == 0
                or completed_lookups == total_lookups
            ):
                _write_progress_checkpoint(
                    rows=enriched_rows,
                    pending_rows=input_rows[index + 1 :],
                    output_path=progress_output_path,
                    finalize=_finalize_same_type_pair_dataframe,
                    progress_label=progress_name,
                    completed_lookups=completed_lookups,
                    total_lookups=total_lookups,
                )

    return pd.DataFrame(enriched_rows)


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
            for value in protein_view.get("query_chain_id", pd.Series(dtype=object)).dropna()
            if str(value).strip()
        }
    )
    ligand_chain_ids = sorted(str(chain_id).strip() for chain_id in ligand_views if str(chain_id).strip())

    if not protein_chain_ids and not ligand_chain_ids:
        return pd.DataFrame(columns=BIAS_TRAINING_DATA_COLUMNS)

    rows: list[dict[str, object]] = []
    protein_iter = protein_chain_ids or [None]
    ligand_iter = ligand_chain_ids or [None]

    for protein_chain_id in protein_iter:
        protein_subset = (
            protein_view[protein_view["query_chain_id"].astype(str) == protein_chain_id].copy()
            if protein_chain_id is not None and "query_chain_id" in protein_view.columns
            else protein_view.iloc[0:0].copy()
        )
        protein_lookup_subset = (
            protein_lookup_view[protein_lookup_view["query_chain_id"].astype(str) == protein_chain_id].copy()
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
    ordered["_pairing_order"] = ordered["pairing_status"].map(_pairing_sort_order()).fillna(99)
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


def _materialize_bias_training_views(
    chain_df: pd.DataFrame,
    ligands_df: pd.DataFrame,
    proteins_df: pd.DataFrame,
    protein_queries: dict[str, str | None],
    ligand_queries: dict[str, str | None],
    boltz_cache_path: Path,
    output_dir: Path,
    components_cif_path: Path | None = None,
) -> tuple[pd.DataFrame, dict[str, pd.DataFrame], pd.DataFrame, dict[str, pd.DataFrame]]:
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
        ligand_view.drop(columns=["query_chain_id"], errors="ignore").to_csv(out_path, index=False)

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
        for value in chain_df.loc[chain_df["ENTITY_TYPE"].astype(str) == "protein", "CHAIN_ID"].astype(str)
        if str(value).strip()
    }
    grouped: dict[str, dict[str, object]] = {}
    for raw_chain_id, query_sequence in protein_queries.items():
        chain_id = str(raw_chain_id).strip().upper()
        if not chain_id or chain_id not in protein_chain_ids:
            continue
        identity_key = f"sequence:{str(query_sequence).strip()}" if query_sequence else f"chain:{chain_id}"
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
                character.isalnum() or character in {"_", "-"} for character in identity_label
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
    frames = [views.get(chain_id, pd.DataFrame()).copy() for chain_id in sorted(chain_ids)]
    frames = [frame for frame in frames if not frame.empty]
    if not frames:
        return pd.DataFrame()
    subset = pd.concat(frames, ignore_index=True)
    return _dedupe_query_group_view(subset)


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
            boltz_cache_path=boltz_cache_path if boltz_cache_path is not None else Path("~/.boltz").expanduser(),
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
            protein_group = component_1 if component_1_type == "protein" else component_2
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

        pair_suffix = (
            f"{_safe_pair_file_token(component_1_label)}__{_safe_pair_file_token(component_2_label)}"
        )
        component_1_view = (
            _subset_group_protein_view(view=protein_view, chain_ids=component_1_chain_ids)
            if component_1_type == "protein"
            else _subset_group_ligand_view(views=ligand_views, chain_ids=component_1_chain_ids)
        )
        component_2_view = (
            _subset_group_protein_view(view=protein_view, chain_ids=component_2_chain_ids)
            if component_2_type == "protein"
            else _subset_group_ligand_view(views=ligand_views, chain_ids=component_2_chain_ids)
        )
        component_1_lookup = (
            _subset_group_protein_view(view=protein_lookup_view, chain_ids=component_1_chain_ids)
            if component_1_type == "protein"
            else _subset_group_ligand_view(views=ligand_lookup_views, chain_ids=component_1_chain_ids)
        )
        component_2_lookup = (
            _subset_group_protein_view(view=protein_lookup_view, chain_ids=component_2_chain_ids)
            if component_2_type == "protein"
            else _subset_group_ligand_view(views=ligand_lookup_views, chain_ids=component_2_chain_ids)
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


def apply_bias_metrics(
    system_df: pd.DataFrame,
    chain_df: pd.DataFrame,
    sys_obj,
    protein_training_data_path: Path | str | None,
    ligand_training_data_path: Path | str | None,
    custom_protein_reference_path: Path | str | None = None,
    custom_ligand_reference_path: Path | str | None = None,
    release_cutoff: str = "2023-06-01",
    bias_chains: set[str] | list[str] | None = None,
    protein_top_n: int = 100,
    protein_similarity_threshold: float = BIAS_PLOT_SEQUENCE_THRESHOLD,
    ligand_similarity_threshold: float = BIAS_PLOT_LIGAND_THRESHOLD,
    boltz_cache_path: Path | str | None = None,
    output_dir: Path | str | None = None,
    logger: logging.Logger | None = None,
    timings: DebugTimingCollector | None = None,
    components_cif_path: Path | str | None = None,
    strict_query_resolution: bool = False,
):
    """Compute lightweight chain-level bias metrics only."""
    if logger is None:
        logger = logging.getLogger(__name__)

    release_policy = parse_bias_release_policy(release_cutoff)
    cutoff = release_policy.cutoff
    protein_training_path = (
        Path(protein_training_data_path) if protein_training_data_path is not None else None
    )
    ligand_training_path = (
        Path(ligand_training_data_path) if ligand_training_data_path is not None else None
    )
    custom_protein_path = (
        Path(custom_protein_reference_path) if custom_protein_reference_path is not None else None
    )
    custom_ligand_path = (
        Path(custom_ligand_reference_path) if custom_ligand_reference_path is not None else None
    )
    proteins_df = _load_protein_training(
        protein_training_path,
        cutoff,
        top_n=protein_top_n,
        custom_reference_path=custom_protein_path,
    )
    ligands_df = _load_ligand_training(
        ligand_training_path,
        cutoff,
        custom_reference_path=custom_ligand_path,
    )
    cache_path = Path(boltz_cache_path).expanduser() if boltz_cache_path else Path("~/.boltz").expanduser()
    components_path = (
        Path(components_cif_path).expanduser() if components_cif_path is not None else None
    )
    logger.info(
        "Bias release policy applied (%s): protein_rows=%d ligand_rows=%d",
        release_policy.label,
        len(proteins_df),
        len(ligands_df),
    )
    selected_chains = {
        str(x).strip().upper()
        for x in (bias_chains or [])
        if str(x).strip()
    }
    if selected_chains:
        logger.info("Bias restricted to selected chains: %s", sorted(selected_chains))

    protein_queries, ligand_queries, unresolved_ccd_by_chain = _extract_system_queries(
        sys_obj,
        ligands_df,
        boltz_cache_path=cache_path,
        components_cif_path=components_path,
    )
    if selected_chains:
        protein_queries = {
            cid: seq
            for cid, seq in protein_queries.items()
            if str(cid).strip().upper() in selected_chains
        }
        ligand_queries = {
            cid: smi
            for cid, smi in ligand_queries.items()
            if str(cid).strip().upper() in selected_chains
        }
        unresolved_ccd_by_chain = {
            cid: ccd_ids
            for cid, ccd_ids in unresolved_ccd_by_chain.items()
            if str(cid).strip().upper() in selected_chains
        }
    if strict_query_resolution and unresolved_ccd_by_chain:
        details = ", ".join(
            f"{chain_id} ({'/'.join(ccd_ids)})"
            for chain_id, ccd_ids in sorted(unresolved_ccd_by_chain.items())
        )
        raise ValueError(
            "Could not resolve SMILES for CCD-backed ligand chain(s): "
            f"{details}. Provide --bias_training_components_cif, use ligand training data "
            "with matching ligand_id values, or prepare the local CCD cache."
        )

    protein_best: dict[str, float | None] = {}
    protein_pairwise_best: dict[str, float | None] = {}
    protein_timer = (
        timings.measure("scores.bias_metrics.protein_similarity", logger=logger)
        if timings is not None
        else nullcontext()
    )
    with protein_timer:
        for chain_id, query_seq in protein_queries.items():
            normalized_chain_id = str(chain_id).strip().upper()
            query_references = _reference_rows_for_query(proteins_df, chain_id)
            protein_best[normalized_chain_id] = (
                _best_protein_hit(query_seq, query_references) if query_seq else None
            )
            protein_pairwise_best[normalized_chain_id] = (
                _best_pairwise_protein_hit(query_seq, query_references)
                if query_seq
                else None
            )

    ligand_best: dict[str, float | None] = {}
    ligand_timer = (
        timings.measure("scores.bias_metrics.ligand_similarity", logger=logger)
        if timings is not None
        else nullcontext()
    )
    with ligand_timer:
        for chain_id, query_smiles in ligand_queries.items():
            query_references = _reference_rows_for_query(ligands_df, chain_id)
            ligand_best[str(chain_id).strip().upper()] = (
                _best_ligand_hit(query_smiles, query_references)
                if query_smiles
                else None
            )

    if "bias_prot_sim_train" not in chain_df.columns:
        chain_df["bias_prot_sim_train"] = pd.NA
    if "bias_prot_sim_train_pairwise" not in chain_df.columns:
        chain_df["bias_prot_sim_train_pairwise"] = pd.NA
    if "bias_lig_sim_train" not in chain_df.columns:
        chain_df["bias_lig_sim_train"] = pd.NA

    for idx, row in chain_df.iterrows():
        entity = str(row.get("ENTITY_TYPE"))
        chain_id = str(row.get("CHAIN_ID")).strip().upper()
        if selected_chains and chain_id not in selected_chains:
            continue
        if entity == "protein":
            chain_df.at[idx, "bias_prot_sim_train"] = protein_best.get(chain_id)
            chain_df.at[idx, "bias_prot_sim_train_pairwise"] = protein_pairwise_best.get(
                chain_id
            )
        elif entity == "ligand":
            sim = ligand_best.get(chain_id)
            mol_id = row.get("ligand_molecule_id")
            if pd.notna(mol_id) and str(mol_id).strip():
                ref = ligands_df[ligands_df.get("ligand_id").astype(str) == _norm_id(mol_id)]
                if not ref.empty:
                    query_smiles = ligand_queries.get(chain_id)
                    if query_smiles and ref["smiles"].astype(str).eq(str(query_smiles)).any():
                        sim = 1.0
                    elif "ecfp_similarity" in ref.columns:
                        ecfp = pd.to_numeric(ref["ecfp_similarity"], errors="coerce").dropna()
                        if len(ecfp):
                            sim = float(ecfp.max())
                    if sim is None:
                        query_smiles = str(ref.iloc[0]["smiles"])
                        sim = _best_ligand_hit(query_smiles, ligands_df)
            chain_df.at[idx, "bias_lig_sim_train"] = sim

    prot_vals = pd.to_numeric(chain_df.get("bias_prot_sim_train"), errors="coerce").dropna()
    prot_pairwise_vals = pd.to_numeric(
        chain_df.get("bias_prot_sim_train_pairwise"), errors="coerce"
    ).dropna()
    lig_vals = pd.to_numeric(chain_df.get("bias_lig_sim_train"), errors="coerce").dropna()

    if len(prot_vals):
        system_df["bias_prot_sim_train_max"] = float(prot_vals.max())
    if len(prot_pairwise_vals):
        system_df["bias_prot_sim_train_pairwise_max"] = float(
            prot_pairwise_vals.max()
        )
    if len(lig_vals):
        system_df["bias_lig_sim_train_max"] = float(lig_vals.max())

    if output_dir is not None:
        ligand_reference_index = _build_ligand_reference_index(ligands_df)
        output_chain_df = chain_df
        if selected_chains:
            output_chain_df = chain_df[
                chain_df["CHAIN_ID"].astype(str).str.strip().str.upper().isin(selected_chains)
            ].copy()
        output_root = Path(output_dir)
        protein_view, ligand_views, protein_lookup_view, ligand_lookup_views = _materialize_bias_training_views(
            chain_df=output_chain_df,
            ligands_df=ligands_df,
            proteins_df=_load_protein_training(
                protein_training_path,
                cutoff,
                top_n=0,
                custom_reference_path=custom_protein_path,
            ),
            protein_queries=protein_queries,
            ligand_queries=ligand_queries,
            boltz_cache_path=cache_path,
            output_dir=output_root,
            components_cif_path=components_path,
        )
        _clear_pair_specific_artifacts(output_root)
        pair_artifacts = _pair_specific_bias_training_datasets(
            chain_df=output_chain_df,
            protein_queries=protein_queries,
            ligand_queries=ligand_queries,
            ligand_reference_df=ligands_df,
            ligand_reference_index=ligand_reference_index,
            protein_view=protein_view,
            ligand_views=ligand_views,
            protein_lookup_view=protein_lookup_view,
            ligand_lookup_views=ligand_lookup_views,
            boltz_cache_path=cache_path,
            components_cif_path=components_path,
            output_dir=output_root,
        )
        bias_training_df = _combined_mixed_pair_bias_training_dataframe(pair_artifacts)
        bias_training_df.to_csv(output_root / "bias_training_data.csv", index=False)
        _write_bias_plot_artifacts(
            bias_training_df=bias_training_df,
            output_dir=output_root,
            file_stem=BIAS_PLOT_FILE_STEM,
            sequence_threshold=protein_similarity_threshold,
            ligand_threshold=ligand_similarity_threshold,
        )
        for pair_artifact in pair_artifacts:
            pair_df = pair_artifact["dataframe"]
            pair_df.to_csv(
                output_root / str(pair_artifact["data_file"]),
                index=False,
            )
            if pair_artifact["pair_type"] == "mixed":
                _write_bias_plot_artifacts(
                    bias_training_df=pair_df,
                    output_dir=output_root,
                    file_stem=str(pair_artifact["plot_file_stem"]),
                    sequence_threshold=protein_similarity_threshold,
                    ligand_threshold=ligand_similarity_threshold,
                )
            else:
                _write_same_type_plot_artifacts(
                    pair_df=pair_df,
                    output_dir=output_root,
                    file_stem=str(pair_artifact["plot_file_stem"]),
                    component_type="protein" if pair_artifact["pair_type"] == "protein_pair" else "ligand",
                    threshold=(
                        protein_similarity_threshold
                        if pair_artifact["pair_type"] == "protein_pair"
                        else ligand_similarity_threshold
                    ),
                )
        _materialize_reference_landscape_summary(
            chain_df=output_chain_df,
            proteins_df=proteins_df,
            ligands_df=ligands_df,
            protein_queries=protein_queries,
            ligand_queries=ligand_queries,
            output_dir=output_root,
        )
        logger.info("Bias training views written to %s", output_dir)

    logger.info(
        "Lightweight bias metrics computed: protein_max=%s ligand_max=%s (protein_top_n=%d)",
        (float(prot_vals.max()) if len(prot_vals) else None),
        (float(lig_vals.max()) if len(lig_vals) else None),
        protein_top_n,
    )

    return system_df, chain_df
