"""Lightweight bias metrics against protein/ligand training references."""

from __future__ import annotations

from contextlib import nullcontext
from datetime import date
import logging
from pathlib import Path
import pickle

from Bio.Align import PairwiseAligner
import gemmi
import pandas as pd
from pandas.errors import EmptyDataError
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

from cofolder.modules.analytics.plots import plot_bias_reference_overlap
from cofolder.modules.utils.timing import DebugTimingCollector

ALIGNER = PairwiseAligner(mode="global")
ALIGNER.match_score = 1.0
ALIGNER.mismatch_score = 0.0
ALIGNER.open_gap_score = 0.0
ALIGNER.extend_gap_score = 0.0
MORGAN_FP = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

PROVENANCE_COLUMNS = [
    "source",
    "dataset_name",
    "source_structure_path",
    "source_reference_path",
]
PROTEIN_BASE_COLUMNS = ["pdb_id", "release_date", "sequence", "sequence_similarity"]
LIGAND_BASE_COLUMNS = ["pdb_id", "release_date", "ligand_id", "smiles", "ecfp_similarity"]
CUSTOM_PROTEIN_REQUIRED_COLUMNS = {"sequence"}
CUSTOM_LIGAND_REQUIRED_COLUMNS = {"smiles"}
REFERENCE_PATH_COLUMNS = ("source_structure_path", "source_reference_path")
BIAS_PROTEIN_VIEW_MIN_SIMILARITY = 25.0
BIAS_LIGAND_VIEW_MIN_SIMILARITY = 0.35
BIAS_PLOT_SEQUENCE_THRESHOLD = 0.25
BIAS_PLOT_LIGAND_THRESHOLD = 0.35
BIAS_PLOT_FILE_STEM = "bias_reference_overlap_scatter"
BIAS_PLOT_SKIPPED_FILE = f"{BIAS_PLOT_FILE_STEM}.skipped.txt"
BIAS_PAIRING_PAIRED = "paired"
BIAS_PAIRING_SYNTHETIC = "synthetic_paired"
BIAS_PAIRING_PROTEIN_ONLY = "protein_only"
BIAS_PAIRING_LIGAND_ONLY = "ligand_only"
BIAS_TRAINING_DATA_COLUMNS = [
    "query_pair_id",
    "query_protein_chain_id",
    "query_ligand_chain_id",
    "reference_key",
    "reference_label",
    "pairing_status",
    "source",
    "dataset_name",
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
    "ecfp_similarity",
    "plot_sequence_similarity",
    "plot_ecfp_similarity",
    "sequence",
    "ligand_id",
    "smiles",
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
            df["smiles"] = df["smiles"].astype(str)
        if "ecfp_similarity" in df.columns:
            df["ecfp_similarity"] = pd.to_numeric(df["ecfp_similarity"], errors="coerce")
    else:
        if "sequence" in df.columns:
            df["sequence"] = df["sequence"].astype(str)
        if "sequence_similarity" in df.columns:
            df["sequence_similarity"] = pd.to_numeric(df["sequence_similarity"], errors="coerce")

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


def _smiles_from_components_cif(ccd_id: str, components_cif_path: Path | None) -> str | None:
    if components_cif_path is None or not components_cif_path.exists():
        return None

    target_id = _norm_id(ccd_id)
    try:
        doc = gemmi.cif.read_file(str(components_cif_path))
    except Exception:
        return None

    descriptor_priority = {
        "SMILES_CANONICAL": 0,
        "SMILES": 1,
    }
    best_rank = None
    best_smiles = None

    for block in doc:
        block_id = _norm_id(block.find_value("_chem_comp.id") or block.name)
        if block_id != target_id:
            continue

        table = block.find(
            "_pdbx_chem_comp_descriptor.",
            ["comp_id", "type", "program", "descriptor"],
        )
        for row in table:
            if len(row) != 4 or _norm_id(row[0]) != target_id:
                continue
            descriptor_type = str(row[1]).strip().upper()
            descriptor = str(row[3]).strip()
            if not descriptor:
                continue
            rank = descriptor_priority.get(descriptor_type, 99)
            if best_rank is None or rank < best_rank:
                best_rank = rank
                best_smiles = descriptor

    return best_smiles


def _load_public_protein_training(path: Path, cutoff: date, top_n: int = 100) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {"pdb_id", "release_date", "sequence"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Protein training data missing required columns: {sorted(missing)}")

    release_mask = df["release_date"].map(lambda x: _is_before_cutoff(x, cutoff)).fillna(False).astype(bool)
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


def _load_public_ligand_training(path: Path, cutoff: date) -> pd.DataFrame:
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

    release_mask = df["release_date"].map(lambda x: _is_before_cutoff(x, cutoff)).fillna(False).astype(bool)
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
    cutoff: date,
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
    cutoff: date,
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


def _sequence_identity_percent(query: str, target: str) -> float:
    if not query or not target:
        return 0.0
    score = ALIGNER.score(query, target)
    denom = max(len(query), len(target))
    if denom == 0:
        return 0.0
    return float(100.0 * float(score) / float(denom))


def _best_protein_hit(query_seq: str, proteins_df: pd.DataFrame) -> float | None:
    if "sequence_similarity" in proteins_df.columns:
        if "sequence" in proteins_df.columns:
            exact_match = proteins_df["sequence"].astype(str) == str(query_seq)
            if exact_match.any():
                return 100.0
        sims = pd.to_numeric(proteins_df["sequence_similarity"], errors="coerce").dropna()
        if len(sims):
            return float(sims.max())

    best_score = None
    for _, row in proteins_df.iterrows():
        score = _sequence_identity_percent(query_seq, str(row["sequence"]))
        if best_score is None or score > best_score:
            best_score = score
    return best_score


def _morgan_fp_from_smiles(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
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


def _protein_similarity_series(query_seq: str, proteins_df: pd.DataFrame) -> pd.Series:
    if proteins_df.empty:
        return pd.Series(dtype=float)
    if "sequence_similarity" in proteins_df.columns:
        sims = pd.to_numeric(proteins_df["sequence_similarity"], errors="coerce")
    else:
        sims = proteins_df["sequence"].apply(lambda seq: _sequence_identity_percent(query_seq, str(seq)))
    if "sequence" in proteins_df.columns:
        exact_match = proteins_df["sequence"].astype(str) == str(query_seq)
        sims = sims.mask(exact_match, 100.0)
    return pd.to_numeric(sims, errors="coerce")


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
    similarities = (
        _ligand_similarity_series(query_value, df, mol_id=mol_id)
        if is_ligand
        else _protein_similarity_series(query_value, df)
    )
    if similarities.dropna().empty:
        return None
    best_idx = similarities.fillna(float("-inf")).idxmax()
    best_row = df.loc[best_idx].copy()
    best_row["best_similarity"] = float(similarities.loc[best_idx])
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
        "nearest_custom_pdb_id": custom_best.get("pdb_id") if custom_best is not None else None,
        "nearest_custom_dataset_name": custom_best.get("dataset_name") if custom_best is not None else None,
        "nearest_custom_similarity": custom_similarity,
        "nearest_overall_source": nearest_overall_source,
        "nearest_overall_pdb_id": overall_best.get("pdb_id") if overall_best is not None else None,
        "nearest_overall_dataset_name": overall_best.get("dataset_name") if overall_best is not None else None,
        "nearest_overall_similarity": (
            float(overall_best["best_similarity"])
            if overall_best is not None and pd.notna(overall_best.get("best_similarity"))
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
            l = seq["ligand"] or {}
            ids = l.get("id")
            if ids is None:
                continue
            if not isinstance(ids, list):
                ids = [ids]

            smiles = l.get("smiles")
            ccd = l.get("ccd")
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
) -> pd.DataFrame:
    prot_rows: list[pd.DataFrame] = []
    for chain_id, query_seq in protein_queries.items():
        if not query_seq:
            continue
        sub = proteins_df.copy()
        sub["sequence_similarity"] = _protein_similarity_series(str(query_seq), sub)
        sub["sequence_similarity"] = pd.to_numeric(sub["sequence_similarity"], errors="coerce").clip(
            lower=0.0,
            upper=100.0,
        )
        sub = sub[sub["sequence_similarity"] >= BIAS_PROTEIN_VIEW_MIN_SIMILARITY].copy()
        sub = sub.sort_values("sequence_similarity", ascending=False).head(100).copy()
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
            by=["query_chain_id", "sequence_similarity", "source", "pdb_id"],
            ascending=[True, False, True, True],
            na_position="last",
        )
    else:
        proteins_out = proteins_df.iloc[0:0].copy()
        proteins_out["sequence_similarity"] = pd.Series(dtype=float)
        proteins_out["query_chain_id"] = pd.Series(dtype=str)

    return _order_protein_training_columns(proteins_out)


def _build_ligand_training_views(
    chain_df: pd.DataFrame,
    ligands_df: pd.DataFrame,
    ligand_queries: dict[str, str | None],
    boltz_cache_path: Path,
    components_cif_path: Path | None = None,
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

        sub = ligands_df.copy()
        sub["ecfp_similarity"] = _ligand_similarity_series(query_smiles, sub, mol_id=mol_id)
        sub = sub.dropna(subset=["ecfp_similarity"]).copy()
        sub = sub.sort_values("ecfp_similarity", ascending=False, na_position="last")
        above = sub[sub["ecfp_similarity"] >= BIAS_LIGAND_VIEW_MIN_SIMILARITY].copy()
        if not above.empty:
            sub = above
        else:
            sub = sub.head(100).copy()
        sub["query_chain_id"] = chain_id
        views[chain_id] = _order_ligand_training_columns(sub)

    return views


def _reference_key_from_row(row: pd.Series, *, is_ligand: bool) -> str:
    source = _canonical_reference_source(row.get("source"), "custom")
    pdb_id = row.get("pdb_id")
    if pd.notna(pdb_id) and str(pdb_id).strip() and source == "public":
        return f"{source}:pdb:{_norm_id(pdb_id)}"

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


def _synthetic_reference_key(protein_reference_key: object, ligand_reference_key: object) -> str:
    protein_component = str(protein_reference_key).strip()
    ligand_component = str(ligand_reference_key).strip()
    return (
        "synthetic:"
        f"{len(protein_component)}:{protein_component}"
        "|"
        f"{len(ligand_component)}:{ligand_component}"
    )


def _synthetic_reference_label(protein_row: pd.Series, ligand_row: pd.Series) -> str:
    protein_label = _reference_label_from_merged_row(
        pd.Series(
            {
                "pdb_id_protein": protein_row.get("pdb_id_protein"),
                "source_structure_path_protein": protein_row.get("source_structure_path_protein"),
                "source_reference_path_protein": protein_row.get("source_reference_path_protein"),
                "dataset_name_protein": protein_row.get("dataset_name_protein"),
            }
        )
    )
    ligand_label = _reference_label_from_merged_row(
        pd.Series(
            {
                "pdb_id_ligand": ligand_row.get("pdb_id_ligand"),
                "source_structure_path_ligand": ligand_row.get("source_structure_path_ligand"),
                "source_reference_path_ligand": ligand_row.get("source_reference_path_ligand"),
                "dataset_name_ligand": ligand_row.get("dataset_name_ligand"),
            }
        )
    )
    return f"{protein_label or 'protein'} x {ligand_label or 'ligand'}"


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
) -> dict[str, object]:
    protein_chain = query_protein_chain_id or None
    ligand_chain = query_ligand_chain_id or None
    protein_reference_key = protein_row.get("reference_key") if protein_row is not None else None
    ligand_reference_key = ligand_row.get("reference_key") if ligand_row is not None else None
    protein_source = protein_row.get("source_protein") if protein_row is not None else None
    ligand_source = ligand_row.get("source_ligand") if ligand_row is not None else None

    sequence_similarity = pd.to_numeric(
        pd.Series([protein_row.get("sequence_similarity_protein") if protein_row is not None else pd.NA]),
        errors="coerce",
    ).iloc[0]
    ecfp_similarity = pd.to_numeric(
        pd.Series([ligand_row.get("ecfp_similarity_ligand") if ligand_row is not None else pd.NA]),
        errors="coerce",
    ).iloc[0]
    plot_sequence_similarity = (
        float(sequence_similarity) / 100.0 if pd.notna(sequence_similarity) else pd.NA
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
        "pdb_id": "synthetic_pair"
        if pairing_status == BIAS_PAIRING_SYNTHETIC
        else _collapse_pair_value(
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
        "ecfp_similarity": float(ecfp_similarity) if pd.notna(ecfp_similarity) else pd.NA,
        "plot_sequence_similarity": plot_sequence_similarity,
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
) -> list[dict[str, object]]:
    left = protein_view.copy()
    if "reference_key" not in left.columns:
        left["reference_key"] = left.apply(lambda row: _reference_key_from_row(row, is_ligand=False), axis=1)
    left = left.rename(columns=lambda column: column if column == "reference_key" else f"{column}_protein")

    right = ligand_view.copy()
    if "reference_key" not in right.columns:
        right["reference_key"] = right.apply(lambda row: _reference_key_from_row(row, is_ligand=True), axis=1)
    right = right.rename(columns=lambda column: column if column == "reference_key" else f"{column}_ligand")

    merged = left.merge(
        right,
        how="outer",
        on="reference_key",
        sort=True,
    )
    if merged.empty:
        return []

    # Public references are PDB-level aggregates in the legacy bias builder.
    # When multiple protein or ligand rows exist for the same PDB, we keep the
    # paired numeric bias rows but avoid presenting sequence/ligand identifiers
    # as if each cross-product row were a validated one-to-one reference pair.
    protein_unique_counts = (
        left.groupby("reference_key")["sequence_protein"].nunique(dropna=True).to_dict()
        if "sequence_protein" in left.columns
        else {}
    )
    ligand_id_unique_counts = (
        right.groupby("reference_key")["ligand_id_ligand"].nunique(dropna=True).to_dict()
        if "ligand_id_ligand" in right.columns
        else {}
    )
    ligand_smiles_unique_counts = (
        right.groupby("reference_key")["smiles_ligand"].nunique(dropna=True).to_dict()
        if "smiles_ligand" in right.columns
        else {}
    )

    rows: list[dict[str, object]] = []
    for _, row in merged.iterrows():
        has_protein = pd.notna(row.get("query_chain_id_protein")) and str(row.get("query_chain_id_protein")).strip()
        has_ligand = pd.notna(row.get("query_chain_id_ligand")) and str(row.get("query_chain_id_ligand")).strip()
        if has_protein and has_ligand:
            pairing_status = BIAS_PAIRING_PAIRED
        elif has_protein:
            pairing_status = BIAS_PAIRING_PROTEIN_ONLY
        else:
            pairing_status = BIAS_PAIRING_LIGAND_ONLY

        rows.append(
            _bias_training_row_from_sources(
                query_protein_chain_id=query_protein_chain_id,
                query_ligand_chain_id=query_ligand_chain_id,
                reference_key=row.get("reference_key"),
                reference_label=_reference_label_from_merged_row(row),
                pairing_status=pairing_status,
                protein_row=row if has_protein else None,
                ligand_row=row if has_ligand else None,
                protein_unique_counts=protein_unique_counts,
                ligand_id_unique_counts=ligand_id_unique_counts,
                ligand_smiles_unique_counts=ligand_smiles_unique_counts,
            )
        )

    left_records = left.to_dict("records")
    right_records = right.to_dict("records")
    for protein_row_dict in left_records:
        protein_row = pd.Series(protein_row_dict)
        for ligand_row_dict in right_records:
            ligand_row = pd.Series(ligand_row_dict)
            if protein_row.get("reference_key") == ligand_row.get("reference_key"):
                continue
            rows.append(
                _bias_training_row_from_sources(
                    query_protein_chain_id=query_protein_chain_id,
                    query_ligand_chain_id=query_ligand_chain_id,
                    reference_key=_synthetic_reference_key(
                        protein_row.get("reference_key"),
                        ligand_row.get("reference_key"),
                    ),
                    reference_label=_synthetic_reference_label(protein_row, ligand_row),
                    pairing_status=BIAS_PAIRING_SYNTHETIC,
                    protein_row=protein_row,
                    ligand_row=ligand_row,
                    protein_unique_counts=protein_unique_counts,
                    ligand_id_unique_counts=ligand_id_unique_counts,
                    ligand_smiles_unique_counts=ligand_smiles_unique_counts,
                )
            )

    return rows


def _build_bias_training_dataset(
    *,
    protein_view: pd.DataFrame,
    ligand_views: dict[str, pd.DataFrame],
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
        for ligand_chain_id in ligand_iter:
            ligand_subset = (
                ligand_views.get(ligand_chain_id, pd.DataFrame()).copy()
                if ligand_chain_id is not None
                else pd.DataFrame()
            )
            rows.extend(
                _build_bias_training_rows(
                    query_protein_chain_id=protein_chain_id,
                    query_ligand_chain_id=ligand_chain_id,
                    protein_view=protein_subset,
                    ligand_view=ligand_subset,
                )
            )

    bias_training_df = pd.DataFrame(rows)
    if bias_training_df.empty:
        return pd.DataFrame(columns=BIAS_TRAINING_DATA_COLUMNS)

    bias_training_df["plot_sequence_similarity"] = pd.to_numeric(
        bias_training_df["plot_sequence_similarity"],
        errors="coerce",
    ).clip(lower=0.0, upper=1.0)
    bias_training_df["plot_ecfp_similarity"] = pd.to_numeric(
        bias_training_df["plot_ecfp_similarity"],
        errors="coerce",
    ).clip(lower=0.0, upper=1.0)
    bias_training_df["_pairing_order"] = bias_training_df["pairing_status"].map(
        {
            BIAS_PAIRING_PAIRED: 0,
            BIAS_PAIRING_SYNTHETIC: 1,
            BIAS_PAIRING_PROTEIN_ONLY: 2,
            BIAS_PAIRING_LIGAND_ONLY: 3,
        }
    )
    bias_training_df = bias_training_df.sort_values(
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

    for column in BIAS_TRAINING_DATA_COLUMNS:
        if column not in bias_training_df.columns:
            bias_training_df[column] = pd.NA
    return bias_training_df.loc[:, BIAS_TRAINING_DATA_COLUMNS].reset_index(drop=True)


def _materialize_bias_training_views(
    chain_df: pd.DataFrame,
    ligands_df: pd.DataFrame,
    proteins_df: pd.DataFrame,
    protein_queries: dict[str, str | None],
    ligand_queries: dict[str, str | None],
    boltz_cache_path: Path,
    output_dir: Path,
    components_cif_path: Path | None = None,
) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    """Write bias training files for debug/inspection under results/bias_train."""
    output_dir.mkdir(parents=True, exist_ok=True)
    for stale_path in output_dir.glob("ligand_training_data*.csv"):
        stale_path.unlink()

    proteins_out = _build_protein_training_view(proteins_df, protein_queries)
    proteins_out.to_csv(output_dir / "protein_training_data.csv", index=False)

    ligand_views = _build_ligand_training_views(
        chain_df=chain_df,
        ligands_df=ligands_df,
        ligand_queries=ligand_queries,
        boltz_cache_path=boltz_cache_path,
        components_cif_path=components_cif_path,
    )

    lig_rows = chain_df[chain_df["ENTITY_TYPE"].astype(str) == "ligand"].copy()
    if not lig_rows.empty and "ligand_molecule_id" not in lig_rows.columns:
        _order_ligand_training_columns(ligands_df).to_csv(
            output_dir / "ligand_training_data.csv",
            index=False,
        )
        return proteins_out, ligand_views

    for chain_id, ligand_view in ligand_views.items():
        out_path = output_dir / f"ligand_training_data_{chain_id}.csv"
        ligand_view.drop(columns=["query_chain_id"], errors="ignore").to_csv(out_path, index=False)

    return proteins_out, ligand_views


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

    cutoff = date.fromisoformat(str(release_cutoff))
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
        "Bias cutoff applied (< %s): protein_rows=%d ligand_rows=%d",
        cutoff.isoformat(),
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
    protein_timer = (
        timings.measure("scores.bias_metrics.protein_similarity", logger=logger)
        if timings is not None
        else nullcontext()
    )
    with protein_timer:
        for chain_id, query_seq in protein_queries.items():
            protein_best[str(chain_id).strip().upper()] = (
                _best_protein_hit(query_seq, proteins_df) if query_seq else None
            )

    ligand_best: dict[str, float | None] = {}
    ligand_timer = (
        timings.measure("scores.bias_metrics.ligand_similarity", logger=logger)
        if timings is not None
        else nullcontext()
    )
    with ligand_timer:
        for chain_id, query_smiles in ligand_queries.items():
            ligand_best[str(chain_id).strip().upper()] = (
                _best_ligand_hit(query_smiles, ligands_df) if query_smiles else None
            )

    if "bias_prot_sim_train" not in chain_df.columns:
        chain_df["bias_prot_sim_train"] = pd.NA
    if "bias_lig_sim_train" not in chain_df.columns:
        chain_df["bias_lig_sim_train"] = pd.NA

    for idx, row in chain_df.iterrows():
        entity = str(row.get("ENTITY_TYPE"))
        chain_id = str(row.get("CHAIN_ID")).strip().upper()
        if selected_chains and chain_id not in selected_chains:
            continue
        if entity == "protein":
            chain_df.at[idx, "bias_prot_sim_train"] = protein_best.get(chain_id)
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
    lig_vals = pd.to_numeric(chain_df.get("bias_lig_sim_train"), errors="coerce").dropna()

    if len(prot_vals):
        system_df["bias_prot_sim_train_max"] = float(prot_vals.max())
    if len(lig_vals):
        system_df["bias_lig_sim_train_max"] = float(lig_vals.max())

    if output_dir is not None:
        output_chain_df = chain_df
        if selected_chains:
            output_chain_df = chain_df[
                chain_df["CHAIN_ID"].astype(str).str.strip().str.upper().isin(selected_chains)
            ].copy()
        output_root = Path(output_dir)
        protein_view, ligand_views = _materialize_bias_training_views(
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
        bias_training_df = _build_bias_training_dataset(
            protein_view=protein_view,
            ligand_views=ligand_views,
        )
        bias_training_df.to_csv(output_root / "bias_training_data.csv", index=False)
        skipped_plot_path = output_root / BIAS_PLOT_SKIPPED_FILE
        plot_paths = plot_bias_reference_overlap(
            bias_training_df,
            output_dir=output_root,
            file_stem=BIAS_PLOT_FILE_STEM,
            sequence_threshold=BIAS_PLOT_SEQUENCE_THRESHOLD,
            ligand_threshold=BIAS_PLOT_LIGAND_THRESHOLD,
        )
        if plot_paths:
            if skipped_plot_path.exists():
                skipped_plot_path.unlink()
        else:
            for suffix in ("png", "pdf"):
                (output_root / f"{BIAS_PLOT_FILE_STEM}.{suffix}").unlink(missing_ok=True)
            skipped_plot_path.write_text(
                (
                    "Bias reference-overlap scatter plot skipped because the merged plotting dataset "
                    "did not contain rows with both protein and ligand similarity axes. "
                    "Bias artifacts remain reference-overlap diagnostics only and do not predict "
                    "affinity, structural confidence, or cofolding success.\n"
                ),
                encoding="utf-8",
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
