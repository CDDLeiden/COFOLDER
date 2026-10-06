"""Reference loading and query resolution for bias analytics."""

from __future__ import annotations

from rdkit import Chem
from pandas.errors import EmptyDataError
from pathlib import Path
from datetime import date
import gemmi
from functools import lru_cache
import pandas as pd
import pickle

from ._common import (
    CUSTOM_LIGAND_REQUIRED_COLUMNS,
    CUSTOM_PROTEIN_REQUIRED_COLUMNS,
    LIGAND_BASE_COLUMNS,
    PROTEIN_BASE_COLUMNS,
    PROVENANCE_COLUMNS,
    _canonical_reference_source,
    _empty_ligand_reference_frame,
    _empty_protein_reference_frame,
    _is_before_cutoff,
    _norm_id,
    _normalize_reference_path_fields,
    _normalize_smiles,
    _path_cache_token,
)


def _finalize_reference_frame(
    df: pd.DataFrame,
    *,
    input_path: Path,
    default_source: str,
    default_dataset_name: str,
    is_ligand: bool,
) -> pd.DataFrame:
    if df.empty and len(df.columns) == 0:
        df = (
            _empty_ligand_reference_frame()
            if is_ligand
            else _empty_protein_reference_frame()
        )

    base_columns = LIGAND_BASE_COLUMNS if is_ligand else PROTEIN_BASE_COLUMNS
    for column in base_columns:
        if column not in df.columns:
            df[column] = pd.NA

    if "source" not in df.columns:
        df["source"] = default_source
    else:
        df["source"] = df["source"].apply(
            lambda value: _canonical_reference_source(value, default_source)
        )

    if "dataset_name" not in df.columns:
        df["dataset_name"] = default_dataset_name
    else:
        df["dataset_name"] = df["dataset_name"].apply(
            lambda value: (
                str(value).strip()
                if pd.notna(value) and str(value).strip()
                else default_dataset_name
            )
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
            df["ecfp_similarity"] = pd.to_numeric(
                df["ecfp_similarity"], errors="coerce"
            )
    else:
        if "sequence" in df.columns:
            df["sequence"] = df["sequence"].astype(str)
        df["sequence_similarity"] = pd.to_numeric(
            df["sequence_similarity"], errors="coerce"
        )
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
        invalid_pairwise = df["sequence_similarity_pairwise"].notna() & method.ne(
            "pairwise_aligner"
        )
        invalid_unavailable = (
            df["sequence_similarity"].isna()
            & df["sequence_similarity_pairwise"].isna()
            & method.ne("unavailable")
        )
        if invalid_mmseqs.any() or invalid_pairwise.any() or invalid_unavailable.any():
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

    return {ccd_id: descriptor for ccd_id, (_, descriptor) in best_by_ccd.items()}


def _smiles_from_components_cif(
    ccd_id: str, components_cif_path: Path | None
) -> str | None:
    if components_cif_path is None or not components_cif_path.exists():
        return None

    target_id = _norm_id(ccd_id)
    return _components_cif_smiles_index(_path_cache_token(components_cif_path)).get(
        target_id
    )


def _load_public_protein_training(
    path: Path, cutoff: date | None, top_n: int = 100
) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {"pdb_id", "release_date", "sequence"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Protein training data missing required columns: {sorted(missing)}"
        )

    release_mask = (
        pd.Series(True, index=df.index)
        if cutoff is None
        else df["release_date"]
        .map(lambda x: _is_before_cutoff(x, cutoff))
        .fillna(False)
        .astype(bool)
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
            out = (
                out.sort_values("sequence_similarity", ascending=False)
                .head(top_n)
                .copy()
            )
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
        release_date = (
            mol.GetProp("release_date") if mol.HasProp("release_date") else None
        )
        ligand_id = mol.GetProp("ligand_id") if mol.HasProp("ligand_id") else None
        source = mol.GetProp("source") if mol.HasProp("source") else None
        dataset_name = (
            mol.GetProp("dataset_name") if mol.HasProp("dataset_name") else None
        )
        source_structure_path = (
            mol.GetProp("source_structure_path")
            if mol.HasProp("source_structure_path")
            else None
        )
        source_reference_path = (
            mol.GetProp("source_reference_path")
            if mol.HasProp("source_reference_path")
            else None
        )
        ecfp_similarity = (
            mol.GetProp("ecfp_similarity") if mol.HasProp("ecfp_similarity") else None
        )
        if smiles:
            rows.append(
                {
                    "pdb_id": str(pdb_id) if pdb_id is not None else None,
                    "release_date": (
                        str(release_date) if release_date is not None else None
                    ),
                    "smiles": str(smiles),
                    "ligand_id": str(ligand_id) if ligand_id is not None else None,
                    "source": str(source) if source is not None else None,
                    "dataset_name": (
                        str(dataset_name) if dataset_name is not None else None
                    ),
                    "source_structure_path": (
                        str(source_structure_path)
                        if source_structure_path is not None
                        else None
                    ),
                    "source_reference_path": (
                        str(source_reference_path)
                        if source_reference_path is not None
                        else None
                    ),
                    "ecfp_similarity": (
                        str(ecfp_similarity) if ecfp_similarity is not None else None
                    ),
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
        raise ValueError(
            f"Ligand training data missing required columns: {sorted(missing)}"
        )

    release_mask = (
        pd.Series(True, index=df.index)
        if cutoff is None
        else df["release_date"]
        .map(lambda x: _is_before_cutoff(x, cutoff))
        .fillna(False)
        .astype(bool)
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


def _concat_reference_frames(
    frames: list[pd.DataFrame], *, is_ligand: bool
) -> pd.DataFrame:
    non_empty_frames = [frame for frame in frames if frame is not None]
    if not non_empty_frames:
        return (
            _empty_ligand_reference_frame()
            if is_ligand
            else _empty_protein_reference_frame()
        )
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
                    protein_by_chain[str(cid).strip()] = (
                        str(sequence) if sequence else None
                    )

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
                        mapped_smiles = _smiles_from_components_cif(
                            ccd_id, components_cif_path
                        )
                    if not mapped_smiles:
                        mapped_smiles = _smiles_from_ccd_cache(ccd_id, boltz_cache_path)
                    if mapped_smiles:
                        break

                ligand_by_chain[chain_id] = mapped_smiles
                if mapped_smiles is None:
                    unresolved_ccd_by_chain[chain_id] = tuple(ccd_ids)

    return protein_by_chain, ligand_by_chain, unresolved_ccd_by_chain


__all__ = [
    "_components_cif_smiles_index",
    "_concat_reference_frames",
    "_extract_system_queries",
    "_finalize_reference_frame",
    "_load_custom_ligand_training",
    "_load_custom_protein_training",
    "_load_ligand_training",
    "_load_protein_training",
    "_load_public_ligand_training",
    "_load_public_protein_training",
    "_read_ligand_training_from_sdf",
    "_smiles_from_ccd_cache",
    "_smiles_from_components_cif",
]
