"""Bias metrics against protein/ligand training references."""

from __future__ import annotations

from datetime import date
from pathlib import Path
import logging

import pandas as pd
from Bio import pairwise2
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem


def _parse_iso_date(value: str) -> date | None:
    try:
        return date.fromisoformat(str(value)[:10])
    except Exception:
        return None


def _is_before_cutoff(value: str, cutoff: date) -> bool:
    d = _parse_iso_date(value)
    return bool(d and d < cutoff)


def _load_protein_training(path: Path, cutoff: date) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {"pdb_id", "release_date", "sequence"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Protein training data missing required columns: {sorted(missing)}")
    out = df[df["release_date"].apply(lambda x: _is_before_cutoff(x, cutoff))].copy()
    out["sequence"] = out["sequence"].astype(str)
    return out


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
        if smiles and pdb_id and release_date:
            rows.append(
                {
                    "pdb_id": str(pdb_id),
                    "release_date": str(release_date),
                    "smiles": str(smiles),
                    "ligand_id": str(ligand_id) if ligand_id is not None else None,
                }
            )
    return pd.DataFrame(rows)


def _load_ligand_training(path: Path, cutoff: date) -> pd.DataFrame:
    if path.suffix.lower() == ".sdf":
        df = _read_ligand_training_from_sdf(path)
    else:
        df = pd.read_csv(path)
    required = {"pdb_id", "release_date", "smiles"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Ligand training data missing required columns: {sorted(missing)}")
    out = df[df["release_date"].apply(lambda x: _is_before_cutoff(x, cutoff))].copy()
    out["smiles"] = out["smiles"].astype(str)
    if "ligand_id" not in out.columns:
        out["ligand_id"] = None
    return out


def _sequence_identity_percent(query: str, target: str) -> float:
    if not query or not target:
        return 0.0
    aln = pairwise2.align.globalxx(query, target, score_only=True)
    denom = max(len(query), len(target))
    if denom == 0:
        return 0.0
    return float(100.0 * float(aln) / float(denom))


def _best_protein_hit(query_seq: str, proteins_df: pd.DataFrame):
    best_score = None
    best_pdb = None
    for _, row in proteins_df.iterrows():
        score = _sequence_identity_percent(query_seq, str(row["sequence"]))
        if best_score is None or score > best_score:
            best_score = score
            best_pdb = str(row["pdb_id"])
    return best_score, best_pdb


def _morgan_fp_from_smiles(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)


def _best_ligand_hit(query_smiles: str, ligands_df: pd.DataFrame):
    query_fp = _morgan_fp_from_smiles(query_smiles)
    if query_fp is None:
        return None, None, None
    best_sim = None
    best_pdb = None
    best_ligand_id = None
    for _, row in ligands_df.iterrows():
        fp = _morgan_fp_from_smiles(str(row["smiles"]))
        if fp is None:
            continue
        sim = float(DataStructs.TanimotoSimilarity(query_fp, fp))
        if best_sim is None or sim > best_sim:
            best_sim = sim
            best_pdb = str(row["pdb_id"])
            best_ligand_id = (
                str(row["ligand_id"]) if pd.notna(row.get("ligand_id")) else None
            )
    return best_sim, best_pdb, best_ligand_id


def _extract_system_queries(sys_obj, ligands_df: pd.DataFrame):
    sequences = sys_obj.find_value(key="sequences") or []
    protein_by_chain = {}
    ligand_by_chain = {}
    ligand_fallback_by_id = {}

    for _, row in ligands_df.dropna(subset=["ligand_id", "smiles"]).iterrows():
        ligand_fallback_by_id[str(row["ligand_id"])] = str(row["smiles"])

    for seq in sequences:
        if not isinstance(seq, dict):
            continue
        if "protein" in seq:
            p = seq["protein"] or {}
            ids = p.get("id")
            if ids is None:
                continue
            if not isinstance(ids, list):
                ids = [ids]
            sequence = p.get("sequence") or p.get("fasta")
            for cid in ids:
                protein_by_chain[str(cid)] = str(sequence) if sequence else None
        if "ligand" in seq:
            l = seq["ligand"] or {}
            ids = l.get("id")
            if ids is None:
                continue
            if not isinstance(ids, list):
                ids = [ids]
            smiles = l.get("smiles")
            ccd = l.get("ccd")
            for cid in ids:
                if smiles:
                    ligand_by_chain[str(cid)] = str(smiles)
                elif ccd:
                    ligand_by_chain[str(cid)] = ligand_fallback_by_id.get(str(ccd))
                else:
                    ligand_by_chain[str(cid)] = None

    return protein_by_chain, ligand_by_chain


def apply_bias_metrics(
    system_df: pd.DataFrame,
    chain_df: pd.DataFrame,
    sys_obj,
    protein_training_data_path: Path | str,
    ligand_training_data_path: Path | str,
    release_cutoff: str = "2023-06-01",
    logger: logging.Logger | None = None,
):
    """Compute chain-level bias metrics and one system-level bias point row."""
    if logger is None:
        logger = logging.getLogger(__name__)

    cutoff = date.fromisoformat(str(release_cutoff))
    proteins_df = _load_protein_training(Path(protein_training_data_path), cutoff)
    ligands_df = _load_ligand_training(Path(ligand_training_data_path), cutoff)

    protein_queries, ligand_queries = _extract_system_queries(sys_obj, ligands_df)

    protein_best = {}
    for chain_id, query_seq in protein_queries.items():
        if not query_seq:
            protein_best[chain_id] = (None, None)
            continue
        protein_best[chain_id] = _best_protein_hit(query_seq, proteins_df)

    ligand_best = {}
    for chain_id, query_smiles in ligand_queries.items():
        if not query_smiles:
            ligand_best[chain_id] = (None, None, None)
            continue
        ligand_best[chain_id] = _best_ligand_hit(query_smiles, ligands_df)

    for idx, row in chain_df.iterrows():
        entity = str(row.get("ENTITY_TYPE"))
        chain_id = str(row.get("CHAIN_ID"))
        if entity == "protein":
            chain_df.at[idx, "protein_sequence_identity_train"] = protein_best.get(chain_id, (None, None))[0]
        elif entity == "ligand":
            chain_df.at[idx, "ligand_fingerprint_similarity_train"] = ligand_best.get(chain_id, (None, None, None))[0]

    protein_vals = pd.to_numeric(chain_df["protein_sequence_identity_train"], errors="coerce").dropna()
    ligand_vals = pd.to_numeric(chain_df["ligand_fingerprint_similarity_train"], errors="coerce").dropna()
    protein_max = float(protein_vals.max()) if len(protein_vals) else None
    ligand_max = float(ligand_vals.max()) if len(ligand_vals) else None

    if protein_max is not None:
        system_df["protein_sequence_identity_train_max"] = protein_max
    if ligand_max is not None:
        system_df["ligand_fingerprint_similarity_train_max"] = ligand_max

    # one system-level point only, independent from repeats/diffusion
    protein_top_pdb = None
    ligand_top_pdb = None
    ligand_top_id = None

    for _, (score, pdb_id) in protein_best.items():
        if score is not None and protein_max is not None and abs(score - protein_max) < 1e-12:
            protein_top_pdb = pdb_id
            break
    for _, (score, pdb_id, ligand_id) in ligand_best.items():
        if score is not None and ligand_max is not None and abs(score - ligand_max) < 1e-12:
            ligand_top_pdb = pdb_id
            ligand_top_id = ligand_id
            break

    system_name = (
        str(system_df["model_name"].iloc[0]) if "model_name" in system_df.columns and len(system_df) else "system"
    )
    bias_df = pd.DataFrame(
        [
            {
                "system_name": system_name,
                "x_ligand_ecfp_similarity": ligand_max,
                "y_protein_sequence_identity": protein_max,
                "x_ligand_top_pdb": ligand_top_pdb,
                "y_protein_top_pdb": protein_top_pdb,
                "x_ligand_top_ligand_id": ligand_top_id,
                "ligand_bias_threshold": 0.35,
                "protein_bias_threshold": 25.0,
                "bias_flag": (
                    "high_bias"
                    if ligand_max is not None and protein_max is not None and ligand_max > 0.35 and protein_max > 25.0
                    else "low_bias"
                ),
            }
        ]
    )

    logger.info(
        "Bias metrics computed: protein_max=%s ligand_max=%s", protein_max, ligand_max
    )
    return system_df, chain_df, bias_df
