"""Lightweight bias metrics against protein/ligand training references."""

from __future__ import annotations

from datetime import date
from pathlib import Path
import logging
import pickle

import pandas as pd
from Bio.Align import PairwiseAligner
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

ALIGNER = PairwiseAligner(mode="global")
ALIGNER.match_score = 1.0
ALIGNER.mismatch_score = 0.0
ALIGNER.open_gap_score = 0.0
ALIGNER.extend_gap_score = 0.0
MORGAN_FP = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def _order_ligand_training_columns(df: pd.DataFrame) -> pd.DataFrame:
    preferred = [
        "query_chain_id",
        "pdb_id",
        "release_date",
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


def _load_protein_training(path: Path, cutoff: date, top_n: int = 100) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {"pdb_id", "release_date", "sequence"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Protein training data missing required columns: {sorted(missing)}")

    release_mask = df["release_date"].map(lambda x: _is_before_cutoff(x, cutoff)).fillna(False).astype(bool)
    out = df.loc[release_mask, :].copy()
    if "sequence" not in out.columns:
        # Defensive fallback for edge-case empty-frame indexing behavior.
        out = out.reindex(columns=df.columns)
    out["sequence"] = out["sequence"].astype(str)

    # Lightweight mode: only keep top-N protein candidates.
    if top_n > 0 and len(out) > top_n:
        if "sequence_similarity" in out.columns:
            out["sequence_similarity"] = pd.to_numeric(out["sequence_similarity"], errors="coerce")
            out = out.sort_values("sequence_similarity", ascending=False).head(top_n).copy()
        else:
            out = out.head(top_n).copy()

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

    release_mask = df["release_date"].map(lambda x: _is_before_cutoff(x, cutoff)).fillna(False).astype(bool)
    out = df.loc[release_mask, :].copy()
    out["smiles"] = out["smiles"].astype(str)
    if "ligand_id" not in out.columns:
        out["ligand_id"] = None
    else:
        out["ligand_id"] = out["ligand_id"].apply(
            lambda x: _norm_id(x) if pd.notna(x) and str(x).strip() else None
        )
    return out


def _sequence_identity_percent(query: str, target: str) -> float:
    if not query or not target:
        return 0.0
    score = ALIGNER.score(query, target)
    denom = max(len(query), len(target))
    if denom == 0:
        return 0.0
    return float(100.0 * float(score) / float(denom))


def _best_protein_hit(query_seq: str, proteins_df: pd.DataFrame) -> float | None:
    # If training data was generated by build_bias_training_data.py, reuse MMseqs
    # identity directly so validate() matches the standalone pipeline.
    if "sequence_similarity" in proteins_df.columns:
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


def _extract_system_queries(
    sys_obj,
    ligands_df: pd.DataFrame,
    boltz_cache_path: Path,
):
    sequences = sys_obj.find_value(key="sequences") or []
    protein_by_chain = {}
    ligand_by_chain = {}
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
                elif ccd_ids:
                    mapped_smiles = None
                    for ccd_id in ccd_ids:
                        mapped_smiles = ligand_fallback_by_id.get(ccd_id)
                        if not mapped_smiles:
                            mapped_smiles = _smiles_from_ccd_cache(ccd_id, boltz_cache_path)
                        if mapped_smiles:
                            break
                    ligand_by_chain[chain_id] = mapped_smiles
                else:
                    ligand_by_chain[chain_id] = None

    return protein_by_chain, ligand_by_chain


def _materialize_bias_training_views(
    chain_df: pd.DataFrame,
    ligands_df: pd.DataFrame,
    proteins_df: pd.DataFrame,
    protein_queries: dict[str, str | None],
    ligand_queries: dict[str, str | None],
    output_dir: Path,
) -> None:
    """Write bias training files for debug/inspection under results/bias_train."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save full protein training entities above threshold for any protein query.
    prot_rows: list[pd.DataFrame] = []
    for chain_id, query_seq in protein_queries.items():
        if not query_seq:
            continue
        sub = proteins_df.copy()
        if "sequence_similarity" in sub.columns:
            sub["sequence_similarity"] = pd.to_numeric(sub["sequence_similarity"], errors="coerce")
        else:
            sub["sequence_similarity"] = sub["sequence"].apply(
                lambda s: _sequence_identity_percent(str(query_seq), str(s))
            )
        # Keep all entities above threshold and cap at 100.
        sub["sequence_similarity"] = sub["sequence_similarity"].clip(lower=0.0, upper=100.0)
        sub = sub[sub["sequence_similarity"] >= 25.0].copy()
        sub = sub.sort_values("sequence_similarity", ascending=False).head(100).copy()
        sub["query_chain_id"] = str(chain_id).strip()
        prot_rows.append(sub)

    if prot_rows:
        proteins_out = pd.concat(prot_rows, ignore_index=True)
        proteins_out = proteins_out.drop_duplicates(
            subset=["pdb_id", "release_date", "sequence", "query_chain_id"]
        )
        proteins_out = proteins_out.sort_values(
            by=["sequence_similarity", "pdb_id"], ascending=[False, True]
        )
    else:
        proteins_out = proteins_df.iloc[0:0].copy()
        proteins_out["sequence_similarity"] = pd.Series(dtype=float)
        proteins_out["query_chain_id"] = pd.Series(dtype=str)

    proteins_out = _order_protein_training_columns(proteins_out)
    proteins_out.to_csv(output_dir / "protein_training_data.csv", index=False)

    # Save one ligand training file per unique ligand_molecule_id.
    lig_rows = chain_df[chain_df["ENTITY_TYPE"].astype(str) == "ligand"].copy()
    if lig_rows.empty:
        return

    if "ligand_molecule_id" not in lig_rows.columns:
        _order_ligand_training_columns(ligands_df).to_csv(
            output_dir / "ligand_training_data.csv",
            index=False,
        )
        return

    ordered = (
        lig_rows[["CHAIN_ID", "ligand_molecule_id"]]
        .dropna()
        .drop_duplicates()
        .drop_duplicates(subset=["ligand_molecule_id"], keep="first")
    )
    for _, row in ordered.iterrows():
        chain_id = str(row["CHAIN_ID"]).strip()
        mol_id = str(row["ligand_molecule_id"]).strip()
        if not chain_id or not mol_id:
            continue

        query_smiles = ligand_queries.get(chain_id)
        out_name = f"ligand_training_data_{chain_id}.csv"
        out_path = output_dir / out_name

        # Prefer chain-specific training rows produced by build_bias_training_data.py.
        if "query_chain_id" in ligands_df.columns:
            sub = ligands_df[
                ligands_df["query_chain_id"].astype(str).str.strip() == chain_id
            ].copy()
            if not sub.empty:
                sub["ecfp_similarity"] = pd.to_numeric(sub.get("ecfp_similarity"), errors="coerce")
                sub = sub.dropna(subset=["ecfp_similarity"]).copy()
                sub = sub[sub["ecfp_similarity"] >= 0.35].copy()
                sub = sub.sort_values("ecfp_similarity", ascending=False, na_position="last")
                _order_ligand_training_columns(sub).to_csv(out_path, index=False)
                continue

        if not query_smiles:
            # Fallback: derive chain-specific ligand query from molecule/CCD ID.
            ref = ligands_df[ligands_df.get("ligand_id").astype(str) == _norm_id(mol_id)]
            if not ref.empty:
                query_smiles = str(ref.iloc[0]["smiles"])

        if not query_smiles:
            # If still unavailable, keep full filtered ligand set for debugging.
            _order_ligand_training_columns(ligands_df).to_csv(out_path, index=False)
            continue

        qfp = _morgan_fp_from_smiles(query_smiles)
        sub = ligands_df.copy()
        sims = []
        for _, lrow in sub.iterrows():
            fp = _morgan_fp_from_smiles(str(lrow["smiles"]))
            if qfp is None or fp is None:
                sims.append(pd.NA)
            else:
                sims.append(float(DataStructs.TanimotoSimilarity(qfp, fp)))
        sub["ecfp_similarity"] = pd.to_numeric(pd.Series(sims), errors="coerce")
        sub = sub.dropna(subset=["ecfp_similarity"]).copy()
        sub = sub[sub["ecfp_similarity"] >= 0.35].copy()
        sub = sub.sort_values("ecfp_similarity", ascending=False, na_position="last")
        _order_ligand_training_columns(sub).to_csv(out_path, index=False)


def apply_bias_metrics(
    system_df: pd.DataFrame,
    chain_df: pd.DataFrame,
    sys_obj,
    protein_training_data_path: Path | str,
    ligand_training_data_path: Path | str,
    release_cutoff: str = "2023-06-01",
    bias_chains: set[str] | list[str] | None = None,
    protein_top_n: int = 100,
    boltz_cache_path: Path | str | None = None,
    output_dir: Path | str | None = None,
    logger: logging.Logger | None = None,
):
    """Compute lightweight chain-level bias metrics only.

    Outputs in `chain_df`:
    - bias_prot_sim_train (protein rows)
    - bias_lig_sim_train (ligand rows)
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    cutoff = date.fromisoformat(str(release_cutoff))
    proteins_df = _load_protein_training(Path(protein_training_data_path), cutoff, top_n=protein_top_n)
    ligands_df = _load_ligand_training(Path(ligand_training_data_path), cutoff)
    cache_path = Path(boltz_cache_path).expanduser() if boltz_cache_path else Path("~/.boltz").expanduser()
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

    protein_queries, ligand_queries = _extract_system_queries(
        sys_obj,
        ligands_df,
        boltz_cache_path=cache_path,
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

    protein_best: dict[str, float | None] = {}
    for chain_id, query_seq in protein_queries.items():
        protein_best[str(chain_id).strip().upper()] = (
            _best_protein_hit(query_seq, proteins_df) if query_seq else None
        )

    ligand_best: dict[str, float | None] = {}
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
            if sim is None:
                # Fallback: derive query from molecule ID when available.
                mol_id = row.get("ligand_molecule_id")
                if pd.notna(mol_id) and str(mol_id).strip():
                    ref = ligands_df[ligands_df.get("ligand_id").astype(str) == _norm_id(mol_id)]
                    if not ref.empty:
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
        _materialize_bias_training_views(
            chain_df=output_chain_df,
            ligands_df=ligands_df,
            proteins_df=_load_protein_training(Path(protein_training_data_path), cutoff, top_n=0),
            protein_queries=protein_queries,
            ligand_queries=ligand_queries,
            output_dir=Path(output_dir),
        )
        logger.info("Bias training views written to %s", output_dir)

    logger.info(
        "Lightweight bias metrics computed: protein_max=%s ligand_max=%s (protein_top_n=%d)",
        (float(prot_vals.max()) if len(prot_vals) else None),
        (float(lig_vals.max()) if len(lig_vals) else None),
        protein_top_n,
    )

    return system_df, chain_df
