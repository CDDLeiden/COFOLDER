"""Model reproduction and bias metrics for validation workflows."""

from __future__ import annotations

from pathlib import Path
import json
import logging
import numpy as np

from Bio.PDB import MMCIFParser, PDBParser
from rdkit import Chem

from cofolder.modules.analytics.align import (
    _align_structures_on_protein_ca,
    _rmsd,
)
from cofolder.modules.analytics.sucos import compute_sucos

SYSTEM_REPRODUCTION_COLUMNS = (
    "ligand_rmsd_ref_mean",
    "protein_rmsd_ref_mean",
    "ifp_similarity_ref_mean",
    "pocket_coverage_ref",
    "pocket_coverage_ref_mean",
    "ligand_pose_overlap_ref",
    "sucos_ref_mean",
    "protein_sequence_identity_train_max",
    "ligand_fingerprint_similarity_train_max",
)

CHAIN_REPRODUCTION_COLUMNS = (
    "ligand_rmsd_ref",
    "protein_rmsd_ref",
    "ifp_similarity_ref",
    "pocket_coverage_ref",
    "ligand_pose_overlap_ref",
    "sucos_ref",
    "sucos_shape_ref",
    "sucos_feature_ref",
    "protein_sequence_identity_train",
    "ligand_fingerprint_similarity_train",
)

DEFAULT_REPRODUCTION_METRICS = {
    "protein_rmsd",
    "ligand_rmsd",
    "sucos",
    "pocket_coverage",
}


# ---------------------------------------------------------------------------
# Structure / molecule loaders
# ---------------------------------------------------------------------------


def _load_structure(path: Path):
    suffix = path.suffix.lower()
    if suffix in {".cif", ".mmcif"}:
        parser = MMCIFParser(QUIET=True)
    elif suffix == ".pdb":
        parser = PDBParser(QUIET=True)
    else:
        raise ValueError(f"Unsupported structure format: {path}")
    return parser.get_structure(path.stem, str(path))


def _collect_matched_backbone_coords(ref_struct, mob_struct, chain_id):
    ref_coords = []
    mob_coords = []

    ref_model = next(iter(ref_struct))
    mob_model = next(iter(mob_struct))
    if chain_id not in ref_model or chain_id not in mob_model:
        return np.array(ref_coords), np.array(mob_coords)

    ref_chain = ref_model[chain_id]
    mob_chain = mob_model[chain_id]

    ref_res = {res.id: res for res in ref_chain if res.id[0] == " "}
    mob_res = {res.id: res for res in mob_chain if res.id[0] == " "}

    for rid in sorted(set(ref_res.keys()) & set(mob_res.keys())):
        r_ref = ref_res[rid]
        r_mob = mob_res[rid]
        for atom_name in ("N", "CA", "C", "O"):
            if atom_name in r_ref and atom_name in r_mob:
                ref_coords.append(r_ref[atom_name].coord)
                mob_coords.append(r_mob[atom_name].coord)

    return np.array(ref_coords, dtype=float), np.array(mob_coords, dtype=float)


def _collect_ligand_residue_heavy_atoms(structure):
    ligands = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == " ":
                    continue
                if residue.resname in {"HOH", "WAT", "DOD"}:
                    continue

                atoms = []
                elements = []
                for atom in residue:
                    element = (getattr(atom, "element", "") or "").upper()
                    if element == "H":
                        continue
                    atoms.append(atom)
                    elements.append(element)

                if not atoms:
                    continue

                ligands.append(
                    {
                        "resname": residue.resname,
                        "chain_id": chain.id,
                        "residue_id": residue.id,
                        "elements": tuple(elements),
                        "elements_sorted": tuple(sorted(elements)),
                        "coords": np.array([a.coord for a in atoms], dtype=float),
                    }
                )
    return ligands


def _collect_predicted_ligand_chain_heavy_atoms(structure, chain_id):
    atoms = []
    elements = []
    for model in structure:
        if chain_id not in model:
            continue
        for residue in model[chain_id]:
            if residue.id[0] == " ":
                continue
            if residue.resname in {"HOH", "WAT", "DOD"}:
                continue
            for atom in residue:
                element = (getattr(atom, "element", "") or "").upper()
                if element == "H":
                    continue
                atoms.append(atom)
                elements.append(element)

    if not atoms:
        return None, None

    return np.array([a.coord for a in atoms], dtype=float), tuple(sorted(elements))


# ---------------------------------------------------------------------------
# Ligand matching helpers
# ---------------------------------------------------------------------------


def _greedy_bipartite_rmsd(coords_a: np.ndarray, coords_b: np.ndarray) -> float:
    if len(coords_a) != len(coords_b):
        raise ValueError("Greedy bipartite RMSD requires equal atom counts.")

    d2 = ((coords_a[:, None, :] - coords_b[None, :, :]) ** 2).sum(axis=2)
    d2 = d2.copy()
    total = 0.0
    n = len(coords_a)

    for _ in range(n):
        i, j = np.unravel_index(np.argmin(d2), d2.shape)
        total += float(d2[i, j])
        d2[i, :] = np.inf
        d2[:, j] = np.inf

    return float(np.sqrt(total / n))


def _find_exact_ligand_match(pred_coords, pred_signature, reference_ligands):
    candidates = [
        lig
        for lig in reference_ligands
        if lig["elements_sorted"] == pred_signature and len(lig["coords"]) == len(pred_coords)
    ]
    if not candidates:
        return None, None

    best = None
    best_rmsd = None
    for candidate in candidates:
        rmsd_val = _greedy_bipartite_rmsd(pred_coords, candidate["coords"])
        if best_rmsd is None or rmsd_val < best_rmsd:
            best = candidate
            best_rmsd = rmsd_val
    return best, best_rmsd


def _to_rdkit_mol(elements: tuple[str, ...], coords: np.ndarray):
    if len(elements) != len(coords):
        return None

    lines = []
    for i, (element, xyz) in enumerate(zip(elements, coords), 1):
        x, y, z = map(float, xyz)
        # PDB HETATM columns
        lines.append(
            f"HETATM{i:5d}  X{i % 100:>2d} LIG A   1    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00          {element:>2}"
        )
    lines.append("END")

    pdb_block = "\n".join(lines) + "\n"
    mol = Chem.MolFromPDBBlock(
        pdb_block,
        removeHs=False,
        sanitize=True,
        proximityBonding=True,
    )

    if mol is None:
        mol = Chem.MolFromPDBBlock(
            pdb_block,
            removeHs=False,
            sanitize=False,
            proximityBonding=True,
        )
        if mol is not None:
            try:
                Chem.SanitizeMol(mol)
            except Exception:
                return None
    return mol


# ---------------------------------------------------------------------------
# Reproduction metrics
# ---------------------------------------------------------------------------


def _ensure_schema(system_df, chain_df):
    for column in SYSTEM_REPRODUCTION_COLUMNS:
        if column not in system_df.columns:
            system_df[column] = None

    for column in CHAIN_REPRODUCTION_COLUMNS:
        if column not in chain_df.columns:
            chain_df[column] = None

    return system_df, chain_df


def _build_aligned_predicted_structures(chain_df, structures_dir: Path, reference_structure, logger=None):
    predicted = []
    for cif_name in chain_df["cif_file"].dropna().unique():
        cif_path = structures_dir / cif_name
        if not cif_path.exists():
            if logger is not None:
                logger.warning("Predicted structure missing for reproduction metric: %s", cif_path)
            continue
        predicted.append((cif_name, _load_structure(cif_path)))

    if not predicted:
        return {}

    structures = [("__reference__", reference_structure)] + predicted
    aligned, _, _ = _align_structures_on_protein_ca(structures, save_dir=None)

    return {name: struct for name, struct in aligned.items() if name != "__reference__"}


def _infer_receptor_chain_id(chain_df, cif_name):
    rows = chain_df[(chain_df["cif_file"] == cif_name) & (chain_df["ENTITY_TYPE"] != "ligand")]
    if rows.empty:
        return None
    return str(rows.iloc[0]["CHAIN_ID"])


def _parse_ifp_vector(value):
    if value is None:
        return None
    if isinstance(value, list):
        return [int(v) for v in value]
    if isinstance(value, str):
        s = value.strip()
        if not s:
            return None
        try:
            parsed = json.loads(s)
            if isinstance(parsed, list):
                return [int(v) for v in parsed]
        except json.JSONDecodeError:
            pass
    return None


def _compute_reference_ifp_vector(reference_structure, receptor_chain_id, reference_ligand_coords, residue_order, cutoff=5.0):
    ref_model = next(iter(reference_structure))
    if receptor_chain_id not in ref_model:
        return None

    chain = ref_model[receptor_chain_id]
    residue_lookup = {res.id: res for res in chain if res.id[0] == " "}

    bits = []
    for rid in residue_order:
        res = residue_lookup.get(rid)
        if res is None:
            bits.append(0)
            continue

        on = 0
        for atom in res:
            element = (getattr(atom, "element", "") or "").upper()
            if element == "H":
                continue
            dists = np.linalg.norm(reference_ligand_coords - atom.coord, axis=1)
            if float(np.min(dists)) <= cutoff:
                on = 1
                break
        bits.append(on)

    return bits


def _compute_chain_reference_metrics(
    chain_df,
    aligned_predicted,
    reference_structure,
    reproduction_metrics,
    logger=None,
):
    """Compute per-row reproduction metrics against a reference structure."""
    reference_ligands = _collect_ligand_residue_heavy_atoms(reference_structure)

    rmsd_cache = {}
    warned_missing_ligand_match = set()
    warned_ifp_mismatch = set()

    for idx, row in chain_df.iterrows():
        cif_name = row["cif_file"]
        chain_id = str(row["CHAIN_ID"])
        entity_type = row["ENTITY_TYPE"]

        cache_key = (cif_name, chain_id, entity_type)
        if cache_key in rmsd_cache:
            cached = rmsd_cache[cache_key]
            for k, v in cached.items():
                chain_df.at[idx, k] = v
            continue

        pred = aligned_predicted.get(cif_name)
        if pred is None:
            continue

        metrics = {
            "ligand_rmsd_ref": None,
            "protein_rmsd_ref": None,
            "sucos_shape_ref": None,
            "sucos_feature_ref": None,
            "sucos_ref": None,
            "ligand_pose_overlap_ref": None,
            "pocket_coverage_ref": None,
        }

        if entity_type == "ligand":
            pred_coords, pred_signature = _collect_predicted_ligand_chain_heavy_atoms(pred, chain_id)
            if pred_coords is not None:
                match, ligand_rmsd = _find_exact_ligand_match(
                    pred_coords=pred_coords,
                    pred_signature=pred_signature,
                    reference_ligands=reference_ligands,
                )
                if "ligand_rmsd" in reproduction_metrics:
                    metrics["ligand_rmsd_ref"] = ligand_rmsd

                if match is None:
                    if logger is not None:
                        warn_key = (cif_name, chain_id)
                        if warn_key not in warned_missing_ligand_match:
                            warned_missing_ligand_match.add(warn_key)
                            logger.warning(
                                "No exact heavy-atom ligand match found for predicted ligand "
                                "(cif=%s, chain=%s, heavy_atoms=%d).",
                                cif_name,
                                chain_id,
                                len(pred_coords),
                            )
                else:
                    if "sucos" in reproduction_metrics:
                        pred_mol = _to_rdkit_mol(match["elements"], pred_coords)
                        ref_mol = _to_rdkit_mol(match["elements"], match["coords"])
                        if pred_mol is not None and ref_mol is not None:
                            sucos = compute_sucos(ref_mol=ref_mol, query_mol=pred_mol)
                            metrics["sucos_shape_ref"] = sucos.shape
                            metrics["sucos_feature_ref"] = sucos.feature
                            metrics["sucos_ref"] = sucos.score
                            metrics["ligand_pose_overlap_ref"] = sucos.score
                        elif logger is not None:
                            logger.warning(
                                "Could not build RDKit ligand molecules for SuCOS "
                                "(cif=%s, chain=%s).",
                                cif_name,
                                chain_id,
                            )

                    if "pocket_coverage" in reproduction_metrics:
                        pred_ifp = _parse_ifp_vector(row.get("ifp_distance"))
                        receptor_chain_id = _infer_receptor_chain_id(chain_df, cif_name)

                        if pred_ifp is not None and receptor_chain_id is not None:
                            pred_model = next(iter(pred))
                            if receptor_chain_id in pred_model:
                                residue_order = [
                                    res.id
                                    for res in pred_model[receptor_chain_id]
                                    if res.id[0] == " "
                                ]
                                ref_ifp = _compute_reference_ifp_vector(
                                    reference_structure=reference_structure,
                                    receptor_chain_id=receptor_chain_id,
                                    reference_ligand_coords=match["coords"],
                                    residue_order=residue_order,
                                    cutoff=5.0,
                                )
                                if ref_ifp is not None:
                                    n = min(len(pred_ifp), len(ref_ifp))
                                    if len(pred_ifp) != len(ref_ifp) and logger is not None:
                                        warn_key = (cif_name, chain_id, "ifp_len")
                                        if warn_key not in warned_ifp_mismatch:
                                            warned_ifp_mismatch.add(warn_key)
                                            logger.warning(
                                                "IFP length mismatch for pocket coverage "
                                                "(cif=%s, chain=%s): pred=%d ref=%d; using min=%d.",
                                                cif_name,
                                                chain_id,
                                                len(pred_ifp),
                                                len(ref_ifp),
                                                n,
                                            )

                                    pred_on = {i for i, v in enumerate(pred_ifp[:n]) if int(v) == 1}
                                    ref_on = {i for i, v in enumerate(ref_ifp[:n]) if int(v) == 1}

                                    if len(ref_on) == 0:
                                        if logger is not None:
                                            logger.warning(
                                                "Reference IFP has no active bits for pocket coverage "
                                                "(cif=%s, chain=%s).",
                                                cif_name,
                                                chain_id,
                                            )
                                        metrics["pocket_coverage_ref"] = None
                                    else:
                                        metrics["pocket_coverage_ref"] = float(
                                            len(pred_on & ref_on) / len(ref_on)
                                        )
                        elif logger is not None:
                            logger.warning(
                                "Skipping pocket coverage (missing ifp_distance or receptor chain) "
                                "for (cif=%s, chain=%s).",
                                cif_name,
                                chain_id,
                            )

        else:
            ref_coords, pred_coords = _collect_matched_backbone_coords(
                ref_struct=reference_structure,
                mob_struct=pred,
                chain_id=chain_id,
            )
            n = min(len(ref_coords), len(pred_coords))
            if n > 0 and "protein_rmsd" in reproduction_metrics:
                metrics["protein_rmsd_ref"] = float(_rmsd(
                    coords_a=ref_coords[:n],
                    coords_b=pred_coords[:n],
                ))

        for key, value in metrics.items():
            chain_df.at[idx, key] = value

        rmsd_cache[cache_key] = metrics

    return chain_df


def _aggregate_system_reference_metrics(system_df, chain_df):
    for idx, row in system_df.iterrows():
        model_rows = chain_df[
            (chain_df["repeat"] == row["repeat"]) &
            (chain_df["diffusion_sample"] == row["diffusion_sample"])
        ]

        ligand_rmsd_vals = model_rows.loc[
            model_rows["ENTITY_TYPE"] == "ligand", "ligand_rmsd_ref"
        ].dropna()
        protein_rmsd_vals = model_rows.loc[
            model_rows["ENTITY_TYPE"] != "ligand", "protein_rmsd_ref"
        ].dropna()

        sucos_vals = model_rows.loc[
            model_rows["ENTITY_TYPE"] == "ligand", "sucos_ref"
        ].dropna()
        pocket_vals = model_rows.loc[
            model_rows["ENTITY_TYPE"] == "ligand", "pocket_coverage_ref"
        ].dropna()

        if len(ligand_rmsd_vals) > 0:
            system_df.at[idx, "ligand_rmsd_ref_mean"] = float(ligand_rmsd_vals.astype(float).mean())
        if len(protein_rmsd_vals) > 0:
            system_df.at[idx, "protein_rmsd_ref_mean"] = float(protein_rmsd_vals.astype(float).mean())

        if len(sucos_vals) > 0:
            mean_sucos = float(sucos_vals.astype(float).mean())
            system_df.at[idx, "sucos_ref_mean"] = mean_sucos
            system_df.at[idx, "ligand_pose_overlap_ref"] = mean_sucos

        if len(pocket_vals) > 0:
            mean_pocket = float(pocket_vals.astype(float).mean())
            system_df.at[idx, "pocket_coverage_ref_mean"] = mean_pocket
            system_df.at[idx, "pocket_coverage_ref"] = mean_pocket

    return system_df


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def scaffold_reproduction_metrics(
    system_df,
    chain_df,
    reference_path: Path | None,
    wrk_dir: Path | str,
    reproduction_metrics: list[str] | None = None,
    logger: logging.Logger | None = None,
):
    """Add reproduction schema and compute reference-based reproduction metrics."""
    system_df, chain_df = _ensure_schema(system_df=system_df, chain_df=chain_df)

    enabled_metrics = set(reproduction_metrics or DEFAULT_REPRODUCTION_METRICS)

    if reference_path is None:
        if logger is not None:
            logger.info(
                "Reproduction schema initialized without reference structure. "
                "Metrics are placeholders and not computed."
            )
        return system_df, chain_df

    reference_structure = _load_structure(Path(reference_path))
    structures_dir = Path(wrk_dir) / "results" / "structures"
    aligned_predicted = _build_aligned_predicted_structures(
        chain_df=chain_df,
        structures_dir=structures_dir,
        reference_structure=reference_structure,
        logger=logger,
    )

    chain_df = _compute_chain_reference_metrics(
        chain_df=chain_df,
        aligned_predicted=aligned_predicted,
        reference_structure=reference_structure,
        reproduction_metrics=enabled_metrics,
        logger=logger,
    )
    system_df = _aggregate_system_reference_metrics(system_df=system_df, chain_df=chain_df)

    if logger is not None:
        logger.info(
            "Computed reproduction metrics using reference structure %s with metrics=%s",
            reference_path,
            sorted(enabled_metrics),
        )

    return system_df, chain_df
