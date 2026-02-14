"""Model reproduction and bias metrics for validation workflows."""

from pathlib import Path
import logging
import numpy as np

from Bio.PDB import MMCIFParser, PDBParser

from cofolder.modules.analytics.align import (
    _align_structures_on_protein_ca,
    _rmsd,
)

SYSTEM_REPRODUCTION_COLUMNS = (
    "ligand_rmsd_ref_mean",
    "protein_rmsd_ref_mean",
    "ifp_similarity_ref_mean",
    "pocket_coverage_ref",
    "ligand_pose_overlap_ref",
    "protein_sequence_identity_train_max",
    "ligand_fingerprint_similarity_train_max",
)

CHAIN_REPRODUCTION_COLUMNS = (
    "ligand_rmsd_ref",
    "protein_rmsd_ref",
    "ifp_similarity_ref",
    "pocket_coverage_ref",
    "ligand_pose_overlap_ref",
    "protein_sequence_identity_train",
    "ligand_fingerprint_similarity_train",
)


def _load_structure(path: Path):
    suffix = path.suffix.lower()
    if suffix in {".cif", ".mmcif"}:
        parser = MMCIFParser(QUIET=True)
    elif suffix == ".pdb":
        parser = PDBParser(QUIET=True)
    else:
        raise ValueError(f"Unsupported structure format: {path}")
    return parser.get_structure(path.stem, str(path))


def _collect_backbone_atoms(structure, chain_id):
    atoms = []
    for model in structure:
        if chain_id not in model:
            continue
        for residue in model[chain_id]:
            if residue.id[0] != " ":
                continue
            for atom_name in ("N", "CA", "C", "O"):
                if atom_name in residue:
                    atoms.append(residue[atom_name])
    return atoms


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
        lig for lig in reference_ligands
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

    aligned_predicted = {
        name: struct
        for name, struct in aligned.items()
        if name != "__reference__"
    }
    return aligned_predicted


def _compute_chain_reference_rmsd(chain_df, aligned_predicted, reference_structure, logger=None):
    """Compute per-row RMSD against a reference structure."""
    reference_ligands = _collect_ligand_residue_heavy_atoms(reference_structure)
    rmsd_cache = {}
    warned_missing_ligand_match = set()

    for idx, row in chain_df.iterrows():
        cif_name = row["cif_file"]
        chain_id = str(row["CHAIN_ID"])
        entity_type = row["ENTITY_TYPE"]

        cache_key = (cif_name, chain_id, entity_type)
        if cache_key in rmsd_cache:
            ligand_rmsd, protein_rmsd = rmsd_cache[cache_key]
            chain_df.at[idx, "ligand_rmsd_ref"] = ligand_rmsd
            chain_df.at[idx, "protein_rmsd_ref"] = protein_rmsd
            continue

        pred = aligned_predicted.get(cif_name)
        if pred is None:
            continue

        ligand_rmsd = None
        protein_rmsd = None

        if entity_type == "ligand":
            pred_coords, pred_signature = _collect_predicted_ligand_chain_heavy_atoms(pred, chain_id)
            if pred_coords is not None:
                _, ligand_rmsd = _find_exact_ligand_match(
                    pred_coords=pred_coords,
                    pred_signature=pred_signature,
                    reference_ligands=reference_ligands,
                )
                if ligand_rmsd is None and logger is not None:
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
                chain_df.at[idx, "ligand_rmsd_ref"] = ligand_rmsd
        else:
            ref_coords, pred_coords = _collect_matched_backbone_coords(
                ref_struct=reference_structure,
                mob_struct=pred,
                chain_id=chain_id,
            )
            n = min(len(ref_coords), len(pred_coords))
            if n > 0:
                protein_rmsd = float(_rmsd(
                    coords_a=ref_coords[:n],
                    coords_b=pred_coords[:n],
                ))
                chain_df.at[idx, "protein_rmsd_ref"] = protein_rmsd

        rmsd_cache[cache_key] = (ligand_rmsd, protein_rmsd)

    return chain_df


def _aggregate_system_reference_rmsd(system_df, chain_df):
    for idx, row in system_df.iterrows():
        model_rows = chain_df[
            (chain_df["repeat"] == row["repeat"]) &
            (chain_df["diffusion_sample"] == row["diffusion_sample"])
        ]

        ligand_vals = model_rows.loc[
            model_rows["ENTITY_TYPE"] == "ligand", "ligand_rmsd_ref"
        ].dropna()
        protein_vals = model_rows.loc[
            model_rows["ENTITY_TYPE"] != "ligand", "protein_rmsd_ref"
        ].dropna()

        if len(ligand_vals) > 0:
            system_df.at[idx, "ligand_rmsd_ref_mean"] = float(ligand_vals.astype(float).mean())
        if len(protein_vals) > 0:
            system_df.at[idx, "protein_rmsd_ref_mean"] = float(protein_vals.astype(float).mean())

    return system_df


def scaffold_reproduction_metrics(
    system_df,
    chain_df,
    reference_path: Path | None,
    wrk_dir: Path | str,
    logger: logging.Logger | None = None,
):
    """Add reproduction schema and compute phase-2 RMSD metrics when possible."""
    system_df, chain_df = _ensure_schema(system_df=system_df, chain_df=chain_df)

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

    chain_df = _compute_chain_reference_rmsd(
        chain_df=chain_df,
        aligned_predicted=aligned_predicted,
        reference_structure=reference_structure,
        logger=logger,
    )
    system_df = _aggregate_system_reference_rmsd(system_df=system_df, chain_df=chain_df)

    if logger is not None:
        logger.info(
            "Computed reference RMSD metrics for reproduction scaffold using %s",
            reference_path,
        )

    return system_df, chain_df
