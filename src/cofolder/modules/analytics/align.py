from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, Superimposer, MMCIFIO


# ==============================================================
# STRUCTURE HELPERS
# ==============================================================

def _save_aligned_structures(aligned_structs, out_dir):
    """
    Save aligned structures as CIF files.

    Parameters
    ----------
    aligned_structs : dict[name -> structure]
    out_dir : Path
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    io = MMCIFIO()

    for name, struct in aligned_structs.items():
        out_path = out_dir / name
        io.set_structure(struct)
        io.save(str(out_path))


def _load_structures(cif_paths):
    parser = MMCIFParser(QUIET=True)
    structures = []
    for p in cif_paths:
        s = parser.get_structure(p.stem, str(p))
        structures.append((p.name, s))
    return structures


def _get_ca_atoms(structure):
    """Collect all protein Cα atoms across all chains."""
    cas = []
    for model in structure:
        for chain in model:
            for res in chain:
                if res.id[0] != " ":
                    continue
                if "CA" in res:
                    cas.append(res["CA"])
    return cas


def _collect_chain_atoms(structure, chain_id):
    atoms = []
    for model in structure:
        if chain_id in model:
            for res in model[chain_id]:
                for atom in res:
                    atoms.append(atom)
    return atoms


def _align_structures_on_protein_ca(structures, save_dir=None):
    """
    Align all structures to the first structure using protein Cα atoms.

    Parameters
    ----------
    structures : list[(name, structure)]
    save_dir : Path | None
        If provided, aligned CIFs are written here.

    Returns
    -------
    aligned_structs : dict[name -> structure]
    ref_struct
    ref_name
    """
    ref_name, ref_struct = structures[0]
    ref_ca = _get_ca_atoms(ref_struct)

    aligned = {ref_name: ref_struct}
    sup = Superimposer()

    for name, struct in structures[1:]:
        mob_ca = _get_ca_atoms(struct)
        n = min(len(ref_ca), len(mob_ca))
        if n == 0:
            continue

        sup.set_atoms(ref_ca[:n], mob_ca[:n])
        sup.apply(struct.get_atoms())
        aligned[name] = struct

    # --- save if requested ---
    if save_dir is not None:
        _save_aligned_structures(aligned, save_dir)

    return aligned, ref_struct, ref_name


def _rmsd(coords_a, coords_b):
    diff = coords_a - coords_b
    return np.sqrt((diff * diff).sum(axis=1).mean())


# ==============================================================
# RMSD COMPUTATION
# ==============================================================

def _compute_chain_rmsd(chain_df, aligned_structs, ref_struct, ref_name):
    """
    Returns:
        dict[(CHAIN_ID, ENTITY_TYPE)] -> list[RMSD]
    """
    rmsd_map = defaultdict(list)

    non_ligand_df = chain_df[chain_df["ENTITY_TYPE"] != "ligand"]

    for (chain_id, entity_type), group in non_ligand_df.groupby(["CHAIN_ID", "ENTITY_TYPE"]):
        ref_atoms = _collect_chain_atoms(ref_struct, chain_id)
        if not ref_atoms:
            continue

        ref_coords = np.array([a.coord for a in ref_atoms])

        for cif_name, struct in aligned_structs.items():
            if cif_name == ref_name:
                continue

            mob_atoms = _collect_chain_atoms(struct, chain_id)
            if not mob_atoms:
                continue

            n = min(len(ref_coords), len(mob_atoms))
            mob_coords = np.array([a.coord for a in mob_atoms[:n]])

            val = _rmsd(ref_coords[:n], mob_coords)
            rmsd_map[(chain_id, entity_type)].append(val)

    return rmsd_map


def _compute_ligand_rmsd(chain_df, aligned_structs, ref_struct, ref_name):
    """
    Returns:
        dict[ligand_molecule_id] -> list[RMSD]
    """
    rmsd_map = defaultdict(list)
    ligand_df = chain_df[chain_df["ENTITY_TYPE"] == "ligand"]

    for ligand_id, group in ligand_df.groupby("ligand_molecule_id"):
        # determine reference chain_id for this ligand
        ref_row = group[group["cif_file"] == ref_name]
        if ref_row.empty:
            ref_row = group.iloc[[0]]

        ref_chain_id = ref_row.iloc[0]["CHAIN_ID"]
        ref_atoms = _collect_chain_atoms(ref_struct, ref_chain_id)
        if not ref_atoms:
            continue

        ref_coords = np.array([a.coord for a in ref_atoms])

        for cif_name, struct in aligned_structs.items():
            if cif_name == ref_name:
                continue

            row = group[group["cif_file"] == cif_name]
            if row.empty:
                continue

            chain_id = row.iloc[0]["CHAIN_ID"]
            mob_atoms = _collect_chain_atoms(struct, chain_id)
            if not mob_atoms:
                continue

            n = min(len(ref_coords), len(mob_atoms))
            mob_coords = np.array([a.coord for a in mob_atoms[:n]])

            val = _rmsd(ref_coords[:n], mob_coords)
            rmsd_map[ligand_id].append(val)

    return rmsd_map