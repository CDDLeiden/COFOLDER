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

def _compute_chain_rmsd(chain_df, aligned_structs, wrk_dir):
    """
    Returns:
        dict[(CHAIN_ID, ENTITY_TYPE)] -> list[RMSD]

    Computes all-against-all RMSD for each chain.
    Saves full NxN RMSD matrices to Path(wrk_dir)/results/matrices/rmsd_matrix_{CHAINID}.csv
    """
    rmsd_map = defaultdict(list)
    matrices_dir = Path(wrk_dir) / "results" / "matrices"
    matrices_dir.mkdir(parents=True, exist_ok=True)

    non_ligand_df = chain_df[chain_df["ENTITY_TYPE"] != "ligand"]

    for (chain_id, entity_type), group in non_ligand_df.groupby(["CHAIN_ID", "ENTITY_TYPE"]):
        # Collect all poses and their coordinates
        pose_names = []
        coords_list = []

        for idx, row in group.iterrows():
            cif_name = row["cif_file"]
            struct = aligned_structs.get(cif_name)
            if struct is None:
                continue

            atoms = _collect_chain_atoms(struct, chain_id)
            if not atoms:
                continue

            pose_names.append(cif_name)
            coords_list.append(np.array([a.coord for a in atoms]))

        n = len(coords_list)
        if n < 2:
            continue

        # Compute NxN RMSD matrix
        mat = np.zeros((n, n), dtype=float)
        for i in range(n):
            for j in range(i, n):
                a = coords_list[i]
                b = coords_list[j]
                m = min(len(a), len(b))
                val = _rmsd(a[:m], b[:m])
                mat[i, j] = val
                mat[j, i] = val

        # --- save full NxN RMSD matrix ---
        df = pd.DataFrame(mat, index=pose_names, columns=pose_names)
        df.to_csv(matrices_dir / f"rmsd_matrix_chain_{entity_type}_{chain_id}.csv")

        # Flatten off-diagonal values to list for existing pipeline
        off_diag = mat[np.triu_indices(n, k=1)]
        rmsd_map[(chain_id, entity_type)] = off_diag.tolist()

    return rmsd_map


def _compute_ligand_rmsd(chain_df, aligned_structs, wrk_dir):
    """
    Returns:
        dict[ligand_molecule_id] -> list[RMSD]

    Computes all-against-all RMSD for each ligand.
    Saves full NxN RMSD matrices to Path(wrk_dir)/results/matrices/rmsd_matrix_{ligand_id}.csv
    """
    rmsd_map = defaultdict(list)
    matrices_dir = Path(wrk_dir) / "results" / "matrices"
    matrices_dir.mkdir(parents=True, exist_ok=True)

    ligand_df = chain_df[chain_df["ENTITY_TYPE"] == "ligand"]

    for ligand_id, group in ligand_df.groupby("ligand_molecule_id"):
        # Collect all poses and their coordinates
        pose_names = []
        coords_list = []

        for idx, row in group.iterrows():
            cif_name = row["cif_file"]
            chain_id = row["CHAIN_ID"]
            struct = aligned_structs.get(cif_name)
            if struct is None:
                continue

            atoms = _collect_chain_atoms(struct, chain_id)
            if not atoms:
                continue

            pose_names.append(cif_name)
            coords_list.append(np.array([a.coord for a in atoms]))

        n = len(coords_list)
        if n < 2:
            continue

        # Compute NxN RMSD matrix
        mat = np.zeros((n, n), dtype=float)
        for i in range(n):
            for j in range(i, n):
                a = coords_list[i]
                b = coords_list[j]
                m = min(len(a), len(b))
                val = _rmsd(a[:m], b[:m])
                mat[i, j] = val
                mat[j, i] = val

        # --- save full NxN RMSD matrix ---
        df = pd.DataFrame(mat, index=pose_names, columns=pose_names)
        df.to_csv(matrices_dir / f"rmsd_matrix_ligand_{ligand_id}.csv")

        # Flatten off-diagonal values to list for existing pipeline
        off_diag = mat[np.triu_indices(n, k=1)]
        rmsd_map[ligand_id] = off_diag.tolist()

    return rmsd_map