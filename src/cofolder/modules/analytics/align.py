from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser, Superimposer, MMCIFIO


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
    structures = []
    for p in cif_paths:
        suffix = p.suffix.lower()
        if suffix in {".cif", ".mmcif"}:
            parser = MMCIFParser(QUIET=True)
        elif suffix == ".pdb":
            parser = PDBParser(QUIET=True)
        else:
            raise ValueError(f"Unsupported structure format: {p}")

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


def _get_matched_ca_atoms(ref_struct, mob_struct):
    """Collect matched CA atoms by chain ID and residue ID."""
    ref_atoms = []
    mob_atoms = []

    ref_model = next(iter(ref_struct))
    mob_model = next(iter(mob_struct))

    common_chains = sorted(set(c.id for c in ref_model) & set(c.id for c in mob_model))

    for chain_id in common_chains:
        ref_chain = ref_model[chain_id]
        mob_chain = mob_model[chain_id]

        ref_res = {
            res.id: res
            for res in ref_chain
            if res.id[0] == " " and "CA" in res
        }
        mob_res = {
            res.id: res
            for res in mob_chain
            if res.id[0] == " " and "CA" in res
        }

        common_res = sorted(set(ref_res.keys()) & set(mob_res.keys()))
        for rid in common_res:
            ref_atoms.append(ref_res[rid]["CA"])
            mob_atoms.append(mob_res[rid]["CA"])

    return ref_atoms, mob_atoms


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
        matched_ref_ca, matched_mob_ca = _get_matched_ca_atoms(ref_struct, struct)

        if matched_ref_ca and matched_mob_ca:
            align_ref = matched_ref_ca
            align_mob = matched_mob_ca
        else:
            # Fallback for edge cases where residue IDs do not overlap.
            mob_ca = _get_ca_atoms(struct)
            n = min(len(ref_ca), len(mob_ca))
            align_ref = ref_ca[:n]
            align_mob = mob_ca[:n]

        n = min(len(align_ref), len(align_mob))
        if n == 0:
            continue

        sup.set_atoms(align_ref[:n], align_mob[:n])
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
