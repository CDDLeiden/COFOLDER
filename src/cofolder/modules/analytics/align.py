from pathlib import Path
from collections import defaultdict
import logging
import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner
from Bio.PDB import MMCIFParser, PDBParser, Superimposer, MMCIFIO
from Bio.SeqUtils import seq1


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
    chain_mappings = []

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
        if common_res:
            chain_mappings.append(
                {
                    "reference_chain_id": chain_id,
                    "mobile_chain_id": chain_id,
                    "matched_residue_pairs": [(rid, rid) for rid in common_res],
                    "reference_sequence_length": len(ref_res),
                    "mobile_sequence_length": len(mob_res),
                }
            )
        for rid in common_res:
            ref_atoms.append(ref_res[rid]["CA"])
            mob_atoms.append(mob_res[rid]["CA"])

    metadata = {
        "strategy": "exact_residue_id",
        "pair_count": len(ref_atoms),
        "chain_mappings": chain_mappings,
    }
    return ref_atoms, mob_atoms, metadata


def _collect_protein_chain_data(structure):
    chains = []
    model = next(iter(structure))
    for chain in model:
        residues = [
            res
            for res in chain
            if res.id[0] == " " and "CA" in res
        ]
        if not residues:
            continue
        sequence = "".join(
            seq1(str(res.get_resname()).strip().upper(), undef_code="X")
            for res in residues
        )
        chains.append(
            {
                "chain_id": chain.id,
                "residues": residues,
                "sequence": sequence,
            }
        )
    return chains


def _build_aligned_index_pairs(alignment):
    pairs = []
    coordinates = alignment.coordinates
    for idx in range(coordinates.shape[1] - 1):
        ref_start = int(coordinates[0, idx])
        ref_end = int(coordinates[0, idx + 1])
        mob_start = int(coordinates[1, idx])
        mob_end = int(coordinates[1, idx + 1])
        span = min(ref_end - ref_start, mob_end - mob_start)
        for offset in range(max(span, 0)):
            pairs.append((ref_start + offset, mob_start + offset))
    return pairs


def _get_sequence_fallback_ca_atoms(ref_struct, mob_struct):
    ref_chains = _collect_protein_chain_data(ref_struct)
    mob_chains = _collect_protein_chain_data(mob_struct)
    if not ref_chains or not mob_chains:
        return [], [], None

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5

    candidate_pairs = []
    for ref_chain in ref_chains:
        for mob_chain in mob_chains:
            alignment = aligner.align(
                ref_chain["sequence"],
                mob_chain["sequence"],
            )[0]
            index_pairs = _build_aligned_index_pairs(alignment)
            if not index_pairs:
                continue
            candidate_pairs.append(
                {
                    "score": float(alignment.score),
                    "pair_count": len(index_pairs),
                    "reference_chain": ref_chain,
                    "mobile_chain": mob_chain,
                    "index_pairs": index_pairs,
                }
            )

    if not candidate_pairs:
        return [], [], None

    candidate_pairs.sort(
        key=lambda item: (item["score"], item["pair_count"]),
        reverse=True,
    )

    used_ref_chains = set()
    used_mob_chains = set()
    ref_atoms = []
    mob_atoms = []
    chain_mappings = []

    for candidate in candidate_pairs:
        ref_chain_id = candidate["reference_chain"]["chain_id"]
        mob_chain_id = candidate["mobile_chain"]["chain_id"]
        if ref_chain_id in used_ref_chains or mob_chain_id in used_mob_chains:
            continue

        residue_pairs = []
        for ref_idx, mob_idx in candidate["index_pairs"]:
            ref_residue = candidate["reference_chain"]["residues"][ref_idx]
            mob_residue = candidate["mobile_chain"]["residues"][mob_idx]
            ref_atoms.append(ref_residue["CA"])
            mob_atoms.append(mob_residue["CA"])
            residue_pairs.append((ref_residue.id, mob_residue.id))

        if not residue_pairs:
            continue

        used_ref_chains.add(ref_chain_id)
        used_mob_chains.add(mob_chain_id)
        chain_mappings.append(
            {
                "reference_chain_id": ref_chain_id,
                "mobile_chain_id": mob_chain_id,
                "matched_residue_pairs": residue_pairs,
                "reference_sequence_length": len(candidate["reference_chain"]["residues"]),
                "mobile_sequence_length": len(candidate["mobile_chain"]["residues"]),
                "alignment_score": candidate["score"],
            }
        )

    if not ref_atoms:
        return [], [], None

    metadata = {
        "strategy": "sequence_fallback",
        "pair_count": len(ref_atoms),
        "chain_mappings": chain_mappings,
    }
    return ref_atoms, mob_atoms, metadata


def _collect_chain_atoms(structure, chain_id):
    atoms = []
    for model in structure:
        if chain_id in model:
            for res in model[chain_id]:
                for atom in res:
                    atoms.append(atom)
    return atoms


def _align_structures_on_protein_ca(
    structures,
    save_dir=None,
    logger: logging.Logger | None = None,
    return_alignment_details: bool = False,
):
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
    alignment_details = {}
    sup = Superimposer()

    for name, struct in structures[1:]:
        matched_ref_ca, matched_mob_ca, metadata = _get_matched_ca_atoms(ref_struct, struct)

        if matched_ref_ca and matched_mob_ca:
            align_ref = matched_ref_ca
            align_mob = matched_mob_ca
        else:
            matched_ref_ca, matched_mob_ca, metadata = _get_sequence_fallback_ca_atoms(
                ref_struct,
                struct,
            )
            if matched_ref_ca and matched_mob_ca:
                align_ref = matched_ref_ca
                align_mob = matched_mob_ca
            else:
                # Last-resort fallback for structures without usable chain or sequence overlap.
                mob_ca = _get_ca_atoms(struct)
                n = min(len(ref_ca), len(mob_ca))
                align_ref = ref_ca[:n]
                align_mob = mob_ca[:n]
                metadata = {
                    "strategy": "naive_first_n",
                    "pair_count": n,
                    "chain_mappings": [],
                    "reference_sequence_length": len(ref_ca),
                    "mobile_sequence_length": len(mob_ca),
                }

        n = min(len(align_ref), len(align_mob))
        if n == 0:
            if logger is not None:
                logger.warning(
                    "Could not align predicted structure %s to reference %s: no usable protein C-alpha pairs found.",
                    name,
                    ref_name,
                )
            continue

        sup.set_atoms(align_ref[:n], align_mob[:n])
        sup.apply(struct.get_atoms())
        aligned[name] = struct
        metadata["pair_count"] = n
        metadata["post_superposition_rms"] = float(sup.rms)
        alignment_details[name] = metadata

        if logger is not None:
            log_message = (
                "Aligned predicted structure %s to reference %s using %s (%d CA pairs, RMS %.3f A)."
            )
            log_args = (
                name,
                ref_name,
                metadata["strategy"],
                n,
                float(sup.rms),
            )
            if metadata["strategy"] == "naive_first_n":
                logger.warning(log_message, *log_args)
            else:
                logger.info(log_message, *log_args)

            chain_summaries = []
            for chain_mapping in metadata.get("chain_mappings", []):
                residue_pairs = chain_mapping.get("matched_residue_pairs", [])
                if residue_pairs:
                    first_pair = residue_pairs[0]
                    last_pair = residue_pairs[-1]
                    residue_span = (
                        f"ref {chain_mapping['reference_chain_id']}:{first_pair[0]}->{last_pair[0]} | "
                        f"mob {chain_mapping['mobile_chain_id']}:{first_pair[1]}->{last_pair[1]}"
                    )
                else:
                    residue_span = "no residue span"
                chain_summaries.append(
                    {
                        "reference_chain_id": chain_mapping["reference_chain_id"],
                        "mobile_chain_id": chain_mapping["mobile_chain_id"],
                        "reference_sequence_length": chain_mapping.get("reference_sequence_length"),
                        "mobile_sequence_length": chain_mapping.get("mobile_sequence_length"),
                        "matched_residue_pairs": len(residue_pairs),
                        "alignment_score": chain_mapping.get("alignment_score"),
                        "residue_span": residue_span,
                    }
                )
            logger.debug(
                "Alignment details for %s: strategy=%s chain_mappings=%s",
                name,
                metadata["strategy"],
                chain_summaries,
            )

    # --- save if requested ---
    if save_dir is not None:
        _save_aligned_structures(aligned, save_dir)

    if return_alignment_details:
        return aligned, ref_struct, ref_name, alignment_details
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
