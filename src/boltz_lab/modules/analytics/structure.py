import gemmi
import numpy as np
from pathlib import Path
import json

import logging

logger = logging.getLogger(__name__)

def distance_interaction_fingerprint(
    cif_path: Path,
    receptor_chain: str,
    ligand_chain: str,
    cutoff: float = 5.0,
    to_json: bool = True
) -> np.ndarray | str:
    """
    Return a residue-level interaction fingerprint for receptor_chain.

    One bit per receptor residue:
      1 = at least one heavy-atom contact with ligand_chain within cutoff
      0 = no contact

    Parameters
    ----------
    cif_path : Path
        Path to CIF file.
    receptor_chain : str
        Receptor protein chain ID.
    ligand_chain : str
        Ligand chain ID.
    cutoff : float
        Distance cutoff in Angstroms for contact.
    to_json : bool
        If True, return JSON string instead of np.ndarray (safe for DataFrame/CSV).

    Returns
    -------
    np.ndarray or str
        Interaction fingerprint, either as np.ndarray or JSON string.
    """
    structure = gemmi.read_structure(str(cif_path))
    model = structure[0]

    try:
        chain_rec = model[receptor_chain]
        chain_lig = model[ligand_chain]
    except KeyError:
        logger.warning(
            "Missing receptor (%s) or ligand (%s) chain in CIF: %s",
            receptor_chain, ligand_chain, cif_path
        )
        return json.dumps([]) if to_json else np.array([], dtype=int)

    rec_residues = list(chain_rec)
    lig_atoms = [
        atom
        for res in chain_lig
        for atom in res
        if atom.element.name != "H"
    ]

    bitvector = np.zeros(len(rec_residues), dtype=int)

    for i, res in enumerate(rec_residues):
        rec_atoms = [atom for atom in res if atom.element.name != "H"]

        for ra in rec_atoms:
            for la in lig_atoms:
                if ra.pos.dist(la.pos) <= cutoff:
                    bitvector[i] = 1
                    break
            if bitvector[i]:
                break

    if to_json:
        return json.dumps(bitvector.tolist())
    return bitvector