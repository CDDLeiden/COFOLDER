"""Isolated ProLIF process entrypoint.

This module is deliberately the only analytics module that imports ProLIF and
MDAnalysis at module scope.  The parent process invokes it only for an explicit
``prolif`` taxonomy request.
"""

from __future__ import annotations

import json
import sys

import MDAnalysis as mda
import prolif as plf


def _field(value, *names, default=None):
    for name in names:
        candidate = getattr(value, name, None)
        if candidate is not None:
            return candidate
    return default


def _residue_payload(residue_id):
    chain = str(_field(residue_id, "chain", "chain_id", default="")).strip()
    number = int(_field(residue_id, "number", "resid", default=0))
    name = str(_field(residue_id, "name", "resname", default="")).strip().upper() or None
    return chain, number, name


def _selection(chain_id, residue_number=None):
    chain_clause = f"(chainID {chain_id} or segid {chain_id})"
    if residue_number is None:
        return chain_clause
    return f"{chain_clause} and resid {int(residue_number)}"


def main() -> None:
    request = json.load(sys.stdin)
    universe = mda.Universe(request["structure_path"])
    ligand_spec = request.get("ligand") or {}
    ligand_chain = ligand_spec.get("chain_id")
    ligand_number = ligand_spec.get("residue_number")

    if ligand_chain is None:
        candidates = universe.select_atoms("not protein and not resname HOH WAT DOD").residues
        if len(candidates) != 1:
            raise ValueError("ambiguous_reference_ligand" if len(candidates) else "ligand_not_found")
        ligand_atoms = candidates[0].atoms
        ligand_chain = str(ligand_atoms.chainIDs[0] or ligand_atoms.segids[0])
        ligand_number = int(candidates[0].resid)
    else:
        ligand_atoms = universe.select_atoms(_selection(ligand_chain, ligand_number))
        if len(ligand_atoms.residues) != 1:
            raise ValueError("ambiguous_ligand_selector" if len(ligand_atoms.residues) else "ligand_not_found")
        ligand_number = int(ligand_atoms.residues[0].resid)

    receptor_chains = request.get("receptor_chains")
    if receptor_chains:
        clauses = " or ".join(_selection(chain) for chain in receptor_chains)
        protein_atoms = universe.select_atoms(f"protein and ({clauses})")
    else:
        protein_atoms = universe.select_atoms("protein")
        receptor_chains = sorted(
            (set(protein_atoms.chainIDs) | set(protein_atoms.segids)) - {""}
        )
    if not len(protein_atoms):
        raise ValueError("receptor_chain_not_found")

    ligand_molecule = plf.Molecule.from_mda(ligand_atoms)
    protein_molecule = plf.Molecule.from_mda(protein_atoms)
    fingerprint = plf.Fingerprint(interactions=request["interactions"])
    generated = fingerprint.generate(ligand_molecule, protein_molecule, metadata=True)

    interactions = []
    for pair, interaction_map in generated.items():
        if not isinstance(pair, tuple) or len(pair) != 2:
            continue
        _, protein_residue = pair
        chain, number, _ = _residue_payload(protein_residue)
        for interaction_name, occurrences in interaction_map.items():
            if occurrences:
                interactions.append(f"{chain}:{number}:{interaction_name}")

    residue = ligand_atoms.residues[0]
    ligand_name = str(residue.resname).strip().upper() or None
    json.dump(
        {
            "ligand": {
                "chain_id": str(ligand_chain),
                "residue_number": int(ligand_number),
                "insertion_code": str(getattr(residue, "icode", "") or "").strip(),
                "residue_name": ligand_name,
            },
            "receptor_chains": sorted(str(chain) for chain in receptor_chains),
            "interactions": sorted(set(interactions)),
        },
        sys.stdout,
    )


if __name__ == "__main__":
    main()
