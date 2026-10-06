"""Isolated ProLIF process entrypoint.

This module is deliberately the only analytics module that imports ProLIF and
MDAnalysis at module scope. The parent process invokes it only for an explicit
``prolif`` taxonomy request.
"""

from __future__ import annotations

import json
import math
import string
import sys
from contextlib import ExitStack
from pathlib import Path
from tempfile import TemporaryDirectory

import gemmi
import MDAnalysis as mda
import prolif as plf

_PDB_CHAIN_IDS = string.ascii_uppercase + string.ascii_lowercase + string.digits
_INTERACTION_ROLES = {
    "Hydrophobic": ("hydrophobic", "hydrophobic"),
    "HBAcceptor": ("acceptor", "donor"),
    "HBDonor": ("donor", "acceptor"),
    "PiStacking": ("pi", "pi"),
    "Anionic": ("anion", "cation"),
    "Cationic": ("cation", "anion"),
    "CationPi": ("cation", "pi"),
    "PiCation": ("pi", "cation"),
    "VdWContact": ("van_der_waals", "van_der_waals"),
}


def _field(value, *names, default=None):
    for name in names:
        candidate = getattr(value, name, None)
        if candidate is not None:
            return candidate
    return default


def _residue_payload(residue_id):
    chain = str(_field(residue_id, "chain", "chain_id", default="")).strip()
    number = int(_field(residue_id, "number", "resid", default=0))
    name = (
        str(_field(residue_id, "name", "resname", default="")).strip().upper() or None
    )
    return chain, number, name


def _selection(chain_id, residue_number=None, insertion_code=""):
    chain_clause = f"(chainID {chain_id} or segid {chain_id})"
    if residue_number is None:
        return chain_clause
    selection = f"{chain_clause} and resid {int(residue_number)}"
    if str(insertion_code).strip():
        selection += f" and icode {str(insertion_code).strip()}"
    return selection


def _prepare_topology(path: Path, stack: ExitStack):
    """Return an MDAnalysis-readable path and original/temporary chain maps."""

    if path.suffix.lower() not in {".cif", ".mmcif"}:
        return path, {}, {}

    structure = gemmi.read_structure(str(path))
    if not len(structure):
        raise ValueError("structure_contains_no_model")
    chain_names = []
    for model in structure:
        for chain in model:
            if chain.name not in chain_names:
                chain_names.append(chain.name)
    if len(chain_names) > len(_PDB_CHAIN_IDS):
        raise ValueError("too_many_chains_for_temporary_pdb")
    original_to_alias = {
        name: _PDB_CHAIN_IDS[index] for index, name in enumerate(chain_names)
    }
    alias_to_original = {alias: name for name, alias in original_to_alias.items()}
    for model in structure:
        for chain in model:
            chain.name = original_to_alias[chain.name]

    temporary_dir = Path(
        stack.enter_context(TemporaryDirectory(prefix="cofolder-prolif-"))
    )
    pdb_path = temporary_dir / "structure.pdb"
    structure.write_pdb(str(pdb_path))
    return pdb_path, original_to_alias, alias_to_original


def _atom_payload(atom, alias_to_original):
    chain = str(_field(atom, "chainID", "segid", default="")).strip()
    chain = alias_to_original.get(chain, chain)
    element = str(_field(atom, "element", "type", default="")).strip() or None
    serial = _field(atom, "id", default=None)
    return {
        "chain_id": chain,
        "residue_number": int(atom.resid),
        "insertion_code": str(_field(atom, "icode", default="") or "").strip(),
        "residue_name": str(atom.resname).strip().upper() or None,
        "atom_name": str(atom.name).strip(),
        "element": element,
        "atom_serial": int(serial) if serial is not None else None,
        "source_index": int(atom.index),
    }


def _geometry_payload(metadata):
    measurements = []
    for name, value in metadata.items():
        if name in {"indices", "parent_indices"} or not isinstance(value, (int, float)):
            continue
        numeric = float(value)
        if not math.isfinite(numeric):
            continue
        lowered = str(name).lower()
        unit = (
            "degree"
            if "angle" in lowered
            else "angstrom" if "distance" in lowered else "unitless"
        )
        measurements.append({"name": str(name), "value": numeric, "unit": unit})
    return sorted(measurements, key=lambda item: item["name"])


def _event_sort_key(event):
    def atom_key(atom):
        return (
            atom["chain_id"],
            atom["residue_number"],
            atom["insertion_code"],
            atom["atom_name"],
            atom["source_index"],
        )

    return (
        event["interaction_key"],
        tuple(atom_key(atom) for atom in event["protein_atoms"]),
        tuple(atom_key(atom) for atom in event["ligand_atoms"]),
        tuple((item["name"], item["value"]) for item in event["geometry"]),
    )


def main() -> None:
    request = json.load(sys.stdin)
    source_path = Path(request["structure_path"])
    with ExitStack() as stack:
        topology_path, original_to_alias, alias_to_original = _prepare_topology(
            source_path, stack
        )
        universe = mda.Universe(str(topology_path))
        ligand_spec = request.get("ligand") or {}
        original_ligand_chain = ligand_spec.get("chain_id")
        ligand_chain = original_to_alias.get(
            original_ligand_chain, original_ligand_chain
        )
        ligand_number = ligand_spec.get("residue_number")
        ligand_icode = ligand_spec.get("insertion_code") or ""

        if ligand_chain is None:
            candidates = universe.select_atoms(
                "not protein and not resname HOH WAT DOD"
            ).residues
            if len(candidates) != 1:
                raise ValueError(
                    "ambiguous_reference_ligand"
                    if len(candidates)
                    else "ligand_not_found"
                )
            ligand_atoms = candidates[0].atoms
            ligand_chain = str(ligand_atoms.chainIDs[0] or ligand_atoms.segids[0])
            original_ligand_chain = alias_to_original.get(ligand_chain, ligand_chain)
            ligand_number = int(candidates[0].resid)
        else:
            ligand_atoms = universe.select_atoms(
                _selection(ligand_chain, ligand_number, ligand_icode)
            )
            if len(ligand_atoms.residues) != 1:
                raise ValueError(
                    "ambiguous_ligand_selector"
                    if len(ligand_atoms.residues)
                    else "ligand_not_found"
                )
            ligand_number = int(ligand_atoms.residues[0].resid)

        original_receptor_chains = request.get("receptor_chains")
        if original_receptor_chains:
            receptor_chains = [
                original_to_alias.get(chain, chain)
                for chain in original_receptor_chains
            ]
            clauses = " or ".join(_selection(chain) for chain in receptor_chains)
            protein_atoms = universe.select_atoms(f"protein and ({clauses})")
        else:
            protein_atoms = universe.select_atoms("protein")
            receptor_chains = sorted(
                (set(protein_atoms.chainIDs) | set(protein_atoms.segids)) - {""}
            )
            original_receptor_chains = [
                alias_to_original.get(str(chain), str(chain))
                for chain in receptor_chains
            ]
        if not len(protein_atoms):
            raise ValueError("receptor_chain_not_found")

        # Prediction outputs commonly omit explicit hydrogens. MDAnalysis can still
        # infer the molecular graph, but requires ``force`` for that topology.
        ligand_molecule = plf.Molecule.from_mda(ligand_atoms, force=True)
        protein_molecule = plf.Molecule.from_mda(protein_atoms, force=True)
        fingerprint = plf.Fingerprint(interactions=request["interactions"])
        generated = fingerprint.generate(
            ligand_molecule, protein_molecule, metadata=True
        )

        events = []
        for pair, interaction_map in generated.items():
            if not isinstance(pair, tuple) or len(pair) != 2:
                continue
            _, protein_residue = pair
            protein_chain, protein_number, _ = _residue_payload(protein_residue)
            original_protein_chain = alias_to_original.get(protein_chain, protein_chain)
            for interaction_name, occurrences in interaction_map.items():
                ligand_role, protein_role = _INTERACTION_ROLES[str(interaction_name)]
                for metadata in occurrences:
                    parent_indices = metadata.get("parent_indices") or {}
                    ligand_event_atoms = tuple(
                        _atom_payload(ligand_atoms[int(index)], alias_to_original)
                        for index in parent_indices.get("ligand", ())
                    )
                    protein_event_atoms = tuple(
                        _atom_payload(protein_atoms[int(index)], alias_to_original)
                        for index in parent_indices.get("protein", ())
                    )
                    insertion_code = (
                        protein_event_atoms[0]["insertion_code"]
                        if protein_event_atoms
                        else ""
                    )
                    suffix = f"{protein_number}{insertion_code}"
                    events.append(
                        {
                            "interaction_key": (
                                f"{original_protein_chain}:{suffix}:{interaction_name}"
                            ),
                            "interaction_type": str(interaction_name),
                            "ligand_role": ligand_role,
                            "protein_role": protein_role,
                            "ligand_atoms": list(ligand_event_atoms),
                            "protein_atoms": list(protein_event_atoms),
                            "geometry": _geometry_payload(metadata),
                        }
                    )

        events.sort(key=_event_sort_key)
        interactions = sorted({event["interaction_key"] for event in events})
        residue = ligand_atoms.residues[0]
        ligand_name = str(residue.resname).strip().upper() or None
        json.dump(
            {
                "ligand": {
                    "chain_id": str(original_ligand_chain),
                    "residue_number": int(ligand_number),
                    "insertion_code": str(getattr(residue, "icode", "") or "").strip(),
                    "residue_name": ligand_name,
                },
                "receptor_chains": sorted(
                    str(chain) for chain in original_receptor_chains
                ),
                "interactions": interactions,
                "events": events,
            },
            sys.stdout,
        )


if __name__ == "__main__":
    main()
