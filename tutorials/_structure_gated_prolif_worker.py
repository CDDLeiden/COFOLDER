"""Isolated, tutorial-local preparation and ProLIF extraction worker."""

from __future__ import annotations

import json
import math
import shutil
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

import gemmi
import prolif as plf
from rdkit import Chem
from rdkit.Chem import AllChem


def _atom_payload(atom: Chem.Atom) -> dict[str, object]:
    info = atom.GetPDBResidueInfo()
    return {
        "chain_id": info.GetChainId().strip() if info else "",
        "residue_number": info.GetResidueNumber() if info else 0,
        "insertion_code": info.GetInsertionCode().strip() if info else "",
        "residue_name": info.GetResidueName().strip() if info else None,
        "atom_name": info.GetName().strip()
        if info
        else f"{atom.GetSymbol()}{atom.GetIdx()}",
        "element": atom.GetSymbol(),
        "source_index": atom.GetIdx(),
    }


def _write_selected_pdb(
    source: Path,
    destination: Path,
    *,
    ligand_chain: str,
    ligand_number: int,
    protein_only: bool,
) -> None:
    structure = gemmi.read_structure(str(source))
    if not len(structure):
        raise ValueError("structure_contains_no_model")
    selected = gemmi.Structure()
    selected.name = structure.name
    model = gemmi.Model("1")
    for chain in structure[0]:
        output_chain = gemmi.Chain(chain.name)
        for residue in chain:
            is_ligand = (
                chain.name == ligand_chain and residue.seqid.num == ligand_number
            )
            keep = (
                residue.entity_type == gemmi.EntityType.Polymer
                if protein_only
                else is_ligand
            )
            if keep:
                output_chain.add_residue(residue.clone())
        if len(output_chain):
            model.add_chain(output_chain)
    selected.add_model(model)
    selected.write_pdb(str(destination))


def _prepare_ligand(
    ligand_pdb: Path, smiles: str, residue_name: str, chain: str, number: int
) -> tuple[plf.Molecule, dict[str, object]]:
    template = Chem.MolFromSmiles(smiles)
    if template is None:
        raise ValueError("invalid_curated_smiles")
    observed = Chem.MolFromPDBFile(
        str(ligand_pdb), sanitize=False, removeHs=False, proximityBonding=True
    )
    if observed is None:
        raise ValueError("ligand_coordinate_parse_failed")
    if observed.GetNumHeavyAtoms() != template.GetNumHeavyAtoms():
        raise ValueError(
            f"ligand_heavy_atom_count_mismatch:{observed.GetNumHeavyAtoms()}!="
            f"{template.GetNumHeavyAtoms()}"
        )
    prepared = AllChem.AssignBondOrdersFromTemplate(template, observed)
    match = prepared.GetSubstructMatch(template)
    if len(match) != template.GetNumAtoms():
        raise ValueError("ligand_atom_mapping_failed")
    # Coordinates may imply stereo that the curated identity intentionally leaves
    # unspecified (the SB203580 sulfoxide is one example). The curated graph is
    # authoritative for stereo as well as bond orders and formal charges.
    Chem.RemoveStereochemistry(prepared)
    for template_atom, prepared_index in zip(template.GetAtoms(), match):
        prepared.GetAtomWithIdx(prepared_index).SetChiralTag(
            template_atom.GetChiralTag()
        )
    for template_bond in template.GetBonds():
        mapped_bond = prepared.GetBondBetweenAtoms(
            match[template_bond.GetBeginAtomIdx()], match[template_bond.GetEndAtomIdx()]
        )
        mapped_bond.SetStereo(template_bond.GetStereo())
    before = prepared.GetConformer()
    heavy_positions = {
        atom.GetIdx(): tuple(
            float(value) for value in before.GetAtomPosition(atom.GetIdx())
        )
        for atom in prepared.GetAtoms()
        if atom.GetAtomicNum() > 1
    }
    prepared = Chem.AddHs(prepared, addCoords=True)
    for atom in prepared.GetAtoms():
        if atom.GetPDBResidueInfo() is None:
            info = Chem.AtomPDBResidueInfo()
            info.SetResidueName(residue_name[:3])
            info.SetResidueNumber(number)
            info.SetChainId(chain)
            info.SetName(f"H{atom.GetIdx():03}"[-4:])
            atom.SetMonomerInfo(info)
    after = prepared.GetConformer()
    max_displacement = max(
        math.dist(
            coordinates,
            tuple(float(value) for value in after.GetAtomPosition(index)),
        )
        for index, coordinates in heavy_positions.items()
    )
    if max_displacement > 1e-6:
        raise ValueError(f"ligand_heavy_atoms_moved:{max_displacement}")
    canonical = Chem.MolToSmiles(Chem.RemoveHs(prepared), isomericSmiles=True)
    expected = Chem.MolToSmiles(template, isomericSmiles=True)
    if canonical != expected:
        raise ValueError(f"prepared_ligand_identity_mismatch:{canonical}!={expected}")
    mapping = []
    for template_index, source_index in enumerate(match):
        source_atom = prepared.GetAtomWithIdx(source_index)
        payload = _atom_payload(source_atom)
        mapping.append(
            {
                "template_index": template_index,
                "source_index": source_index,
                "source_atom_name": payload["atom_name"],
                "element": source_atom.GetSymbol(),
            }
        )
    return plf.Molecule(prepared), {
        "canonical_isomeric_smiles": canonical,
        "formal_charge": Chem.GetFormalCharge(template),
        "heavy_atom_count": template.GetNumHeavyAtoms(),
        "added_hydrogens": prepared.GetNumAtoms() - prepared.GetNumHeavyAtoms(),
        "max_heavy_atom_displacement_angstrom": max_displacement,
        "atom_mapping": mapping,
    }


def _heavy_coordinates(
    path: Path,
) -> dict[tuple[str, int, str], tuple[float, float, float]]:
    molecule = Chem.MolFromPDBFile(
        str(path), sanitize=False, removeHs=False, proximityBonding=False
    )
    if molecule is None:
        raise ValueError(f"protein_coordinate_parse_failed:{path.name}")
    conformer = molecule.GetConformer()
    result = {}
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() <= 1:
            continue
        info = atom.GetPDBResidueInfo()
        if info is None:
            continue
        key = (
            info.GetChainId().strip(),
            info.GetResidueNumber(),
            info.GetName().strip(),
        )
        result[key] = tuple(
            float(value) for value in conformer.GetAtomPosition(atom.GetIdx())
        )
    return result


def _histidine_states(molecule: Chem.Mol) -> dict[str, str]:
    residues: dict[tuple[str, int], set[str]] = {}
    for atom in molecule.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info and info.GetResidueName().strip() in {"HIS", "HID", "HIE", "HIP"}:
            key = (info.GetChainId().strip(), info.GetResidueNumber())
            residues.setdefault(key, set()).add(info.GetName().strip())
    states = {}
    for (chain, number), names in sorted(residues.items()):
        state = "HIP" if {"HD1", "HE2"} <= names else "HID" if "HD1" in names else "HIE"
        states[f"{chain}:{number}"] = state
    return states


def _prepare_protein(
    protein_pdb: Path, directory: Path, settings: dict
) -> tuple[plf.Molecule, dict]:
    executable = shutil.which("pdb2pqr")
    if executable is None:
        raise ValueError("pdb2pqr_executable_unavailable")
    output_pqr = directory / "protein.pqr"
    output_pdb = directory / "protein-protonated.pdb"
    command = [
        executable,
        f"--ff={settings['force_field']}",
        "--keep-chain",
        "--drop-water",
        f"--with-ph={settings['ph']}",
        f"--pdb-output={output_pdb}",
        str(protein_pdb),
        str(output_pqr),
    ]
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    if completed.returncode:
        raise ValueError(f"pdb2pqr_failed:{completed.stderr[-500:]}")
    original = _heavy_coordinates(protein_pdb)
    protonated = _heavy_coordinates(output_pdb)
    common = sorted(set(original) & set(protonated))
    if not common:
        raise ValueError("protein_heavy_atom_mapping_failed")
    pre_restore_displacement = max(
        math.dist(original[key], protonated[key]) for key in common
    )
    molecule = Chem.MolFromPDBFile(
        str(output_pdb), sanitize=False, removeHs=False, proximityBonding=True
    )
    if molecule is None:
        raise ValueError("prepared_protein_parse_failed")
    editable = Chem.RWMol(molecule)
    removed_spurious_bonds = 0
    for bond in list(editable.GetBonds()):
        first = bond.GetBeginAtom()
        second = bond.GetEndAtom()
        first_info = first.GetPDBResidueInfo()
        second_info = second.GetPDBResidueInfo()
        if first_info is None or second_info is None:
            continue
        first_residue = (first_info.GetChainId(), first_info.GetResidueNumber())
        second_residue = (second_info.GetChainId(), second_info.GetResidueNumber())
        if first_residue == second_residue:
            continue
        names = {first_info.GetName().strip(), second_info.GetName().strip()}
        numbers = {first_info.GetResidueNumber(), second_info.GetResidueNumber()}
        is_peptide = names == {"C", "N"} and max(numbers) - min(numbers) == 1
        is_disulfide = names == {"SG"}
        if not is_peptide and not is_disulfide:
            editable.RemoveBond(first.GetIdx(), second.GetIdx())
            removed_spurious_bonds += 1
    molecule = editable.GetMol()
    try:
        Chem.SanitizeMol(molecule)
    except Exception as exc:
        raise ValueError(f"prepared_protein_parse_failed:{exc}") from exc
    conformer = molecule.GetConformer()
    restored = 0
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() <= 1 or atom.GetPDBResidueInfo() is None:
            continue
        info = atom.GetPDBResidueInfo()
        key = (
            info.GetChainId().strip(),
            info.GetResidueNumber(),
            info.GetName().strip(),
        )
        if key in original:
            conformer.SetAtomPosition(atom.GetIdx(), original[key])
            restored += 1
    restored_coordinates = {
        key: tuple(float(value) for value in conformer.GetAtomPosition(atom.GetIdx()))
        for atom in molecule.GetAtoms()
        if atom.GetAtomicNum() > 1 and atom.GetPDBResidueInfo() is not None
        for info in (atom.GetPDBResidueInfo(),)
        for key in (
            (
                info.GetChainId().strip(),
                info.GetResidueNumber(),
                info.GetName().strip(),
            ),
        )
        if key in original
    }
    max_displacement = max(
        math.dist(original[key], restored_coordinates[key])
        for key in restored_coordinates
    )
    if max_displacement > float(settings["max_heavy_atom_displacement_angstrom"]):
        raise ValueError(f"protein_heavy_atoms_moved:{max_displacement}")
    return plf.Molecule(molecule), {
        "method": "PDB2PQR",
        "force_field": settings["force_field"],
        "ph": float(settings["ph"]),
        "hydrogen_optimization": True,
        "input_heavy_atoms": len(original),
        "mapped_heavy_atoms": len(common),
        "restored_heavy_atoms": restored,
        "added_heavy_atoms": len(set(protonated) - set(original)),
        "removed_spurious_inter_residue_bonds": removed_spurious_bonds,
        "max_pre_restore_displacement_angstrom": pre_restore_displacement,
        "max_heavy_atom_displacement_angstrom": max_displacement,
        "histidine_states": _histidine_states(molecule),
    }


def _geometry(metadata: dict) -> list[dict[str, object]]:
    rows = []
    for name, value in metadata.items():
        if name in {"indices", "parent_indices"} or not isinstance(value, (int, float)):
            continue
        rows.append(
            {
                "name": str(name),
                "value": float(value),
                "unit": "degree" if "angle" in str(name).lower() else "angstrom",
            }
        )
    return sorted(rows, key=lambda row: str(row["name"]))


def _extract_events(
    generated, ligand_molecule: Chem.Mol, protein_molecule: Chem.Mol
) -> list[dict[str, object]]:
    roles = {
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
    events = []
    for (_, protein_residue), interaction_map in generated.items():
        for interaction_name, occurrences in interaction_map.items():
            ligand_role, protein_role = roles[str(interaction_name)]
            normalized = {
                "HBAcceptor": "hb_acceptor",
                "HBDonor": "hb_donor",
                "VdWContact": "vdw_contact",
                "PiStacking": "pi_stacking",
                "CationPi": "cation_pi",
                "PiCation": "pi_cation",
            }.get(str(interaction_name), str(interaction_name).lower())
            for occurrence in occurrences:
                indices = (
                    occurrence.get("parent_indices") or occurrence.get("indices") or {}
                )
                events.append(
                    {
                        "interaction_key": (
                            f"{protein_residue.chain}:{protein_residue.number}:{normalized}"
                        ),
                        "interaction_type": normalized,
                        "ligand_role": ligand_role,
                        "protein_role": protein_role,
                        "ligand_atoms": [
                            _atom_payload(ligand_molecule.GetAtomWithIdx(int(index)))
                            for index in indices.get("ligand", ())
                        ],
                        "protein_atoms": [
                            _atom_payload(protein_molecule.GetAtomWithIdx(int(index)))
                            for index in indices.get("protein", ())
                        ],
                        "geometry": _geometry(occurrence),
                    }
                )
    return events


def main() -> None:
    request = json.load(sys.stdin)
    source = Path(request["structure_path"])
    ligand = request["ligand"]
    with TemporaryDirectory(prefix="cofolder-structure-gated-") as temporary:
        directory = Path(temporary)
        ligand_pdb = directory / "ligand.pdb"
        protein_pdb = directory / "protein.pdb"
        _write_selected_pdb(
            source,
            ligand_pdb,
            ligand_chain=ligand["chain_id"],
            ligand_number=int(ligand["residue_number"]),
            protein_only=False,
        )
        _write_selected_pdb(
            source,
            protein_pdb,
            ligand_chain=ligand["chain_id"],
            ligand_number=int(ligand["residue_number"]),
            protein_only=True,
        )
        ligand_molecule, ligand_report = _prepare_ligand(
            ligand_pdb,
            request["smiles"],
            ligand["residue_name"],
            ligand["chain_id"],
            int(ligand["residue_number"]),
        )
        protein_molecule, protein_report = _prepare_protein(
            protein_pdb, directory, request["preparation"]
        )
        fingerprint = plf.Fingerprint(interactions=request["interactions"])
        generated = fingerprint.generate(
            ligand_molecule, protein_molecule, metadata=True
        )
        events = _extract_events(generated, ligand_molecule, protein_molecule)
        result = {
            "schema_version": 1,
            "preparation": {
                "ligand": ligand_report,
                "protein": protein_report,
                "prolif_version": plf.__version__,
                "rdkit_version": Chem.rdBase.rdkitVersion,
            },
            "events": events,
        }
        json.dump(result, sys.stdout, sort_keys=True)


if __name__ == "__main__":
    main()
