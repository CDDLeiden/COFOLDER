"""Workflow-independent ligand validation, selection, and preparation."""

from __future__ import annotations

import copy
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor

from cofolder.modules.input.config import LigandSelectionError, LigandValidationError
from cofolder.modules.input.system import System, iter_system_chains


@dataclass(frozen=True, slots=True)
class LigandSourceIdentity:
    entity_id: str
    chain_ids: tuple[str, ...]
    source_record_id: str | None = None


@dataclass(frozen=True, slots=True)
class NormalizedLigand:
    source: LigandSourceIdentity
    source_smiles: str
    canonical_smiles: str
    molecule: Chem.Mol


@dataclass(frozen=True, slots=True)
class LigandPreparationCapabilities:
    native_smiles: bool
    conformer_modes: frozenset[Literal["2D", "3D", "sdf"]]


@dataclass(frozen=True, slots=True)
class PreparedLigand:
    ligand: NormalizedLigand
    mode: Literal["native", "2D", "3D", "sdf"]
    molecule: Chem.Mol
    sdf_path: Path | None = None


@dataclass(frozen=True, slots=True)
class LigandTarget:
    entity_id: str
    sequence_index: int
    chain_ids: tuple[str, ...]


def validate_smiles(
    smiles: object,
    *,
    source: LigandSourceIdentity,
    source_path: Path | None = None,
    field_path: tuple[str | int, ...] = (),
) -> NormalizedLigand:
    text = str(smiles).strip() if smiles is not None else ""
    try:
        # Disable RDKit's implicit sanitization so parsing and sanitization each
        # happen exactly once at this shared boundary.
        mol = Chem.MolFromSmiles(text, sanitize=False) if text else None
        if mol is None:
            raise ValueError("RDKit could not parse the SMILES string")
    except Exception as exc:
        raise LigandValidationError(
            f"Invalid SMILES for ligand {source.entity_id}: {exc}",
            source_path=source_path,
            field_path=field_path,
            entity_id=source.entity_id,
            chain_id=source.chain_ids[0] if source.chain_ids else None,
            source_record_id=source.source_record_id,
        ) from exc
    return validate_molecule(
        mol,
        source=source,
        source_path=source_path,
        field_path=field_path,
        source_smiles=text,
    )


def validate_molecule(
    molecule: Chem.Mol,
    *,
    source: LigandSourceIdentity,
    source_path: Path | None = None,
    field_path: tuple[str | int, ...] = (),
    line: int | None = None,
    source_smiles: str | None = None,
) -> NormalizedLigand:
    """Sanitize and normalize an already-parsed molecule at the shared boundary."""
    mol = Chem.Mol(molecule)
    try:
        Chem.SanitizeMol(mol)
        canonical_smiles = Chem.MolToSmiles(
            mol, canonical=True, isomericSmiles=True
        )
    except Exception as exc:
        raise LigandValidationError(
            f"Invalid molecule for ligand {source.entity_id}: {exc}",
            source_path=source_path,
            field_path=field_path,
            line=line,
            entity_id=source.entity_id,
            chain_id=source.chain_ids[0] if source.chain_ids else None,
            source_record_id=source.source_record_id,
        ) from exc
    return NormalizedLigand(
        source=source,
        source_smiles=source_smiles if source_smiles is not None else canonical_smiles,
        canonical_smiles=canonical_smiles,
        molecule=mol,
    )


def prepare_ligand(
    ligand: NormalizedLigand,
    *,
    mode: Literal["native", "2D", "3D", "sdf"],
    sdf_path: Path | None,
    capabilities: LigandPreparationCapabilities,
    work_dir: Path,
) -> PreparedLigand:
    if mode == "native":
        if not capabilities.native_smiles:
            raise LigandValidationError(
                "The selected runner does not support native SMILES input.",
                entity_id=ligand.source.entity_id,
                chain_id=ligand.source.chain_ids[0] if ligand.source.chain_ids else None,
                source_record_id=ligand.source.source_record_id,
            )
        return PreparedLigand(ligand, mode, Chem.Mol(ligand.molecule))
    if mode not in capabilities.conformer_modes:
        raise LigandValidationError(
            f"The selected runner does not support {mode} ligand conformers.",
            entity_id=ligand.source.entity_id,
            chain_id=ligand.source.chain_ids[0] if ligand.source.chain_ids else None,
            source_record_id=ligand.source.source_record_id,
        )
    if mode == "sdf":
        if sdf_path is None or not sdf_path.is_file():
            raise LigandValidationError(
                "SDF conformer mode requires a readable SDF file.",
                source_path=sdf_path,
                entity_id=ligand.source.entity_id,
                chain_id=ligand.source.chain_ids[0]
                if ligand.source.chain_ids
                else None,
                source_record_id=ligand.source.source_record_id,
            )
        try:
            supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
            mol = next((item for item in supplier if item is not None), None)
        except Exception as exc:
            raise LigandValidationError(
                f"Unable to read SDF conformer file {sdf_path}: {exc}",
                source_path=sdf_path,
                entity_id=ligand.source.entity_id,
                chain_id=ligand.source.chain_ids[0]
                if ligand.source.chain_ids
                else None,
                source_record_id=ligand.source.source_record_id,
            ) from exc
        if mol is None:
            raise LigandValidationError(
                f"No valid molecule found in SDF: {sdf_path}",
                source_path=sdf_path,
                entity_id=ligand.source.entity_id,
                chain_id=ligand.source.chain_ids[0]
                if ligand.source.chain_ids
                else None,
                source_record_id=ligand.source.source_record_id,
            )
        sdf_smiles = Chem.MolToSmiles(
            Chem.RemoveHs(mol), canonical=True, isomericSmiles=True
        )
        if sdf_smiles != ligand.canonical_smiles:
            raise LigandValidationError(
                f"SDF molecule does not match ligand {ligand.source.entity_id}.",
                source_path=sdf_path,
                entity_id=ligand.source.entity_id,
                chain_id=ligand.source.chain_ids[0]
                if ligand.source.chain_ids
                else None,
                source_record_id=ligand.source.source_record_id,
            )
        return PreparedLigand(ligand, mode, mol, sdf_path)
    mol = Chem.AddHs(Chem.Mol(ligand.molecule))
    if mode == "2D":
        rdDepictor.Compute2DCoords(mol)
    else:
        if AllChem.EmbedMolecule(mol, AllChem.ETKDGv3()) != 0:
            raise LigandValidationError(
                f"Unable to generate a 3D conformer for ligand {ligand.source.entity_id}."
            )
        AllChem.UFFOptimizeMolecule(mol)
    work_dir.mkdir(parents=True, exist_ok=True)
    destination = work_dir / f"{ligand.source.entity_id.replace(':', '_')}.sdf"
    writer = Chem.SDWriter(str(destination))
    writer.write(mol)
    writer.close()
    return PreparedLigand(ligand, mode, mol, destination)


def resolve_ligand_target(system_obj: System | dict[str, Any], chain_id: str) -> LigandTarget:
    matches = [
        chain
        for chain in iter_system_chains(system_obj)
        if chain.chain_id == str(chain_id).strip()
    ]
    if not matches:
        raise LigandSelectionError(f"Ligand selector references unknown chain {chain_id!r}.")
    if any(chain.entity_type != "ligand" for chain in matches):
        raise LigandSelectionError(f"Chain {chain_id!r} is not a ligand entity.")
    indices = {chain.sequence_index for chain in matches}
    if len(indices) != 1:
        raise LigandSelectionError(
            f"Ligand selector {chain_id!r} resolves to multiple ligand entities."
        )
    index = next(iter(indices))
    chain_ids = tuple(
        chain.chain_id
        for chain in iter_system_chains(system_obj)
        if chain.sequence_index == index
    )
    return LigandTarget(f"entity:{index}", index, chain_ids)


def replace_ligand_smiles(
    system_obj: System | dict[str, Any],
    *,
    target: LigandTarget,
    ligand: NormalizedLigand,
) -> System:
    source = system_obj.system if isinstance(system_obj, System) else system_obj
    result = copy.deepcopy(source)
    sequences = result.get("sequences")
    try:
        payload = sequences[target.sequence_index]["ligand"]
    except (TypeError, KeyError, IndexError) as exc:
        raise LigandSelectionError(
            f"Resolved ligand target {target.entity_id} no longer exists."
        ) from exc
    before = copy.deepcopy(result)
    payload.pop("ccd", None)
    payload.pop("ccd_codes", None)
    payload["smiles"] = ligand.canonical_smiles

    before_payload = before["sequences"][target.sequence_index]["ligand"]
    after_payload = result["sequences"][target.sequence_index]["ligand"]
    for key in ("smiles", "ccd", "ccd_codes"):
        before_payload.pop(key, None)
        after_payload.pop(key, None)
    if before != result:
        raise LigandSelectionError(
            "Ligand replacement changed fields outside the selected ligand chemistry."
        )
    # Restore the intended chemistry after the comparison projection removed it.
    result["sequences"][target.sequence_index]["ligand"]["smiles"] = ligand.canonical_smiles
    return System(system=result)
