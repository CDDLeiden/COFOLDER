"""Runner-aware validation for molecular system YAML inputs."""

from __future__ import annotations

import pickle
import math
from pathlib import Path
from typing import Any

from rdkit import Chem

from cofolder.modules.input.system import SUPPORTED_ENTITY_TYPES, System, SystemChain
from cofolder.modules.runners.contracts import RunnerInputCapabilities


class SystemInputValidationError(ValueError):
    """Raised when a system cannot be faithfully handled by a selected runner."""


_PROTEIN_SIDECHAIN_ATOMS = {
    "A": {"CB"},
    "R": {"CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"},
    "N": {"CB", "CG", "OD1", "ND2"},
    "D": {"CB", "CG", "OD1", "OD2"},
    "C": {"CB", "SG"},
    "Q": {"CB", "CG", "CD", "OE1", "NE2"},
    "E": {"CB", "CG", "CD", "OE1", "OE2"},
    "G": set(),
    "H": {"CB", "CG", "ND1", "CD2", "CE1", "NE2"},
    "I": {"CB", "CG1", "CG2", "CD1"},
    "L": {"CB", "CG", "CD1", "CD2"},
    "K": {"CB", "CG", "CD", "CE", "NZ"},
    "M": {"CB", "CG", "SD", "CE"},
    "F": {"CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "P": {"CB", "CG", "CD"},
    "S": {"CB", "OG"},
    "T": {"CB", "OG1", "CG2"},
    "W": {"CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "Y": {"CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"},
    "V": {"CB", "CG1", "CG2"},
}
_PROTEIN_BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}
_NUCLEIC_BACKBONE_ATOMS = {
    "P",
    "OP1",
    "OP2",
    "OP3",
    "O5'",
    "C5'",
    "C4'",
    "O4'",
    "C3'",
    "O3'",
    "C2'",
    "C1'",
}
_RNA_BASE_ATOMS = {
    "A": {"N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4"},
    "C": {"N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"},
    "G": {"N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4"},
    "U": {"N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"},
}
_DNA_BASE_ATOMS = {
    **_RNA_BASE_ATOMS,
    "T": {"N1", "C2", "O2", "N3", "C4", "O4", "C5", "C7", "C6"},
}


def _root(system_obj: Any) -> dict[str, Any]:
    value = system_obj.system if isinstance(system_obj, System) else system_obj
    if not isinstance(value, dict):
        raise SystemInputValidationError("System YAML root must be a mapping.")
    return value


def _error(runner_name: str, message: str) -> SystemInputValidationError:
    return SystemInputValidationError(
        f"Runner '{runner_name}' cannot accept this system: {message}"
    )


def _validated_chains(
    root: dict[str, Any], runner_name: str, capabilities: RunnerInputCapabilities
) -> list[SystemChain]:
    sequences = root.get("sequences")
    if not isinstance(sequences, list) or not sequences:
        raise _error(runner_name, "'sequences' must be a non-empty list.")

    chains: list[SystemChain] = []
    seen: set[str] = set()
    for index, entry in enumerate(sequences):
        if not isinstance(entry, dict) or len(entry) != 1:
            raise _error(
                runner_name,
                f"sequence entry {index} must contain exactly one entity type.",
            )
        entity_type, payload = next(iter(entry.items()))
        entity_type = str(entity_type).lower()
        if entity_type not in SUPPORTED_ENTITY_TYPES:
            raise _error(
                runner_name,
                f"sequence entry {index} has unknown entity type {entity_type!r}; "
                f"recognized types are {sorted(SUPPORTED_ENTITY_TYPES)}.",
            )
        if entity_type not in capabilities.entity_types:
            raise _error(
                runner_name,
                f"entity type {entity_type!r} is unsupported; supported entity types are "
                f"{sorted(capabilities.entity_types)}.",
            )
        if not isinstance(payload, dict):
            raise _error(runner_name, f"{entity_type} entry {index} must be a mapping.")
        raw_ids = payload.get("id")
        ids = raw_ids if isinstance(raw_ids, list) else [raw_ids]
        if not ids or any(value is None or not str(value).strip() for value in ids):
            raise _error(
                runner_name,
                f"{entity_type} entry {index} requires non-empty chain id(s).",
            )
        if entity_type in {"protein", "dna", "rna"} and not (
            payload.get("sequence") or payload.get("fasta")
        ):
            raise _error(
                runner_name, f"{entity_type} entry {index} requires 'sequence'."
            )
        if entity_type == "ligand" and not any(
            payload.get(key) for key in ("smiles", "ccd", "ccd_codes")
        ):
            raise _error(
                runner_name, f"ligand entry {index} requires 'smiles' or 'ccd'."
            )
        for raw_id in ids:
            chain_id = str(raw_id)
            if chain_id in seen:
                raise _error(
                    runner_name, f"chain id {chain_id!r} is declared more than once."
                )
            seen.add(chain_id)
            chains.append(SystemChain(chain_id, entity_type, payload, index))
    return chains


def _sequence(chain: SystemChain) -> str:
    return str(
        chain.entity_data.get("sequence") or chain.entity_data.get("fasta") or ""
    ).upper()


def _polymer_atom_names(chain: SystemChain, residue_id: int) -> set[str]:
    sequence = _sequence(chain)
    if residue_id < 1 or residue_id > len(sequence):
        return set()
    residue = sequence[residue_id - 1]
    if chain.entity_type == "protein":
        return _PROTEIN_BACKBONE_ATOMS | _PROTEIN_SIDECHAIN_ATOMS.get(residue, set())
    bases = _RNA_BASE_ATOMS if chain.entity_type == "rna" else _DNA_BASE_ATOMS
    atoms = _NUCLEIC_BACKBONE_ATOMS | bases.get(residue, set())
    if chain.entity_type == "rna":
        atoms.add("O2'")
    return atoms


def _smiles_atom_names(smiles: str) -> set[str]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return set()
    mol = Chem.AddHs(mol)
    ranks = Chem.CanonicalRankAtoms(mol)
    names = {
        f"{atom.GetSymbol().upper()}{rank + 1}"
        for atom, rank in zip(mol.GetAtoms(), ranks)
        if atom.GetAtomicNum() != 1
    }
    return names


def _ccd_atom_names(
    chain: SystemChain, residue_id: int, cache_path: str | None
) -> set[str] | None:
    raw_codes = chain.entity_data.get("ccd_codes", chain.entity_data.get("ccd"))
    codes = raw_codes if isinstance(raw_codes, list) else [raw_codes]
    codes = [str(code) for code in codes if code is not None]
    if residue_id < 1 or residue_id > len(codes):
        return set()
    if not cache_path:
        return None
    mol_path = Path(cache_path).expanduser() / "mols" / f"{codes[residue_id - 1]}.pkl"
    if not mol_path.is_file():
        return None
    try:
        with mol_path.open("rb") as handle:
            mol = pickle.load(handle)  # noqa: S301 - runner-owned local cache artifact
    except Exception as exc:
        raise SystemInputValidationError(
            f"Unable to read ligand CCD cache entry {mol_path}: {exc}"
        ) from exc
    names = {atom.GetProp("name") for atom in mol.GetAtoms() if atom.HasProp("name")}
    if not names:
        names = {
            f"{atom.GetSymbol().upper()}{index}"
            for index, atom in enumerate(mol.GetAtoms(), 1)
        }
    return names


def _ligand_atom_names(
    chain: SystemChain, residue_id: int, cache_path: str | None
) -> set[str] | None:
    if (
        chain.entity_data.get("ccd") is not None
        or chain.entity_data.get("ccd_codes") is not None
    ):
        return _ccd_atom_names(chain, residue_id, cache_path)
    if residue_id != 1:
        return set()
    smiles = chain.entity_data.get("smiles")
    return _smiles_atom_names(str(smiles)) if smiles else set()


def _require_chain(
    value: Any, chains: dict[str, SystemChain], runner_name: str, label: str
) -> SystemChain:
    chain_id = str(value)
    if chain_id not in chains:
        raise _error(runner_name, f"{label} references unknown chain {chain_id!r}.")
    return chains[chain_id]


def _require_residue(
    value: Any, chain: SystemChain, runner_name: str, label: str
) -> int:
    if isinstance(value, bool):
        raise _error(
            runner_name, f"{label} residue identifier must be a 1-based integer."
        )
    try:
        residue_id = int(value)
    except (TypeError, ValueError) as exc:
        raise _error(
            runner_name, f"{label} residue identifier must be a 1-based integer."
        ) from exc
    if chain.entity_type != "ligand":
        max_residue = len(_sequence(chain))
    else:
        raw_codes = chain.entity_data.get("ccd_codes", chain.entity_data.get("ccd"))
        max_residue = len(raw_codes) if isinstance(raw_codes, list) else 1
    if residue_id < 1 or residue_id > max_residue:
        raise _error(
            runner_name,
            f"{label} residue {residue_id} is outside chain {chain.chain_id!r} "
            f"(valid range: 1-{max_residue}).",
        )
    return residue_id


def _validate_atom_reference(
    raw: Any,
    chains: dict[str, SystemChain],
    runner_name: str,
    label: str,
    cache_path: str | None,
    check_atom_names: bool,
) -> None:
    if not isinstance(raw, (list, tuple)) or len(raw) != 3:
        raise _error(runner_name, f"{label} must be [chain_id, residue_id, atom_name].")
    chain = _require_chain(raw[0], chains, runner_name, label)
    residue_id = _require_residue(raw[1], chain, runner_name, label)
    atom_name = str(raw[2]).strip()
    if not atom_name:
        raise _error(runner_name, f"{label} atom name must be non-empty.")
    if not check_atom_names:
        return
    names = (
        _ligand_atom_names(chain, residue_id, cache_path)
        if chain.entity_type == "ligand"
        else _polymer_atom_names(chain, residue_id)
    )
    if names is None:
        raise _error(
            runner_name,
            f"{label} atom {atom_name!r} cannot be verified because the CCD for chain "
            f"{chain.chain_id!r} is absent from the configured cache.",
        )
    if atom_name not in names:
        available = ", ".join(sorted(names)) or "<none>"
        raise _error(
            runner_name,
            f"{label} atom {atom_name!r} does not exist at chain {chain.chain_id!r}, "
            f"residue {residue_id}; available atoms: {available}.",
        )


def _validate_token_reference(
    raw: Any,
    chains: dict[str, SystemChain],
    runner_name: str,
    label: str,
    cache_path: str | None,
    check_atom_names: bool,
) -> None:
    if not isinstance(raw, (list, tuple)) or len(raw) != 2:
        raise _error(
            runner_name, f"{label} must be [chain_id, residue_id_or_atom_name]."
        )
    chain = _require_chain(raw[0], chains, runner_name, label)
    if chain.entity_type != "ligand":
        _require_residue(raw[1], chain, runner_name, label)
        return
    atom_name = str(raw[1]).strip()
    if not atom_name:
        raise _error(runner_name, f"{label} ligand atom name must be non-empty.")
    if check_atom_names:
        names = _ligand_atom_names(chain, 1, cache_path)
        if names is None:
            raise _error(
                runner_name,
                f"{label} atom {atom_name!r} cannot be verified because the CCD for chain "
                f"{chain.chain_id!r} is absent from the configured cache.",
            )
        if atom_name not in names:
            raise _error(
                runner_name,
                f"{label} atom {atom_name!r} does not exist on ligand chain {chain.chain_id!r}.",
            )


def _validate_distance(payload: dict[str, Any], runner_name: str, label: str) -> None:
    if "max_distance" not in payload:
        return
    try:
        distance = float(payload["max_distance"])
    except (TypeError, ValueError) as exc:
        raise _error(
            runner_name, f"{label} max_distance must be a positive number."
        ) from exc
    if not math.isfinite(distance) or distance <= 0:
        raise _error(runner_name, f"{label} max_distance must be a positive number.")


def validate_system_input(
    system_obj: Any,
    *,
    runner_name: str,
    capabilities: RunnerInputCapabilities,
    cache_path: str | None = None,
    check_atom_names: bool = True,
) -> None:
    """Validate entity and constraint data before a backend command is launched."""

    root = _root(system_obj)
    chain_list = _validated_chains(root, runner_name, capabilities)
    chains = {chain.chain_id: chain for chain in chain_list}
    constraints = root.get("constraints", [])
    if constraints is None:
        constraints = []
    if not isinstance(constraints, list):
        raise _error(runner_name, "'constraints' must be a list.")

    pocket_count = 0
    for index, item in enumerate(constraints):
        prefix = f"constraint {index}"
        if not isinstance(item, dict) or len(item) != 1:
            raise _error(
                runner_name, f"{prefix} must contain exactly one constraint type."
            )
        constraint_type, payload = next(iter(item.items()))
        constraint_type = str(constraint_type).lower()
        if constraint_type not in capabilities.constraint_types:
            raise _error(
                runner_name,
                f"{prefix} type {constraint_type!r} is unsupported; supported constraint types are "
                f"{sorted(capabilities.constraint_types) or ['none']}.",
            )
        if not isinstance(payload, dict):
            raise _error(runner_name, f"{prefix} payload must be a mapping.")
        if "force" in payload and not capabilities.supports_constraint_force:
            raise _error(
                runner_name,
                f"{prefix} uses 'force', which is unsupported by this runner.",
            )

        if constraint_type == "bond":
            if set(payload) != {"atom1", "atom2"}:
                raise _error(
                    runner_name, f"{prefix} bond requires only 'atom1' and 'atom2'."
                )
            _validate_atom_reference(
                payload["atom1"],
                chains,
                runner_name,
                f"{prefix}.bond.atom1",
                cache_path,
                check_atom_names,
            )
            _validate_atom_reference(
                payload["atom2"],
                chains,
                runner_name,
                f"{prefix}.bond.atom2",
                cache_path,
                check_atom_names,
            )
        elif constraint_type == "pocket":
            pocket_count += 1
            allowed = {"binder", "contacts", "max_distance", "force"}
            if (
                "binder" not in payload
                or "contacts" not in payload
                or set(payload) - allowed
            ):
                raise _error(
                    runner_name, f"{prefix} pocket requires 'binder' and 'contacts'."
                )
            binder = _require_chain(
                payload["binder"], chains, runner_name, f"{prefix}.pocket.binder"
            )
            if binder.entity_type != "ligand":
                raise _error(
                    runner_name,
                    f"{prefix}.pocket.binder must reference a ligand chain.",
                )
            contacts = payload["contacts"]
            if not isinstance(contacts, list) or not contacts:
                raise _error(
                    runner_name, f"{prefix}.pocket.contacts must be a non-empty list."
                )
            for contact_index, contact in enumerate(contacts):
                _validate_token_reference(
                    contact,
                    chains,
                    runner_name,
                    f"{prefix}.pocket.contacts[{contact_index}]",
                    cache_path,
                    check_atom_names,
                )
                if capabilities.pocket_contacts_must_be_polymers:
                    contact_chain = _require_chain(
                        contact[0],
                        chains,
                        runner_name,
                        f"{prefix}.pocket.contacts[{contact_index}]",
                    )
                    if contact_chain.entity_type == "ligand":
                        raise _error(
                            runner_name,
                            f"{prefix}.pocket.contacts[{contact_index}] must reference a polymer chain.",
                        )
            _validate_distance(payload, runner_name, f"{prefix}.pocket")
            required = capabilities.required_pocket_distance
            if (
                required is not None
                and float(payload.get("max_distance", required)) != required
            ):
                raise _error(
                    runner_name,
                    f"{prefix}.pocket max_distance must be {required:g} for this runner.",
                )
        elif constraint_type == "contact":
            allowed = {"token1", "token2", "max_distance", "force"}
            if (
                "token1" not in payload
                or "token2" not in payload
                or set(payload) - allowed
            ):
                raise _error(
                    runner_name, f"{prefix} contact requires 'token1' and 'token2'."
                )
            _validate_token_reference(
                payload["token1"],
                chains,
                runner_name,
                f"{prefix}.contact.token1",
                cache_path,
                check_atom_names,
            )
            _validate_token_reference(
                payload["token2"],
                chains,
                runner_name,
                f"{prefix}.contact.token2",
                cache_path,
                check_atom_names,
            )
            _validate_distance(payload, runner_name, f"{prefix}.contact")

    maximum = capabilities.max_pocket_constraints
    if maximum is not None and pocket_count > maximum:
        raise _error(
            runner_name,
            f"at most {maximum} pocket constraint(s) are supported, found {pocket_count}.",
        )
