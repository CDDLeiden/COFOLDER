"""Prepare paired custom bias references from explicitly selected structures."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
import hashlib
import json
from pathlib import Path
import shutil
import tempfile

import gemmi
import pandas as pd
import yaml
from rdkit import Chem


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _selector_chain(value) -> str:
    return str(value.get("chain_id") if isinstance(value, dict) else value).strip()


def _embedded_smiles(path: Path, residue_name: str) -> str | None:
    if path.suffix.lower() not in {".cif", ".mmcif"}:
        return None
    document = gemmi.cif.read(str(path))
    target = residue_name.strip().upper()
    for block in document:
        table = block.find(
            "_pdbx_chem_comp_descriptor.", ["comp_id", "type", "descriptor"]
        )
        choices: list[tuple[int, str]] = []
        for row in table:
            if str(row[0]).strip().upper() != target or "SMILES" not in str(row[1]).upper():
                continue
            choices.append((0 if "CANONICAL" in str(row[1]).upper() else 1, str(row[2])))
        if choices:
            return sorted(choices)[0][1]
    return None


def _component_smiles(path: Path | None, residue_name: str) -> str | None:
    if path is None:
        return None
    from cofolder.modules.analytics.bias import _components_cif_smiles_index, _path_cache_token

    return _components_cif_smiles_index(_path_cache_token(path)).get(
        residue_name.strip().upper()
    )


def prepare(manifest_path: Path, output_root: Path, *, overwrite: bool = False) -> Path:
    payload = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or payload.get("schema_version") != 1:
        raise ValueError("Custom-complex manifest must use schema_version: 1.")
    dataset = str(payload.get("dataset_name") or "").strip()
    complexes = payload.get("complexes")
    if not dataset or not isinstance(complexes, list) or not complexes:
        raise ValueError("Custom-complex manifest requires dataset_name and complexes.")
    components_path = payload.get("components_cif")
    components = (
        (manifest_path.parent / components_path).resolve()
        if components_path and not Path(components_path).is_absolute()
        else Path(components_path).expanduser().resolve() if components_path else None
    )
    if components is not None and not components.is_file():
        raise ValueError(f"components.cif not found: {components}")
    seen: set[str] = set()
    protein_rows: list[dict[str, object]] = []
    ligand_rows: list[dict[str, object]] = []
    diagnostics: list[dict[str, object]] = []
    for entry in complexes:
        if not isinstance(entry, dict):
            raise ValueError("Each custom complex must be a mapping.")
        complex_id = str(entry.get("complex_id") or "").strip()
        if not complex_id or complex_id in seen:
            raise ValueError(f"Duplicate or empty complex_id: {complex_id!r}")
        seen.add(complex_id)
        raw_path = Path(str(entry.get("structure_path") or ""))
        structure_path = (
            raw_path.expanduser().resolve()
            if raw_path.is_absolute()
            else (manifest_path.parent / raw_path).resolve()
        )
        if not structure_path.is_file() or structure_path.suffix.lower() not in {".pdb", ".cif", ".mmcif"}:
            raise ValueError(f"Structure for {complex_id} is missing or unsupported: {structure_path}")
        structure = gemmi.read_structure(str(structure_path))
        if len(structure) == 0:
            raise ValueError(f"Structure has no model: {structure_path}")
        model = structure[0]
        proteins = entry.get("proteins")
        ligands = entry.get("ligands")
        if not isinstance(proteins, list) or not proteins or not isinstance(ligands, list) or not ligands:
            raise ValueError(f"Complex {complex_id} requires protein and ligand selectors.")
        for selector in proteins:
            chain_id = _selector_chain(selector)
            chain = model.find_chain(chain_id)
            if not chain:
                raise ValueError(f"Protein chain {chain_id!r} not found in {complex_id}.")
            residues = [res.name for res in chain.get_polymer()]
            sequence = gemmi.one_letter_code(residues).replace(" ", "").replace("?", "X")
            if not sequence:
                raise ValueError(f"Protein chain {chain_id!r} has no polymer sequence in {complex_id}.")
            protein_rows.append(
                {
                    "complex_id": complex_id,
                    "pdb_id": complex_id,
                    "sequence": sequence,
                    "source": "custom",
                    "dataset_name": dataset,
                    "source_structure_path": str(structure_path),
                    "source_reference_path": str(manifest_path.resolve()),
                    "source_component_id": chain_id,
                }
            )
        ligand_keys: set[tuple[str, int, str, str]] = set()
        for selector in ligands:
            if not isinstance(selector, dict):
                raise ValueError(f"Ligands in {complex_id} require explicit mapping selectors.")
            chain_id = str(selector.get("chain_id") or "").strip()
            residue_name = str(selector.get("residue_name") or "").strip().upper()
            residue_number = int(selector.get("residue_number"))
            insertion_code = str(selector.get("insertion_code") or "").strip()
            key = (chain_id, residue_number, insertion_code, residue_name)
            if key in ligand_keys:
                raise ValueError(f"Duplicate ligand selector in {complex_id}: {key}")
            ligand_keys.add(key)
            chain = model.find_chain(chain_id)
            matches = [
                residue
                for residue in chain
                if residue.name.upper() == residue_name
                and residue.seqid.num == residue_number
                and str(residue.seqid.icode).strip() == insertion_code
            ] if chain else []
            if len(matches) != 1:
                raise ValueError(f"Ligand selector {key} resolved {len(matches)} residues in {complex_id}.")
            smiles = str(selector.get("smiles") or "").strip()
            smiles = smiles or _embedded_smiles(structure_path, residue_name) or _component_smiles(components, residue_name) or ""
            molecule = Chem.MolFromSmiles(smiles) if smiles else None
            if molecule is None:
                raise ValueError(
                    f"Ligand {key} in {complex_id} needs explicit SMILES or resolvable CCD chemistry."
                )
            ligand_rows.append(
                {
                    "complex_id": complex_id,
                    "pdb_id": complex_id,
                    "ligand_id": residue_name,
                    "smiles": Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True),
                    "source": "custom",
                    "dataset_name": dataset,
                    "source_structure_path": str(structure_path),
                    "source_reference_path": str(manifest_path.resolve()),
                    "source_component_id": f"{chain_id}:{residue_number}{insertion_code}:{residue_name}",
                }
            )
        diagnostics.append(
            {"complex_id": complex_id, "proteins": len(proteins), "ligands": len(ligands), "status": "prepared"}
        )

    destination = output_root.expanduser().resolve()
    if destination.exists() and not overwrite:
        raise FileExistsError(f"Custom bias reference output already exists: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".custom-bias-", dir=destination.parent) as raw:
        staged = Path(raw) / destination.name
        staged.mkdir()
        pd.DataFrame(protein_rows).to_csv(staged / "custom_protein_references.csv", index=False)
        pd.DataFrame(ligand_rows).to_csv(staged / "custom_ligand_references.csv", index=False)
        (staged / "preparation_report.json").write_text(json.dumps(diagnostics, indent=2), encoding="utf-8")
        files = {
            path.name: _sha(path) for path in staged.iterdir() if path.is_file()
        }
        (staged / "manifest.json").write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "kind": "custom_complexes",
                    "dataset_name": dataset,
                    "complex_ids": sorted(seen),
                    "source_manifest_sha256": _sha(manifest_path),
                    "files": files,
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        if destination.exists():
            backup = destination.with_name(f".{destination.name}.previous")
            if backup.exists():
                shutil.rmtree(backup)
            destination.replace(backup)
            try:
                staged.replace(destination)
            except Exception:
                backup.replace(destination)
                raise
            shutil.rmtree(backup)
        else:
            staged.replace(destination)
    return destination


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="cofolder-tools prepare-bias-custom-complexes")
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--output_root", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    destination = prepare(args.manifest.resolve(), args.output_root, overwrite=args.overwrite)
    print(f"[done] custom_bias_reference_path={destination}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
