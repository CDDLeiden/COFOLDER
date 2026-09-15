"""Characterization tests for Boltz ligand-preparation orchestration."""

from __future__ import annotations

import logging
import subprocess
import sys
from pathlib import Path

from rdkit import Chem

from cofolder.modules.entities import ligand
from cofolder.modules.input.command import Command
from cofolder.modules.input.system import System
from cofolder.modules.runners import _ligand_preparation


def _system() -> System:
    return System(
        system={
            "sequences": [
                {"protein": {"id": "A", "sequence": "AC"}},
                {"ligand": {"id": "L1", "smiles": "CCO"}},
                {"ligand": {"id": "L2", "smiles": "CCN"}},
            ]
        }
    )


def test_prepare_ligand_conformers_uses_explicit_cache_and_updates_each_ligand(
    monkeypatch, temp_dir
):
    system = _system()
    cache_path = temp_dir / "cache with spaces"
    generated: list[Path] = []
    converted: list[tuple[str, str, Path]] = []

    def fake_smiles_to_sdf(*, data, output_sdf_path):
        path = Path(output_sdf_path)
        path.write_text(data, encoding="utf-8")
        generated.append(path)

    monkeypatch.setattr(ligand, "smiles_to_sdf", fake_smiles_to_sdf)
    monkeypatch.setattr(ligand, "generate_3d_conformers", lambda path: None)
    monkeypatch.setattr(
        _ligand_preparation,
        "_allocate_ccd_resname",
        lambda smiles, boltz_path: {"CCO": "ETH", "CCN": "ETN"}[smiles],
    )
    monkeypatch.setattr(
        _ligand_preparation.read,
        "read_sdf",
        lambda path: [Chem.MolFromSmiles(Path(path).read_text(encoding="utf-8"))],
    )
    monkeypatch.setattr(
        ligand,
        "mol_to_ccd",
        lambda resname, mol, boltz_path: converted.append(
            (resname, Chem.MolToSmiles(mol), Path(boltz_path))
        ),
    )

    result = _ligand_preparation.prepare_ligand_conformers(
        system,
        cache_path=cache_path,
        wrk_dir=temp_dir,
        conformers="3D",
        logger=logging.getLogger(__name__),
    )

    assert result == {"L1": "ETH", "L2": "ETN"}
    assert converted == [
        ("ETH", "CCO", cache_path),
        ("ETN", "CCN", cache_path),
    ]
    assert [entry["ligand"] for entry in system.system["sequences"][1:]] == [
        {"id": "L1", "ccd": "ETH"},
        {"id": "L2", "ccd": "ETN"},
    ]
    assert generated
    assert all(not path.exists() for path in generated)


def test_prepare_ligand_conformers_preserves_supplied_sdf_order_and_failures(
    monkeypatch, temp_dir
):
    system = _system()
    sdf_path = temp_dir / "ligands.sdf"
    sdf_path.write_text("fixture", encoding="utf-8")
    mols = [Chem.MolFromSmiles("CCO"), Chem.MolFromSmiles("CCN")]
    attempted: list[tuple[str, str]] = []

    monkeypatch.setattr(_ligand_preparation.read, "read_sdf", lambda path: mols)
    monkeypatch.setattr(
        _ligand_preparation,
        "_allocate_ccd_resname",
        lambda smiles, boltz_path: {"CCO": "ETH", "CCN": "ETN"}[smiles],
    )

    def fail_first(resname, mol, boltz_path):
        attempted.append((resname, Chem.MolToSmiles(mol)))
        if resname == "ETH":
            raise RuntimeError("conversion failed")

    monkeypatch.setattr(ligand, "mol_to_ccd", fail_first)

    result = _ligand_preparation.prepare_ligand_conformers(
        system,
        cache_path=temp_dir / "cache",
        wrk_dir=temp_dir,
        conformers="sdf",
        sdf_file=sdf_path,
        logger=logging.getLogger(__name__),
    )

    assert result == {"L1": "ETH", "L2": "ETN"}
    assert attempted == [("ETH", "CCO"), ("ETN", "CCN")]
    assert sdf_path.exists()


def test_legacy_handle_conformers_facade_preserves_command_cache(monkeypatch, temp_dir):
    captured = {}

    def fake_prepare(*args, **kwargs):
        captured.update(kwargs)
        return {"L": "ABCDE"}

    monkeypatch.setattr(_ligand_preparation, "prepare_ligand_conformers", fake_prepare)
    options = Command(options={"options": [{"cache": str(temp_dir / "legacy")}]})

    result = ligand.handle_conformers(
        object(), options, str(temp_dir), conformers="2D", logger=logging.getLogger(__name__)
    )

    assert result == {"L": "ABCDE"}
    assert captured["cache_path"] == temp_dir / "legacy"
    assert captured["conformers"] == "2D"


def test_importing_extracted_boundaries_has_no_execution_side_effects(temp_dir):
    code = """
import sys
def deny(event, args):
    if event in {'subprocess.Popen', 'socket.connect'}:
        raise RuntimeError(event)
sys.addaudithook(deny)
import cofolder.modules.runners._ligand_preparation
import cofolder.modules.analytics.aggregation
"""
    before = tuple(temp_dir.iterdir())
    completed = subprocess.run(
        [sys.executable, "-c", code],
        cwd=temp_dir,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert tuple(temp_dir.iterdir()) == before
