"""Regression checks for the marimo tutorial surface."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
from unittest.mock import patch

from rdkit import Chem
import yaml

from cofolder import cli
from cofolder.modules.input.compound_library import CompoundMember, load_compound_library
from cofolder.modules.input.ligand import LigandTarget

REPO_ROOT = Path(__file__).resolve().parents[1]
TUTORIALS_DIR = REPO_ROOT / "tutorials"
EXAMPLES_DIR = REPO_ROOT / "examples"
sys.path.insert(0, str(TUTORIALS_DIR))

from _workflow_examples import (  # noqa: E402
    build_bias_command,
    build_oracle_command,
    build_screen_command,
    build_validate_command,
    prepare_bias_example,
    with_preflight,
)


def _load_module(path: Path) -> None:
    spec = importlib.util.spec_from_file_location(path.stem, path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)


def test_expected_marimo_tutorial_files_exist() -> None:
    expected = {
        "bias.py",
        "boltz_system_inputs.py",
        "oracle.py",
        "ligand_handling.py",
        "openfold3_system_inputs.py",
        "runners.py",
        "screen.py",
        "validate.py",
    }
    actual = {path.name for path in TUTORIALS_DIR.glob("*.py") if not path.name.startswith("_")}
    assert expected.issubset(actual)


def test_legacy_ipynb_tutorial_files_are_gone() -> None:
    assert list(TUTORIALS_DIR.glob("*.ipynb")) == []


def test_marimo_tutorial_modules_import() -> None:
    for path in sorted(TUTORIALS_DIR.glob("*.py")):
        if path.name.startswith("_"):
            continue
        _load_module(path)


def test_offline_bias_tutorial_executes_end_to_end(tmp_path: Path) -> None:
    inputs = prepare_bias_example(tmp_path / "inputs")
    output = tmp_path / "bias output"
    command = build_bias_command(inputs, output)

    assert cli.main(command[1:]) == 0

    manifest = json.loads((output / "results" / "manifest.json").read_text())
    assert manifest["status"] == "success"
    assert (
        output / "results" / "bias_train" / "reference_landscape_summary.csv"
    ).is_file()


def test_runner_backed_tutorial_commands_pass_no_service_preflight(
    tmp_path: Path,
) -> None:
    commands = [
        build_validate_command(EXAMPLES_DIR, tmp_path / "validate"),
        build_validate_command(
            EXAMPLES_DIR,
            tmp_path / "validate sdf",
            system_name="system_screen.yaml",
            conformer_sdf=True,
        ),
        build_validate_command(
            EXAMPLES_DIR,
            tmp_path / "validate ccd",
            system_name="system_custom_ccd.yaml",
        ),
        build_screen_command(EXAMPLES_DIR, tmp_path / "screen csv"),
        build_screen_command(
            EXAMPLES_DIR, tmp_path / "screen sdf", library_name="ethanol.sdf"
        ),
        build_screen_command(
            EXAMPLES_DIR, tmp_path / "screen mol", library_name="ethanol.mol"
        ),
        build_oracle_command(EXAMPLES_DIR, tmp_path / "oracle"),
    ]

    with patch(
        "cofolder.modules.runners.boltz2_runner.Boltz2Runner.check_availability",
        return_value=(True, None),
    ):
        for command in commands:
            assert cli.main(with_preflight(command)[1:]) == 0
            output = Path(command[command.index("-w") + 1])
            assert not output.exists()


def test_packaged_ligand_examples_are_valid_and_equivalent() -> None:
    system = yaml.safe_load(
        (EXAMPLES_DIR / "system_custom_ccd.yaml").read_text(encoding="utf-8")
    )
    assert system["sequences"][1]["ligand"]["ccd"] == "ET5"

    sdf_supplier = Chem.SDMolSupplier(str(EXAMPLES_DIR / "ethanol.sdf"))
    sdf_molecule = next(molecule for molecule in sdf_supplier if molecule is not None)
    mol_molecule = Chem.MolFromMolFile(str(EXAMPLES_DIR / "ethanol.mol"))
    assert mol_molecule is not None
    assert sdf_molecule.GetProp("ID") == "ET5"
    assert {
        Chem.MolToSmiles(sdf_molecule, isomericSmiles=True),
        Chem.MolToSmiles(mol_molecule, isomericSmiles=True),
    } == {"CCO"}

    target = LigandTarget("entity:1", 1, ("B",))
    for name in ("ethanol.sdf", "ethanol.mol"):
        library = load_compound_library(EXAMPLES_DIR / name, target=target)
        assert len(library.outcomes) == 1
        assert isinstance(library.outcomes[0], CompoundMember)
        assert library.outcomes[0].ligand.canonical_smiles == "CCO"


def test_public_tutorials_do_not_reference_removed_or_fictional_apis() -> None:
    content = "\n".join(
        path.read_text(encoding="utf-8")
        for path in sorted((REPO_ROOT / "docs" / "tutorials").glob("*.md"))
    )
    forbidden = (
        "run_prediction(",
        "calculate_ifp(",
        'load_predictions("',
        'load_reference("',
        "system_df[[",
    )
    assert all(symbol not in content for symbol in forbidden)
