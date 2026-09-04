"""Regression checks for the marimo tutorial surface."""

from __future__ import annotations

import importlib.util
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
TUTORIALS_DIR = REPO_ROOT / "tutorials"


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
