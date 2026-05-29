"""Packaging metadata regression tests."""

from __future__ import annotations

from pathlib import Path
import tomllib


def test_default_dependencies_include_core_runtime_packages():
    pyproject_path = Path(__file__).resolve().parents[1] / "pyproject.toml"
    pyproject = tomllib.loads(pyproject_path.read_text(encoding="utf-8"))
    dependencies = pyproject["project"]["dependencies"]

    assert "pandas" in dependencies
    assert "pyyaml" in dependencies
    assert "biopython" in dependencies
    assert "scipy" in dependencies
    assert "scikit-learn" in dependencies
