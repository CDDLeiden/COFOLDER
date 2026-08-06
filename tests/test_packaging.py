"""Packaging metadata regression tests."""

from __future__ import annotations

from pathlib import Path
import importlib.util
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


def test_tutorial_extra_includes_marimo():
    pyproject_path = Path(__file__).resolve().parents[1] / "pyproject.toml"
    pyproject = tomllib.loads(pyproject_path.read_text(encoding="utf-8"))
    tutorials = pyproject["project"]["optional-dependencies"]["tutorials"]

    assert any(dep.startswith("marimo") for dep in tutorials)


def test_bias_training_builder_is_an_importable_package_module():
    spec = importlib.util.find_spec(
        "cofolder.modules.analytics.build_bias_training_data"
    )

    assert spec is not None
    assert spec.origin is not None
    assert Path(spec.origin).name == "build_bias_training_data.py"
