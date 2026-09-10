"""Packaging metadata regression tests."""

from __future__ import annotations

from pathlib import Path
import importlib.util
import runpy
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


def test_unsupported_ui_is_absent_from_install_metadata():
    repository = Path(__file__).resolve().parents[1]
    pyproject = tomllib.loads(
        (repository / "pyproject.toml").read_text(encoding="utf-8")
    )

    optional_dependencies = pyproject["project"]["optional-dependencies"]
    all_dependencies = list(pyproject["project"]["dependencies"])
    for dependencies in optional_dependencies.values():
        all_dependencies.extend(dependencies)

    assert "ui" not in optional_dependencies
    assert not any(
        dependency.lower().startswith("streamlit")
        for dependency in all_dependencies
    )

    package_discovery = pyproject["tool"]["setuptools"]["packages"]["find"]
    assert set(package_discovery["exclude"]) >= {"cofolder.ui", "cofolder.ui.*"}
    assert pyproject["tool"]["setuptools"]["include-package-data"] is False

    scripts = pyproject["project"]["scripts"]
    assert scripts["cofolder"] == "cofolder.cli:main"
    assert not any("ui" in script.lower() for script in scripts)


def test_bias_training_builder_is_an_importable_package_module():
    spec = importlib.util.find_spec(
        "cofolder.modules.analytics.build_bias_training_data"
    )

    assert spec is not None
    assert spec.origin is not None
    assert Path(spec.origin).name == "build_bias_training_data.py"


def test_package_version_is_consistent():
    repository = Path(__file__).resolve().parents[1]
    pyproject_path = repository / "pyproject.toml"
    pyproject = tomllib.loads(pyproject_path.read_text(encoding="utf-8"))
    package_version = runpy.run_path(repository / "src/cofolder/__init__.py")[
        "__version__"
    ]

    assert "version" not in pyproject["project"]
    assert "version" in pyproject["project"]["dynamic"]
    assert pyproject["tool"]["setuptools"]["dynamic"]["version"] == {
        "attr": "cofolder.__version__"
    }
    assert package_version == "1.0.0"


def test_package_declares_mit_license_metadata():
    repository = Path(__file__).resolve().parents[1]
    pyproject = tomllib.loads(
        (repository / "pyproject.toml").read_text(encoding="utf-8")
    )

    assert pyproject["project"]["license"] == "MIT"
    assert pyproject["project"]["license-files"] == ["LICENSE"]
    assert "setuptools >= 77.0.3" in pyproject["build-system"]["requires"]
