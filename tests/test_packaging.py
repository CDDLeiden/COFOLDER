"""Packaging metadata regression tests."""

from __future__ import annotations

from fnmatch import fnmatch
from pathlib import Path
import importlib.util
import runpy
import tomllib

import yaml


REPOSITORY = Path(__file__).resolve().parents[1]

EXAMPLE_RESOURCE_NAMES = {
    "4HJO.cif",
    "4HJO.pdb",
    "README.md",
    "ifp_clustering_demo.py",
    "ligand_screen.csv",
    "options.yaml",
    "system.yaml",
    "system_covalent.yaml",
    "system_nucleic_acid.yaml",
    "system_screen.yaml",
}

ACCEPTANCE_RESOURCE_NAMES = {
    "ligand_screen.csv",
    "ligand_training_data.csv",
    "options_acceptance.yaml",
    "options_boltz1_acceptance.yaml",
    "options_openfold3_acceptance.yaml",
    "protein_training_data.csv",
    "system.yaml",
    "system_constraints_boltz1.yaml",
    "system_constraints_boltz2.yaml",
    "system_constraints_openfold3.yaml",
    "system_nucleic_acid.yaml",
    "system_screen.yaml",
}


def _load_pyproject() -> dict:
    return tomllib.loads((REPOSITORY / "pyproject.toml").read_text(encoding="utf-8"))


def _manifest_directives() -> set[str]:
    return {
        line.strip()
        for line in (REPOSITORY / "MANIFEST.in").read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    }


def _resource_names(path: Path) -> set[str]:
    return {
        item.name
        for item in path.iterdir()
        if item.is_file() and item.name != "__init__.py"
    }


def _assert_names_match_package_data(names: set[str], patterns: list[str]) -> None:
    assert all(any(fnmatch(name, pattern) for pattern in patterns) for name in names)


def test_default_dependencies_include_core_runtime_packages():
    pyproject = _load_pyproject()
    dependencies = pyproject["project"]["dependencies"]

    assert "pandas" in dependencies
    assert "pyyaml" in dependencies
    assert "biopython" in dependencies
    assert "scipy" in dependencies
    assert "scikit-learn" in dependencies


def test_backend_extras_use_supported_immutable_versions():
    pyproject = _load_pyproject()
    optional_dependencies = pyproject["project"]["optional-dependencies"]

    assert optional_dependencies["boltz1"] == [
        "boltz==1.0.0; python_version < '3.13'"
    ]
    assert optional_dependencies["boltz2"] == [
        "boltz[cuda]>=2.0.0,<3; python_version < '3.13'"
    ]
    assert optional_dependencies["boltz-community"] == [
        "boltz-community[cuda]==2.10.12"
    ]
    assert not any(
        "git+" in dependency
        for dependencies in optional_dependencies.values()
        for dependency in dependencies
    )


def test_tutorial_extra_includes_marimo():
    pyproject = _load_pyproject()
    tutorials = pyproject["project"]["optional-dependencies"]["tutorials"]

    assert any(dep.startswith("marimo") for dep in tutorials)


def test_acceptance_extra_includes_marimo():
    pyproject = _load_pyproject()
    acceptance = pyproject["project"]["optional-dependencies"]["acceptance"]

    assert any(dep.startswith("marimo") for dep in acceptance)


def test_test_extra_is_pytest_focused():
    pyproject = _load_pyproject()
    test = pyproject["project"]["optional-dependencies"]["test"]

    assert any(dep.startswith("pytest") for dep in test)
    assert not any(dep.startswith("marimo") for dep in test)
    assert not any(dep.lower().startswith("streamlit") for dep in test)


def test_development_extra_exposes_repository_tools():
    pyproject = _load_pyproject()
    optional_dependencies = pyproject["project"]["optional-dependencies"]
    development = optional_dependencies["development"]

    assert any(dep.startswith("black") for dep in development)
    assert any(dep.startswith("build") for dep in development)
    assert any(dep.startswith("ruff") for dep in development)
    assert not any(dep.lower().startswith("streamlit") for dep in development)


def test_routine_ci_uses_supported_lanes_without_ui_or_backends():
    workflow = (REPOSITORY / ".github/workflows/quality.yml").read_text(
        encoding="utf-8"
    )
    normalized = workflow.lower()
    jobs = yaml.safe_load(workflow)["jobs"]
    matrix = jobs["supported-tests"]["strategy"]["matrix"]["include"]
    expected_test_jobs = {
        (version, lane)
        for version in ("3.11", "3.12")
        for lane in ("core", "contracts-tutorial", "acceptance")
    }

    assert {(job["python-version"], job["lane"]) for job in matrix} == (
        expected_test_jobs
    )
    for lane in ("core", "contracts-tutorial", "artifact", "acceptance"):
        assert lane in workflow
    assert "streamlit" not in normalized
    assert "[acceptance,boltz" not in normalized
    assert "[acceptance,openfold3" not in normalized


def test_unsupported_ui_is_absent_from_install_metadata():
    pyproject = _load_pyproject()

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
    assert scripts == {
        "cofolder": "cofolder.cli:main",
        "cofolder-tools": "cofolder.tools.cli:main",
    }
    assert not any("ui" in script.lower() for script in scripts)


def test_bias_training_builder_is_an_importable_package_module():
    spec = importlib.util.find_spec(
        "cofolder.modules.analytics.build_bias_training_data"
    )

    assert spec is not None
    assert spec.origin is not None
    assert Path(spec.origin).name == "build_bias_training_data.py"

    tools_spec = importlib.util.find_spec("cofolder.tools.build_bias_training_data")
    assert tools_spec is not None


def test_package_version_is_consistent():
    pyproject = _load_pyproject()
    package_version = runpy.run_path(REPOSITORY / "src/cofolder/__init__.py")[
        "__version__"
    ]

    assert "version" not in pyproject["project"]
    assert "version" in pyproject["project"]["dynamic"]
    assert pyproject["tool"]["setuptools"]["dynamic"]["version"] == {
        "attr": "cofolder.__version__"
    }
    assert package_version == "1.0.0"


def test_package_declares_mit_license_metadata():
    pyproject = _load_pyproject()

    assert pyproject["project"]["license"] == "MIT"
    assert pyproject["project"]["license-files"] == [
        "LICENSE",
        "THIRD_PARTY_SOFTWARE.md",
    ]
    assert "setuptools >= 77.0.3" in pyproject["build-system"]["requires"]


def test_packaged_examples_match_canonical_examples_byte_for_byte():
    canonical_dir = REPOSITORY / "examples"
    packaged_dir = REPOSITORY / "src/cofolder/resources/examples"

    assert importlib.util.find_spec("cofolder.resources.examples") is not None
    assert _resource_names(canonical_dir) == EXAMPLE_RESOURCE_NAMES
    assert _resource_names(packaged_dir) == EXAMPLE_RESOURCE_NAMES
    for name in EXAMPLE_RESOURCE_NAMES:
        assert (packaged_dir / name).read_bytes() == (canonical_dir / name).read_bytes()


def test_package_data_declares_acceptance_and_example_resources():
    package_data = _load_pyproject()["tool"]["setuptools"]["package-data"]
    acceptance_patterns = package_data["cofolder.acceptance"]
    example_patterns = package_data["cofolder.resources.examples"]

    acceptance_dir = REPOSITORY / "src/cofolder/acceptance/data"
    assert _resource_names(acceptance_dir) == ACCEPTANCE_RESOURCE_NAMES
    _assert_names_match_package_data(
        {f"data/{name}" for name in ACCEPTANCE_RESOURCE_NAMES},
        acceptance_patterns,
    )
    _assert_names_match_package_data(EXAMPLE_RESOURCE_NAMES, example_patterns)


def test_sdist_manifest_has_explicit_supported_source_inventory():
    directives = _manifest_directives()
    required = {
        "include pyproject.toml",
        "include README.md",
        "include LICENSE",
        "include THIRD_PARTY_SOFTWARE.md",
        "include mkdocs.yml",
        "recursive-include .github/workflows *.yml",
        "recursive-include src/cofolder *.py",
        "recursive-include src/cofolder/acceptance/data *.csv *.yaml",
        (
            "recursive-include src/cofolder/resources/examples "
            "*.cif *.csv *.md *.pdb *.py *.yaml"
        ),
        "recursive-include tests *.py",
        "recursive-include docs *.md",
        "recursive-include tutorials *.md *.py",
        "recursive-include examples *.cif *.csv *.md *.pdb *.py *.yaml",
        "recursive-include scripts *.md *.py *.sh",
    }

    assert required <= directives
    assert not any(line.startswith("recursive-include . ") for line in directives)


def test_sdist_manifest_excludes_unsupported_and_generated_material():
    directives = _manifest_directives()
    required_exclusions = {
        "exclude tests/ui_development_check.py",
        "prune src/cofolder/ui",
        "prune legacy",
        "prune .streamlit",
        "prune docs/project-knowledge",
        "global-exclude __pycache__",
        "global-exclude *.py[cod]",
        "global-exclude .openfold3-cache",
        "global-exclude __marimo__",
        "global-exclude *.ckpt",
        "global-exclude *.pt",
        "global-exclude *.safetensors",
        "global-exclude *.onnx",
    }

    assert required_exclusions <= directives
