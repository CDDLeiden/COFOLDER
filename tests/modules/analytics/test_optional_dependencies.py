"""Regression tests for the optional analysis dependency boundary."""

from __future__ import annotations

import subprocess
import sys

import pandas as pd
import pytest

from cofolder.modules.analytics import ifp_clustering, plots, reference_ifp, stats, structure
from cofolder.modules.utils import write
from cofolder.modules.utils import _optional_dependencies


OPTIONAL_ROOTS = {
    "scipy",
    "sklearn",
    "prolif",
    "MDAnalysis",
    "pdb2pqr",
    "matplotlib",
    "seaborn",
}
BASE_IMPORT_BLOCKED_ROOTS = OPTIONAL_ROOTS | {"boltz"}


def _block_optional_imports(monkeypatch, blocked_roots: set[str]) -> None:
    original_import_module = _optional_dependencies.importlib.import_module

    def guarded_import_module(name: str, package: str | None = None):
        if name.split(".", maxsplit=1)[0] in blocked_roots:
            raise ModuleNotFoundError(
                f"No module named '{name}'",
                name=name,
            )
        return original_import_module(name, package)

    monkeypatch.setattr(
        _optional_dependencies.importlib,
        "import_module",
        guarded_import_module,
    )


def test_base_safe_import_graph_does_not_load_analysis_dependencies():
    script = f"""
import importlib.abc
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

blocked = {BASE_IMPORT_BLOCKED_ROOTS!r}

class BlockAnalysis(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.', 1)[0] in blocked:
            raise ModuleNotFoundError(f"blocked optional analysis import: {{fullname}}", name=fullname)
        return None

sys.meta_path.insert(0, BlockAnalysis())

import cofolder
import cofolder.cli
import cofolder.recipes.bias
import cofolder.recipes.oracle
import cofolder.recipes.screen
import cofolder.recipes.validate
import cofolder.modules.runners.boltz_runner
from cofolder.modules.analytics.stats import affinity_to_pic50_and_ic50
from cofolder.modules.utils import write
from cofolder.resources.examples import copy_examples

assert affinity_to_pic50_and_ic50(6.0) == (0.0, 1.0)
assert not blocked.intersection(sys.modules)
assert callable(copy_examples)

with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    write.write_csv(pd.DataFrame({{"value": [1]}}), str(root / "table.csv"))
    write.write_yaml(SimpleNamespace(system={{"version": 1}}), root / "input.yaml")
    copy_examples(root / "examples")
    assert (root / "table.csv").is_file()
    assert (root / "input.yaml").is_file()
    assert (root / "examples" / "system.yaml").is_file()

for argv in (["--help"], ["--version"], *([name, "--help"] for name in ("bias", "validate", "screen", "oracle"))):
    try:
        cofolder.cli.main(argv)
    except SystemExit as exc:
        assert exc.code == 0

assert not blocked.intersection(sys.modules)
"""
    completed = subprocess.run(
        [sys.executable, "-c", script],
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr


@pytest.mark.parametrize(
    ("blocked_roots", "operation", "feature", "dependency"),
    [
        (
            {"scipy"},
            lambda: ifp_clustering.cluster_binary_ifps(
                [[1, 0]], ["member"], similarity_threshold=0.5
            ),
            "Interaction-fingerprint clustering",
            "scipy",
        ),
        (
            {"scipy"},
            lambda: stats.calculate_affinity_correlations(
                pd.DataFrame({"pred": [1.0, 2.0], "exp": [1.0, 2.0]}),
                "pred",
                "exp",
            ),
            "Affinity correlation statistics",
            "scipy",
        ),
        (
            {"sklearn"},
            lambda: stats.calculate_affinity_correlations(
                pd.DataFrame({"pred": [1.0, 2.0], "exp": [1.0, 2.0]}),
                "pred",
                "exp",
            ),
            "Affinity correlation statistics",
            "scikit-learn",
        ),
        (
            {"matplotlib"},
            lambda: plots.plot_reference_overlap_scatter(
                pd.DataFrame(
                    {
                        "x": [0.8],
                        "y": [0.9],
                        "query_1": ["A"],
                        "query_2": ["B"],
                    }
                ),
                output_dir="unused",
                file_stem="unused",
                x_col="x",
                y_col="y",
                x_threshold=0.0,
                y_threshold=0.0,
                x_label="x",
                y_label="y",
                title="title",
                query_1_col="query_1",
                query_2_col="query_2",
            ),
            "Analysis plotting",
            "matplotlib",
        ),
        (
            {"seaborn"},
            lambda: plots.plot_reference_overlap_scatter(
                pd.DataFrame(
                    {
                        "x": [0.8],
                        "y": [0.9],
                        "query_1": ["A"],
                        "query_2": ["B"],
                    }
                ),
                output_dir="unused",
                file_stem="unused",
                x_col="x",
                y_col="y",
                x_threshold=0.0,
                y_threshold=0.0,
                x_label="x",
                y_label="y",
                title="title",
                query_1_col="query_1",
                query_2_col="query_2",
            ),
            "Analysis plotting",
            "seaborn",
        ),
        (
            {"matplotlib"},
            lambda: write.save_plot("unused.png"),
            "Plot saving",
            "matplotlib",
        ),
        (
            {"MDAnalysis"},
            structure._load_prolif_dependencies,
            "ProLIF interaction processing",
            "MDAnalysis",
        ),
        (
            {"prolif"},
            structure._load_prolif_dependencies,
            "ProLIF interaction processing",
            "prolif",
        ),
        (
            {"pdb2pqr"},
            lambda: structure._run_pdb2pqr(["input.pdb", "output.pdb"]),
            "PDB2PQR interaction preparation",
            "pdb2pqr",
        ),
        (
            {"prolif"},
            lambda: reference_ifp._extract_prolif(None, None, None, None),
            "ProLIF interaction fingerprints",
            "prolif",
        ),
    ],
)
def test_missing_analysis_dependency_has_feature_specific_guidance(
    monkeypatch,
    blocked_roots,
    operation,
    feature,
    dependency,
):
    _block_optional_imports(monkeypatch, blocked_roots)

    with pytest.raises(ImportError) as error:
        operation()

    message = str(error.value)
    assert feature in message
    assert dependency in message
    assert 'pip install "cofolder[analysis]"' in message
