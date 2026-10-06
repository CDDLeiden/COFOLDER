"""Checks for the repository-owned supported-test partition."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


REPOSITORY = Path(__file__).resolve().parents[1]
RUNNER_PATH = REPOSITORY / "scripts/run_test_lane.py"


def _load_runner():
    spec = importlib.util.spec_from_file_location("run_test_lane", RUNNER_PATH)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_supported_tests_are_partitioned_exactly_once():
    runner = _load_runner()
    discovered = runner.discover_supported_tests(REPOSITORY)
    classified = runner.classify_tests(discovered)
    assigned = [path for paths in classified.values() for path in paths]

    assert set(classified) == set(runner.LANE_NAMES)
    assert len(assigned) == len(set(assigned)) == len(discovered)
    assert set(assigned) == set(discovered)


def test_lane_ownership_matches_release_boundaries():
    runner = _load_runner()
    classified = runner.classify_tests(runner.discover_supported_tests(REPOSITORY))

    assert Path("tests/test_cli.py") in classified["core"]
    assert (
        Path("tests/modules/contracts/test_public_contracts.py")
        in classified["contracts-tutorial"]
    )
    assert Path("tests/test_tutorials.py") in classified["contracts-tutorial"]
    assert (
        Path("tests/test_structure_gated_oracle_tutorial.py")
        in classified["contracts-tutorial"]
    )
    assert Path("tests/test_packaging.py") in classified["artifact"]
    assert Path("tests/acceptance/test_shared.py") in classified["acceptance"]


@pytest.mark.parametrize(
    "paths",
    [
        [Path("tests/unknown/test_example.py")],
        [Path("tests/test_packaging.py"), Path("tests/unknown/test_example.py")],
    ],
)
def test_unclassified_test_modules_fail_the_partition(paths):
    runner = _load_runner()

    with pytest.raises(runner.LaneConfigurationError, match="exactly one lane"):
        runner.classify_tests(paths)


def test_overlapping_test_modules_fail_the_partition():
    runner = _load_runner()
    runner.LANE_RULES["duplicate-artifact"] = lambda path: path == Path(
        "tests/test_packaging.py"
    )

    with pytest.raises(runner.LaneConfigurationError, match="duplicate-artifact"):
        runner.classify_tests([Path("tests/test_packaging.py")])
