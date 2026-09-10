"""Tests for the boltz1 runner."""

import json
from importlib.metadata import PackageNotFoundError
from unittest.mock import patch

import pandas as pd
import pytest
import yaml

from cofolder.modules.contracts import (
    BackendVersionStatus,
    RepeatSeedProvenance,
    RunnerBackendIdentity,
    SeedOrigin,
)
from cofolder.modules.input import OptionsValidationError
from cofolder.modules.runners.boltz1_runner import Boltz1Runner
from cofolder.modules.runners.contracts import RunnerExecutionRequest


class _MockSystem:
    def __init__(self):
        self._data = {
            "sequences": [
                {"protein": {"id": "A", "fasta": "MKRAAT"}},
                {"ligand": {"id": "B", "smiles": "CCO"}},
            ],
            "properties": [{"affinity": {"binder": "B"}}],
        }

    def find_value(self, key=None, path=None):
        if key == "sequences":
            return self._data["sequences"]
        if key == "properties":
            return self._data["properties"]
        return None


def _typed_options(temp_dir):
    options_path = temp_dir / "typed_options.yaml"
    options_path.write_text(
        "version: 1\nruntime:\n  cache_path: ~/.boltz\n  diffusion_samples: 1\nrunner: {}\n",
        encoding="utf-8",
    )
    return Boltz1Runner().load_options(options_path)


def test_boltz1_runner_requires_exact_100_package_line():
    runner = Boltz1Runner()

    def _fake_version(name):
        if name == "boltz":
            return "2.2.1"
        raise PackageNotFoundError

    with patch(
        "cofolder.modules.runners.base.metadata.version",
        side_effect=_fake_version,
    ):
        available, message = runner.check_availability()

    assert available is False
    assert "requires boltz==1.0.0" in message
    assert "cofolder[boltz1]" in message


def test_boltz1_runner_accepts_exact_100_package_line():
    runner = Boltz1Runner()

    def _fake_version(name):
        if name == "boltz":
            return "1.0.0"
        raise PackageNotFoundError

    with patch(
        "cofolder.modules.runners.base.metadata.version",
        side_effect=_fake_version,
    ):
        available, message = runner.check_availability()

    assert available is True
    assert message is None
    assert runner.capabilities == {"confidence_metrics"}


def test_boltz1_runner_rejects_python_313_before_package_discovery():
    runner = Boltz1Runner()

    with (
        patch("cofolder.modules.runners.base.sys.version_info", (3, 13)),
        patch("cofolder.modules.runners.base.metadata.version") as version,
    ):
        available, message = runner.check_availability()

    assert available is False
    assert "does not support Python 3.13" in message
    assert "Python 3.11 or 3.12" in message
    version.assert_not_called()


def test_boltz1_runner_rejects_workflow_owned_model_flag(temp_dir):
    options_path = temp_dir / "options.yaml"
    options_path.write_text(
        yaml.safe_dump(
            {
                "version": 1,
                "runtime": {"cache_path": "~/.boltz"},
                "runner": {"model": "boltz2"},
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(OptionsValidationError, match="unknown key 'model'"):
        Boltz1Runner().load_options(options_path)


def test_boltz1_runner_marks_affinity_groups_unsupported_without_affinity_payload(
    monkeypatch,
    temp_dir,
):
    runner = Boltz1Runner()
    repeat_dir = temp_dir / "repeat_1"
    request = RunnerExecutionRequest(
        runner_name="boltz1",
        system_name="system",
        system_path=temp_dir / "system.yaml",
        system_obj=_MockSystem(),
        options_path=temp_dir / "options.yaml",
        options_obj=_typed_options(temp_dir),
        repeat=1,
        seed=123,
        seed_provenance=RepeatSeedProvenance(1, 123, 123, 123, 123, SeedOrigin.USER_SPECIFIED),
        backend_identity=RunnerBackendIdentity(
            "boltz1", "boltz", "1.0.0", BackendVersionStatus.DETECTED, "1.0.0"
        ),
        repeat_dir=repeat_dir,
        raw_dir=temp_dir,
        logger=None,
        timings=None,
        label_prefix="repeat_1",
    )

    def _fake_run_boltz(cmd, check=True, timings=None, label_prefix=None):
        prediction_dir = repeat_dir / "boltz_results_system" / "predictions" / "system"
        prediction_dir.mkdir(parents=True, exist_ok=True)
        (prediction_dir / "system_model_0.cif").write_text("data_test", encoding="utf-8")
        (prediction_dir / "confidence_system_model_0.json").write_text(
            json.dumps(
                {
                    "ptm": 0.8,
                    "iptm": 0.7,
                    "confidence_score": 0.9,
                    "chains_ptm": {"0": 0.85, "1": 0.65},
                    "pair_chains_iptm": {"0": {"1": 0.55}, "1": {"0": 0.55}},
                }
            ),
            encoding="utf-8",
        )

    monkeypatch.setattr("cofolder.modules.runners.boltz_runner.run_boltz", _fake_run_boltz)

    result = runner.run(request)

    chain_df = pd.read_csv(result.chain_metrics_path)
    assert result.metric_outcomes["confidence_metrics"].state == "computed"
    assert result.metric_outcomes["affinity_metrics"].state == "unsupported"
    assert result.metric_outcomes["affinity_metrics_ext"].state == "unsupported"
    assert {
        "affinity_pred_value",
        "affinity_probability_binary",
        "pIC50",
        "IC50_M",
        "pIC50_kcal_per_mol",
    }.isdisjoint(chain_df.columns)


def test_boltz1_runner_rejects_unexpected_affinity_payload(
    monkeypatch,
    temp_dir,
):
    runner = Boltz1Runner()
    repeat_dir = temp_dir / "repeat_1"
    request = RunnerExecutionRequest(
        runner_name="boltz1",
        system_name="system",
        system_path=temp_dir / "system.yaml",
        system_obj=_MockSystem(),
        options_path=temp_dir / "options.yaml",
        options_obj=_typed_options(temp_dir),
        repeat=1,
        seed=123,
        seed_provenance=RepeatSeedProvenance(1, 123, 123, 123, 123, SeedOrigin.USER_SPECIFIED),
        backend_identity=RunnerBackendIdentity(
            "boltz1", "boltz", "1.0.0", BackendVersionStatus.DETECTED, "1.0.0"
        ),
        repeat_dir=repeat_dir,
        raw_dir=temp_dir,
        logger=None,
        timings=None,
        label_prefix="repeat_1",
    )

    def _fake_run_boltz(cmd, check=True, timings=None, label_prefix=None):
        prediction_dir = repeat_dir / "boltz_results_system" / "predictions" / "system"
        prediction_dir.mkdir(parents=True, exist_ok=True)
        (prediction_dir / "system_model_0.cif").write_text("data_test", encoding="utf-8")
        (prediction_dir / "confidence_system_model_0.json").write_text(
            json.dumps(
                {
                    "ptm": 0.8,
                    "iptm": 0.7,
                    "confidence_score": 0.9,
                    "chains_ptm": {"0": 0.85, "1": 0.65},
                    "pair_chains_iptm": {"0": {"1": 0.55}, "1": {"0": 0.55}},
                }
            ),
            encoding="utf-8",
        )
        (prediction_dir / "affinity_system.json").write_text(
            json.dumps(
                {
                    "affinity_pred_value": 6.1,
                    "affinity_probability_binary": 0.8,
                }
            ),
            encoding="utf-8",
        )

    monkeypatch.setattr("cofolder.modules.runners.boltz_runner.run_boltz", _fake_run_boltz)

    with pytest.raises(ValueError, match="does not declare affinity capabilities"):
        runner.run(request)
