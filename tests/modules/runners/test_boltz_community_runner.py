"""Tests for the boltz-community runner."""

import json
from importlib.metadata import PackageNotFoundError
from unittest.mock import patch

import pandas as pd

from cofolder.modules.contracts import (
    BackendVersionStatus,
    RepeatSeedProvenance,
    RunnerBackendIdentity,
    SeedOrigin,
)
from cofolder.modules.input.command import Command
from cofolder.modules.runners.boltz_community_runner import BoltzCommunityRunner
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


def test_boltz_community_runner_accepts_installed_package_line():
    runner = BoltzCommunityRunner()

    def _fake_version(name):
        if name == "boltz-community":
            return "0.6.0"
        raise PackageNotFoundError

    with patch(
        "cofolder.modules.runners.base.metadata.version",
        side_effect=_fake_version,
    ):
        available, message = runner.check_availability()

    assert available is True
    assert message is None
    assert runner.capabilities == {
        "confidence_metrics",
        "affinity_metrics",
        "affinity_metrics_ext",
    }


def test_boltz_community_runner_uses_same_normalized_bundle_as_boltz(monkeypatch, temp_dir):
    runner = BoltzCommunityRunner()
    repeat_dir = temp_dir / "repeat_1"
    request = RunnerExecutionRequest(
        runner_name="boltz-community",
        system_name="system",
        system_path=temp_dir / "system.yaml",
        system_obj=_MockSystem(),
        options_path=temp_dir / "options.yaml",
        options_obj=Command(options={"options": [{"diffusion_samples": 1}, {"cache": "~/.boltz"}]}),
        repeat=1,
        seed=123,
        seed_provenance=RepeatSeedProvenance(1, 123, 123, 123, 123, SeedOrigin.USER_SPECIFIED),
        backend_identity=RunnerBackendIdentity(
            "boltz-community",
            "boltz-community",
            None,
            BackendVersionStatus.UNAVAILABLE,
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

    result = runner.run(request)

    assert result.runner_name == "boltz-community"
    assert result.system_metrics_path.exists()
    assert result.chain_metrics_path.exists()
    assert result.manifest_path.exists()
    assert (result.structures_dir / "1_system_model_0.cif").exists()

    system_df = pd.read_csv(result.system_metrics_path)
    chain_df = pd.read_csv(result.chain_metrics_path)
    assert {"ptm", "iptm", "confidence_score"}.issubset(system_df.columns)
    assert {"chains_ptm", "affinity_pred_value", "affinity_probability_binary", "pIC50"}.issubset(
        chain_df.columns
    )
    assert result.metric_outcomes["confidence_metrics"].state == "computed"
    assert result.metric_outcomes["affinity_metrics"].state == "computed"
    assert result.metric_outcomes["affinity_metrics_ext"].state == "computed"
    manifest = json.loads(result.manifest_path.read_text(encoding="utf-8"))
    assert manifest["backend"]["version"] is None
    assert manifest["backend"]["version_status"] == "unavailable"
    assert manifest["backend"]["raw_version"] is None


def test_boltz_community_runner_rejects_environment_with_boltz_installed():
    runner = BoltzCommunityRunner()

    with patch(
        "cofolder.modules.runners.base.metadata.version",
        side_effect=lambda name: "1.0" if name in {"boltz", "boltz-community"} else None,
    ):
        available, message = runner.check_availability()

    assert available is False
    assert "same environment" in message
    assert "boltz" in message
