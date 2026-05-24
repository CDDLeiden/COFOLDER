"""Tests for the boltz1 runner."""

import json
from importlib.metadata import PackageNotFoundError
from unittest.mock import patch

import pandas as pd
import pytest
import yaml

from cofolder.modules.input.command import Command
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


def test_boltz1_runner_load_options_removes_model_flag(temp_dir):
    options_path = temp_dir / "options.yaml"
    options_path.write_text(
        yaml.safe_dump({"options": [{"cache": "~/.boltz"}, {"model": "boltz2"}]}),
        encoding="utf-8",
    )

    command = Boltz1Runner().load_options(options_path)

    assert command.find_value(key="model") is None


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
        options_obj=Command(options={"options": [{"diffusion_samples": 1}, {"cache": "~/.boltz"}]}),
        repeat=1,
        seed=123,
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
        options_obj=Command(options={"options": [{"diffusion_samples": 1}, {"cache": "~/.boltz"}]}),
        repeat=1,
        seed=123,
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
