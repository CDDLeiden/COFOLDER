"""Tests for Boltz runner timing parsing."""

import json
from pathlib import Path

import pandas as pd

from cofolder.modules.input.command import Command
from cofolder.modules.runners.base import RunnerRequest
from cofolder.modules.runners.boltz_runner import BoltzRunner
from cofolder.modules.runners.boltz_runner import _parse_boltz_stage_timings


def test_parse_boltz_stage_timings_with_msa_and_affinity():
    timed_lines = [
        (0.5, "Calling MSA server for target system with 1 sequences"),
        (1.0, "SUBMIT:   0%|          | 0/150 [elapsed: 00:00 remaining: ?]"),
        (4.0, "COMPLETE: 100%|██████████| 150/150 [elapsed: 00:03 remaining: 00:00]"),
        (7.0, "Running affinity prediction for 1 input."),
    ]

    out = _parse_boltz_stage_timings(timed_lines, total_elapsed=15.0)

    assert out["boltz.msa"] == 3.5
    assert out["boltz.affinity_prediction"] == 8.0


def test_parse_boltz_stage_timings_without_msa_or_affinity():
    timed_lines = [
        (0.2, "Checking input data."),
        (3.0, "Predicting structure for 1 input."),
    ]

    out = _parse_boltz_stage_timings(timed_lines, total_elapsed=10.0)

    assert out == {}


def test_parse_boltz_stage_timings_skips_affinity_when_outputs_already_exist():
    timed_lines = [
        (0.4, "Predicting property: affinity"),
        (0.5, "Found some existing affinity predictions (1), skipping and running only the missing ones, if any. If you wish to override these existing affinity predictions, please set the --override flag."),
        (0.6, "Found existing affinity predictions for all inputs, skipping."),
    ]

    out = _parse_boltz_stage_timings(timed_lines, total_elapsed=5.0)

    assert out == {}


def test_parse_boltz_stage_timings_ignores_incomplete_msa_window():
    timed_lines = [
        (0.5, "Calling MSA server for target system with 1 sequences"),
        (2.0, "SUBMIT:   0%|          | 0/150 [elapsed: 00:00 remaining: ?]"),
        (4.0, "Running affinity prediction for 1 input."),
    ]

    out = _parse_boltz_stage_timings(timed_lines, total_elapsed=9.0)

    assert "boltz.msa" not in out
    assert out["boltz.affinity_prediction"] == 5.0


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


def test_boltz_runner_writes_canonical_bundle(monkeypatch, temp_dir):
    runner = BoltzRunner()
    repeat_dir = temp_dir / "repeat_1"
    request = RunnerRequest(
        runner_name="boltz",
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

    result = runner.run(request)

    assert result.system_metrics_path.exists()
    assert result.chain_metrics_path.exists()
    assert result.manifest_path.exists()
    assert (result.structures_dir / "1_system_model_0.cif").exists()

    system_df = pd.read_csv(result.system_metrics_path)
    chain_df = pd.read_csv(result.chain_metrics_path)
    assert {"ptm", "iptm", "confidence_score"}.issubset(system_df.columns)
    assert {"chains_ptm", "affinity_pred_value", "affinity_probability_binary", "pIC50"}.issubset(chain_df.columns)
