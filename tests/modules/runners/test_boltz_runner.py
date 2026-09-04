"""Tests for Boltz-2 runner timing parsing."""

import json
import pickle
from pathlib import Path
from unittest.mock import patch
from importlib.metadata import PackageNotFoundError

import pandas as pd
import pytest
import yaml
from rdkit import Chem

from cofolder.modules.input.command import Command
from cofolder.modules.input.system import System
from cofolder.modules.input.validation import SystemInputValidationError
from cofolder.modules.runners.contracts import RunnerExecutionRequest
from cofolder.modules.runners.boltz_runner import _parse_boltz_stage_timings
from cofolder.modules.runners.boltz2_runner import Boltz2Runner


def _prepared_constrained_system(atom_name: str) -> System:
    return System(
        system={
            "sequences": [
                {"protein": {"id": "A", "sequence": "AC"}},
                {"ligand": {"id": "L", "ccd": "NEW"}},
            ],
            "constraints": [
                {"bond": {"atom1": ["A", 2, "SG"], "atom2": ["L", 1, atom_name]}}
            ],
        }
    )


def test_conformer_preparation_preserves_compatible_constraint(monkeypatch, temp_dir):
    cache = temp_dir / "cache"
    mols = cache / "mols"
    mols.mkdir(parents=True)
    mol = Chem.MolFromSmiles("CCO")
    for index, atom in enumerate(mol.GetAtoms(), 1):
        atom.SetProp("name", f"{atom.GetSymbol().upper()}{index}")
    with (mols / "NEW.pkl").open("wb") as handle:
        pickle.dump(mol, handle)

    system = _prepared_constrained_system("O3")
    original_constraints = json.loads(json.dumps(system.system["constraints"]))
    runner = Boltz2Runner()
    options = Command(options={"options": [{"cache": str(cache)}]})
    monkeypatch.setattr(
        "cofolder.modules.entities.ligand.handle_conformers",
        lambda **kwargs: {"L": "NEW"},
    )

    prepared = runner.prepare_system(system, options, temp_dir, "3D", None, None)
    runner.validate_system(prepared.system_obj, options, check_atom_names=True)

    assert prepared.system_obj.system["constraints"] == original_constraints


def test_post_conformer_validation_rejects_missing_ligand_atom(temp_dir):
    cache = temp_dir / "cache"
    mols = cache / "mols"
    mols.mkdir(parents=True)
    mol = Chem.MolFromSmiles("CCO")
    for index, atom in enumerate(mol.GetAtoms(), 1):
        atom.SetProp("name", f"{atom.GetSymbol().upper()}{index}")
    with (mols / "NEW.pkl").open("wb") as handle:
        pickle.dump(mol, handle)

    with pytest.raises(SystemInputValidationError, match="atom 'N99'.*does not exist"):
        Boltz2Runner().validate_system(
            _prepared_constrained_system("N99"),
            Command(options={"options": [{"cache": str(cache)}]}),
            check_atom_names=True,
        )


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


def test_msa_reuse_settings_include_generation_options_and_backend_version(monkeypatch):
    runner = Boltz2Runner()
    options = Command(
        options={
            "options": [
                {"msa_server_url": "https://msa.example.test"},
                {"msa_pairing_strategy": "complete"},
                {"max_msa_seqs": 4096},
                {"msa_server_password": "must-not-be-persisted"},
                {"diffusion_samples": 2},
            ]
        }
    )
    monkeypatch.setattr(runner, "get_distribution_version", lambda name: "2.2.1")

    assert runner.msa_reuse_settings(options) == {
        "msa_server_url": "https://msa.example.test",
        "msa_pairing_strategy": "complete",
        "max_msa_seqs": 4096,
        "backend_version": "2.2.1",
    }


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
    runner = Boltz2Runner()
    repeat_dir = temp_dir / "repeat_1"
    request = RunnerExecutionRequest(
        runner_name="boltz2",
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
    assert result.metric_outcomes["confidence_metrics"].state == "computed"
    assert result.metric_outcomes["affinity_metrics"].state == "computed"
    assert result.metric_outcomes["affinity_metrics_ext"].state == "computed"


def test_boltz_runner_reports_mixed_metric_outcomes_without_affinity_payload(monkeypatch, temp_dir):
    runner = Boltz2Runner()
    repeat_dir = temp_dir / "repeat_1"
    request = RunnerExecutionRequest(
        runner_name="boltz2",
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

    assert result.metric_outcomes["confidence_metrics"].state == "computed"
    assert result.metric_outcomes["affinity_metrics"].state == "missing"
    assert result.metric_outcomes["affinity_metrics"].required_columns == (
        "chain_metrics.affinity_pred_value",
        "chain_metrics.affinity_probability_binary",
    )
    assert result.metric_outcomes["affinity_metrics_ext"].state == "missing"
    assert result.metric_outcomes["affinity_metrics_ext"].required_columns == (
        "chain_metrics.pIC50",
        "chain_metrics.IC50_M",
        "chain_metrics.pIC50_kcal_per_mol",
    )


def test_boltz2_runner_rejects_environment_with_boltz_community_installed():
    runner = Boltz2Runner()

    with patch(
        "cofolder.modules.runners.base.metadata.version",
        side_effect=lambda name: "1.0" if name in {"boltz", "boltz-community"} else None,
    ):
        available, message = runner.check_availability()

    assert available is False
    assert "same environment" in message
    assert "boltz-community" in message


def test_boltz2_runner_rejects_boltz1_package_line():
    runner = Boltz2Runner()

    def _fake_version(name):
        if name == "boltz":
            return "1.0.0"
        raise PackageNotFoundError

    with patch(
        "cofolder.modules.runners.base.metadata.version",
        side_effect=_fake_version,
    ):
        available, message = runner.check_availability()

    assert available is False
    assert "Boltz-2 package line" in message
    assert "cofolder[boltz2]" in message


def test_boltz2_runner_load_options_forces_boltz2_model(temp_dir):
    options_path = temp_dir / "options.yaml"
    options_path.write_text(
        yaml.safe_dump({"options": [{"cache": "~/.boltz"}, {"diffusion_samples": 1}]}),
        encoding="utf-8",
    )

    command = Boltz2Runner().load_options(options_path)

    assert command.find_value(key="model") == "boltz2"
