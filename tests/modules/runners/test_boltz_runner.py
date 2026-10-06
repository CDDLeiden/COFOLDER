"""Tests for Boltz-2 runner timing parsing."""

import json
import pickle
from importlib.metadata import PackageNotFoundError
from unittest.mock import patch

import pandas as pd
import pytest
import yaml
from rdkit import Chem

from cofolder.modules.contracts import (
    BackendVersionStatus,
    RepeatSeedProvenance,
    RunnerBackendIdentity,
    RunnerProvenanceError,
    SeedOrigin,
)
from cofolder.modules.input.config import RunnerOptions, RunnerRuntimeOptions
from cofolder.modules.input.system import System
from cofolder.modules.input.validation import SystemInputValidationError
from cofolder.modules.runners.boltz2_runner import Boltz2Runner
from cofolder.modules.runners.boltz_runner import _parse_boltz_stage_timings
from cofolder.modules.runners.contracts import RunnerExecutionRequest
from cofolder.modules.runners.validators import validate_runner_bundle


def _request_provenance(
    *,
    version: str | None = "2.1.0",
    status: BackendVersionStatus = BackendVersionStatus.DETECTED,
    raw_version: str | None = "2.1.0",
) -> dict:
    return {
        "seed_provenance": RepeatSeedProvenance(
            repeat_id=1,
            requested_base_seed=123,
            resolved_base_seed=123,
            derived_seed=123,
            effective_seed=123,
            origin=SeedOrigin.USER_SPECIFIED,
        ),
        "backend_identity": RunnerBackendIdentity(
            runner_name="boltz2",
            backend_name="boltz",
            version=version,
            version_status=status,
            raw_version=raw_version,
        ),
    }


def _typed_options(temp_dir, *, runner_options=None, diffusion_samples=1):
    options_path = temp_dir / "typed_options.yaml"
    options_path.write_text(
        yaml.safe_dump(
            {
                "version": 1,
                "runtime": {
                    "cache_path": str(temp_dir / "cache"),
                    "diffusion_samples": diffusion_samples,
                },
                "runner": runner_options or {},
            }
        ),
        encoding="utf-8",
    )
    return Boltz2Runner().load_options(options_path)


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
    options = _typed_options(temp_dir)
    captured = {}

    def fake_prepare(system_obj, **kwargs):
        captured["system_obj"] = system_obj
        captured.update(kwargs)
        return {"L": "NEW"}

    monkeypatch.setattr(
        "cofolder.modules.runners.boltz_runner.prepare_ligand_conformers",
        fake_prepare,
    )

    prepared = runner.prepare_system(system, options, temp_dir, "3D", None, None)
    runner.validate_system(prepared.system_obj, options, check_atom_names=True)

    assert captured["system_obj"] is system
    assert captured["cache_path"] == temp_dir / "cache"
    assert captured["conformers"] == "3D"
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
            RunnerOptions(
                version=1,
                runtime=RunnerRuntimeOptions(cache_path=cache),
                runner={},
            ),
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


def test_msa_reuse_settings_include_generation_options_and_backend_version(
    monkeypatch, temp_dir
):
    runner = Boltz2Runner()
    options = _typed_options(
        temp_dir,
        diffusion_samples=2,
        runner_options={
            "msa_server_url": "https://msa.example.test",
            "msa_pairing_strategy": "complete",
            "max_msa_seqs": 4096,
            "msa_server_password": "must-not-be-persisted",
        },
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
        options_obj=_typed_options(temp_dir),
        repeat=1,
        seed=123,
        **_request_provenance(),
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
    manifest = json.loads(result.manifest_path.read_text(encoding="utf-8"))
    assert manifest["backend"] == {
        "runner_name": "boltz2",
        "backend_name": "boltz",
        "version": "2.1.0",
        "version_status": "detected",
        "raw_version": "2.1.0",
        "detail": None,
    }
    assert manifest["seed"]["effective_seed"] == 123
    assert system_df.loc[0, "effective_seed"] == 123
    assert result.sample_records[0]["backend_name"] == "boltz"

    manifest["seed"]["effective_seed"] = 124
    result.manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(RunnerProvenanceError, match="manifest seed provenance"):
        validate_runner_bundle(result, requested_metric_groups={"confidence_metrics"})


def test_boltz_runner_reports_mixed_metric_outcomes_without_affinity_payload(monkeypatch, temp_dir):
    runner = Boltz2Runner()
    repeat_dir = temp_dir / "repeat_1"
    request = RunnerExecutionRequest(
        runner_name="boltz2",
        system_name="system",
        system_path=temp_dir / "system.yaml",
        system_obj=_MockSystem(),
        options_path=temp_dir / "options.yaml",
        options_obj=_typed_options(temp_dir),
        repeat=1,
        seed=123,
        **_request_provenance(
            version=None,
            status=BackendVersionStatus.UNPARSEABLE,
            raw_version="not-a-version",
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
    manifest = json.loads(result.manifest_path.read_text(encoding="utf-8"))
    assert manifest["backend"]["version"] is None
    assert manifest["backend"]["version_status"] == "unparseable"
    assert manifest["backend"]["raw_version"] == "not-a-version"


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


@pytest.mark.parametrize("installed_version", ["2.0.0", "2.2.1", "2.99.0"])
def test_boltz2_runner_accepts_supported_major_line(installed_version):
    runner = Boltz2Runner()

    def _fake_version(name):
        if name == "boltz":
            return installed_version
        raise PackageNotFoundError

    with patch(
        "cofolder.modules.runners.base.metadata.version",
        side_effect=_fake_version,
    ):
        available, message = runner.check_availability()

    assert available is True
    assert message is None


def test_boltz2_runner_rejects_next_major_package_line():
    runner = Boltz2Runner()

    def _fake_version(name):
        if name == "boltz":
            return "3.0.0"
        raise PackageNotFoundError

    with patch(
        "cofolder.modules.runners.base.metadata.version",
        side_effect=_fake_version,
    ):
        available, message = runner.check_availability()

    assert available is False
    assert "Boltz-2 package line" in message
    assert "cofolder[boltz2]" in message


def test_boltz2_runner_rejects_python_313_before_package_discovery():
    runner = Boltz2Runner()

    with (
        patch("cofolder.modules.runners.base.sys.version_info", (3, 13)),
        patch("cofolder.modules.runners.base.metadata.version") as version,
    ):
        available, message = runner.check_availability()

    assert available is False
    assert "does not support Python 3.13" in message
    assert "Python 3.11 or 3.12" in message
    version.assert_not_called()


def test_boltz2_runner_load_options_keeps_model_runner_owned(temp_dir):
    options_path = temp_dir / "options.yaml"
    options_path.write_text(
        yaml.safe_dump(
            {
                "version": 1,
                "runtime": {"cache_path": "~/.boltz", "diffusion_samples": 1},
                "runner": {},
            }
        ),
        encoding="utf-8",
    )

    options = Boltz2Runner().load_options(options_path)

    assert options.find_value(key="model") is None
    assert Boltz2Runner.model_name == "boltz2"
