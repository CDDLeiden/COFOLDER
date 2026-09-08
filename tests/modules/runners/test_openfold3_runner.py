"""Tests for the OpenFold3 runner."""

from __future__ import annotations

import json
from importlib.metadata import PackageNotFoundError
from pathlib import Path
from subprocess import CompletedProcess
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
from cofolder.modules.input.system import System
from cofolder.modules.runners.contracts import RunnerExecutionRequest
from cofolder.modules.runners.openfold3_runner import (
    OpenFold3Runner,
    check_openfold3_setup_ready,
    run_openfold3,
)
from cofolder.modules.runners.validators import validate_runner_bundle


def _make_system(*, smiles: str | None = "CCO", ccd: str | None = "ETH") -> System:
    ligand: dict[str, object] = {"id": "B"}
    if smiles is not None:
        ligand["smiles"] = smiles
    if ccd is not None:
        ligand["ccd"] = ccd

    return System(
        system={
            "sequences": [
                {
                    "protein": {
                        "id": "A",
                        "sequence": "MKRAAT",
                    }
                },
                {"ligand": ligand},
            ],
            "properties": [{"affinity": {"binder": "B"}}],
        }
    )


def test_openfold3_builds_nucleic_acids_and_translates_pocket_constraint():
    system = System(
        system={
            "sequences": [
                {"protein": {"id": "A", "sequence": "AC"}},
                {"dna": {"id": "D", "sequence": "AT"}},
                {"rna": {"id": ["R", "S"], "sequence": "GU"}},
                {"ligand": {"id": "L", "smiles": "CCO"}},
            ],
            "constraints": [
                {
                    "pocket": {
                        "binder": "L",
                        "contacts": [["A", 1], ["D", 2]],
                        "max_distance": 7.5,
                    }
                }
            ],
        }
    )

    payload = OpenFold3Runner()._build_query_payload("mixed", system)
    query = payload["queries"]["mixed"]

    assert query["chains"] == [
        {"molecule_type": "protein", "chain_ids": ["A"], "sequence": "AC"},
        {"molecule_type": "dna", "chain_ids": ["D"], "sequence": "AT"},
        {"molecule_type": "rna", "chain_ids": ["R", "S"], "sequence": "GU"},
        {"molecule_type": "ligand", "chain_ids": ["L"], "smiles": "CCO"},
    ]
    assert query["pocket_constraint"] == {
        "ligand_chain_id": "L",
        "pocket_residues": [["A", 1], ["D", 2]],
        "max_distance": 7.5,
    }
    assert OpenFold3Runner._chain_order(system) == ["A", "D", "R", "S", "L"]


def test_openfold3_payload_never_silently_omits_unsupported_constraints():
    system = _make_system(smiles="CCO", ccd=None)
    system.system["constraints"] = [
        {"bond": {"atom1": ["A", 1, "N"], "atom2": ["B", 1, "C1"]}}
    ]

    with pytest.raises(ValueError, match="bond.*unsupported"):
        OpenFold3Runner().validate_system(system, {}, check_atom_names=False)


def _write_options_yaml(temp_dir: Path, *, samples_per_seed: int = 2) -> Path:
    options_path = temp_dir / "openfold3-options.yaml"
    options_path.write_text(
        yaml.safe_dump(
            {
                "version": 1,
                "runtime": {
                    "cache_path": "/tmp/openfold3-cache",
                    "diffusion_samples": samples_per_seed,
                    "extra_args": ["--use-msa-server=False"],
                },
                "runner": {},
            }
        ),
        encoding="utf-8",
    )
    return options_path


def _make_request(temp_dir: Path, *, system_obj: System, samples_per_seed: int = 2) -> RunnerExecutionRequest:
    runner = OpenFold3Runner()
    options_path = _write_options_yaml(temp_dir, samples_per_seed=samples_per_seed)
    options_obj = runner.load_options(options_path)
    repeat_dir = temp_dir / "repeat_1"
    return RunnerExecutionRequest(
        runner_name="openfold3",
        system_name="system",
        system_path=temp_dir / "system.yaml",
        system_obj=system_obj,
        options_path=options_path,
        options_obj=options_obj,
        repeat=1,
        seed=123,
        seed_provenance=RepeatSeedProvenance(
            repeat_id=1,
            requested_base_seed=123,
            resolved_base_seed=123,
            derived_seed=123,
            effective_seed=123,
            origin=SeedOrigin.USER_SPECIFIED,
        ),
        backend_identity=RunnerBackendIdentity(
            runner_name="openfold3",
            backend_name="openfold3",
            version="0.3.0",
            version_status=BackendVersionStatus.DETECTED,
            raw_version="0.3.0",
        ),
        repeat_dir=repeat_dir,
        raw_dir=temp_dir,
        logger=None,
    )


def _write_openfold3_sample_outputs(
    output_dir: Path,
    *,
    system_name: str,
    seed: int,
    pair_map_style: str = "nested",
) -> None:
    seed_dir = output_dir / system_name / f"seed_{seed}"
    seed_dir.mkdir(parents=True, exist_ok=True)

    for sample_number in (1, 2):
        prefix = f"{system_name}_seed_{seed}_sample_{sample_number}"
        (seed_dir / f"{prefix}_model.cif").write_text("data_test", encoding="utf-8")
        (seed_dir / f"{prefix}_confidences.json").write_text(
            json.dumps(
                {
                    "plddt": [91.2, 88.4],
                    "pae": [[0.1, 0.2], [0.2, 0.1]],
                    "pde": [[0.3, 0.4], [0.4, 0.3]],
                }
            ),
            encoding="utf-8",
        )
        if pair_map_style == "nested":
            chain_pair_iptm = {
                "A": {"A": 0.91, "B": 0.61},
                "B": {"A": 0.61, "B": 0.88},
            }
            bespoke_iptm = {
                "A": {"B": 0.58},
                "B": {"A": 0.58},
            }
        elif pair_map_style == "tuple_keys":
            chain_pair_iptm = {
                "(A, B)": 0.61,
            }
            bespoke_iptm = {
                "(A, B)": 0.58,
            }
        else:
            raise ValueError(f"Unsupported pair_map_style {pair_map_style!r}")

        (seed_dir / f"{prefix}_confidences_aggregated.json").write_text(
            json.dumps(
                {
                    "ptm": 0.80,
                    "iptm": 0.70,
                    "avg_plddt": 0.79,
                    "gpde": 0.31,
                    "disorder": 0.14,
                    "has_clash": 0.0,
                    "sample_ranking_score": 0.83,
                    "chain_ptm": {"A": 0.85, "B": 0.65},
                    "chain_pair_iptm": chain_pair_iptm,
                    "bespoke_iptm": bespoke_iptm,
                }
            ),
            encoding="utf-8",
        )
    (seed_dir / "timing.json").write_text(json.dumps({"total_seconds": 42.0}), encoding="utf-8")


def test_run_openfold3_uses_documented_cli_flags(monkeypatch, temp_dir):
    options = OpenFold3Runner().load_options(_write_options_yaml(temp_dir, samples_per_seed=3))
    query_json_path = temp_dir / "query.json"
    runner_yaml_path = temp_dir / "runner.yaml"
    output_dir = temp_dir / "output"
    query_json_path.write_text("{}", encoding="utf-8")
    runner_yaml_path.write_text("experiment_settings:\n  seeds:\n    - 123\n", encoding="utf-8")

    captured: dict[str, object] = {}

    def _fake_subprocess_run(cmd, check, capture_output, text):
        captured["cmd"] = cmd
        return CompletedProcess(args=cmd, returncode=0, stdout="", stderr="")

    monkeypatch.setattr(
        "cofolder.modules.runners.openfold3_runner.subprocess.run",
        _fake_subprocess_run,
    )

    run_openfold3(
        query_json_path=query_json_path,
        runner_yaml_path=runner_yaml_path,
        output_dir=output_dir,
        diffusion_samples=3,
        options=options,
    )

    command = captured["cmd"]
    assert "--query-json=" + str(query_json_path) in command
    assert "--output-dir=" + str(output_dir) in command
    assert "--num-diffusion-samples=3" in command
    assert "--runner-yaml=" + str(runner_yaml_path) in command
    assert all("--config_json" not in token for token in command)
    assert all("--seed" not in token for token in command)


def test_openfold3_runner_reports_install_hint_when_missing():
    runner = OpenFold3Runner()

    with patch(
        "cofolder.modules.runners.base.metadata.version",
        side_effect=PackageNotFoundError,
    ):
        available, message = runner.check_availability()

    assert available is False
    assert 'python -m pip install -e ".[openfold3]"' in message


def test_openfold3_runner_check_availability_only_requires_package_install(monkeypatch):
    runner = OpenFold3Runner()

    monkeypatch.setattr(
        "cofolder.modules.runners.base.metadata.version",
        lambda distribution_name: "1.0.0",
    )

    available, message = runner.check_availability()

    assert available is True
    assert message is None


def test_openfold3_runner_reports_setup_script_when_cache_is_unprepared(monkeypatch, temp_dir):
    cache_root = temp_dir / "openfold3-cache"
    cache_root.mkdir(parents=True, exist_ok=True)
    available, message = check_openfold3_setup_ready({"OPENFOLD_CACHE": str(cache_root)})

    assert available is False
    assert "scripts/setup_openfold3.sh" in message
    assert "ckpt_root" in message


def test_check_openfold3_setup_ready_accepts_existing_checkpoint_root(temp_dir):
    cache_root = temp_dir / "openfold3-cache"
    checkpoint_root = cache_root / "checkpoints"
    checkpoint_root.mkdir(parents=True, exist_ok=True)
    (cache_root / "ckpt_root").write_text(str(checkpoint_root), encoding="utf-8")

    available, message = check_openfold3_setup_ready({"OPENFOLD_CACHE": str(cache_root)})

    assert available is True
    assert message is None


def test_openfold3_runner_prefers_ccd_over_smiles_in_query_json(monkeypatch, temp_dir):
    runner = OpenFold3Runner()
    request = _make_request(temp_dir, system_obj=_make_system(smiles="CCO", ccd="ETH"))

    def _fake_run_openfold3(*, query_json_path, runner_yaml_path, output_dir, diffusion_samples, options, timings, label_prefix):
        query_payload = json.loads(query_json_path.read_text(encoding="utf-8"))
        assert query_payload == {
            "queries": {
                "system": {
                    "chains": [
                        {
                            "molecule_type": "protein",
                            "chain_ids": ["A"],
                            "sequence": "MKRAAT",
                        },
                        {
                            "molecule_type": "ligand",
                            "chain_ids": ["B"],
                            "ccd_codes": ["ETH"],
                        },
                    ]
                }
            }
        }
        runner_settings = yaml.safe_load(runner_yaml_path.read_text(encoding="utf-8"))
        assert runner_settings["experiment_settings"]["seeds"] == [123]
        assert diffusion_samples == 2
        _write_openfold3_sample_outputs(output_dir, system_name="system", seed=123)

    monkeypatch.setattr(
        "cofolder.modules.runners.openfold3_runner.run_openfold3",
        _fake_run_openfold3,
    )

    result = runner.run(request)
    bundle = validate_runner_bundle(
        result,
        requested_metric_groups={"confidence_metrics"},
    )

    system_df = pd.read_csv(result.system_metrics_path)
    chain_df = pd.read_csv(result.chain_metrics_path)
    manifest = json.loads(result.manifest_path.read_text(encoding="utf-8"))

    assert {"ptm", "iptm", "avg_plddt", "gpde", "disorder", "has_clash", "sample_ranking_score"}.issubset(
        system_df.columns
    )
    assert "confidence_score" not in system_df.columns
    assert {"chain_ptm", "chain_pair_iptm_A_B", "bespoke_iptm_A_B"}.issubset(chain_df.columns)
    assert result.metric_outcomes["confidence_metrics"].state == "computed"
    assert result.metric_outcomes["confidence_metrics"].required_artifacts == (
        "companion_artifacts.plddt",
        "companion_artifacts.pae",
        "companion_artifacts.pde",
    )
    assert result.metric_outcomes["affinity_metrics"].state == "unsupported"
    assert result.metric_outcomes["affinity_metrics_ext"].state == "unsupported"
    assert [artifact.label for artifact in result.companion_artifacts] == ["plddt", "pae", "pde"]
    assert (result.normalized_dir / "artifacts" / "plddt" / "system_seed_123_sample_1_plddt.json").exists()
    assert (result.normalized_dir / "artifacts" / "pae" / "system_seed_123_sample_2_pae.json").exists()
    assert bundle.metric_outcomes["confidence_metrics"].state == "computed"
    assert manifest["runtime_context"]["timing_path"].endswith("/seed_123/timing.json")
    assert manifest["backend"]["version"] == "0.3.0"
    assert manifest["seed"]["effective_seed"] == 123
    assert manifest["seed"]["origin"] == "user_specified"
    assert "companion_artifacts" in manifest
    assert all(entry["label"] != "timing" for entry in manifest["companion_artifacts"])


def test_openfold3_runner_uses_ccd_codes_when_smiles_is_absent(monkeypatch, temp_dir):
    runner = OpenFold3Runner()
    request = _make_request(temp_dir, system_obj=_make_system(smiles=None, ccd="ATP"), samples_per_seed=1)

    def _fake_run_openfold3(*, query_json_path, runner_yaml_path, output_dir, diffusion_samples, options, timings, label_prefix):
        query_payload = json.loads(query_json_path.read_text(encoding="utf-8"))
        assert query_payload["queries"]["system"]["chains"][1] == {
            "molecule_type": "ligand",
            "chain_ids": ["B"],
            "ccd_codes": ["ATP"],
        }
        _write_openfold3_sample_outputs(output_dir, system_name="system", seed=123)

    monkeypatch.setattr(
        "cofolder.modules.runners.openfold3_runner.run_openfold3",
        _fake_run_openfold3,
    )

    result = runner.run(request)

    assert result.runtime.model_name == "openfold3"
    assert result.runtime.cache_path == "/tmp/openfold3-cache"
    assert result.diffusion_samples == 2


def test_openfold3_runner_accepts_tuple_key_pair_confidence_maps(
    monkeypatch,
    temp_dir,
):
    runner = OpenFold3Runner()
    request = _make_request(temp_dir, system_obj=_make_system(smiles="CCO", ccd="ETH"))

    def _fake_run_openfold3(*, query_json_path, runner_yaml_path, output_dir, diffusion_samples, options, timings, label_prefix):
        _write_openfold3_sample_outputs(
            output_dir,
            system_name="system",
            seed=123,
            pair_map_style="tuple_keys",
        )

    monkeypatch.setattr(
        "cofolder.modules.runners.openfold3_runner.run_openfold3",
        _fake_run_openfold3,
    )

    result = runner.run(request)
    bundle = validate_runner_bundle(
        result,
        requested_metric_groups={"confidence_metrics"},
    )
    chain_df = pd.read_csv(result.chain_metrics_path)

    assert result.metric_outcomes["confidence_metrics"].state == "computed"
    assert "chain_pair_iptm_A_B" in chain_df.columns
    assert "bespoke_iptm_A_B" in chain_df.columns
    assert chain_df["chain_pair_iptm_A_B"].dropna().tolist() == [0.61, 0.61]
    assert chain_df["bespoke_iptm_A_B"].dropna().tolist() == [0.58, 0.58]
    assert bundle.metric_outcomes["confidence_metrics"].state == "computed"


def test_openfold3_runner_marks_confidence_failed_when_any_sample_is_incomplete(
    monkeypatch,
    temp_dir,
):
    runner = OpenFold3Runner()
    request = _make_request(temp_dir, system_obj=_make_system(smiles="CCO", ccd="ETH"))

    def _fake_run_openfold3(*, query_json_path, runner_yaml_path, output_dir, diffusion_samples, options, timings, label_prefix):
        _write_openfold3_sample_outputs(output_dir, system_name="system", seed=123)
        incomplete_path = (
            output_dir
            / "system"
            / "seed_123"
            / "system_seed_123_sample_2_confidences.json"
        )
        incomplete_path.unlink()

    monkeypatch.setattr(
        "cofolder.modules.runners.openfold3_runner.run_openfold3",
        _fake_run_openfold3,
    )

    result = runner.run(request)

    assert result.metric_outcomes["confidence_metrics"].state == "failed"
    assert "sample 2 is missing confidences.json" in (result.metric_outcomes["confidence_metrics"].message or "")
