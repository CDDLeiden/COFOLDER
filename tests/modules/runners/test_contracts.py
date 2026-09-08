"""Tests for the shared runner contract surface."""

from __future__ import annotations

from importlib import metadata
from pathlib import Path

import pytest

from cofolder.modules.contracts import (
    BackendVersionStatus,
    RepeatSeedProvenance,
    RunnerBackendIdentity,
    SeedOrigin,
)
from cofolder.modules.runners.base import BaseRunner
from cofolder.modules.runners.boltz1_runner import Boltz1Runner
from cofolder.modules.runners.boltz2_runner import Boltz2Runner
from cofolder.modules.runners.boltz_community_runner import BoltzCommunityRunner
from cofolder.modules.runners.contracts import (
    RunnerCompanionArtifact,
    RunnerExecutionRequest,
    RunnerExecutionResult,
    RunnerMetricOutcome,
    RunnerNormalizedBundle,
    RunnerPreparationResult,
    RunnerRuntime,
    merge_runner_runtime,
)
from cofolder.modules.runners.openfold3_runner import OpenFold3Runner


class _SimpleRunner(BaseRunner):
    name = "simple"
    capabilities = {"confidence_metrics"}

    def check_availability(self) -> tuple[bool, str | None]:
        return True, None

    def load_options(self, options_path: Path):
        return {"options_path": str(options_path)}

    def run(self, request: RunnerExecutionRequest) -> RunnerExecutionResult:
        return RunnerExecutionResult(
            runner_name=self.name,
            raw_output_dir=request.repeat_dir,
            normalized_dir=request.repeat_dir / "normalized",
            structures_dir=request.repeat_dir / "normalized" / "structures",
            system_metrics_path=request.repeat_dir / "normalized" / "system_metrics.csv",
            chain_metrics_path=request.repeat_dir / "normalized" / "chain_metrics.csv",
            manifest_path=request.repeat_dir / "normalized" / "manifest.json",
        )


def test_base_runner_prepare_system_defaults_to_noop():
    runner = _SimpleRunner()
    system_obj = object()
    options_obj = {"hello": "world"}

    preparation = runner.prepare_system(
        system_obj=system_obj,
        options_obj=options_obj,
        wrk_dir=Path("/tmp/wrk"),
        conformers=None,
        sdf_file=None,
        logger=None,
    )

    assert isinstance(preparation, RunnerPreparationResult)
    assert preparation.system_obj is system_obj
    assert preparation.options_obj is options_obj
    assert preparation.warnings == []
    assert preparation.runtime == RunnerRuntime()


def test_backend_detection_normalizes_distribution_version(monkeypatch):
    monkeypatch.setattr(
        "cofolder.modules.runners.base.metadata.version", lambda name: "1.0RC1"
    )
    identity = _SimpleRunner().detect_backend_identity()
    assert identity.version == "1.0rc1"
    assert identity.version_status == BackendVersionStatus.DETECTED
    assert identity.raw_version == "1.0RC1"


@pytest.mark.parametrize(
    ("runner", "backend_name", "distribution"),
    [
        (Boltz1Runner(), "boltz", "boltz"),
        (Boltz2Runner(), "boltz", "boltz"),
        (BoltzCommunityRunner(), "boltz-community", "boltz-community"),
        (OpenFold3Runner(), "openfold3", "openfold3"),
    ],
)
def test_selectable_runner_backend_declarations(runner, backend_name, distribution):
    assert runner.backend_name == backend_name
    assert runner.backend_distribution == distribution


def test_backend_detection_records_unavailable_without_raising(monkeypatch):
    def missing(name):
        raise metadata.PackageNotFoundError(name)

    monkeypatch.setattr("cofolder.modules.runners.base.metadata.version", missing)
    identity = _SimpleRunner().detect_backend_identity()
    assert identity.version is None
    assert identity.raw_version is None
    assert identity.version_status == BackendVersionStatus.UNAVAILABLE


def test_backend_detection_sanitizes_unexpected_lookup_error(monkeypatch):
    def broken(name):
        raise RuntimeError("secret filesystem detail")

    monkeypatch.setattr("cofolder.modules.runners.base.metadata.version", broken)
    identity = _SimpleRunner().detect_backend_identity()
    assert identity.version_status == BackendVersionStatus.UNAVAILABLE
    assert identity.detail == "Version lookup failed (RuntimeError)."
    assert "secret" not in identity.detail


@pytest.mark.parametrize("raw_version", ["", "not a version"])
def test_backend_detection_preserves_unparseable_value(monkeypatch, raw_version):
    monkeypatch.setattr(
        "cofolder.modules.runners.base.metadata.version", lambda name: raw_version
    )
    identity = _SimpleRunner().detect_backend_identity()
    assert identity.version is None
    assert identity.raw_version == raw_version
    assert identity.version_status == BackendVersionStatus.UNPARSEABLE


def test_runtime_metadata_is_typed_and_explicit():
    runtime = RunnerRuntime(
        cache_path="~/.boltz",
        diffusion_samples=4,
        model_name="boltz2",
    )

    assert runtime.cache_path == "~/.boltz"
    assert runtime.diffusion_samples == 4
    assert runtime.model_name == "boltz2"


def test_execution_request_carries_typed_runtime_metadata():
    runtime = RunnerRuntime(diffusion_samples=2)
    request = RunnerExecutionRequest(
        runner_name="simple",
        system_name="system",
        system_path=Path("system.yaml"),
        system_obj=object(),
        options_path=Path("options.yaml"),
        options_obj={},
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
            runner_name="simple",
            backend_name="simple",
            version="1.0",
            version_status=BackendVersionStatus.DETECTED,
            raw_version="1.0",
        ),
        repeat_dir=Path("repeat_1"),
        raw_dir=Path("raw"),
        logger=None,
        runtime=runtime,
    )

    assert request.runtime is runtime
    assert request.runtime.diffusion_samples == 2


def test_runner_preparation_result_no_longer_accepts_legacy_runtime_context():
    with pytest.raises(TypeError, match="unexpected keyword argument 'runtime_context'"):
        RunnerPreparationResult(
            system_obj=object(),
            options_obj={},
            runtime_context={"diffusion_samples": 3},
        )


def test_runner_execution_result_no_longer_accepts_legacy_runtime_context():
    with pytest.raises(TypeError, match="unexpected keyword argument 'runtime_context'"):
        RunnerExecutionResult(
            runner_name="simple",
            raw_output_dir=Path("raw"),
            normalized_dir=Path("normalized"),
            structures_dir=Path("normalized/structures"),
            system_metrics_path=Path("normalized/system_metrics.csv"),
            chain_metrics_path=Path("normalized/chain_metrics.csv"),
            manifest_path=Path("normalized/manifest.json"),
            runtime_context={"diffusion_samples": 4},
        )


def test_merge_runner_runtime_preserves_later_explicit_single_sample_value():
    merged = merge_runner_runtime(
        RunnerRuntime(
            diffusion_samples=9,
            cache_path="/prep/cache",
            model_name="prep-model",
        ),
        RunnerRuntime(
            diffusion_samples=1,
            model_name="run-model",
        ),
    )

    assert merged.diffusion_samples == 1
    assert merged.cache_path == "/prep/cache"
    assert merged.model_name == "run-model"


def test_execution_result_exposes_authoritative_normalized_bundle():
    result = RunnerExecutionResult(
        runner_name="simple",
        raw_output_dir=Path("raw"),
        normalized_dir=Path("normalized"),
        structures_dir=Path("normalized/structures"),
        system_metrics_path=Path("normalized/system_metrics.csv"),
        chain_metrics_path=Path("normalized/chain_metrics.csv"),
        manifest_path=Path("normalized/manifest.json"),
        capabilities={"confidence_metrics"},
        runtime=RunnerRuntime(cache_path="~/.boltz", diffusion_samples=2, model_name="boltz2"),
        sample_records=[{"repeat": 1, "diffusion_sample": 0, "cif_file": "1_system_model_0.cif"}],
        metric_outcomes={
            "confidence_metrics": RunnerMetricOutcome(
                state="computed",
                required_artifacts=("companion_artifacts.pae",),
            ),
        },
        companion_artifacts=[
            RunnerCompanionArtifact(
                label="pae",
                relative_path="artifacts/pae.npy",
                kind="matrix",
            )
        ],
    )

    assert isinstance(result.normalized_bundle, RunnerNormalizedBundle)
    assert result.normalized_bundle.runner_name == "simple"
    assert result.normalized_bundle.capabilities == {"confidence_metrics"}
    assert result.metric_outcomes["confidence_metrics"].state == "computed"
    assert result.metric_outcomes["confidence_metrics"].required_artifacts == (
        "companion_artifacts.pae",
    )
    assert result.companion_artifacts[0].label == "pae"
    assert result.normalized_bundle.companion_artifacts[0].relative_path == "artifacts/pae.npy"
