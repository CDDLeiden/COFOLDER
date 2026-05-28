"""Tests for normalized runner bundle validation."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from cofolder.modules.runners.contracts import (
    RunnerCompanionArtifact,
    RunnerExecutionResult,
    RunnerMetricOutcome,
    RunnerRuntime,
)
from cofolder.modules.runners.validators import (
    RunnerBundleValidationError,
    validate_runner_bundle,
)


def _write_bundle_files(
    repeat_dir: Path,
    *,
    runner_name: str = "simple",
    include_confidence: bool = True,
    include_affinity: bool = True,
    capabilities: set[str] | None = None,
    system_metric_overrides: dict[str, object] | None = None,
    chain_metric_overrides: dict[str, object] | None = None,
    companion_artifacts: list[RunnerCompanionArtifact] | None = None,
) -> tuple[Path, Path, Path, Path]:
    normalized_dir = repeat_dir / "normalized"
    structures_dir = normalized_dir / "structures"
    structures_dir.mkdir(parents=True, exist_ok=True)
    (structures_dir / "1_system_model_0.cif").write_text("data", encoding="utf-8")

    system_row = {
        "cif_file": "1_system_model_0.cif",
        "model_name": "system",
        "repeat": 1,
        "diffusion_sample": 0,
    }
    if include_confidence:
        system_row.update(
            {
                "ptm": 0.8,
                "iptm": 0.7,
                "confidence_score": 0.9,
            }
        )

    chain_row = {
        "conf_chain_id": 0,
        "cif_file": "1_system_model_0.cif",
        "model_name": "system",
        "repeat": 1,
        "diffusion_sample": 0,
    }
    if include_confidence:
        chain_row["chains_ptm"] = 0.85
    if include_affinity:
        chain_row.update(
            {
                "affinity_pred_value": 6.1,
                "affinity_probability_binary": 0.8,
                "pIC50": 5.0,
                "IC50_M": 1e-6,
                "pIC50_kcal_per_mol": 7.0,
            }
        )
    if system_metric_overrides:
        system_row.update(system_metric_overrides)
    if chain_metric_overrides:
        chain_row.update(chain_metric_overrides)

    system_metrics_path = normalized_dir / "system_metrics.csv"
    chain_metrics_path = normalized_dir / "chain_metrics.csv"
    manifest_path = normalized_dir / "manifest.json"

    pd.DataFrame([system_row]).to_csv(system_metrics_path, index=False)
    pd.DataFrame([chain_row]).to_csv(chain_metrics_path, index=False)
    manifest_path.write_text(
        json.dumps(
            {
                "runner": runner_name,
                "capabilities": sorted(capabilities or []),
                "companion_artifacts": [
                    {
                        "label": artifact.label,
                        "relative_path": artifact.relative_path,
                        **({"kind": artifact.kind} if artifact.kind is not None else {}),
                        **(
                            {"description": artifact.description}
                            if artifact.description is not None
                            else {}
                        ),
                    }
                    for artifact in (companion_artifacts or [])
                ],
            }
        ),
        encoding="utf-8",
    )

    for artifact in companion_artifacts or []:
        artifact_path = normalized_dir / artifact.relative_path
        artifact_path.parent.mkdir(parents=True, exist_ok=True)
        artifact_path.write_text("artifact", encoding="utf-8")

    return normalized_dir, structures_dir, system_metrics_path, chain_metrics_path


def test_validate_runner_bundle_accepts_valid_explicit_outcomes(temp_dir):
    repeat_dir = temp_dir / "repeat_1"
    normalized_dir, structures_dir, system_metrics_path, chain_metrics_path = _write_bundle_files(
        repeat_dir,
        capabilities={"confidence_metrics", "affinity_metrics", "affinity_metrics_ext"},
    )
    manifest_path = normalized_dir / "manifest.json"

    result = RunnerExecutionResult(
        runner_name="simple",
        raw_output_dir=repeat_dir,
        normalized_dir=normalized_dir,
        structures_dir=structures_dir,
        system_metrics_path=system_metrics_path,
        chain_metrics_path=chain_metrics_path,
        manifest_path=manifest_path,
        capabilities={"confidence_metrics", "affinity_metrics", "affinity_metrics_ext"},
        runtime=RunnerRuntime(diffusion_samples=1),
        sample_records=[{"repeat": 1, "diffusion_sample": 0, "cif_file": "1_system_model_0.cif"}],
        metric_outcomes={
            "confidence_metrics": RunnerMetricOutcome(state="computed"),
            "affinity_metrics": RunnerMetricOutcome(state="computed"),
            "affinity_metrics_ext": RunnerMetricOutcome(state="computed"),
        },
    )

    bundle = validate_runner_bundle(
        result,
        requested_metric_groups={"confidence_metrics", "affinity_metrics", "affinity_metrics_ext"},
    )

    assert bundle.metric_outcomes["confidence_metrics"].state == "computed"
    assert bundle.metric_outcomes["affinity_metrics"].state == "computed"
    assert bundle.metric_outcomes["confidence_metrics"].required_columns == ()
    assert bundle.metric_outcomes["confidence_metrics"].required_artifacts == ()


def test_validate_runner_bundle_preserves_runner_authored_confidence_payload_metadata(temp_dir):
    repeat_dir = temp_dir / "repeat_1"
    artifacts = [
        RunnerCompanionArtifact(label="pae", relative_path="artifacts/pae.npy"),
    ]
    normalized_dir, structures_dir, system_metrics_path, chain_metrics_path = _write_bundle_files(
        repeat_dir,
        capabilities={"confidence_metrics"},
        system_metric_overrides={"sample_ranking_score": 0.83},
        companion_artifacts=artifacts,
    )
    manifest_path = normalized_dir / "manifest.json"

    result = RunnerExecutionResult(
        runner_name="simple",
        raw_output_dir=repeat_dir,
        normalized_dir=normalized_dir,
        structures_dir=structures_dir,
        system_metrics_path=system_metrics_path,
        chain_metrics_path=chain_metrics_path,
        manifest_path=manifest_path,
        capabilities={"confidence_metrics"},
        runtime=RunnerRuntime(diffusion_samples=1),
        metric_outcomes={
            "confidence_metrics": RunnerMetricOutcome(
                state="computed",
                required_columns=("system_metrics.sample_ranking_score",),
                required_artifacts=("companion_artifacts.pae",),
            ),
        },
        companion_artifacts=artifacts,
    )

    bundle = validate_runner_bundle(
        result,
        requested_metric_groups={"confidence_metrics"},
    )

    assert bundle.metric_outcomes["confidence_metrics"].state == "computed"
    assert bundle.metric_outcomes["confidence_metrics"].required_columns == (
        "system_metrics.sample_ranking_score",
    )
    assert bundle.metric_outcomes["confidence_metrics"].required_artifacts == (
        "companion_artifacts.pae",
    )


def test_validate_runner_bundle_rejects_runner_authored_computed_payload_when_declared_entries_are_missing(
    temp_dir,
):
    repeat_dir = temp_dir / "repeat_1"
    normalized_dir, structures_dir, system_metrics_path, chain_metrics_path = _write_bundle_files(
        repeat_dir,
        capabilities={"confidence_metrics"},
    )
    manifest_path = normalized_dir / "manifest.json"

    result = RunnerExecutionResult(
        runner_name="simple",
        raw_output_dir=repeat_dir,
        normalized_dir=normalized_dir,
        structures_dir=structures_dir,
        system_metrics_path=system_metrics_path,
        chain_metrics_path=chain_metrics_path,
        manifest_path=manifest_path,
        capabilities={"confidence_metrics"},
        runtime=RunnerRuntime(diffusion_samples=1),
        metric_outcomes={
            "confidence_metrics": RunnerMetricOutcome(
                state="computed",
                required_columns=("system_metrics.sample_ranking_score",),
                required_artifacts=("companion_artifacts.pae",),
            ),
        },
    )

    with pytest.raises(
        RunnerBundleValidationError,
        match="runner-declared normalized payload entries are missing",
    ):
        validate_runner_bundle(
            result,
            requested_metric_groups={"confidence_metrics"},
        )


def test_validate_runner_bundle_accepts_runner_authored_confidence_surface(temp_dir):
    repeat_dir = temp_dir / "repeat_1"
    artifacts = [
        RunnerCompanionArtifact(label="plddt", relative_path="artifacts/plddt.npy"),
        RunnerCompanionArtifact(label="pae", relative_path="artifacts/pae.npy"),
        RunnerCompanionArtifact(label="pde", relative_path="artifacts/pde.npy"),
        RunnerCompanionArtifact(label="timing", relative_path="runtime/timing.json"),
    ]
    normalized_dir, structures_dir, system_metrics_path, chain_metrics_path = _write_bundle_files(
        repeat_dir,
        runner_name="openfold3",
        capabilities={"confidence_metrics"},
        system_metric_overrides={
            "avg_plddt": 0.79,
            "gpde": 0.31,
            "disorder": 0.14,
            "has_clash": 0.0,
            "sample_ranking_score": 0.83,
        },
        chain_metric_overrides={
            "chain_ptm": 0.85,
            "chain_pair_iptm_0_0": 0.72,
            "bespoke_iptm_0_0": 0.69,
        },
        companion_artifacts=artifacts,
    )
    system_df = pd.read_csv(system_metrics_path).drop(columns=["ptm", "iptm", "confidence_score"])
    chain_df = pd.read_csv(chain_metrics_path).drop(columns=["chains_ptm"])
    system_df.to_csv(system_metrics_path, index=False)
    chain_df.to_csv(chain_metrics_path, index=False)
    manifest_path = normalized_dir / "manifest.json"

    result = RunnerExecutionResult(
        runner_name="openfold3",
        raw_output_dir=repeat_dir,
        normalized_dir=normalized_dir,
        structures_dir=structures_dir,
        system_metrics_path=system_metrics_path,
        chain_metrics_path=chain_metrics_path,
        manifest_path=manifest_path,
        capabilities={"confidence_metrics"},
        runtime=RunnerRuntime(diffusion_samples=1),
        metric_outcomes={
            "confidence_metrics": RunnerMetricOutcome(
                state="computed",
                required_columns=(
                    "system_metrics.avg_plddt",
                    "system_metrics.gpde",
                    "system_metrics.disorder",
                    "system_metrics.has_clash",
                    "system_metrics.sample_ranking_score",
                    "chain_metrics.chain_ptm",
                    "chain_metrics.chain_pair_iptm_0_0",
                    "chain_metrics.bespoke_iptm_0_0",
                ),
                required_artifacts=(
                    "companion_artifacts.plddt",
                    "companion_artifacts.pae",
                    "companion_artifacts.pde",
                ),
            ),
        },
        companion_artifacts=artifacts,
    )

    bundle = validate_runner_bundle(
        result,
        requested_metric_groups={"confidence_metrics"},
    )

    assert bundle.metric_outcomes["confidence_metrics"].required_columns == (
        "system_metrics.avg_plddt",
        "system_metrics.gpde",
        "system_metrics.disorder",
        "system_metrics.has_clash",
        "system_metrics.sample_ranking_score",
        "chain_metrics.chain_ptm",
        "chain_metrics.chain_pair_iptm_0_0",
        "chain_metrics.bespoke_iptm_0_0",
    )
    assert bundle.metric_outcomes["confidence_metrics"].required_artifacts == (
        "companion_artifacts.plddt",
        "companion_artifacts.pae",
        "companion_artifacts.pde",
    )


def test_validate_runner_bundle_uses_runner_authored_confidence_missing_details(temp_dir):
    repeat_dir = temp_dir / "repeat_1"
    normalized_dir, structures_dir, system_metrics_path, chain_metrics_path = _write_bundle_files(
        repeat_dir,
        capabilities={"confidence_metrics"},
    )
    manifest_path = normalized_dir / "manifest.json"

    result = RunnerExecutionResult(
        runner_name="simple",
        raw_output_dir=repeat_dir,
        normalized_dir=normalized_dir,
        structures_dir=structures_dir,
        system_metrics_path=system_metrics_path,
        chain_metrics_path=chain_metrics_path,
        manifest_path=manifest_path,
        capabilities={"confidence_metrics"},
        runtime=RunnerRuntime(diffusion_samples=1),
        metric_outcomes={
            "confidence_metrics": RunnerMetricOutcome(
                state="missing",
                required_columns=("system_metrics.avg_plddt",),
                required_artifacts=("companion_artifacts.pae",),
                message="OpenFold3 confidence outputs were incomplete.",
            ),
        },
    )

    with pytest.raises(
        RunnerBundleValidationError,
        match="system_metrics.avg_plddt.*companion_artifacts.pae",
    ):
        validate_runner_bundle(
            result,
            requested_metric_groups={"confidence_metrics"},
        )


def test_validate_runner_bundle_rejects_computed_affinity_without_required_payload(temp_dir):
    repeat_dir = temp_dir / "repeat_1"
    normalized_dir, structures_dir, system_metrics_path, chain_metrics_path = _write_bundle_files(
        repeat_dir,
        include_affinity=False,
        capabilities={"affinity_metrics"},
    )
    manifest_path = normalized_dir / "manifest.json"

    result = RunnerExecutionResult(
        runner_name="simple",
        raw_output_dir=repeat_dir,
        normalized_dir=normalized_dir,
        structures_dir=structures_dir,
        system_metrics_path=system_metrics_path,
        chain_metrics_path=chain_metrics_path,
        manifest_path=manifest_path,
        capabilities={"affinity_metrics"},
        runtime=RunnerRuntime(diffusion_samples=1),
        metric_outcomes={
            "affinity_metrics": RunnerMetricOutcome(state="computed"),
        },
    )

    with pytest.raises(RunnerBundleValidationError, match="marked 'computed'"):
        validate_runner_bundle(
            result,
            requested_metric_groups={"affinity_metrics"},
        )


def test_validate_runner_bundle_rejects_supported_group_marked_unsupported(temp_dir):
    repeat_dir = temp_dir / "repeat_1"
    normalized_dir, structures_dir, system_metrics_path, chain_metrics_path = _write_bundle_files(
        repeat_dir,
        capabilities={"confidence_metrics"},
    )
    manifest_path = normalized_dir / "manifest.json"

    result = RunnerExecutionResult(
        runner_name="simple",
        raw_output_dir=repeat_dir,
        normalized_dir=normalized_dir,
        structures_dir=structures_dir,
        system_metrics_path=system_metrics_path,
        chain_metrics_path=chain_metrics_path,
        manifest_path=manifest_path,
        capabilities={"confidence_metrics"},
        runtime=RunnerRuntime(diffusion_samples=1),
        metric_outcomes={
            "confidence_metrics": RunnerMetricOutcome(state="unsupported"),
        },
    )

    with pytest.raises(RunnerBundleValidationError, match="cannot be marked 'unsupported'"):
        validate_runner_bundle(
            result,
            requested_metric_groups={"confidence_metrics"},
        )


def test_validate_runner_bundle_rejects_unsupported_confidence_when_canonical_payload_is_present(
    temp_dir,
):
    repeat_dir = temp_dir / "repeat_1"
    normalized_dir, structures_dir, system_metrics_path, chain_metrics_path = _write_bundle_files(
        repeat_dir,
        capabilities=set(),
    )
    manifest_path = normalized_dir / "manifest.json"
    manifest_path.write_text(
        json.dumps({"runner": "simple", "capabilities": [], "companion_artifacts": []}),
        encoding="utf-8",
    )

    result = RunnerExecutionResult(
        runner_name="simple",
        raw_output_dir=repeat_dir,
        normalized_dir=normalized_dir,
        structures_dir=structures_dir,
        system_metrics_path=system_metrics_path,
        chain_metrics_path=chain_metrics_path,
        manifest_path=manifest_path,
        capabilities=set(),
        runtime=RunnerRuntime(diffusion_samples=1),
        metric_outcomes={
            "confidence_metrics": RunnerMetricOutcome(state="unsupported"),
        },
    )

    with pytest.raises(
        RunnerBundleValidationError,
        match="cannot be marked 'unsupported' while confidence payload is present or declared",
    ):
        validate_runner_bundle(
            result,
            requested_metric_groups=set(),
        )


def test_validate_runner_bundle_rejects_unsupported_confidence_when_backend_native_payload_is_present(
    temp_dir,
):
    repeat_dir = temp_dir / "repeat_1"
    artifacts = [
        RunnerCompanionArtifact(label="plddt", relative_path="artifacts/plddt.json"),
    ]
    normalized_dir, structures_dir, system_metrics_path, chain_metrics_path = _write_bundle_files(
        repeat_dir,
        capabilities=set(),
        system_metric_overrides={
            "sample_ranking_score": 0.83,
            "avg_plddt": 0.79,
        },
        chain_metric_overrides={
            "chain_ptm": 0.85,
            "chain_pair_iptm_A_B": 0.61,
        },
        companion_artifacts=artifacts,
    )
    system_df = pd.read_csv(system_metrics_path).drop(columns=["ptm", "iptm", "confidence_score"])
    chain_df = pd.read_csv(chain_metrics_path).drop(columns=["chains_ptm"])
    system_df.to_csv(system_metrics_path, index=False)
    chain_df.to_csv(chain_metrics_path, index=False)
    manifest_path = normalized_dir / "manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "runner": "simple",
                "capabilities": [],
                "companion_artifacts": [
                    {"label": "plddt", "relative_path": "artifacts/plddt.json"}
                ],
            }
        ),
        encoding="utf-8",
    )

    result = RunnerExecutionResult(
        runner_name="simple",
        raw_output_dir=repeat_dir,
        normalized_dir=normalized_dir,
        structures_dir=structures_dir,
        system_metrics_path=system_metrics_path,
        chain_metrics_path=chain_metrics_path,
        manifest_path=manifest_path,
        capabilities=set(),
        runtime=RunnerRuntime(diffusion_samples=1),
        metric_outcomes={
            "confidence_metrics": RunnerMetricOutcome(state="unsupported"),
        },
        companion_artifacts=artifacts,
    )

    with pytest.raises(
        RunnerBundleValidationError,
        match="cannot be marked 'unsupported' while confidence payload is present or declared",
    ):
        validate_runner_bundle(
            result,
            requested_metric_groups=set(),
        )


def test_validate_runner_bundle_rejects_unsupported_group_when_payload_columns_exist(temp_dir):
    repeat_dir = temp_dir / "repeat_1"
    normalized_dir, structures_dir, system_metrics_path, chain_metrics_path = _write_bundle_files(
        repeat_dir,
        capabilities={"confidence_metrics"},
    )
    manifest_path = normalized_dir / "manifest.json"

    result = RunnerExecutionResult(
        runner_name="simple",
        raw_output_dir=repeat_dir,
        normalized_dir=normalized_dir,
        structures_dir=structures_dir,
        system_metrics_path=system_metrics_path,
        chain_metrics_path=chain_metrics_path,
        manifest_path=manifest_path,
        capabilities={"confidence_metrics"},
        runtime=RunnerRuntime(diffusion_samples=1),
        metric_outcomes={
            "confidence_metrics": RunnerMetricOutcome(state="computed"),
            "affinity_metrics": RunnerMetricOutcome(state="unsupported"),
            "affinity_metrics_ext": RunnerMetricOutcome(state="unsupported"),
        },
    )

    with pytest.raises(RunnerBundleValidationError, match="while normalized payload is present"):
        validate_runner_bundle(
            result,
            requested_metric_groups={"confidence_metrics"},
        )


def test_validate_runner_bundle_requires_explicit_outcome_for_requested_group(temp_dir):
    repeat_dir = temp_dir / "repeat_1"
    normalized_dir, structures_dir, system_metrics_path, chain_metrics_path = _write_bundle_files(
        repeat_dir,
        include_confidence=False,
        include_affinity=False,
    )
    manifest_path = normalized_dir / "manifest.json"
    manifest_path.write_text(
        json.dumps({"runner": "simple", "capabilities": []}),
        encoding="utf-8",
    )

    result = RunnerExecutionResult(
        runner_name="simple",
        raw_output_dir=repeat_dir,
        normalized_dir=normalized_dir,
        structures_dir=structures_dir,
        system_metrics_path=system_metrics_path,
        chain_metrics_path=chain_metrics_path,
        manifest_path=manifest_path,
        capabilities=set(),
        runtime=RunnerRuntime(diffusion_samples=1),
    )

    with pytest.raises(RunnerBundleValidationError, match="omits an explicit outcome"):
        validate_runner_bundle(
            result,
            requested_metric_groups={"confidence_metrics"},
        )


def test_validate_runner_bundle_rejects_missing_structure_referenced_by_csv(temp_dir):
    repeat_dir = temp_dir / "repeat_1"
    normalized_dir, structures_dir, system_metrics_path, chain_metrics_path = _write_bundle_files(
        repeat_dir,
        capabilities={"confidence_metrics"},
    )
    manifest_path = normalized_dir / "manifest.json"
    structure_path = structures_dir / "1_system_model_0.cif"
    structure_path.unlink()
    (structures_dir / "different_model_0.cif").write_text("data", encoding="utf-8")

    result = RunnerExecutionResult(
        runner_name="simple",
        raw_output_dir=repeat_dir,
        normalized_dir=normalized_dir,
        structures_dir=structures_dir,
        system_metrics_path=system_metrics_path,
        chain_metrics_path=chain_metrics_path,
        manifest_path=manifest_path,
        capabilities={"confidence_metrics"},
        runtime=RunnerRuntime(diffusion_samples=1),
    )

    with pytest.raises(RunnerBundleValidationError, match="reference missing structure files"):
        validate_runner_bundle(
            result,
            requested_metric_groups={"confidence_metrics"},
        )


def test_validate_runner_bundle_rejects_invalid_unrequested_explicit_outcome(temp_dir):
    repeat_dir = temp_dir / "repeat_1"
    normalized_dir, structures_dir, system_metrics_path, chain_metrics_path = _write_bundle_files(
        repeat_dir,
        capabilities={"confidence_metrics", "affinity_metrics"},
    )
    manifest_path = normalized_dir / "manifest.json"

    result = RunnerExecutionResult(
        runner_name="simple",
        raw_output_dir=repeat_dir,
        normalized_dir=normalized_dir,
        structures_dir=structures_dir,
        system_metrics_path=system_metrics_path,
        chain_metrics_path=chain_metrics_path,
        manifest_path=manifest_path,
        capabilities={"confidence_metrics", "affinity_metrics"},
        runtime=RunnerRuntime(diffusion_samples=1),
        metric_outcomes={
            "confidence_metrics": RunnerMetricOutcome(state="computed"),
            "affinity_metrics": RunnerMetricOutcome(state="unsupported"),
        },
    )

    with pytest.raises(RunnerBundleValidationError, match="cannot be marked 'unsupported'"):
        validate_runner_bundle(
            result,
            requested_metric_groups={"confidence_metrics"},
        )


def test_validate_runner_bundle_preserves_explicit_backend_native_score_columns(temp_dir):
    repeat_dir = temp_dir / "repeat_1"
    normalized_dir, structures_dir, system_metrics_path, chain_metrics_path = _write_bundle_files(
        repeat_dir,
        capabilities={"confidence_metrics"},
        system_metric_overrides={
            "sample_ranking_score": 0.83,
            "avg_plddt": 0.79,
            "has_clash": 0.0,
        },
        chain_metric_overrides={"chain_pair_iptm_A_B": 0.61},
    )
    manifest_path = normalized_dir / "manifest.json"

    result = RunnerExecutionResult(
        runner_name="simple",
        raw_output_dir=repeat_dir,
        normalized_dir=normalized_dir,
        structures_dir=structures_dir,
        system_metrics_path=system_metrics_path,
        chain_metrics_path=chain_metrics_path,
        manifest_path=manifest_path,
        capabilities={"confidence_metrics"},
        runtime=RunnerRuntime(diffusion_samples=1),
        metric_outcomes={
            "confidence_metrics": RunnerMetricOutcome(state="computed"),
        },
    )

    bundle = validate_runner_bundle(
        result,
        requested_metric_groups={"confidence_metrics"},
    )
    system_df = pd.read_csv(bundle.system_metrics_path)
    chain_df = pd.read_csv(bundle.chain_metrics_path)

    assert {"confidence_score", "sample_ranking_score", "avg_plddt", "has_clash"}.issubset(
        system_df.columns
    )
    assert "chain_pair_iptm_A_B" in chain_df.columns


def test_validate_runner_bundle_does_not_treat_timing_as_openfold3_confidence_artifact(
    temp_dir,
):
    repeat_dir = temp_dir / "repeat_1"
    artifacts = [
        RunnerCompanionArtifact(label="plddt", relative_path="artifacts/plddt.npy"),
        RunnerCompanionArtifact(label="pae", relative_path="artifacts/pae.npy"),
        RunnerCompanionArtifact(label="pde", relative_path="artifacts/pde.npy"),
        RunnerCompanionArtifact(label="timing", relative_path="runtime/timing.json"),
    ]
    normalized_dir, structures_dir, system_metrics_path, chain_metrics_path = _write_bundle_files(
        repeat_dir,
        runner_name="openfold3",
        capabilities={"confidence_metrics"},
        system_metric_overrides={
            "avg_plddt": 0.79,
            "gpde": 0.31,
            "disorder": 0.14,
            "has_clash": 0.0,
            "sample_ranking_score": 0.83,
        },
        chain_metric_overrides={
            "chain_ptm": 0.85,
            "chain_pair_iptm_0_0": 0.72,
            "bespoke_iptm_0_0": 0.69,
        },
        companion_artifacts=artifacts,
    )
    system_df = pd.read_csv(system_metrics_path).drop(columns=["confidence_score"])
    chain_df = pd.read_csv(chain_metrics_path).drop(columns=["chains_ptm"])
    system_df.to_csv(system_metrics_path, index=False)
    chain_df.to_csv(chain_metrics_path, index=False)
    manifest_path = normalized_dir / "manifest.json"

    result = RunnerExecutionResult(
        runner_name="openfold3",
        raw_output_dir=repeat_dir,
        normalized_dir=normalized_dir,
        structures_dir=structures_dir,
        system_metrics_path=system_metrics_path,
        chain_metrics_path=chain_metrics_path,
        manifest_path=manifest_path,
        capabilities={"confidence_metrics"},
        runtime=RunnerRuntime(diffusion_samples=1),
        metric_outcomes={
            "confidence_metrics": RunnerMetricOutcome(
                state="computed",
                required_artifacts=(
                    "companion_artifacts.plddt",
                    "companion_artifacts.pae",
                    "companion_artifacts.pde",
                ),
            ),
        },
        companion_artifacts=artifacts,
    )

    bundle = validate_runner_bundle(
        result,
        requested_metric_groups={"confidence_metrics"},
    )

    assert bundle.metric_outcomes["confidence_metrics"].required_artifacts == (
        "companion_artifacts.plddt",
        "companion_artifacts.pae",
        "companion_artifacts.pde",
    )


def test_validate_runner_bundle_accepts_companion_artifacts(temp_dir):
    repeat_dir = temp_dir / "repeat_1"
    artifacts = [
        RunnerCompanionArtifact(
            label="pae",
            relative_path="artifacts/pae.npy",
            kind="matrix",
            description="Predicted aligned error matrix.",
        ),
        RunnerCompanionArtifact(
            label="plddt",
            relative_path="artifacts/plddt.npy",
            kind="per_atom",
        ),
    ]
    normalized_dir, structures_dir, system_metrics_path, chain_metrics_path = _write_bundle_files(
        repeat_dir,
        capabilities={"confidence_metrics"},
        companion_artifacts=artifacts,
    )
    manifest_path = normalized_dir / "manifest.json"

    result = RunnerExecutionResult(
        runner_name="simple",
        raw_output_dir=repeat_dir,
        normalized_dir=normalized_dir,
        structures_dir=structures_dir,
        system_metrics_path=system_metrics_path,
        chain_metrics_path=chain_metrics_path,
        manifest_path=manifest_path,
        capabilities={"confidence_metrics"},
        runtime=RunnerRuntime(diffusion_samples=1),
        metric_outcomes={
            "confidence_metrics": RunnerMetricOutcome(state="computed"),
        },
        companion_artifacts=artifacts,
    )

    bundle = validate_runner_bundle(
        result,
        requested_metric_groups={"confidence_metrics"},
    )

    assert [artifact.label for artifact in bundle.companion_artifacts] == ["pae", "plddt"]


def test_validate_runner_bundle_accepts_companion_artifacts_when_manifest_order_differs(
    temp_dir,
):
    repeat_dir = temp_dir / "repeat_1"
    artifacts = [
        RunnerCompanionArtifact(label="pae", relative_path="artifacts/pae.npy"),
        RunnerCompanionArtifact(label="plddt", relative_path="artifacts/plddt.npy"),
    ]
    normalized_dir, structures_dir, system_metrics_path, chain_metrics_path = _write_bundle_files(
        repeat_dir,
        capabilities={"confidence_metrics"},
        companion_artifacts=artifacts,
    )
    manifest_path = normalized_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["companion_artifacts"] = list(reversed(manifest["companion_artifacts"]))
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    result = RunnerExecutionResult(
        runner_name="simple",
        raw_output_dir=repeat_dir,
        normalized_dir=normalized_dir,
        structures_dir=structures_dir,
        system_metrics_path=system_metrics_path,
        chain_metrics_path=chain_metrics_path,
        manifest_path=manifest_path,
        capabilities={"confidence_metrics"},
        runtime=RunnerRuntime(diffusion_samples=1),
        metric_outcomes={
            "confidence_metrics": RunnerMetricOutcome(state="computed"),
        },
        companion_artifacts=artifacts,
    )

    bundle = validate_runner_bundle(
        result,
        requested_metric_groups={"confidence_metrics"},
    )

    assert [artifact.label for artifact in bundle.companion_artifacts] == ["pae", "plddt"]


def test_validate_runner_bundle_rejects_missing_companion_artifact(temp_dir):
    repeat_dir = temp_dir / "repeat_1"
    artifacts = [
        RunnerCompanionArtifact(
            label="pae",
            relative_path="artifacts/pae.npy",
            kind="matrix",
        )
    ]
    normalized_dir, structures_dir, system_metrics_path, chain_metrics_path = _write_bundle_files(
        repeat_dir,
        capabilities={"confidence_metrics"},
        companion_artifacts=artifacts,
    )
    (normalized_dir / "artifacts" / "pae.npy").unlink()
    manifest_path = normalized_dir / "manifest.json"

    result = RunnerExecutionResult(
        runner_name="simple",
        raw_output_dir=repeat_dir,
        normalized_dir=normalized_dir,
        structures_dir=structures_dir,
        system_metrics_path=system_metrics_path,
        chain_metrics_path=chain_metrics_path,
        manifest_path=manifest_path,
        capabilities={"confidence_metrics"},
        runtime=RunnerRuntime(diffusion_samples=1),
        metric_outcomes={
            "confidence_metrics": RunnerMetricOutcome(state="computed"),
        },
        companion_artifacts=artifacts,
    )

    with pytest.raises(RunnerBundleValidationError, match="missing path"):
        validate_runner_bundle(
            result,
            requested_metric_groups={"confidence_metrics"},
        )
