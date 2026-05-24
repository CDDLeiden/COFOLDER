"""Tests for normalized runner bundle validation."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from cofolder.modules.runners.contracts import (
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
    include_confidence: bool = True,
    include_affinity: bool = True,
    capabilities: set[str] | None = None,
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

    system_metrics_path = normalized_dir / "system_metrics.csv"
    chain_metrics_path = normalized_dir / "chain_metrics.csv"
    manifest_path = normalized_dir / "manifest.json"

    pd.DataFrame([system_row]).to_csv(system_metrics_path, index=False)
    pd.DataFrame([chain_row]).to_csv(chain_metrics_path, index=False)
    manifest_path.write_text(
        json.dumps(
            {
                "runner": "simple",
                "capabilities": sorted(capabilities or []),
            }
        ),
        encoding="utf-8",
    )

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
    assert bundle.metric_outcomes["confidence_metrics"].required_columns == (
        "system_metrics.ptm",
        "system_metrics.iptm",
        "system_metrics.confidence_score",
        "chain_metrics.chains_ptm",
    )


def test_validate_runner_bundle_rejects_computed_without_required_payload(temp_dir):
    repeat_dir = temp_dir / "repeat_1"
    normalized_dir, structures_dir, system_metrics_path, chain_metrics_path = _write_bundle_files(
        repeat_dir,
        include_confidence=False,
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
        },
    )

    with pytest.raises(RunnerBundleValidationError, match="marked 'computed'"):
        validate_runner_bundle(
            result,
            requested_metric_groups={"confidence_metrics"},
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

    with pytest.raises(RunnerBundleValidationError, match="while normalized payload columns are present"):
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
