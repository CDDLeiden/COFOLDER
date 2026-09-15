"""Characterization tests for the private workflow extraction seams."""

from __future__ import annotations

import json
import logging
import subprocess
import sys
from dataclasses import replace

import pandas as pd
import pytest

from cofolder.modules.contracts import (
    BackendVersionStatus,
    FailureStage,
    OutputIdentity,
    RepeatSeedProvenance,
    RunnerBackendIdentity,
    SeedAdjustment,
    SeedOrigin,
    SeedPlan,
    SeedResolutionError,
    WorkflowKind,
)
from cofolder.modules.utils.timing import DebugTimingCollector
from cofolder.recipes._diagnostics import FailureStageTracker
from cofolder.recipes._execution import adjust_seed_plan, build_execution_plan
from cofolder.recipes._results import (
    existing_directory_artifacts,
    standard_evidence,
    write_failure_bundle,
    write_frame_bundle,
)


def _identity() -> OutputIdentity:
    return OutputIdentity(
        workflow=WorkflowKind.VALIDATE,
        run_id="run-1",
        system_id="system-1",
        runner_id="runner-1",
    )


class _PlanningRunner:
    def resolve_effective_seed(self, seed):
        return replace(seed, effective_seed=seed.effective_seed + 1)

    def detect_backend_identity(self):
        return RunnerBackendIdentity(
            "runner-1",
            "backend-1",
            "1.2.3",
            BackendVersionStatus.DETECTED,
        )

    def execution_models(self, _options):
        from cofolder.modules.runners import RunnerModelSlot

        return (RunnerModelSlot("model-a", 0), RunnerModelSlot("model-b", 1))


def test_screen_execution_plan_preserves_historical_unvalidated_adjustment():
    plan = build_execution_plan(
        _PlanningRunner(),
        "runner-1",
        object(),
        repeats=1,
        seed=7,
        logger=logging.getLogger(__name__),
    )

    assert plan.backend.version == "1.2.3"
    assert plan.seed_plan.repeats[0].effective_seed == 8
    assert [
        (item.repeat_id, item.model_id, item.sample_id, item.effective_seed)
        for item in plan.executions
    ] == [(1, "model-a", 0, 8), (1, "model-b", 1, 8)]


def test_validate_seed_adjustment_remains_strict():
    original = RepeatSeedProvenance(
        1,
        7,
        7,
        7,
        7,
        SeedOrigin.USER_SPECIFIED,
    )
    plan = SeedPlan(7, 7, SeedOrigin.USER_SPECIFIED, (original,))
    with pytest.raises(SeedResolutionError, match="requires an explicit"):
        adjust_seed_plan(plan, _PlanningRunner())

    adjusted = replace(
        original,
        effective_seed=8,
        adjustment=SeedAdjustment.BACKEND_ADJUSTED,
        adjustment_reason="Backend reserves seed 7.",
    )
    runner = _PlanningRunner()
    runner.resolve_effective_seed = lambda _seed: adjusted
    assert adjust_seed_plan(plan, runner).repeats == (adjusted,)


def test_failure_stage_tracker_retains_the_failing_stage():
    tracker = FailureStageTracker(
        DebugTimingCollector(),
        logging.getLogger(__name__),
    )
    with pytest.raises(RuntimeError, match="broken"):
        with tracker.measure("runner.prepare_system"):
            raise RuntimeError("broken")
    assert tracker.current is FailureStage.PREPARATION


def test_shared_result_helpers_preserve_terminal_bundle_semantics(temp_dir):
    reference = temp_dir / "reference.cif"
    reference.write_text("data_reference\n", encoding="utf-8")
    (temp_dir / "structures").mkdir()
    evidence = standard_evidence(
        reference_path=reference,
        pocket_coverage_reference="A:1-4",
    )
    artifacts = existing_directory_artifacts(
        temp_dir,
        (("predicted_structures", "structures"), ("missing", "missing")),
    )
    assert [item.kind for item in evidence] == [
        "reference_structure",
        "custom_pocket",
    ]
    assert [item.label for item in artifacts] == ["predicted_structures"]

    output_dir = temp_dir / "success"
    write_frame_bundle(
        pd.DataFrame(
            [{"model_name": "model-a", "repeat": 1, "diffusion_sample": 0, "ptm": 0.8}]
        ),
        pd.DataFrame(),
        output_dir=output_dir,
        identity=_identity(),
        evidence=evidence,
        requested_metrics={"confidence_metrics"},
        artifacts=artifacts,
    )
    manifest = json.loads((output_dir / "manifest.json").read_text())
    assert manifest["status"] == "success"
    assert manifest["requested_metrics"] == ["confidence_metrics"]
    assert [item["label"] for item in manifest["artifacts"]] == [
        "predicted_structures"
    ]

    failed_dir = temp_dir / "failed"
    failure = write_failure_bundle(
        output_dir=failed_dir,
        identity=_identity(),
        exc=ValueError("invalid request"),
        stage=FailureStage.INPUT_VALIDATION,
        error_code="validate_input_validation_failed",
    )
    failed_manifest = json.loads((failed_dir / "manifest.json").read_text())
    records = [json.loads(line) for line in (failed_dir / "records.jsonl").read_text().splitlines()]
    assert failed_manifest["status"] == "failed"
    assert failure.error_code == "validate_input_validation_failed"
    assert records[0]["stage"] == "input_validation"
    assert records[0]["message"] == "invalid request"


def test_importing_private_workflow_helpers_has_no_execution_side_effects(temp_dir):
    code = """
import sys
def deny(event, args):
    if event in {'subprocess.Popen', 'socket.connect'}:
        raise RuntimeError(event)
sys.addaudithook(deny)
import cofolder.recipes._execution
import cofolder.recipes._diagnostics
import cofolder.recipes._results
import cofolder.recipes._screen_postprocess
"""
    before = tuple(temp_dir.iterdir())
    completed = subprocess.run(
        [sys.executable, "-c", code],
        cwd=temp_dir,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert tuple(temp_dir.iterdir()) == before
