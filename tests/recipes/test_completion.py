"""Tests for top-level workflow completion reporting."""

import logging
from dataclasses import replace
from pathlib import Path

import pytest

from cofolder import cli
from cofolder.modules.contracts import (
    ExecutionRecord,
    ExecutionStatus,
    FailureStage,
    OutputIdentity,
    PublicManifest,
    PublicOutputBundle,
    RecordKind,
    StructuredExecutionError,
    SuccessRecord,
    WorkflowFailureRecord,
    WorkflowKind,
    make_envelope,
    write_public_bundle,
)
from cofolder.recipes._completion import report_completion


def _identity(workflow: WorkflowKind, *, compound_id: str | None = None):
    values = {
        "workflow": workflow,
        "run_id": "completion-test",
        "system_id": "system",
        "compound_id": compound_id,
    }
    if workflow is not WorkflowKind.BIAS:
        values["runner_id"] = "test-runner"
    return OutputIdentity(
        **values,
    )


def _write_bundle(
    root: Path,
    workflow: WorkflowKind,
    status: str,
    records: tuple,
) -> None:
    write_public_bundle(
        PublicOutputBundle(
            manifest=PublicManifest(
                schema_version="1.0.0",
                identity=_identity(workflow),
                status=status,
            ),
            records=records,
        ),
        root / "results",
    )


class _Recipe:
    def __init__(self, wrk_dir: Path):
        self.wrk_dir = wrk_dir
        self.logger = logging.getLogger("cofolder.test.completion")


@pytest.mark.parametrize("workflow", [WorkflowKind.VALIDATE, WorkflowKind.BIAS])
def test_success_summary_uses_metric_table(workflow, temp_dir, caplog):
    identity = _identity(workflow)

    class Recipe(_Recipe):
        @report_completion(workflow)
        def run(self):
            _write_bundle(
                self.wrk_dir,
                workflow,
                "success",
                (SuccessRecord(make_envelope(RecordKind.SUCCESS, identity, "success")),),
            )
            return "result"

    with caplog.at_level(logging.INFO, logger="cofolder.test.completion"):
        assert Recipe(temp_dir).run() == "result"

    assert f"workflow={workflow.value} status=success total=1 success=1 failed=0" in caplog.text
    assert f"primary_table={(temp_dir / 'results/metrics.csv').resolve()}" in caplog.text
    assert f"manifest={(temp_dir / 'results/manifest.json').resolve()}" in caplog.text


def test_partial_screen_summary_counts_execution_outcomes(temp_dir, caplog):
    success_identity = _identity(WorkflowKind.SCREEN, compound_id="success")
    success_identity = replace(
        success_identity,
        execution_directory="success run",
        repeat_id=1,
        model_id="model_0",
    )
    failed_identity = _identity(WorkflowKind.SCREEN, compound_id="failed")
    failed_identity = replace(
        failed_identity,
        execution_directory="failed run",
        repeat_id=1,
        model_id="model_0",
    )
    records = (
        ExecutionRecord(
            envelope=make_envelope(RecordKind.EXECUTION, success_identity, "execution"),
            status=ExecutionStatus.SUCCESS,
            execution_directory="success run",
        ),
        ExecutionRecord(
            envelope=make_envelope(RecordKind.EXECUTION, failed_identity, "execution"),
            status=ExecutionStatus.FAILED,
            execution_directory="failed run",
            error=StructuredExecutionError(
                stage=FailureStage.BACKEND_EXECUTION,
                exception_type="RuntimeError",
                error_code="backend_failed",
                message="backend failed",
            ),
        ),
    )

    class Recipe(_Recipe):
        @report_completion(WorkflowKind.SCREEN)
        def run(self):
            _write_bundle(self.wrk_dir, WorkflowKind.SCREEN, "partial", records)

    with caplog.at_level(logging.INFO, logger="cofolder.test.completion"):
        Recipe(temp_dir).run()

    assert "workflow=screen status=partial total=2 success=1 failed=1" in caplog.text
    assert f"primary_table={(temp_dir / 'results/executions.csv').resolve()}" in caplog.text


def test_failed_oracle_summary_reports_unavailable_value(temp_dir, caplog):
    identity = _identity(WorkflowKind.ORACLE, compound_id="sha256:candidate")
    failure = WorkflowFailureRecord(
        envelope=make_envelope(
            RecordKind.FAILURE, identity, FailureStage.ANALYTICS.value, "missing_metric"
        ),
        stage=FailureStage.ANALYTICS,
        exception_type="ValueError",
        error_code="missing_metric",
        message="metric unavailable",
    )

    class Recipe(_Recipe):
        output_metric = "affinity_pred_value"
        aggregate = "mean"
        score_components = {}

        @report_completion(WorkflowKind.ORACLE)
        def run(self):
            _write_bundle(
                self.wrk_dir, WorkflowKind.ORACLE, "failed", (failure,)
            )
            raise RuntimeError("metric unavailable")

    with caplog.at_level(logging.INFO, logger="cofolder.test.completion"):
        with pytest.raises(RuntimeError, match="metric unavailable"):
            Recipe(temp_dir).run()

    assert "workflow=oracle status=failed total=1 success=0 failed=1" in caplog.text
    assert f"primary_table={(temp_dir / 'results/failures.csv').resolve()}" in caplog.text
    assert f"manifest={(temp_dir / 'results/manifest.json').resolve()}" in caplog.text
    assert "Oracle result | metric=affinity_pred_value aggregate=mean value=unavailable" in caplog.text


def test_successful_oracle_summary_reports_scalar(temp_dir, caplog):
    identity = _identity(WorkflowKind.ORACLE, compound_id="sha256:candidate")

    class Recipe(_Recipe):
        output_metric = "affinity_pred_value"
        aggregate = "median"
        score_components = {}

        @report_completion(WorkflowKind.ORACLE)
        def run(self):
            _write_bundle(
                self.wrk_dir,
                WorkflowKind.ORACLE,
                "success",
                (SuccessRecord(make_envelope(RecordKind.SUCCESS, identity, "success")),),
            )
            return 0.625

    with caplog.at_level(logging.INFO, logger="cofolder.test.completion"):
        assert Recipe(temp_dir).run() == 0.625

    assert "workflow=oracle status=success total=1 success=1 failed=0" in caplog.text
    assert f"primary_table={(temp_dir / 'results/metrics.csv').resolve()}" in caplog.text
    assert f"manifest={(temp_dir / 'results/manifest.json').resolve()}" in caplog.text
    assert "Oracle result | metric=affinity_pred_value aggregate=median value=0.625" in caplog.text


def test_nested_recipe_reports_only_outer_completion(temp_dir, caplog):
    inner_root = temp_dir / "inner"
    outer_identity = _identity(WorkflowKind.SCREEN, compound_id="compound")
    outer_identity = replace(
        outer_identity,
        execution_directory="compound",
        repeat_id=1,
        model_id="model_0",
    )

    class Inner(_Recipe):
        @report_completion(WorkflowKind.VALIDATE)
        def run(self):
            identity = _identity(WorkflowKind.VALIDATE)
            _write_bundle(
                self.wrk_dir,
                WorkflowKind.VALIDATE,
                "success",
                (SuccessRecord(make_envelope(RecordKind.SUCCESS, identity, "success")),),
            )

    class Outer(_Recipe):
        @report_completion(WorkflowKind.SCREEN)
        def run(self):
            Inner(inner_root).run()
            record = ExecutionRecord(
                envelope=make_envelope(
                    RecordKind.EXECUTION, outer_identity, "execution"
                ),
                status=ExecutionStatus.SUCCESS,
                execution_directory="compound",
            )
            _write_bundle(
                self.wrk_dir, WorkflowKind.SCREEN, "success", (record,)
            )

    with caplog.at_level(logging.INFO, logger="cofolder.test.completion"):
        Outer(temp_dir).run()

    assert caplog.text.count("Workflow summary |") == 1
    assert "workflow=screen" in caplog.text
    assert "workflow=validate" not in caplog.text


def test_cli_emits_shared_summary_on_stdout(
    sample_system_yaml, sample_options_yaml, temp_dir, capsys, monkeypatch
):
    class Recipe(_Recipe):
        def __init__(self, **kwargs):
            super().__init__(Path(kwargs["wrk_dir"]))

        @report_completion(WorkflowKind.VALIDATE)
        def run(self):
            identity = _identity(WorkflowKind.VALIDATE)
            _write_bundle(
                self.wrk_dir,
                WorkflowKind.VALIDATE,
                "success",
                (SuccessRecord(make_envelope(RecordKind.SUCCESS, identity, "success")),),
            )

    root = logging.getLogger()
    handlers = list(root.handlers)
    level = root.level
    try:
        root.handlers.clear()
        monkeypatch.setattr(cli, "_load_recipe_class", lambda *_args: Recipe)
        result = cli.main(
            [
                "validate",
                "-s",
                str(sample_system_yaml),
                "-o",
                str(sample_options_yaml),
                "-w",
                str(temp_dir / "path with spaces"),
            ]
        )
        captured = capsys.readouterr()
    finally:
        root.handlers[:] = handlers
        root.setLevel(level)

    assert result == 0
    assert captured.err == ""
    assert "INFO | Workflow summary | workflow=validate status=success" in captured.out
    assert "primary_table=" in captured.out
    assert "path with spaces/results/metrics.csv" in captured.out
    assert "Validation pipeline completed." not in captured.out
