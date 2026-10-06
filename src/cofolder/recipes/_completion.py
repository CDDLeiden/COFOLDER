"""Internal terminal-summary reporting for top-level recipe executions."""

from __future__ import annotations

import json
import logging
from contextvars import ContextVar
from functools import wraps
from pathlib import Path
from typing import Any, Callable, TypeVar

from cofolder.modules.contracts import (
    ExecutionRecord,
    ExecutionStatus,
    SuccessRecord,
    WorkflowFailureRecord,
    WorkflowKind,
    read_public_records,
)

_T = TypeVar("_T")
_completion_depth: ContextVar[int] = ContextVar("cofolder_completion_depth", default=0)


def report_completion(workflow: WorkflowKind) -> Callable[[Callable[..., _T]], Callable[..., _T]]:
    """Report one terminal summary for the outermost public recipe call."""

    def decorate(method: Callable[..., _T]) -> Callable[..., _T]:
        @wraps(method)
        def wrapped(recipe: Any, *args: Any, **kwargs: Any) -> _T:
            depth = _completion_depth.get()
            token = _completion_depth.set(depth + 1)
            result: Any = None
            try:
                result = method(recipe, *args, **kwargs)
            except Exception:
                if depth == 0:
                    _log_completion(recipe, workflow, result=None)
                raise
            else:
                if depth == 0:
                    _log_completion(recipe, workflow, result=result)
                return result
            finally:
                _completion_depth.reset(token)

        return wrapped

    return decorate


def _log_completion(recipe: Any, workflow: WorkflowKind, *, result: Any) -> None:
    """Read canonical outputs and log a summary without masking workflow results."""

    logger = getattr(recipe, "logger", logging.getLogger(f"cofolder.{workflow.value}"))
    output_dir = (Path(recipe.wrk_dir) / "results").resolve()
    manifest_path = output_dir / "manifest.json"
    records_path = output_dir / "records.jsonl"
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        records = read_public_records(records_path)
        status = str(manifest["status"])

        executions = [record for record in records if isinstance(record, ExecutionRecord)]
        if executions:
            successes = sum(
                record.status is ExecutionStatus.SUCCESS for record in executions
            )
            failures = len(executions) - successes
            total = len(executions)
        else:
            successes = sum(isinstance(record, SuccessRecord) for record in records)
            failures = sum(
                isinstance(record, WorkflowFailureRecord) for record in records
            )
            total = successes + failures

        primary_name = (
            "failures.csv"
            if status == "failed"
            else "executions.csv"
            if workflow is WorkflowKind.SCREEN
            else "metrics.csv"
        )
        level = logging.ERROR if status == "failed" else logging.INFO
        logger.log(
            level,
            "Workflow summary | workflow=%s status=%s total=%d success=%d failed=%d",
            workflow.value,
            status,
            total,
            successes,
            failures,
        )
        logger.log(
            level,
            "Results | output_dir=%s primary_table=%s manifest=%s",
            output_dir,
            output_dir / primary_name,
            manifest_path,
        )
        if workflow is WorkflowKind.ORACLE:
            metric = getattr(recipe, "output_metric", None)
            if not metric:
                metric = "composite" if getattr(recipe, "score_components", None) else "custom"
            value = result if status != "failed" and result is not None else "unavailable"
            logger.log(
                level,
                "Oracle result | metric=%s aggregate=%s value=%s",
                metric,
                getattr(recipe, "aggregate", "first"),
                value,
            )
    except Exception as exc:
        logger.warning(
            "Completion summary unavailable | workflow=%s output_dir=%s reason=%s",
            workflow.value,
            output_dir,
            exc,
        )
