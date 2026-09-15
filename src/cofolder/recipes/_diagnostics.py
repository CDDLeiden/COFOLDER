"""Workflow-local diagnostic stage tracking and typed failure construction."""

from __future__ import annotations

import logging
from contextlib import contextmanager
from typing import Iterator, Mapping

from cofolder.modules.contracts import (
    FailureStage,
    OutputIdentity,
    WorkflowFailureRecord,
    failure_from_exception,
)
from cofolder.modules.utils.timing import DebugTimingCollector


def failure_stage_for_label(label: str) -> FailureStage:
    """Map Validate timing labels to their stable public failure stages."""
    if label.startswith("runner.prepare"):
        return FailureStage.PREPARATION
    if label.startswith("runner.bundle"):
        return FailureStage.OUTPUT_VALIDATION
    if label.startswith(("structures.gather", "results.merge", "results.add_chain")):
        return FailureStage.GATHER
    if label.startswith(("scores.", "structure.")):
        return FailureStage.ANALYTICS
    if label.startswith("results.write"):
        return FailureStage.SERIALIZATION
    return FailureStage.INPUT_VALIDATION


class FailureStageTracker:
    """Track the active public failure stage alongside debug-only timings."""

    def __init__(
        self,
        timings: DebugTimingCollector,
        logger: logging.Logger,
    ) -> None:
        self.timings = timings
        self.logger = logger
        self.current = FailureStage.INPUT_VALIDATION

    @contextmanager
    def measure(self, label: str) -> Iterator[None]:
        previous = self.current
        self.current = failure_stage_for_label(label)
        try:
            with self.timings.measure(label, logger=self.logger):
                yield
        except Exception:
            raise
        else:
            self.current = previous


def workflow_failure(
    exc: BaseException,
    *,
    identity: OutputIdentity,
    stage: FailureStage,
    error_code: str | None = None,
    retryable: bool = False,
    details: Mapping[str, object] | None = None,
) -> WorkflowFailureRecord:
    """Construct a typed workflow failure at the recipe boundary."""
    return failure_from_exception(
        exc,
        identity=identity,
        stage=stage,
        error_code=error_code,
        retryable=retryable,
        details=details,
    )

