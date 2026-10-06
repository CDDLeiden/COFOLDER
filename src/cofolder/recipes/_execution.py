"""Shared, side-effect-free construction of workflow execution provenance."""

from __future__ import annotations

import logging
from dataclasses import replace
from typing import Any, Callable

from cofolder.modules.contracts import (
    BackendVersionStatus,
    RepeatSeedProvenance,
    RunnerBackendIdentity,
    SeedAdjustment,
    SeedPlan,
    SeedResolutionError,
)
from cofolder.modules.runners import (
    PlannedExecution,
    RunnerExecutionPlan,
    RunnerModelSlot,
)
from cofolder.modules.utils.helpers import resolve_seed_plan


def validate_effective_seed(
    original: RepeatSeedProvenance,
    adjusted: RepeatSeedProvenance,
) -> RepeatSeedProvenance:
    """Validate a runner's adjustment without losing seed provenance."""
    if not isinstance(adjusted, RepeatSeedProvenance):
        raise SeedResolutionError(
            "Runner seed adjustment must return RepeatSeedProvenance."
        )
    immutable_fields = (
        "repeat_id",
        "requested_base_seed",
        "resolved_base_seed",
        "derived_seed",
        "origin",
    )
    if any(
        getattr(original, field) != getattr(adjusted, field)
        for field in immutable_fields
    ):
        raise SeedResolutionError(
            "Runner seed adjustment changed immutable seed provenance."
        )
    effective = adjusted.effective_seed
    if (
        isinstance(effective, bool)
        or not isinstance(effective, int)
        or not (0 <= effective <= 2**32 - 1)
    ):
        raise SeedResolutionError(
            "Runner effective seed must be an integer between 0 and 4294967295."
        )
    changed = effective != original.derived_seed
    if changed and (
        adjusted.adjustment != SeedAdjustment.BACKEND_ADJUSTED
        or not adjusted.adjustment_reason
    ):
        raise SeedResolutionError(
            "A backend-adjusted seed requires an explicit adjustment status and reason."
        )
    if not changed and (
        adjusted.adjustment != SeedAdjustment.UNCHANGED
        or adjusted.adjustment_reason is not None
    ):
        raise SeedResolutionError(
            "An unchanged seed must use adjustment='unchanged' without a reason."
        )
    return adjusted


def adjust_seed_plan(
    seed_plan: SeedPlan,
    runner: Any,
    *,
    validate_adjustments: bool = True,
) -> SeedPlan:
    """Apply a runner's optional seed adjustment with caller-selected validation."""
    resolver = getattr(runner, "resolve_effective_seed", None)
    adjusted: list[RepeatSeedProvenance] = []
    for original in seed_plan.repeats:
        candidate = resolver(original) if callable(resolver) else original
        adjusted.append(
            validate_effective_seed(original, candidate)
            if validate_adjustments
            else candidate
        )
    return replace(seed_plan, repeats=tuple(adjusted))


def detect_backend_identity(
    runner: Any,
    runner_name: str,
    *,
    unavailable_detail: str,
) -> RunnerBackendIdentity:
    """Return typed backend identity, including an explicit unavailable result."""
    detect = getattr(runner, "detect_backend_identity", None)
    identity = detect() if callable(detect) else None
    return identity or RunnerBackendIdentity(
        runner_name=runner_name,
        backend_name=runner_name,
        version=None,
        version_status=BackendVersionStatus.UNAVAILABLE,
        detail=unavailable_detail,
    )


def build_execution_plan(
    runner: Any,
    runner_name: str,
    options_obj: Any,
    *,
    repeats: int,
    seed: int | None,
    logger: logging.Logger,
) -> RunnerExecutionPlan:
    """Resolve immutable backend, seed, repeat, model, and sample axes once."""
    seed_plan = adjust_seed_plan(
        resolve_seed_plan(repeats, seed, logger),
        runner,
        validate_adjustments=False,
    )
    backend = detect_backend_identity(
        runner,
        runner_name,
        unavailable_detail="Runner does not expose backend version detection.",
    )
    describe: Callable[[Any], Any] | None = getattr(runner, "execution_models", None)
    models = (
        tuple(describe(options_obj))
        if callable(describe)
        else (RunnerModelSlot(model_id=runner_name, sample_id=0),)
    )
    if not models:
        raise ValueError("Runner execution model plan must not be empty.")
    executions = tuple(
        PlannedExecution(
            repeat_id=repeat.repeat_id,
            model_id=model.model_id,
            sample_id=model.sample_id,
            effective_seed=repeat.effective_seed,
        )
        for repeat in seed_plan.repeats
        for model in models
    )
    return RunnerExecutionPlan(backend, seed_plan, executions)
