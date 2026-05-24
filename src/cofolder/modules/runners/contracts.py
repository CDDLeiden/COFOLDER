from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal
from typing import Any
from typing import Mapping


@dataclass(slots=True)
class RunnerRuntime:
    """Shared runtime metadata that recipes may safely consume."""

    diffusion_samples: int = 1
    cache_path: str | None = None
    model_name: str | None = None


RunnerMetricState = Literal["computed", "unsupported", "missing", "failed"]


@dataclass(slots=True)
class RunnerMetricOutcome:
    """Explicit per-metric-group state produced at the runner boundary."""

    state: RunnerMetricState
    required_columns: tuple[str, ...] = ()
    message: str | None = None


@dataclass(slots=True)
class RunnerNormalizedBundle:
    """Authoritative normalized output bundle consumed by shared validation."""

    runner_name: str
    raw_output_dir: Path
    normalized_dir: Path
    structures_dir: Path
    system_metrics_path: Path
    chain_metrics_path: Path
    manifest_path: Path
    diffusion_samples: int = 1
    capabilities: set[str] = field(default_factory=set)
    runtime: RunnerRuntime = field(default_factory=RunnerRuntime)
    sample_records: list[dict[str, Any]] = field(default_factory=list)
    metric_outcomes: dict[str, RunnerMetricOutcome] = field(default_factory=dict)


def merge_runner_runtime(*values: RunnerRuntime | None) -> RunnerRuntime:
    """Merge typed runtime metadata without falling back to ad hoc dict lookups."""

    merged = RunnerRuntime()
    for value in values:
        if value is None:
            continue
        if value.cache_path is not None:
            merged.cache_path = value.cache_path
        if value.model_name is not None:
            merged.model_name = value.model_name
        # Authoritative typed runtimes are merged in precedence order, so an
        # explicit single-sample value from a later runtime must still win.
        merged.diffusion_samples = value.diffusion_samples
    return merged


def copy_metric_outcomes(
    metric_outcomes: Mapping[str, RunnerMetricOutcome] | None,
) -> dict[str, RunnerMetricOutcome]:
    """Copy explicit metric outcomes into a mutable normalized-bundle mapping."""

    copied: dict[str, RunnerMetricOutcome] = {}
    if metric_outcomes is None:
        return copied

    for group_name, outcome in metric_outcomes.items():
        if not isinstance(outcome, RunnerMetricOutcome):
            raise TypeError(
                "metric_outcomes values must be RunnerMetricOutcome instances, "
                f"got {type(outcome)!r} for group {group_name!r}."
            )
        copied[str(group_name)] = RunnerMetricOutcome(
            state=outcome.state,
            required_columns=tuple(outcome.required_columns),
            message=outcome.message,
        )
    return copied


@dataclass(slots=True, init=False)
class RunnerPreparationResult:
    system_obj: Any
    options_obj: Any
    warnings: list[str] = field(default_factory=list)
    runtime: RunnerRuntime = field(default_factory=RunnerRuntime)

    def __init__(
        self,
        system_obj: Any,
        options_obj: Any,
        warnings: list[str] | None = None,
        runtime: RunnerRuntime | None = None,
    ) -> None:
        self.system_obj = system_obj
        self.options_obj = options_obj
        self.warnings = list(warnings or [])
        self.runtime = runtime or RunnerRuntime()


@dataclass(slots=True)
class RunnerExecutionRequest:
    runner_name: str
    system_name: str
    system_path: Path
    system_obj: Any
    options_path: Path
    options_obj: Any
    repeat: int
    seed: int
    repeat_dir: Path
    raw_dir: Path
    logger: logging.Logger | None
    timings: Any | None = None
    label_prefix: str | None = None
    runtime: RunnerRuntime = field(default_factory=RunnerRuntime)


@dataclass(slots=True, init=False)
class RunnerExecutionResult:
    runner_name: str
    raw_output_dir: Path
    normalized_dir: Path
    structures_dir: Path
    system_metrics_path: Path
    chain_metrics_path: Path
    manifest_path: Path
    diffusion_samples: int = 1
    capabilities: set[str] = field(default_factory=set)
    warnings: list[str] = field(default_factory=list)
    runtime: RunnerRuntime = field(default_factory=RunnerRuntime)
    sample_records: list[dict[str, Any]] = field(default_factory=list)
    normalized_bundle: RunnerNormalizedBundle = field(init=False)

    def __init__(
        self,
        runner_name: str,
        raw_output_dir: Path,
        normalized_dir: Path,
        structures_dir: Path,
        system_metrics_path: Path,
        chain_metrics_path: Path,
        manifest_path: Path,
        diffusion_samples: int = 1,
        capabilities: set[str] | None = None,
        warnings: list[str] | None = None,
        runtime: RunnerRuntime | None = None,
        sample_records: list[dict[str, Any]] | None = None,
        metric_outcomes: Mapping[str, RunnerMetricOutcome] | None = None,
    ) -> None:
        self.runner_name = runner_name
        self.raw_output_dir = raw_output_dir
        self.normalized_dir = normalized_dir
        self.structures_dir = structures_dir
        self.system_metrics_path = system_metrics_path
        self.chain_metrics_path = chain_metrics_path
        self.manifest_path = manifest_path
        self.diffusion_samples = int(diffusion_samples)
        self.capabilities = set(capabilities or set())
        self.warnings = list(warnings or [])
        self.runtime = runtime or RunnerRuntime(diffusion_samples=self.diffusion_samples)
        self.sample_records = list(sample_records or [])
        self.normalized_bundle = RunnerNormalizedBundle(
            runner_name=self.runner_name,
            raw_output_dir=self.raw_output_dir,
            normalized_dir=self.normalized_dir,
            structures_dir=self.structures_dir,
            system_metrics_path=self.system_metrics_path,
            chain_metrics_path=self.chain_metrics_path,
            manifest_path=self.manifest_path,
            diffusion_samples=self.diffusion_samples,
            capabilities=set(self.capabilities),
            runtime=self.runtime,
            sample_records=list(self.sample_records),
            metric_outcomes=copy_metric_outcomes(metric_outcomes),
        )

    @property
    def metric_outcomes(self) -> dict[str, RunnerMetricOutcome]:
        return self.normalized_bundle.metric_outcomes
