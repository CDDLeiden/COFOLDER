"""Runner discovery and public import surface.

New code should author against the explicit contract types imported from
`cofolder.modules.runners.contracts`. Legacy alias names remain exportable for
import compatibility, but they are not the recommended surface for new runner
implementations or contributor docs.
"""

from __future__ import annotations

import importlib
import pkgutil

from cofolder.modules.runners.base import BaseRunner, Runner, RunnerPreparation, RunnerRequest, RunnerResult
from cofolder.modules.runners.contracts import (
    RunnerCompanionArtifact,
    RunnerExecutionRequest,
    RunnerExecutionResult,
    RunnerMetricOutcome,
    RunnerInputCapabilities,
    RunnerNormalizedBundle,
    RunnerPreparationResult,
    RunnerRuntime,
    merge_runner_runtime,
)
from cofolder.modules.runners.validators import RunnerBundleValidationError, validate_runner_bundle


def discover_runners() -> dict[str, Runner]:
    runners: dict[str, Runner] = {}
    for module_info in pkgutil.iter_modules(__path__):
        if module_info.name.startswith("_") or module_info.name == "base":
            continue

        module = importlib.import_module(f"{__name__}.{module_info.name}")
        runner = getattr(module, "RUNNER", None)
        if runner is None:
            continue
        runners[runner.name] = runner
    return dict(sorted(runners.items()))


def list_runner_names() -> list[str]:
    return sorted(discover_runners())


def get_runner(name: str) -> Runner:
    runners = discover_runners()
    try:
        return runners[name]
    except KeyError as exc:
        available = ", ".join(sorted(runners)) or "<none>"
        raise ValueError(f"Unknown runner '{name}'. Available runners: {available}") from exc


__all__ = [
    "Runner",
    "BaseRunner",
    "RunnerCompanionArtifact",
    "RunnerExecutionRequest",
    "RunnerExecutionResult",
    "RunnerMetricOutcome",
    "RunnerInputCapabilities",
    "RunnerNormalizedBundle",
    "RunnerPreparation",
    "RunnerPreparationResult",
    "RunnerRequest",
    "RunnerResult",
    "RunnerRuntime",
    "RunnerBundleValidationError",
    "discover_runners",
    "get_runner",
    "list_runner_names",
    "merge_runner_runtime",
    "validate_runner_bundle",
]
