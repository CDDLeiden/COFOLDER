"""Runner discovery and public import surface."""

from __future__ import annotations

import importlib
import pkgutil

from cofolder.modules.runners.base import (
    BaseRunner,
    Runner,
)
from cofolder.modules.runners.contracts import (
    RUNNER_PROVENANCE_COLUMNS,
    PlannedExecution,
    RepeatSeedProvenance,
    RunnerBackendIdentity,
    RunnerChainIdentity,
    RunnerCompanionArtifact,
    RunnerExecutionPlan,
    RunnerExecutionRequest,
    RunnerExecutionResult,
    RunnerInputCapabilities,
    RunnerMetricOutcome,
    RunnerModelSlot,
    RunnerNormalizedBundle,
    RunnerPreparationResult,
    RunnerRuntime,
    attach_runner_provenance,
    attach_sample_provenance,
    backend_manifest_value,
    build_runner_chain_identities,
    build_runner_public_records,
    merge_runner_runtime,
    runner_provenance_values,
    seed_manifest_value,
)
from cofolder.modules.runners.validators import (
    RunnerBundleValidationError,
    validate_runner_bundle,
)


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
        raise ValueError(
            f"Unknown runner '{name}'. Available runners: {available}"
        ) from exc


__all__ = [
    "RUNNER_PROVENANCE_COLUMNS",
    "BaseRunner",
    "PlannedExecution",
    "RepeatSeedProvenance",
    "Runner",
    "RunnerBackendIdentity",
    "RunnerBundleValidationError",
    "RunnerChainIdentity",
    "RunnerCompanionArtifact",
    "RunnerExecutionPlan",
    "RunnerExecutionRequest",
    "RunnerExecutionResult",
    "RunnerInputCapabilities",
    "RunnerMetricOutcome",
    "RunnerModelSlot",
    "RunnerNormalizedBundle",
    "RunnerPreparationResult",
    "RunnerRuntime",
    "attach_runner_provenance",
    "attach_sample_provenance",
    "backend_manifest_value",
    "build_runner_chain_identities",
    "build_runner_public_records",
    "discover_runners",
    "get_runner",
    "list_runner_names",
    "merge_runner_runtime",
    "runner_provenance_values",
    "seed_manifest_value",
    "validate_runner_bundle",
]
