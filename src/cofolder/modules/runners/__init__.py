from __future__ import annotations

import importlib
import pkgutil
from typing import Any

from cofolder.modules.runners.base import BaseRunner, Runner, RunnerPreparation, RunnerRequest, RunnerResult


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
    "RunnerPreparation",
    "RunnerRequest",
    "RunnerResult",
    "discover_runners",
    "get_runner",
    "list_runner_names",
]
