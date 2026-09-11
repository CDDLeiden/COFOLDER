"""Internal helpers for loading optional analysis dependencies on demand."""

from __future__ import annotations

import importlib
from types import ModuleType


ANALYSIS_INSTALL_GUIDANCE = 'pip install "cofolder[analysis]"'


def require_analysis_dependency(
    module_name: str,
    *,
    feature: str,
    dependency_name: str | None = None,
) -> ModuleType:
    """Import one analysis dependency or raise feature-specific guidance."""
    display_name = dependency_name or module_name.split(".", maxsplit=1)[0]
    try:
        return importlib.import_module(module_name)
    except (ImportError, ModuleNotFoundError) as exc:
        raise ImportError(
            f"{feature} requires the optional analysis dependency "
            f"'{display_name}'. Install the COFOLDER analysis extra with "
            f"`{ANALYSIS_INSTALL_GUIDANCE}`."
        ) from exc
