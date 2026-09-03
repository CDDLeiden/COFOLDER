from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from importlib import metadata
from pathlib import Path
from typing import Any, Protocol

from cofolder.modules.runners.contracts import (
    RunnerExecutionRequest,
    RunnerExecutionResult,
    RunnerPreparationResult,
    RunnerRuntime,
    RunnerInputCapabilities,
)

RunnerPreparation = RunnerPreparationResult
RunnerRequest = RunnerExecutionRequest
RunnerResult = RunnerExecutionResult


class BaseRunner(ABC):
    name: str
    capabilities: set[str]
    input_capabilities = RunnerInputCapabilities(
        entity_types=frozenset(),
        constraint_types=frozenset(),
    )

    @staticmethod
    def _distribution_installed(distribution_name: str) -> bool:
        try:
            metadata.version(distribution_name)
        except metadata.PackageNotFoundError:
            return False
        return True

    @staticmethod
    def get_distribution_version(distribution_name: str) -> str | None:
        try:
            return metadata.version(distribution_name)
        except metadata.PackageNotFoundError:
            return None

    @staticmethod
    def check_distribution_available(
        distribution_name: str,
        missing_message: str,
        conflicting_distributions: tuple[str, ...] = (),
        conflict_message: str | None = None,
    ) -> tuple[bool, str | None]:
        if not BaseRunner._distribution_installed(distribution_name):
            return False, missing_message

        installed_conflicts = [
            name for name in conflicting_distributions if BaseRunner._distribution_installed(name)
        ]
        if installed_conflicts:
            if conflict_message is not None:
                return False, conflict_message
            conflicts = ", ".join(sorted(installed_conflicts))
            return False, (
                f"Runner '{distribution_name}' cannot be used because conflicting runner "
                f"distributions are installed in the same environment: {conflicts}."
            )
        return True, None

    @abstractmethod
    def check_availability(self) -> tuple[bool, str | None]:
        """Return whether the runner dependencies are available."""

    def ensure_available(self) -> None:
        available, message = self.check_availability()
        if available:
            return
        raise RuntimeError(message or f"Runner '{self.name}' is not available.")

    def validate_system(
        self,
        system_obj: Any,
        options_obj: Any,
        *,
        check_atom_names: bool = True,
    ) -> None:
        """Validate the shared YAML contract against this runner's capabilities."""

        from cofolder.modules.input.validation import validate_system_input

        cache_path = None
        find_value = getattr(options_obj, "find_value", None)
        if callable(find_value):
            cache_path = find_value(key="cache") or find_value(key="cache_path")
        validate_system_input(
            system_obj,
            runner_name=self.name,
            capabilities=self.input_capabilities,
            cache_path=cache_path,
            check_atom_names=check_atom_names,
        )

    def prepare_system(
        self,
        system_obj: Any,
        options_obj: Any,
        wrk_dir: Path,
        conformers: str | None,
        sdf_file: Path | None,
        logger: logging.Logger | None,
    ) -> RunnerPreparationResult:
        """Default no-op preparation path for simple runners."""

        return RunnerPreparationResult(
            system_obj=system_obj,
            options_obj=options_obj,
        )

    @abstractmethod
    def load_options(self, options_path: Path):
        """Load any runner-specific options from disk."""

    @abstractmethod
    def run(self, request: RunnerExecutionRequest) -> RunnerExecutionResult:
        """Execute the backend and return canonical normalized outputs."""


class Runner(Protocol):
    name: str
    capabilities: set[str]
    input_capabilities: RunnerInputCapabilities

    def check_availability(self) -> tuple[bool, str | None]:
        ...

    def ensure_available(self) -> None:
        ...

    def validate_system(
        self,
        system_obj: Any,
        options_obj: Any,
        *,
        check_atom_names: bool = True,
    ) -> None:
        ...

    def load_options(self, options_path: Path):
        ...

    def prepare_system(
        self,
        system_obj: Any,
        options_obj: Any,
        wrk_dir: Path,
        conformers: str | None,
        sdf_file: Path | None,
        logger: logging.Logger | None,
    ) -> RunnerPreparationResult:
        ...

    def run(self, request: RunnerExecutionRequest) -> RunnerExecutionResult:
        ...
