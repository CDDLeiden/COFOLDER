from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from importlib import metadata
from pathlib import Path
from typing import Any, Protocol


@dataclass
class RunnerPreparation:
    system_obj: Any
    options_obj: Any
    warnings: list[str] = field(default_factory=list)
    runtime_context: dict[str, Any] = field(default_factory=dict)


@dataclass
class RunnerRequest:
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
    logger: logging.Logger
    timings: Any | None = None
    label_prefix: str | None = None


@dataclass
class RunnerResult:
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
    runtime_context: dict[str, Any] = field(default_factory=dict)
    sample_records: list[dict[str, Any]] = field(default_factory=list)


class BaseRunner(ABC):
    name: str
    capabilities: set[str]

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


class Runner(Protocol):
    name: str
    capabilities: set[str]

    def check_availability(self) -> tuple[bool, str | None]:
        ...

    def ensure_available(self) -> None:
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
        logger: logging.Logger,
    ) -> RunnerPreparation:
        ...

    def run(self, request: RunnerRequest) -> RunnerResult:
        ...
