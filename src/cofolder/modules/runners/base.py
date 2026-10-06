from __future__ import annotations

import logging
import sys
from abc import ABC, abstractmethod
from importlib import metadata
from pathlib import Path
from typing import Any, Protocol

from packaging.version import InvalidVersion, Version

from cofolder.modules.contracts import (
    BackendVersionStatus,
    RepeatSeedProvenance,
    RunnerBackendIdentity,
)
from cofolder.modules.input.config import RunnerOptionsSchema, load_runner_options
from cofolder.modules.input.ligand import LigandPreparationCapabilities
from cofolder.modules.runners.contracts import (
    RunnerExecutionRequest,
    RunnerExecutionResult,
    RunnerInputCapabilities,
    RunnerModelSlot,
    RunnerPreparationResult,
)

class BaseRunner(ABC):
    name: str
    backend_name: str = ""
    backend_distribution: str = ""
    capabilities: set[str]
    input_capabilities = RunnerInputCapabilities(
        entity_types=frozenset(),
        constraint_types=frozenset(),
    )
    supports_msa_reuse = False
    options_schema: RunnerOptionsSchema
    ligand_preparation_capabilities = LigandPreparationCapabilities(
        native_smiles=True, conformer_modes=frozenset()
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

    def detect_backend_identity(self) -> RunnerBackendIdentity:
        """Detect backend package identity without making execution depend on it."""
        backend_name = self.backend_name or self.name
        distribution = self.backend_distribution or backend_name
        try:
            raw_version = metadata.version(distribution)
        except metadata.PackageNotFoundError:
            return RunnerBackendIdentity(
                runner_name=self.name,
                backend_name=backend_name,
                version=None,
                version_status=BackendVersionStatus.UNAVAILABLE,
                detail=f"Distribution {distribution!r} was not found.",
            )
        except Exception as exc:  # noqa: BLE001 - detection must remain non-fatal
            return RunnerBackendIdentity(
                runner_name=self.name,
                backend_name=backend_name,
                version=None,
                version_status=BackendVersionStatus.UNAVAILABLE,
                detail=f"Version lookup failed ({type(exc).__name__}).",
            )

        raw_version = str(raw_version)
        try:
            normalized = str(Version(raw_version))
        except InvalidVersion:
            return RunnerBackendIdentity(
                runner_name=self.name,
                backend_name=backend_name,
                version=None,
                version_status=BackendVersionStatus.UNPARSEABLE,
                raw_version=raw_version,
                detail="Distribution metadata did not contain a valid PEP 440 version.",
            )
        return RunnerBackendIdentity(
            runner_name=self.name,
            backend_name=backend_name,
            version=normalized,
            version_status=BackendVersionStatus.DETECTED,
            raw_version=raw_version,
        )

    def resolve_effective_seed(
        self, seed: RepeatSeedProvenance
    ) -> RepeatSeedProvenance:
        """Return the backend-facing seed; runners may explicitly adjust it."""
        return seed

    def execution_models(self, options_obj: Any) -> tuple[RunnerModelSlot, ...]:
        """Describe the model/sample axis before backend execution."""
        samples = getattr(options_obj, "diffusion_samples", None)
        find_value = getattr(options_obj, "find_value", None)
        if samples is None and callable(find_value):
            samples = find_value(key="diffusion_samples")
        samples = int(samples or 1)
        if samples < 1:
            raise ValueError("Runner diffusion_samples must be positive.")
        model_id = str(getattr(self, "model_name", None) or self.name)
        return tuple(
            RunnerModelSlot(model_id=model_id, sample_id=index)
            for index in range(samples)
        )

    @staticmethod
    def _check_python_compatibility(
        *, runner_name: str, maximum_exclusive: tuple[int, int]
    ) -> tuple[bool, str | None]:
        """Reject Python versions outside a backend's supported range."""
        current = sys.version_info[:2]
        if current < maximum_exclusive:
            return True, None

        current_label = ".".join(str(part) for part in current)
        maximum_label = ".".join(str(part) for part in maximum_exclusive)
        return False, (
            f"The selected '{runner_name}' runner does not support Python "
            f"{current_label}. It requires Python below {maximum_label}; use a "
            "Python 3.11 or 3.12 environment."
        )

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
            name
            for name in conflicting_distributions
            if BaseRunner._distribution_installed(name)
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
        source_path: Path | None = None,
        requirements: Any = None,
    ):
        """Validate the shared YAML contract against this runner's capabilities."""

        from cofolder.modules.input.validation import validate_system_input

        cache_path = None
        find_value = getattr(options_obj, "find_value", None)
        if callable(find_value):
            cache_path = find_value(key="cache") or find_value(key="cache_path")
        return validate_system_input(
            system_obj,
            source_path=source_path,
            runner_name=self.name,
            capabilities=self.input_capabilities,
            requirements=requirements,
            cache_path=cache_path,
            check_atom_names=check_atom_names,
        )

    def _load_typed_options(self, options_path: Path):
        return load_runner_options(options_path, schema=self.options_schema)

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

    def inject_reusable_msas(
        self,
        system_obj: Any,
        cache_dir: Path,
        *,
        settings: dict[str, Any] | None = None,
    ) -> int:
        """Inject runner-compatible cached MSAs, returning the entity count."""

        return 0

    def msa_reuse_settings(self, options_obj: Any) -> dict[str, Any]:
        """Return non-secret settings that affect reusable MSA generation."""

        return {}

    def capture_reusable_msas(
        self,
        system_obj: Any,
        *,
        generated_dir: Path,
        cache_dir: Path,
        settings: dict[str, Any] | None = None,
    ) -> int:
        """Persist runner-generated MSAs, returning the captured entity count."""

        return 0

    @abstractmethod
    def load_options(self, options_path: Path):
        """Load any runner-specific options from disk."""

    @abstractmethod
    def run(self, request: RunnerExecutionRequest) -> RunnerExecutionResult:
        """Execute the backend and return canonical normalized outputs."""


class Runner(Protocol):
    name: str
    backend_name: str
    backend_distribution: str
    capabilities: set[str]
    input_capabilities: RunnerInputCapabilities
    supports_msa_reuse: bool
    options_schema: RunnerOptionsSchema
    ligand_preparation_capabilities: LigandPreparationCapabilities

    def check_availability(self) -> tuple[bool, str | None]: ...

    def ensure_available(self) -> None: ...

    def detect_backend_identity(self) -> RunnerBackendIdentity: ...

    def resolve_effective_seed(
        self, seed: RepeatSeedProvenance
    ) -> RepeatSeedProvenance: ...

    def execution_models(self, options_obj: Any) -> tuple[RunnerModelSlot, ...]: ...

    def validate_system(
        self,
        system_obj: Any,
        options_obj: Any,
        *,
        check_atom_names: bool = True,
        source_path: Path | None = None,
        requirements: Any = None,
    ): ...

    def load_options(self, options_path: Path): ...

    def prepare_system(
        self,
        system_obj: Any,
        options_obj: Any,
        wrk_dir: Path,
        conformers: str | None,
        sdf_file: Path | None,
        logger: logging.Logger | None,
    ) -> RunnerPreparationResult: ...

    def inject_reusable_msas(
        self,
        system_obj: Any,
        cache_dir: Path,
        *,
        settings: dict[str, Any] | None = None,
    ) -> int: ...

    def msa_reuse_settings(self, options_obj: Any) -> dict[str, Any]: ...

    def capture_reusable_msas(
        self,
        system_obj: Any,
        *,
        generated_dir: Path,
        cache_dir: Path,
        settings: dict[str, Any] | None = None,
    ) -> int: ...

    def run(self, request: RunnerExecutionRequest) -> RunnerExecutionResult: ...
