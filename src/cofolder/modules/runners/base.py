from __future__ import annotations

import logging
from dataclasses import dataclass, field
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


class Runner(Protocol):
    name: str
    capabilities: set[str]

    def is_available(self) -> tuple[bool, str | None]:
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
