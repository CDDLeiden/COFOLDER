"""Typed construction of Boltz backend command arguments."""

from __future__ import annotations

from pathlib import Path

from cofolder.modules.input.config import BOLTZ_RUNNER_FIELDS, RunnerOptions


def _append_option(argv: list[str], name: str, value: object) -> None:
    """Append one validated option without conflating numbers and booleans."""

    if value is None or value is False:
        return
    argv.append(f"--{name}")
    if value is not True:
        argv.append(str(value))


def build_boltz_command(
    *,
    options: RunnerOptions,
    system_path: Path | str,
    output_dir: Path | str,
    seed: int,
    model_name: str | None,
    use_msa_server: bool,
) -> list[str]:
    """Build deterministic ``boltz predict`` argv from validated options."""

    argv = [
        "boltz",
        "predict",
        str(system_path),
        "--out_dir",
        str(output_dir),
        "--seed",
        str(seed),
    ]
    if use_msa_server:
        argv.append("--use_msa_server")

    _append_option(
        argv,
        "cache",
        str(options.runtime.cache_path) if options.runtime.cache_path else None,
    )
    _append_option(argv, "diffusion_samples", options.runtime.diffusion_samples)

    # Schema order is the canonical order, so equivalent YAML documents produce
    # identical commands regardless of their source mapping order.
    for name in BOLTZ_RUNNER_FIELDS:
        if name in options.runner:
            _append_option(argv, name, options.runner[name])

    _append_option(argv, "model", model_name)
    return argv
