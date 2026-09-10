#!/usr/bin/env python3
"""Run one exhaustive, non-overlapping COFOLDER pytest lane."""

from __future__ import annotations

import argparse
from collections.abc import Callable, Sequence
from pathlib import Path
import subprocess
import sys


REPOSITORY = Path(__file__).resolve().parents[1]
LANE_NAMES = ("core", "contracts-tutorial", "artifact", "acceptance")


def _under(path: Path, directory: str) -> bool:
    return path == Path(directory) or Path(directory) in path.parents


def _is_core(path: Path) -> bool:
    core_directories = (
        "tests/modules/analytics",
        "tests/modules/entities",
        "tests/modules/input",
        "tests/modules/runners",
        "tests/modules/utils",
        "tests/recipes",
        "tests/tools",
    )
    return path in {Path("tests/test_cli.py"), Path("tests/test_test_lanes.py")} or any(
        _under(path, directory) for directory in core_directories
    )


LANE_RULES: dict[str, Callable[[Path], bool]] = {
    "core": _is_core,
    "contracts-tutorial": lambda path: (
        _under(path, "tests/modules/contracts")
        or path == Path("tests/test_tutorials.py")
    ),
    "artifact": lambda path: path == Path("tests/test_packaging.py"),
    "acceptance": lambda path: _under(path, "tests/acceptance"),
}


class LaneConfigurationError(RuntimeError):
    """Raised when the supported test suite is not partitioned exactly once."""


def discover_supported_tests(repository: Path = REPOSITORY) -> tuple[Path, ...]:
    """Return repository-relative supported pytest modules."""
    tests_directory = repository / "tests"
    return tuple(
        path.relative_to(repository)
        for path in sorted(tests_directory.rglob("test_*.py"))
        if "__pycache__" not in path.parts
    )


def classify_tests(
    tests: Sequence[Path],
) -> dict[str, tuple[Path, ...]]:
    """Partition tests into lanes, rejecting missing or overlapping ownership."""
    classified: dict[str, list[Path]] = {name: [] for name in LANE_NAMES}
    errors: list[str] = []

    for path in tests:
        matches = [name for name, rule in LANE_RULES.items() if rule(path)]
        if len(matches) != 1:
            ownership = ", ".join(matches) if matches else "none"
            errors.append(f"{path.as_posix()}: {ownership}")
            continue
        classified[matches[0]].append(path)

    if errors:
        details = "\n".join(f"  - {error}" for error in errors)
        raise LaneConfigurationError(
            "Every supported test module must belong to exactly one lane:\n" + details
        )

    return {name: tuple(paths) for name, paths in classified.items()}


def lane_tests(
    lane: str, repository: Path = REPOSITORY
) -> tuple[Path, ...]:
    """Return the test modules for a named lane or the complete supported suite."""
    classified = classify_tests(discover_supported_tests(repository))
    if lane == "all":
        return tuple(
            sorted(path for paths in classified.values() for path in paths)
        )
    return classified[lane]


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("lane", choices=(*LANE_NAMES, "all"))
    parser.add_argument(
        "pytest_args",
        nargs=argparse.REMAINDER,
        help="Additional pytest arguments, optionally following '--'.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    pytest_args = list(args.pytest_args)
    if pytest_args[:1] == ["--"]:
        pytest_args.pop(0)

    tests = lane_tests(args.lane)
    print(f"COFOLDER test lane '{args.lane}': {len(tests)} modules", flush=True)
    command = [
        sys.executable,
        "-m",
        "pytest",
        "-q",
        *(path.as_posix() for path in tests),
        *pytest_args,
    ]
    return subprocess.run(command, cwd=REPOSITORY, check=False).returncode


if __name__ == "__main__":
    raise SystemExit(main())
