#!/usr/bin/env python3
"""Build and smoke-test the COFOLDER sdist-to-wheel release path."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
import os
from pathlib import Path, PurePosixPath
import subprocess
import sys
import tarfile
import tempfile
import venv
import zipfile


REPOSITORY = Path(__file__).resolve().parents[1]
EXPECTED_VERSION = "1.0.0"
TOOL_COMMANDS = (
    "setup-openfold3",
    "setup-boltz2-cache",
    "populate-ccd-cache",
    "install-mmseqs",
    "fetch-bias-training-data",
    "build-bias-training-data",
)


def _run(command: Sequence[str], *, cwd: Path) -> None:
    print("+ " + " ".join(command), flush=True)
    subprocess.run(command, cwd=cwd, check=True)


def _single_match(directory: Path, pattern: str) -> Path:
    matches = sorted(directory.glob(pattern))
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected exactly one {pattern!r} in {directory}, found {len(matches)}"
        )
    return matches[0]


def _assert_safe_sdist_members(members: Sequence[tarfile.TarInfo]) -> None:
    for member in members:
        path = PurePosixPath(member.name)
        if path.is_absolute() or ".." in path.parts or member.issym() or member.islnk():
            raise RuntimeError(f"Unsafe sdist member: {member.name}")


def _verify_sdist(sdist: Path, source_root: Path) -> str:
    with tarfile.open(sdist, "r:gz") as archive:
        members = archive.getmembers()
        _assert_safe_sdist_members(members)
        names = {member.name for member in members if member.isfile()}

    roots = {PurePosixPath(name).parts[0] for name in names}
    if len(roots) != 1:
        raise RuntimeError(f"Expected one sdist root, found: {sorted(roots)}")
    root = roots.pop()
    required = {
        f"{root}/LICENSE",
        f"{root}/THIRD_PARTY_SOFTWARE.md",
        f"{root}/README.md",
        f"{root}/pyproject.toml",
        f"{root}/.github/workflows/quality.yml",
        f"{root}/scripts/run_test_lane.py",
        f"{root}/scripts/verify_release_artifacts.py",
        f"{root}/tests/fixtures/structural_parity_s6.json",
    }
    required.update(
        f"{root}/{path.relative_to(source_root).as_posix()}"
        for path in (source_root / "tests").rglob("test_*.py")
        if "__pycache__" not in path.parts
        and path.name != "ui_development_check.py"
    )
    missing = required - names
    if missing:
        raise RuntimeError(f"Missing required sdist files: {sorted(missing)}")

    forbidden = (
        f"{root}/src/cofolder/ui/",
        f"{root}/tests/ui_development_check.py",
        f"{root}/legacy/",
        f"{root}/.streamlit/",
    )
    leaked = sorted(name for name in names if any(item in name for item in forbidden))
    if leaked:
        raise RuntimeError(f"Unsupported files leaked into sdist: {leaked}")
    return root


def _extract_sdist(sdist: Path, destination: Path) -> Path:
    with tarfile.open(sdist, "r:gz") as archive:
        members = archive.getmembers()
        _assert_safe_sdist_members(members)
        archive.extractall(destination)
    directories = [path for path in destination.iterdir() if path.is_dir()]
    if len(directories) != 1:
        raise RuntimeError(
            f"Expected one extracted source directory, found {len(directories)}"
        )
    return directories[0]


def _verify_wheel(wheel: Path) -> None:
    with zipfile.ZipFile(wheel) as archive:
        names = set(archive.namelist())

    required_suffixes = {
        ".dist-info/licenses/LICENSE",
        ".dist-info/licenses/THIRD_PARTY_SOFTWARE.md",
    }
    if "cofolder/__init__.py" not in names:
        raise RuntimeError("Wheel does not contain the cofolder package")
    for suffix in required_suffixes:
        if not any(name.endswith(suffix) for name in names):
            raise RuntimeError(f"Wheel is missing {suffix}")
    if not any(name.startswith("cofolder/acceptance/data/") for name in names):
        raise RuntimeError("Wheel is missing acceptance resources")
    if not any(name.startswith("cofolder/resources/examples/") for name in names):
        raise RuntimeError("Wheel is missing runnable examples")
    leaked = sorted(
        name for name in names if name.startswith("cofolder/ui/") or "streamlit" in name.lower()
    )
    if leaked:
        raise RuntimeError(f"Unsupported UI files leaked into wheel: {leaked}")


def _venv_python(environment: Path) -> Path:
    directory = "Scripts" if os.name == "nt" else "bin"
    executable = "python.exe" if os.name == "nt" else "python"
    return environment / directory / executable


def _installed_smoke(wheel: Path, workspace: Path) -> None:
    environment = workspace / "installed"
    venv.EnvBuilder(with_pip=True, system_site_packages=True).create(environment)
    python = _venv_python(environment)
    bin_directory = python.parent
    cofolder = bin_directory / ("cofolder.exe" if os.name == "nt" else "cofolder")
    tools = bin_directory / ("cofolder-tools.exe" if os.name == "nt" else "cofolder-tools")
    outside_checkout = workspace / "outside-checkout"
    outside_checkout.mkdir()

    _run(
        [
            str(python),
            "-m",
            "pip",
            "install",
            "--force-reinstall",
            "--no-deps",
            str(wheel),
        ],
        cwd=outside_checkout,
    )
    metadata_check = "\n".join(
        [
            "from importlib.metadata import metadata, version",
            "from pathlib import Path",
            "import cofolder, sys",
            f"assert version('cofolder') == '{EXPECTED_VERSION}'",
            f"assert cofolder.__version__ == '{EXPECTED_VERSION}'",
            "assert metadata('cofolder')['License-Expression'] == 'MIT'",
            "assert Path(cofolder.__file__).resolve().is_relative_to(Path(sys.prefix).resolve())",
        ]
    )
    _run([str(python), "-c", metadata_check], cwd=outside_checkout)
    _run([str(cofolder), "--version"], cwd=outside_checkout)
    _run([str(python), "-m", "cofolder", "--version"], cwd=outside_checkout)
    _run([str(tools), "--help"], cwd=outside_checkout)
    for command in TOOL_COMMANDS:
        _run([str(tools), command, "--help"], cwd=outside_checkout)
    examples = outside_checkout / "copied examples"
    _run([str(tools), "copy-examples", str(examples)], cwd=outside_checkout)
    if not (examples / "system.yaml").is_file():
        raise RuntimeError("Installed example-copy smoke did not produce system.yaml")


def verify(repository: Path, workspace: Path) -> None:
    sdist_directory = workspace / "sdist"
    extracted_directory = workspace / "extracted"
    wheel_directory = workspace / "wheel"
    sdist_directory.mkdir(parents=True)
    extracted_directory.mkdir()
    wheel_directory.mkdir()

    _run(
        [sys.executable, "-m", "build", "--sdist", "--outdir", str(sdist_directory)],
        cwd=repository,
    )
    sdist = _single_match(sdist_directory, "*.tar.gz")
    _verify_sdist(sdist, repository)
    extracted_source = _extract_sdist(sdist, extracted_directory)
    _run(
        [
            sys.executable,
            "-m",
            "build",
            "--wheel",
            "--outdir",
            str(wheel_directory),
        ],
        cwd=extracted_source,
    )
    wheel = _single_match(wheel_directory, "*.whl")
    _verify_wheel(wheel)
    _installed_smoke(wheel, workspace)
    print(f"Verified {sdist.name} -> {wheel.name}", flush=True)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workspace",
        type=Path,
        help="Keep build and installation artifacts in this empty directory.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if args.workspace is not None:
        workspace = args.workspace.resolve()
        if workspace.exists() and any(workspace.iterdir()):
            raise SystemExit(f"Workspace must be empty: {workspace}")
        workspace.mkdir(parents=True, exist_ok=True)
        verify(REPOSITORY, workspace)
    else:
        with tempfile.TemporaryDirectory(prefix="cofolder-artifacts-") as temporary:
            verify(REPOSITORY, Path(temporary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
