"""Prepare the OpenFold3 cache and checkpoints through upstream setup prompts."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
import os
from pathlib import Path
import re
import select
import shutil
import subprocess
import sys


def _print_help() -> None:
    parser = argparse.ArgumentParser(
        prog="cofolder-tools setup-openfold3",
        description=__doc__,
        epilog=(
            "Configuration uses OPENFOLD_CACHE, OPENFOLD3_PARAMETER_CHOICE, "
            "OPENFOLD3_PARAMETER_NAME, OPENFOLD3_FORCE_DOWNLOAD_PARAMETERS, and "
            "OPENFOLD3_RUN_INTEGRATION_TESTS. Arguments after `--` are forwarded "
            "to setup_openfold."
        ),
    )
    parser.print_help()


def _validated_environment() -> tuple[str, str, str, str, str]:
    cache_root = os.environ.get("OPENFOLD_CACHE", str(Path.home() / ".openfold3"))
    choice = os.environ.get("OPENFOLD3_PARAMETER_CHOICE", "1")
    parameter_name = os.environ.get("OPENFOLD3_PARAMETER_NAME", "")
    force_download = os.environ.get("OPENFOLD3_FORCE_DOWNLOAD_PARAMETERS", "no")
    run_integration_tests = os.environ.get("OPENFOLD3_RUN_INTEGRATION_TESTS", "no")
    if choice not in {"1", "2", "3"}:
        raise ValueError("OPENFOLD3_PARAMETER_CHOICE must be one of: 1, 2, 3")
    if run_integration_tests not in {"yes", "no"}:
        raise ValueError("OPENFOLD3_RUN_INTEGRATION_TESTS must be one of: yes, no")
    if force_download not in {"yes", "no"}:
        raise ValueError("OPENFOLD3_FORCE_DOWNLOAD_PARAMETERS must be one of: yes, no")
    return cache_root, choice, parameter_name, force_download, run_integration_tests


def _emit_status(key: str, choice: str, force_download: str, integration: str) -> None:
    if key == "parameter_choice":
        if choice == "1":
            message = (
                "Proceeding with the default OpenFold3 checkpoint setup. "
                "Download and extraction may take a while."
            )
        elif choice == "2":
            message = (
                "Proceeding with all published OpenFold3 checkpoints. "
                "Download and extraction may take a while."
            )
        else:
            message = "Proceeding with a specific OpenFold3 checkpoint selection."
        print(f"\n[cofolder] {message}", flush=True)
    elif key == "integration_tests":
        print(
            f"\n[cofolder] Upstream integration tests selection: {integration}.",
            flush=True,
        )
    elif key == "force_download_parameters":
        print(
            f"\n[cofolder] Upstream force-download selection: {force_download}.",
            flush=True,
        )


def _run_setup(command: list[str], responses: dict[str, str]) -> int:
    patterns = [
        (
            re.compile(r"Please specify the OpenFold cache directory \(default: .*?\):\s*$"),
            responses["cache_root"] + "\n",
            "cache_dir",
        ),
        (
            re.compile(r"Please specify the directory for parameter download \(default: .*?\):\s*$"),
            responses["cache_root"] + "\n",
            "parameter_dir",
        ),
        (
            re.compile(r"Enter your choice \(1/2/3, default: [^)]+\):\s*$"),
            responses["choice"] + "\n",
            "parameter_choice",
        ),
        (
            re.compile(
                r"Force re-download parameters even if they already exist\? "
                r"\((?:yes/no|y/n)(?:, default: [^)]+)?\)\s*$",
                re.IGNORECASE,
            ),
            responses["force_download"] + "\n",
            "force_download_parameters",
        ),
        (
            re.compile(r"Run integration tests\? \((?:yes/no|y/n)\)\s*$", re.IGNORECASE),
            responses["integration"] + "\n",
            "integration_tests",
        ),
    ]
    if responses["choice"] == "3" and responses["parameter_name"]:
        patterns.append(
            (
                re.compile(r"Enter .*parameter.*:\s*$", re.IGNORECASE),
                responses["parameter_name"] + "\n",
                "parameter_name",
            )
        )

    master_fd, slave_fd = os.openpty()
    try:
        process_environment = dict(os.environ)
        process_environment.update(
            {
                "OPENFOLD_CACHE": responses["cache_root"],
                "OPENFOLD3_PARAMETER_CHOICE": responses["choice"],
                "OPENFOLD3_PARAMETER_NAME": responses["parameter_name"],
                "OPENFOLD3_FORCE_DOWNLOAD_PARAMETERS": responses["force_download"],
                "OPENFOLD3_RUN_INTEGRATION_TESTS": responses["integration"],
            }
        )
        process = subprocess.Popen(
            command,
            stdin=slave_fd,
            stdout=slave_fd,
            stderr=slave_fd,
            text=False,
            close_fds=True,
            env=process_environment,
        )
    finally:
        os.close(slave_fd)

    buffer = ""
    answered: set[str] = set()
    stdin_open = True
    try:
        while True:
            read_fds = [master_fd]
            if process.poll() is None and stdin_open:
                read_fds.append(sys.stdin.fileno())
            ready, _, _ = select.select(read_fds, [], [], 0.1)

            if master_fd in ready:
                try:
                    chunk = os.read(master_fd, 4096)
                except OSError:
                    chunk = b""
                if chunk:
                    output = chunk.decode(errors="replace")
                    sys.stdout.write(output)
                    sys.stdout.flush()
                    buffer = (buffer + output)[-8000:]
                    for pattern, response, key in patterns:
                        if key in answered or not pattern.search(buffer):
                            continue
                        os.write(master_fd, response.encode())
                        _emit_status(
                            key,
                            responses["choice"],
                            responses["force_download"],
                            responses["integration"],
                        )
                        answered.add(key)
                elif process.poll() is not None:
                    break

            if sys.stdin.fileno() in ready:
                forwarded = os.read(sys.stdin.fileno(), 4096)
                if forwarded:
                    os.write(master_fd, forwarded)
                else:
                    stdin_open = False

            if process.poll() is not None and master_fd not in ready:
                break
        return process.wait()
    finally:
        os.close(master_fd)


def main(argv: Sequence[str] | None = None) -> int:
    upstream_args = list(argv or ())
    if upstream_args in (["-h"], ["--help"]):
        _print_help()
        return 0
    if upstream_args[:1] == ["--"]:
        upstream_args.pop(0)

    try:
        cache_root, choice, parameter_name, force_download, integration = (
            _validated_environment()
        )
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    if shutil.which("setup_openfold") is None:
        print("setup_openfold was not found in PATH.", file=sys.stderr)
        print(
            "Install the OpenFold3 extra first, for example: "
            'python -m pip install ".[openfold3]"',
            file=sys.stderr,
        )
        return 1

    Path(cache_root).expanduser().mkdir(parents=True, exist_ok=True)
    print(f"Using OPENFOLD_CACHE={cache_root}")
    print(f"OpenFold3 parameter download choice: {choice}")
    print(f"Force OpenFold3 parameter download: {force_download}")
    print(f"Run OpenFold3 integration tests: {integration}")
    return _run_setup(
        ["setup_openfold", *upstream_args],
        {
            "cache_root": cache_root,
            "choice": choice,
            "parameter_name": parameter_name,
            "force_download": force_download,
            "integration": integration,
        },
    )


if __name__ == "__main__":
    raise SystemExit(main())
