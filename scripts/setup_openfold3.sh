#!/usr/bin/env bash

set -euo pipefail

OPENFOLD_CACHE="${OPENFOLD_CACHE:-$HOME/.openfold3}"
OPENFOLD3_PARAMETER_CHOICE="${OPENFOLD3_PARAMETER_CHOICE:-1}"
OPENFOLD3_PARAMETER_NAME="${OPENFOLD3_PARAMETER_NAME:-}"
OPENFOLD3_FORCE_DOWNLOAD_PARAMETERS="${OPENFOLD3_FORCE_DOWNLOAD_PARAMETERS:-no}"
OPENFOLD3_RUN_INTEGRATION_TESTS="${OPENFOLD3_RUN_INTEGRATION_TESTS:-no}"
export OPENFOLD_CACHE
export OPENFOLD3_PARAMETER_CHOICE
export OPENFOLD3_PARAMETER_NAME
export OPENFOLD3_FORCE_DOWNLOAD_PARAMETERS
export OPENFOLD3_RUN_INTEGRATION_TESTS

if [[ "$OPENFOLD3_PARAMETER_CHOICE" != "1" && "$OPENFOLD3_PARAMETER_CHOICE" != "2" && "$OPENFOLD3_PARAMETER_CHOICE" != "3" ]]; then
  echo "OPENFOLD3_PARAMETER_CHOICE must be one of: 1, 2, 3" >&2
  exit 1
fi

if [[ "$OPENFOLD3_RUN_INTEGRATION_TESTS" != "yes" && "$OPENFOLD3_RUN_INTEGRATION_TESTS" != "no" ]]; then
  echo "OPENFOLD3_RUN_INTEGRATION_TESTS must be one of: yes, no" >&2
  exit 1
fi

if [[ "$OPENFOLD3_FORCE_DOWNLOAD_PARAMETERS" != "yes" && "$OPENFOLD3_FORCE_DOWNLOAD_PARAMETERS" != "no" ]]; then
  echo "OPENFOLD3_FORCE_DOWNLOAD_PARAMETERS must be one of: yes, no" >&2
  exit 1
fi

if ! command -v setup_openfold >/dev/null 2>&1; then
  echo "setup_openfold was not found in PATH." >&2
  echo "Install the OpenFold3 extra first, for example: python -m pip install -e \".[openfold3]\"" >&2
  exit 1
fi

mkdir -p "$OPENFOLD_CACHE"
echo "Using OPENFOLD_CACHE=$OPENFOLD_CACHE"
echo "OpenFold3 parameter download choice: $OPENFOLD3_PARAMETER_CHOICE"
echo "Force OpenFold3 parameter download: $OPENFOLD3_FORCE_DOWNLOAD_PARAMETERS"
echo "Run OpenFold3 integration tests: $OPENFOLD3_RUN_INTEGRATION_TESTS"

python /dev/fd/3 "$@" 3<<'PY'
from __future__ import annotations

import os
import re
import select
import subprocess
import sys


cache_root = os.environ["OPENFOLD_CACHE"]
choice = os.environ["OPENFOLD3_PARAMETER_CHOICE"]
parameter_name = os.environ.get("OPENFOLD3_PARAMETER_NAME", "")
force_download_parameters = os.environ["OPENFOLD3_FORCE_DOWNLOAD_PARAMETERS"]
run_integration_tests = os.environ["OPENFOLD3_RUN_INTEGRATION_TESTS"]
command = ["setup_openfold", *sys.argv[1:]]

prompt_patterns: list[tuple[re.Pattern[str], str | None, str]] = [
    (
        re.compile(r"Please specify the OpenFold cache directory \(default: .*?\):\s*$"),
        cache_root + "\n",
        "cache_dir",
    ),
    (
        re.compile(r"Please specify the directory for parameter download \(default: .*?\):\s*$"),
        cache_root + "\n",
        "parameter_dir",
    ),
    (
        re.compile(r"Enter your choice \(1/2/3, default: [^)]+\):\s*$"),
        choice + "\n",
        "parameter_choice",
    ),
    (
        re.compile(
            r"Force re-download parameters even if they already exist\? "
            r"\((?:yes/no|y/n)(?:, default: [^)]+)?\)\s*$",
            re.IGNORECASE,
        ),
        force_download_parameters + "\n",
        "force_download_parameters",
    ),
    (
        re.compile(r"Run integration tests\? \((?:yes/no|y/n)\)\s*$", re.IGNORECASE),
        run_integration_tests + "\n",
        "integration_tests",
    ),
]
if choice == "3" and parameter_name:
    prompt_patterns.append(
        (
            re.compile(r"Enter .*parameter.*:\s*$", re.IGNORECASE),
            parameter_name + "\n",
            "parameter_name",
        )
    )

master_fd, slave_fd = os.openpty()
try:
    process = subprocess.Popen(
        command,
        stdin=slave_fd,
        stdout=slave_fd,
        stderr=slave_fd,
        text=False,
        close_fds=True,
    )
finally:
    os.close(slave_fd)

buffer = ""
answered: set[str] = set()


def emit_wrapper_status(key: str) -> None:
    if key == "parameter_choice":
        if choice == "1":
            sys.stdout.write(
                "\n[cofolder] Proceeding with the default OpenFold3 checkpoint setup. "
                "Download and extraction may take a while.\n"
            )
        elif choice == "2":
            sys.stdout.write(
                "\n[cofolder] Proceeding with all published OpenFold3 checkpoints. "
                "Download and extraction may take a while.\n"
            )
        else:
            sys.stdout.write(
                "\n[cofolder] Proceeding with a specific OpenFold3 checkpoint selection.\n"
            )
        sys.stdout.flush()
    elif key == "integration_tests":
        sys.stdout.write(
            f"\n[cofolder] Upstream integration tests selection: {run_integration_tests}.\n"
        )
        sys.stdout.flush()
    elif key == "force_download_parameters":
        sys.stdout.write(
            "\n[cofolder] Upstream force-download selection: "
            f"{force_download_parameters}.\n"
        )
        sys.stdout.flush()

stdin_open = True
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
            text = chunk.decode(errors="replace")
            sys.stdout.write(text)
            sys.stdout.flush()
            buffer = (buffer + text)[-8000:]
            for pattern, response, key in prompt_patterns:
                if key in answered:
                    continue
                if pattern.search(buffer):
                    if response is not None:
                        os.write(master_fd, response.encode())
                    emit_wrapper_status(key)
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

return_code = process.wait()
os.close(master_fd)
raise SystemExit(return_code)
PY
