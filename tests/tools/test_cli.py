from __future__ import annotations

from pathlib import Path
import subprocess
import sys

import pytest

from cofolder.tools.cli import COMMANDS


REPOSITORY = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize("command", [None, *COMMANDS])
def test_help_paths_are_side_effect_free(command, temp_dir):
    argv = [sys.executable, "-m", "cofolder.tools"]
    if command is not None:
        argv.append(command)
    argv.append("--help")
    before = list(temp_dir.iterdir())

    result = subprocess.run(
        argv,
        cwd=temp_dir,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "usage:" in result.stdout
    assert list(temp_dir.iterdir()) == before


def test_copy_examples_command_reports_destination(temp_dir):
    destination = temp_dir / "copied examples"
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "cofolder.tools",
            "copy-examples",
            str(destination),
        ],
        cwd=REPOSITORY,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert f"copied_examples={destination.resolve()}" in result.stdout
    assert (destination / "system.yaml").is_file()


@pytest.mark.parametrize(
    ("command", "expected"),
    [
        (["bash", "scripts/install_mmseqs_vendor.sh", "--help"], "install-mmseqs"),
        ([sys.executable, "scripts/fetch_bias_training_data.py", "--help"], "output_root"),
    ],
)
def test_source_wrappers_remain_invokable(command, expected):
    result = subprocess.run(
        command,
        cwd=REPOSITORY,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert expected in result.stdout
