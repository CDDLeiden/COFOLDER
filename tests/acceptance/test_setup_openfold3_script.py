from __future__ import annotations

import os
import subprocess
from pathlib import Path
import sys


def _write_fake_setup_openfold(
    bin_dir: Path,
    output_path: Path,
    args_path: Path | None = None,
) -> Path:
    script_path = bin_dir / "setup_openfold"
    script_path.write_text(
        "\n".join(
            [
                "#!/usr/bin/env bash",
                "set -euo pipefail",
                *(
                    [f'printf "%s\\n" "$*" > "{args_path}"']
                    if args_path is not None
                    else []
                ),
                'printf "Please specify the OpenFold cache directory (default: /home/remco/.openfold3): "',
                "read -r cache_dir",
                'printf "Please specify the directory for parameter download (default: /home/remco/.openfold3): "',
                "read -r parameter_dir",
                'printf "Enter your choice (1/2/3, default: 1): "',
                "read -r parameter_choice",
                'printf "Force re-download parameters even if they already exist? (yes/no, default: no) "',
                "read -r force_download",
                'printf "Run integration tests? (yes/no) "',
                "read -r integration_choice",
                f'printf "%s\\n%s\\n%s\\n%s\\n%s\\n" "$cache_dir" "$parameter_dir" "$parameter_choice" "$force_download" "$integration_choice" > "{output_path}"',
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    script_path.chmod(0o755)
    return script_path


def test_setup_openfold3_script_answers_standard_prompts_and_forwards_later_input(temp_dir):
    fake_bin = temp_dir / "bin"
    fake_bin.mkdir()
    output_path = temp_dir / "captured.txt"
    _write_fake_setup_openfold(fake_bin, output_path)

    cache_root = temp_dir / "cache-root"
    env = dict(os.environ)
    env["PATH"] = f"{fake_bin}:{env['PATH']}"
    env["OPENFOLD_CACHE"] = str(cache_root)

    result = subprocess.run(
        ["bash", "scripts/setup_openfold3.sh"],
        cwd=Path(__file__).resolve().parents[2],
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "Proceeding with the default OpenFold3 checkpoint setup" in result.stdout
    assert "Upstream force-download selection: no." in result.stdout
    assert "Upstream integration tests selection: no." in result.stdout
    assert output_path.read_text(encoding="utf-8").splitlines() == [
        str(cache_root),
        str(cache_root),
        "1",
        "no",
        "no",
    ]


def test_installed_setup_openfold3_command_answers_standard_prompts(temp_dir):
    fake_bin = temp_dir / "bin"
    fake_bin.mkdir()
    output_path = temp_dir / "captured.txt"
    args_path = temp_dir / "args.txt"
    _write_fake_setup_openfold(fake_bin, output_path, args_path)

    cache_root = temp_dir / "installed-cache-root"
    env = dict(os.environ)
    env["PATH"] = f"{fake_bin}:{env['PATH']}"
    env["OPENFOLD_CACHE"] = str(cache_root)

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "cofolder.tools",
            "setup-openfold3",
            "--",
            "--upstream-option",
            "value with spaces",
        ],
        cwd=Path(__file__).resolve().parents[2],
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "Proceeding with the default OpenFold3 checkpoint setup" in result.stdout
    assert output_path.read_text(encoding="utf-8").splitlines() == [
        str(cache_root),
        str(cache_root),
        "1",
        "no",
        "no",
    ]
    assert args_path.read_text(encoding="utf-8").strip() == (
        "--upstream-option value with spaces"
    )


def test_setup_openfold3_script_respects_parameter_choice_override(temp_dir):
    fake_bin = temp_dir / "bin"
    fake_bin.mkdir()
    output_path = temp_dir / "captured.txt"
    _write_fake_setup_openfold(fake_bin, output_path)

    cache_root = temp_dir / "cache-root"
    env = dict(os.environ)
    env["PATH"] = f"{fake_bin}:{env['PATH']}"
    env["OPENFOLD_CACHE"] = str(cache_root)
    env["OPENFOLD3_PARAMETER_CHOICE"] = "2"

    result = subprocess.run(
        ["bash", "scripts/setup_openfold3.sh"],
        cwd=Path(__file__).resolve().parents[2],
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "Proceeding with all published OpenFold3 checkpoints" in result.stdout
    lines = output_path.read_text(encoding="utf-8").splitlines()
    assert lines[2] == "2"
    assert lines[3] == "no"
    assert lines[4] == "no"


def test_setup_openfold3_script_respects_integration_test_override(temp_dir):
    fake_bin = temp_dir / "bin"
    fake_bin.mkdir()
    output_path = temp_dir / "captured.txt"
    _write_fake_setup_openfold(fake_bin, output_path)

    cache_root = temp_dir / "cache-root"
    env = dict(os.environ)
    env["PATH"] = f"{fake_bin}:{env['PATH']}"
    env["OPENFOLD_CACHE"] = str(cache_root)
    env["OPENFOLD3_RUN_INTEGRATION_TESTS"] = "yes"

    result = subprocess.run(
        ["bash", "scripts/setup_openfold3.sh"],
        cwd=Path(__file__).resolve().parents[2],
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "Upstream integration tests selection: yes." in result.stdout
    assert output_path.read_text(encoding="utf-8").splitlines()[4] == "yes"


def test_setup_openfold3_script_respects_force_download_override(temp_dir):
    fake_bin = temp_dir / "bin"
    fake_bin.mkdir()
    output_path = temp_dir / "captured.txt"
    _write_fake_setup_openfold(fake_bin, output_path)

    cache_root = temp_dir / "cache-root"
    env = dict(os.environ)
    env["PATH"] = f"{fake_bin}:{env['PATH']}"
    env["OPENFOLD_CACHE"] = str(cache_root)
    env["OPENFOLD3_FORCE_DOWNLOAD_PARAMETERS"] = "yes"

    result = subprocess.run(
        ["bash", "scripts/setup_openfold3.sh"],
        cwd=Path(__file__).resolve().parents[2],
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "Upstream force-download selection: yes." in result.stdout
    assert output_path.read_text(encoding="utf-8").splitlines()[3] == "yes"
