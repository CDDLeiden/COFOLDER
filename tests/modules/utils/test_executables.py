"""Tests for shared executable discovery."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from cofolder.modules.analytics import bias_training
from cofolder.tools import build_bias_training_data
from cofolder.modules.utils.executables import resolve_mmseqs_executable
from cofolder.tools import fetch_bias_training_data


def _executable(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    path.chmod(0o755)
    return path


def test_precedence_preserves_every_lookup_route(temp_dir):
    home = temp_dir / "home"
    checkout = temp_dir / "checkout"
    explicit = _executable(temp_dir / "explicit/mmseqs")
    environment = _executable(temp_dir / "environment/mmseqs")
    user_vendor = _executable(home / ".cofolder/vendor/mmseqs/bin/mmseqs")
    source_bin = _executable(checkout / "vendor/mmseqs/bin/mmseqs")
    source_root = _executable(checkout / "vendor/mmseqs/mmseqs")
    path_binary = _executable(temp_dir / "path/mmseqs")

    def which(command: str) -> str | None:
        return str(path_binary) if command == "mmseqs" else None

    env = {"COFOLDER_MMSEQS_BIN": str(environment)}

    assert resolve_mmseqs_executable(
        str(explicit), environ=env, home=home, source_root=checkout, which=which
    ).source == "explicit argument"
    assert resolve_mmseqs_executable(
        environ=env, home=home, source_root=checkout, which=which
    ).source == "COFOLDER_MMSEQS_BIN"

    env = {}
    assert resolve_mmseqs_executable(
        environ=env, home=home, source_root=checkout, which=which
    ).source == "user vendor"
    user_vendor.unlink()
    assert resolve_mmseqs_executable(
        environ=env, home=home, source_root=checkout, which=which
    ).source == "source vendor bin"
    source_bin.unlink()
    assert resolve_mmseqs_executable(
        environ=env, home=home, source_root=checkout, which=which
    ).source == "source vendor root"
    source_root.unlink()
    assert resolve_mmseqs_executable(
        environ=env, home=home, source_root=checkout, which=which
    ).source == "PATH"


def test_explicit_command_name_has_precedence_over_environment(temp_dir):
    explicit = _executable(temp_dir / "commands/mmseqs-explicit")
    environment = _executable(temp_dir / "environment/mmseqs")
    resolution = resolve_mmseqs_executable(
        "mmseqs-explicit",
        environ={"COFOLDER_MMSEQS_BIN": str(environment)},
        home=temp_dir / "home",
        source_root=None,
        which=lambda command: str(explicit) if command == "mmseqs-explicit" else None,
    )

    assert resolution.path == str(explicit.resolve())
    assert resolution.source == "explicit argument"


def test_rejected_candidates_fall_through_with_diagnostics(temp_dir, caplog):
    missing = temp_dir / "missing/mmseqs"
    directory = temp_dir / "directory/mmseqs"
    directory.mkdir(parents=True)
    not_executable = (
        temp_dir / "not-executable/.cofolder/vendor/mmseqs/bin/mmseqs"
    )
    not_executable.parent.mkdir(parents=True)
    not_executable.write_text("binary", encoding="utf-8")
    selected = _executable(temp_dir / "path with spaces/mmseqs")

    with caplog.at_level(logging.DEBUG, logger="cofolder.modules.utils.executables"):
        resolution = resolve_mmseqs_executable(
            str(missing),
            environ={"COFOLDER_MMSEQS_BIN": str(directory)},
            home=temp_dir / "not-executable",
            source_root=None,
            which=lambda command: str(selected) if command == "mmseqs" else None,
        )

    assert resolution.path == str(selected.resolve())
    assert [candidate.outcome for candidate in resolution.candidates] == [
        "missing",
        "not a regular file",
        "not executable",
        "selected",
    ]
    assert "Rejected MMseqs candidate" in caplog.text
    assert "Selected MMseqs executable from PATH" in caplog.text
    assert "path with spaces" in resolution.diagnostic_text()


def test_duplicate_candidates_are_considered_once(temp_dir):
    missing = temp_dir / "same/mmseqs"
    selected = _executable(temp_dir / "selected/mmseqs")
    resolution = resolve_mmseqs_executable(
        str(missing),
        environ={"COFOLDER_MMSEQS_BIN": str(missing)},
        home=temp_dir / "home",
        source_root=None,
        which=lambda command: str(selected) if command == "mmseqs" else None,
    )

    assert sum(candidate.value == str(missing) for candidate in resolution.candidates) == 1


def test_executable_symlink_and_tilde_expansion(monkeypatch, temp_dir):
    target = _executable(temp_dir / "real/mmseqs")
    link = temp_dir / "home/bin/mmseqs"
    link.parent.mkdir(parents=True)
    link.symlink_to(target)
    monkeypatch.setenv("HOME", str(temp_dir / "home"))

    resolution = resolve_mmseqs_executable(
        "~/bin/mmseqs", environ={}, home=temp_dir / "unused", source_root=None
    )

    assert resolution.path == str(target.resolve())


def test_installed_layout_does_not_probe_source_vendor(temp_dir):
    resolution = resolve_mmseqs_executable(
        environ={},
        home=temp_dir / "home",
        source_root=None,
        which=lambda _command: None,
    )

    assert resolution.path is None
    assert all(
        not candidate.source.startswith("source vendor")
        for candidate in resolution.candidates
    )
    assert "unresolved command" in resolution.diagnostic_text()


@pytest.mark.parametrize(
    "resolver",
    [
        bias_training._resolve_mmseqs_bin,
        build_bias_training_data._resolve_mmseqs_bin,
        fetch_bias_training_data._resolve_mmseqs_bin,
    ],
)
def test_existing_resolution_seams_share_the_environment_candidate(
    resolver, monkeypatch, temp_dir
):
    executable = _executable(temp_dir / "shared path/mmseqs")
    monkeypatch.setenv("COFOLDER_MMSEQS_BIN", str(executable))

    assert resolver() == str(executable.resolve())


def test_bias_training_launcher_preserves_space_containing_executable_as_one_argument(
    monkeypatch, temp_dir
):
    executable = _executable(temp_dir / "space containing path/mmseqs")
    commands: list[list[str]] = []
    monkeypatch.setattr(
        bias_training,
        "_resolve_mmseqs_bin",
        lambda: str(executable),
    )
    monkeypatch.setattr(
        bias_training,
        "_run_bias_training_subprocess",
        lambda command: (commands.append(command) or (0, "", [])),
    )

    bias_training.run_build_bias_training_data(
        system_path=temp_dir / "system.yaml",
        components_cif=temp_dir / "components.cif",
        output_protein_csv=temp_dir / "protein.csv",
        output_ligand_csv=temp_dir / "ligand.csv",
        release_cutoff="2023-06-01",
    )

    index = commands[0].index("--mmseqs_bin")
    assert commands[0][index + 1] == str(executable)


def test_missing_required_executable_reports_rejection_reasons(
    monkeypatch, temp_dir
):
    candidate = temp_dir / "not-executable/mmseqs"
    candidate.parent.mkdir(parents=True)
    candidate.write_text("binary", encoding="utf-8")
    monkeypatch.setenv("COFOLDER_MMSEQS_BIN", str(candidate))
    monkeypatch.setattr(
        fetch_bias_training_data, "_resolve_mmseqs_bin", lambda _value: None
    )

    with pytest.raises(RuntimeError, match="not executable"):
        fetch_bias_training_data._run_mmseqs_databases(
            "PDB", temp_dir / "db", temp_dir / "tmp"
        )


def test_optional_resolution_returns_unavailable_without_execution(temp_dir):
    calls: list[str] = []

    resolution = resolve_mmseqs_executable(
        environ={},
        home=temp_dir / "home",
        source_root=None,
        which=lambda command: calls.append(command) or None,
    )

    assert calls == ["mmseqs"]
    assert resolution.path is None
