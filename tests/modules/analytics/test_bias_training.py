"""Tests for bias-training builder orchestration/backoff."""

from __future__ import annotations

import logging

from cofolder.modules.analytics.bias_training import (
    _split_bias_build_phase_timings,
    run_build_bias_training_data,
)
from cofolder.modules.utils.timing import DebugTimingCollector


def _threshold_from_cmd(cmd: list[str]) -> float:
    idx = cmd.index("--ligand_similarity_threshold")
    return float(cmd[idx + 1])


def _assert_packaged_builder_command(cmd: list[str]) -> None:
    assert cmd[:4] == [
        cmd[0],
        "-u",
        "-m",
        "cofolder.modules.analytics.build_bias_training_data",
    ]
    assert "scripts/build_bias_training_data.py" not in " ".join(cmd)


def test_backoff_when_pre_cutoff_hits_are_zero(monkeypatch, temp_dir):
    calls: list[list[str]] = []

    def _fake_run(cmd):
        calls.append(list(cmd))
        if len(calls) == 1:
            return (
                0,
                "[info] selected_pdb_entries_pre_cutoff=0\n",
                [(0.2, "[info] selected_pdb_entries_pre_cutoff=0")],
            )
        return (
            0,
            (
                "[info] selected_pdb_entries_pre_cutoff=3\n"
                "[info] using protein query sequence from system.yaml\n"
                "[info] protein MMseqs hits kept: 5\n"
            ),
            [
                (0.3, "[info] selected_pdb_entries_pre_cutoff=3"),
                (0.5, "[info] using protein query sequence from system.yaml"),
                (1.1, "[info] protein MMseqs hits kept: 5"),
            ],
        )

    monkeypatch.setattr(
        "cofolder.modules.analytics.bias_training._run_bias_training_subprocess",
        _fake_run,
    )

    run_build_bias_training_data(
        system_path=temp_dir / "system.yaml",
        components_cif=temp_dir / "components.cif",
        output_protein_csv=temp_dir / "protein_training_data.csv",
        output_ligand_csv=temp_dir / "ligand_training_data.csv",
        release_cutoff="2023-06-01",
        ligand_similarity_threshold=0.35,
    )

    assert len(calls) == 2
    _assert_packaged_builder_command(calls[0])
    assert _threshold_from_cmd(calls[0]) == 0.35
    assert _threshold_from_cmd(calls[1]) == 0.3


def test_backoff_exhaustion_writes_header_only_training_csvs(monkeypatch, temp_dir):
    calls: list[list[str]] = []

    def _always_zero_pre_cutoff(cmd):
        calls.append(list(cmd))
        return (
            0,
            "[info] selected_pdb_entries_pre_cutoff=0\n",
            [(0.2, "[info] selected_pdb_entries_pre_cutoff=0")],
        )

    monkeypatch.setattr(
        "cofolder.modules.analytics.bias_training._run_bias_training_subprocess",
        _always_zero_pre_cutoff,
    )

    protein_out = temp_dir / "protein_training_data.csv"
    ligand_out = temp_dir / "ligand_training_data.csv"

    run_build_bias_training_data(
        system_path=temp_dir / "system.yaml",
        components_cif=temp_dir / "components.cif",
        output_protein_csv=protein_out,
        output_ligand_csv=ligand_out,
        release_cutoff="2023-06-01",
        ligand_similarity_threshold=0.05,
    )

    assert len(calls) == 2
    assert _threshold_from_cmd(calls[0]) == 0.05
    assert _threshold_from_cmd(calls[1]) == 0.0

    protein_header = protein_out.read_text(encoding="utf-8").strip().splitlines()[0]
    ligand_header = ligand_out.read_text(encoding="utf-8").strip().splitlines()[0]

    assert protein_header == "query_chain_id,pdb_id,release_date,sequence_similarity,sequence"
    assert ligand_header == "query_chain_id,pdb_id,release_date,ligand_id,ecfp_similarity,smiles"


def test_split_bias_build_phase_timings_detects_protein_boundary():
    ligand_elapsed, protein_elapsed = _split_bias_build_phase_timings(
        [
            (0.4, "[info] selected_ligand_chains=B"),
            (1.7, "[done] chain_ligand_csv[B]=/tmp/ligand.csv rows=25"),
            (2.0, "[info] using protein query sequence from system.yaml"),
            (5.4, "[info] protein MMseqs hits kept: 100"),
        ]
    )

    assert ligand_elapsed == 2.0
    assert protein_elapsed == 3.4


def test_split_bias_build_phase_timings_defaults_to_ligand_when_no_boundary():
    ligand_elapsed, protein_elapsed = _split_bias_build_phase_timings(
        [
            (0.3, "[info] selected_ligand_chains=B"),
            (1.2, "[info] selected_pdb_entries_pre_cutoff=0"),
        ]
    )

    assert ligand_elapsed == 1.2
    assert protein_elapsed == 0.0


def test_records_split_build_timings_when_debug_enabled(monkeypatch, temp_dir):
    def _fake_run(cmd):
        return (
            0,
            (
                "[done] chain_ligand_csv[B]=/tmp/ligand.csv rows=25\n"
                "[info] using protein query sequence from system.yaml\n"
                "[info] protein MMseqs hits kept: 100\n"
            ),
            [
                (1.5, "[done] chain_ligand_csv[B]=/tmp/ligand.csv rows=25"),
                (2.0, "[info] using protein query sequence from system.yaml"),
                (6.5, "[info] protein MMseqs hits kept: 100"),
            ],
        )

    monkeypatch.setattr(
        "cofolder.modules.analytics.bias_training._run_bias_training_subprocess",
        _fake_run,
    )

    logger = logging.getLogger("cofolder.test.bias_training")
    logger.setLevel(logging.DEBUG)
    collector = DebugTimingCollector(logger=logger)

    run_build_bias_training_data(
        system_path=temp_dir / "system.yaml",
        components_cif=temp_dir / "components.cif",
        output_protein_csv=temp_dir / "protein_training_data.csv",
        output_ligand_csv=temp_dir / "ligand_training_data.csv",
        release_cutoff="2023-06-01",
        ligand_similarity_threshold=0.35,
        timings=collector,
    )

    assert ("bias.training_data.build.ligand", 2.0) in collector._entries
    assert ("bias.training_data.build.protein", 4.5) in collector._entries
