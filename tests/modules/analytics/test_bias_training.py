"""Tests for bias-training builder orchestration/backoff."""

from __future__ import annotations

import subprocess

from cofolder.modules.analytics.bias_training import run_build_bias_training_data


def _threshold_from_cmd(cmd: list[str]) -> float:
    idx = cmd.index("--ligand_similarity_threshold")
    return float(cmd[idx + 1])


def test_backoff_when_pre_cutoff_hits_are_zero(monkeypatch, temp_dir):
    calls: list[list[str]] = []

    def _fake_run(cmd, check=False, capture_output=True, text=True):
        calls.append(list(cmd))
        if len(calls) == 1:
            return subprocess.CompletedProcess(
                cmd,
                0,
                stdout="[info] selected_pdb_entries_pre_cutoff=0\n",
                stderr="",
            )
        return subprocess.CompletedProcess(
            cmd,
            0,
            stdout="[info] selected_pdb_entries_pre_cutoff=3\n",
            stderr="",
        )

    monkeypatch.setattr("cofolder.modules.analytics.bias_training.subprocess.run", _fake_run)

    run_build_bias_training_data(
        system_path=temp_dir / "system.yaml",
        components_cif=temp_dir / "components.cif",
        output_protein_csv=temp_dir / "protein_training_data.csv",
        output_ligand_csv=temp_dir / "ligand_training_data.csv",
        release_cutoff="2023-06-01",
        ligand_similarity_threshold=0.35,
    )

    assert len(calls) == 2
    assert _threshold_from_cmd(calls[0]) == 0.35
    assert _threshold_from_cmd(calls[1]) == 0.3


def test_backoff_exhaustion_writes_header_only_training_csvs(monkeypatch, temp_dir):
    calls: list[list[str]] = []

    def _always_zero_pre_cutoff(cmd, check=False, capture_output=True, text=True):
        calls.append(list(cmd))
        return subprocess.CompletedProcess(
            cmd,
            0,
            stdout="[info] selected_pdb_entries_pre_cutoff=0\n",
            stderr="",
        )

    monkeypatch.setattr(
        "cofolder.modules.analytics.bias_training.subprocess.run",
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
