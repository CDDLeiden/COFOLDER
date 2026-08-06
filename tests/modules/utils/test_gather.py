"""Tests for cofolder.modules.utils.gather alignment orchestration."""

from pathlib import Path

import pandas as pd

from cofolder.modules.utils import gather


def _make_system_df():
    return pd.DataFrame(
        [
            {"idx": 0, "cif_file": "1_system_model_0.cif", "model_name": "system", "repeat": 1, "diffusion_sample": 0},
            {"idx": 1, "cif_file": "1_system_model_1.cif", "model_name": "system", "repeat": 1, "diffusion_sample": 1},
            {"idx": 2, "cif_file": "2_system_model_0.cif", "model_name": "system", "repeat": 2, "diffusion_sample": 0},
        ]
    )


def _make_chain_df():
    return pd.DataFrame(
        [
            {"idx": 0, "CHAIN_ID": "A", "ENTITY_TYPE": "protein", "ligand_molecule_id": "protein_A", "cif_file": "1_system_model_0.cif", "model_name": "system", "repeat": 1, "diffusion_sample": 0},
            {"idx": 1, "CHAIN_ID": "L", "ENTITY_TYPE": "ligand", "ligand_molecule_id": "LIG", "cif_file": "1_system_model_0.cif", "model_name": "system", "repeat": 1, "diffusion_sample": 0},
            {"idx": 2, "CHAIN_ID": "A", "ENTITY_TYPE": "protein", "ligand_molecule_id": "protein_A", "cif_file": "1_system_model_1.cif", "model_name": "system", "repeat": 1, "diffusion_sample": 1},
            {"idx": 3, "CHAIN_ID": "L", "ENTITY_TYPE": "ligand", "ligand_molecule_id": "LIG", "cif_file": "1_system_model_1.cif", "model_name": "system", "repeat": 1, "diffusion_sample": 1},
            {"idx": 4, "CHAIN_ID": "A", "ENTITY_TYPE": "protein", "ligand_molecule_id": "protein_A", "cif_file": "2_system_model_0.cif", "model_name": "system", "repeat": 2, "diffusion_sample": 0},
            {"idx": 5, "CHAIN_ID": "L", "ENTITY_TYPE": "ligand", "ligand_molecule_id": "LIG", "cif_file": "2_system_model_0.cif", "model_name": "system", "repeat": 2, "diffusion_sample": 0},
        ]
    )


def test_gather_robustness_aligns_all_runs_to_reference_when_provided(monkeypatch, temp_dir):
    captured = {}

    def _fake_load_structures(cif_paths):
        captured["cif_paths"] = [Path(p).name for p in cif_paths]
        return [(Path(p).name, object()) for p in cif_paths]

    def _fake_align_structures_on_protein_ca(structures, save_dir=None):
        captured["aligned_names"] = [name for name, _ in structures]
        captured["save_dir"] = save_dir
        return ({name: struct for name, struct in structures}, structures[0][1], structures[0][0])

    monkeypatch.setattr("cofolder.modules.utils.gather.align._load_structures", _fake_load_structures)
    monkeypatch.setattr(
        "cofolder.modules.utils.gather.align._align_structures_on_protein_ca",
        _fake_align_structures_on_protein_ca,
    )
    monkeypatch.setattr("cofolder.modules.utils.gather.align._compute_chain_rmsd", lambda chain_df, aligned_structs, wrk_dir: {})
    monkeypatch.setattr("cofolder.modules.utils.gather.align._compute_ligand_rmsd", lambda chain_df, aligned_structs, wrk_dir: {})

    reference_path = temp_dir / "reference.pdb"
    reference_path.write_text("HEADER TEST\n", encoding="utf-8")

    gather.gather_robustness_results(
        system_df=_make_system_df(),
        chain_df=_make_chain_df(),
        wrk_dir=temp_dir,
        reference_path=reference_path,
    )

    assert captured["cif_paths"] == [
        "reference.pdb",
        "1_system_model_0.cif",
        "1_system_model_1.cif",
        "2_system_model_0.cif",
    ]
    assert captured["aligned_names"] == captured["cif_paths"]
    assert captured["save_dir"] == temp_dir / "results" / "structures_aligned"


def test_gather_robustness_aligns_all_runs_to_first_predicted_when_no_reference(monkeypatch, temp_dir):
    captured = {}

    def _fake_load_structures(cif_paths):
        captured["cif_paths"] = [Path(p).name for p in cif_paths]
        return [(Path(p).name, object()) for p in cif_paths]

    def _fake_align_structures_on_protein_ca(structures, save_dir=None):
        captured["aligned_names"] = [name for name, _ in structures]
        return ({name: struct for name, struct in structures}, structures[0][1], structures[0][0])

    monkeypatch.setattr("cofolder.modules.utils.gather.align._load_structures", _fake_load_structures)
    monkeypatch.setattr(
        "cofolder.modules.utils.gather.align._align_structures_on_protein_ca",
        _fake_align_structures_on_protein_ca,
    )
    monkeypatch.setattr("cofolder.modules.utils.gather.align._compute_chain_rmsd", lambda chain_df, aligned_structs, wrk_dir: {})
    monkeypatch.setattr("cofolder.modules.utils.gather.align._compute_ligand_rmsd", lambda chain_df, aligned_structs, wrk_dir: {})

    gather.gather_robustness_results(
        system_df=_make_system_df(),
        chain_df=_make_chain_df(),
        wrk_dir=temp_dir,
        reference_path=None,
    )

    assert captured["cif_paths"] == [
        "1_system_model_0.cif",
        "1_system_model_1.cif",
        "2_system_model_0.cif",
    ]
    assert captured["aligned_names"] == captured["cif_paths"]
    assert captured["aligned_names"][0] == "1_system_model_0.cif"
