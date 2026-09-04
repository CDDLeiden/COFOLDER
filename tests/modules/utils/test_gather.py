"""Tests for cofolder.modules.utils.gather alignment orchestration."""

from pathlib import Path
import warnings

import pandas as pd

from cofolder.modules.utils import gather
from cofolder.modules.input.system import System


def test_bitstring_similarity_with_one_pair_has_no_runtime_warning(temp_dir):
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        result = gather.assess_bitstring_similarity(
            [[1, 0], [1, 1]],
            prefix="ifp_distance",
            id="ligand",
            wrk_dir=temp_dir,
        )

    assert result == {
        "ifp_distance_mean": 0.5,
        "ifp_distance_std": None,
    }
    assert (
        temp_dir
        / "results"
        / "matrices"
        / "similarity_matrix_ifp_distance_ligand.csv"
    ).exists()


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


def test_add_chain_info_preserves_protein_nucleic_acid_ligand_order():
    system = System(
        system={
            "sequences": [
                {"protein": {"id": "A", "sequence": "AC"}},
                {"dna": {"id": ["D", "E"], "sequence": "AT"}},
                {"rna": {"id": "R", "sequence": "GU"}},
                {"ligand": {"id": "L", "smiles": "CCO"}},
            ]
        }
    )
    frame = pd.DataFrame({"conf_chain_id": range(5), "score": [1, 2, 3, 4, 5]})

    result = gather.add_chain_info(frame, system)

    assert result["CHAIN_ID"].tolist() == ["A", "D", "E", "R", "L"]
    assert result["ENTITY_TYPE"].tolist() == ["protein", "dna", "dna", "rna", "ligand"]
    assert result["ligand_molecule_id"].tolist() == [
        "protein_A", "dna_D", "dna_E", "rna_R", "CCO"
    ]
