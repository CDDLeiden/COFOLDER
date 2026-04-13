"""Tests for bias metric analytics."""

from datetime import date

import pandas as pd

from cofolder.modules.analytics.bias import _load_ligand_training, apply_bias_metrics


class _MockSystem:
    def __init__(self, sequences):
        self._sequences = sequences

    def find_value(self, key=None, path=None):
        if key == "sequences":
            return self._sequences
        return None


def test_apply_bias_metrics_basic(temp_dir):
    protein_ref = temp_dir / "protein_training.csv"
    ligand_ref = temp_dir / "ligand_training.csv"

    protein_ref.write_text(
        "pdb_id,release_date,sequence\n"
        "1ABC,2022-01-01,MAAA\n"
        "2XYZ,2024-01-01,QQQQ\n"
    )
    ligand_ref.write_text(
        "pdb_id,release_date,ligand_id,smiles\n"
        "1ABC,2022-01-01,LIG,CCO\n"
        "2XYZ,2024-01-01,LIG,CCCC\n"
    )

    system_df = pd.DataFrame(
        [{"model_name": "system", "repeat": 1, "diffusion_sample": 0}]
    )
    chain_df = pd.DataFrame(
        [
            {"CHAIN_ID": "A", "ENTITY_TYPE": "protein", "ligand_molecule_id": "protein_A"},
            {"CHAIN_ID": "E", "ENTITY_TYPE": "ligand", "ligand_molecule_id": "LIG"},
        ]
    )
    sys_obj = _MockSystem(
        sequences=[
            {"protein": {"id": "A", "sequence": "MAAA"}},
            {"ligand": {"id": "E", "smiles": "CCO"}},
        ]
    )

    out_system, out_chain = apply_bias_metrics(
        system_df=system_df,
        chain_df=chain_df,
        sys_obj=sys_obj,
        protein_training_data_path=protein_ref,
        ligand_training_data_path=ligand_ref,
        release_cutoff="2023-06-01",
    )

    prot_val = float(out_chain.loc[out_chain["ENTITY_TYPE"] == "protein", "bias_prot_sim_train"].iloc[0])
    lig_val = float(out_chain.loc[out_chain["ENTITY_TYPE"] == "ligand", "bias_lig_sim_train"].iloc[0])

    assert prot_val == 100.0
    assert lig_val == 1.0
    assert float(out_system["bias_prot_sim_train_max"].iloc[0]) == 100.0
    assert float(out_system["bias_lig_sim_train_max"].iloc[0]) == 1.0


def test_apply_bias_metrics_ccd_multi_ligands_and_cutoff(temp_dir):
    protein_ref = temp_dir / "protein_training.csv"
    ligand_ref = temp_dir / "ligand_training.csv"

    protein_ref.write_text(
        "pdb_id,release_date,sequence,sequence_similarity\n"
        "1ABC,2022-01-01,MAAA,99\n"
        "2XYZ,2024-01-01,QQQQ,100\n"
    )
    # ccd-like ligand_id stored in lower case to verify case-insensitive mapping.
    ligand_ref.write_text(
        "pdb_id,release_date,ligand_id,smiles,ecfp_similarity\n"
        "1ABC,2022-01-01,edo,CCO,0.6\n"
        "2XYZ,2024-01-01,edo,CCCC,0.9\n"
    )

    system_df = pd.DataFrame(
        [{"model_name": "system", "repeat": 1, "diffusion_sample": 0}]
    )
    chain_df = pd.DataFrame(
        [
            {"CHAIN_ID": "A", "ENTITY_TYPE": "protein"},
            {"CHAIN_ID": "B", "ENTITY_TYPE": "ligand"},
            {"CHAIN_ID": "C", "ENTITY_TYPE": "ligand"},
        ]
    )
    sys_obj = _MockSystem(
        sequences=[
            {"protein": {"id": "A", "sequence": "MAAA"}},
            {"ligand": {"id": ["B", "C"], "ccd": "EDO"}},
        ]
    )

    out_system, out_chain = apply_bias_metrics(
        system_df=system_df,
        chain_df=chain_df,
        sys_obj=sys_obj,
        protein_training_data_path=protein_ref,
        ligand_training_data_path=ligand_ref,
        release_cutoff="2023-06-01",
    )

    # Cutoff excludes 2024 rows, so matches must come from 1ABC only.
    prot_val = float(out_chain.loc[out_chain["CHAIN_ID"] == "A", "bias_prot_sim_train"].iloc[0])
    lig_b = float(out_chain.loc[out_chain["CHAIN_ID"] == "B", "bias_lig_sim_train"].iloc[0])
    lig_c = float(out_chain.loc[out_chain["CHAIN_ID"] == "C", "bias_lig_sim_train"].iloc[0])

    assert prot_val == 100.0
    assert lig_b == 1.0
    assert lig_c == 1.0
    assert float(out_system["bias_prot_sim_train_max"].iloc[0]) == 100.0
    assert float(out_system["bias_lig_sim_train_max"].iloc[0]) == 1.0


def test_apply_bias_metrics_prefers_precomputed_similarity_columns(temp_dir):
    protein_ref = temp_dir / "protein_training.csv"
    ligand_ref = temp_dir / "ligand_training.csv"

    protein_ref.write_text(
        "pdb_id,release_date,sequence,sequence_similarity\n"
        "1AAA,2022-01-01,MMMM,31.0\n"
        "1BBB,2022-01-01,QQQQ,42.5\n"
    )
    ligand_ref.write_text(
        "pdb_id,release_date,ligand_id,smiles,ecfp_similarity\n"
        "1AAA,2022-01-01,LIG,CCO,0.36\n"
        "1BBB,2022-01-01,LIG,CCC,0.77\n"
    )

    system_df = pd.DataFrame(
        [{"model_name": "system", "repeat": 1, "diffusion_sample": 0}]
    )
    chain_df = pd.DataFrame(
        [
            {"CHAIN_ID": "A", "ENTITY_TYPE": "protein"},
            {"CHAIN_ID": "B", "ENTITY_TYPE": "ligand", "ligand_molecule_id": "LIG"},
        ]
    )
    sys_obj = _MockSystem(
        sequences=[
            {"protein": {"id": "A", "sequence": "UNRELATEDSEQ"}},
            {"ligand": {"id": "B", "smiles": "N#N"}},
        ]
    )

    out_system, out_chain = apply_bias_metrics(
        system_df=system_df,
        chain_df=chain_df,
        sys_obj=sys_obj,
        protein_training_data_path=protein_ref,
        ligand_training_data_path=ligand_ref,
        release_cutoff="2023-06-01",
    )

    prot_val = float(out_chain.loc[out_chain["CHAIN_ID"] == "A", "bias_prot_sim_train"].iloc[0])
    lig_val = float(out_chain.loc[out_chain["CHAIN_ID"] == "B", "bias_lig_sim_train"].iloc[0])
    assert prot_val == 42.5
    assert lig_val == 0.77
    assert float(out_system["bias_prot_sim_train_max"].iloc[0]) == 42.5
    assert float(out_system["bias_lig_sim_train_max"].iloc[0]) == 0.77


def test_load_ligand_training_tolerates_empty_csv(temp_dir):
    ligand_ref = temp_dir / "ligand_training.csv"
    ligand_ref.write_text("")

    out = _load_ligand_training(ligand_ref, date(2023, 6, 1))

    assert out.empty
    assert {"pdb_id", "release_date", "smiles", "ligand_id"}.issubset(set(out.columns))


def test_apply_bias_metrics_with_empty_ligand_csv_no_crash(temp_dir):
    protein_ref = temp_dir / "protein_training.csv"
    ligand_ref = temp_dir / "ligand_training.csv"

    protein_ref.write_text(
        "pdb_id,release_date,sequence\n"
        "1ABC,2022-01-01,MAAA\n"
    )
    ligand_ref.write_text("")

    system_df = pd.DataFrame(
        [{"model_name": "system", "repeat": 1, "diffusion_sample": 0}]
    )
    chain_df = pd.DataFrame(
        [
            {"CHAIN_ID": "A", "ENTITY_TYPE": "protein", "ligand_molecule_id": "protein_A"},
            {"CHAIN_ID": "B", "ENTITY_TYPE": "ligand", "ligand_molecule_id": "LIG"},
        ]
    )
    sys_obj = _MockSystem(
        sequences=[
            {"protein": {"id": "A", "sequence": "MAAA"}},
            {"ligand": {"id": "B", "smiles": "CCO"}},
        ]
    )

    out_system, out_chain = apply_bias_metrics(
        system_df=system_df,
        chain_df=chain_df,
        sys_obj=sys_obj,
        protein_training_data_path=protein_ref,
        ligand_training_data_path=ligand_ref,
        release_cutoff="2023-06-01",
    )

    prot_val = float(out_chain.loc[out_chain["CHAIN_ID"] == "A", "bias_prot_sim_train"].iloc[0])
    lig_val = out_chain.loc[out_chain["CHAIN_ID"] == "B", "bias_lig_sim_train"].iloc[0]

    assert prot_val == 100.0
    assert pd.isna(lig_val)
    assert float(out_system["bias_prot_sim_train_max"].iloc[0]) == 100.0
