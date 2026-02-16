"""Tests for bias metric analytics."""

import pandas as pd

from cofolder.modules.analytics.bias import apply_bias_metrics


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

    out_system, out_chain, bias_df = apply_bias_metrics(
        system_df=system_df,
        chain_df=chain_df,
        sys_obj=sys_obj,
        protein_training_data_path=protein_ref,
        ligand_training_data_path=ligand_ref,
        release_cutoff="2023-06-01",
    )

    prot_val = float(out_chain.loc[out_chain["ENTITY_TYPE"] == "protein", "protein_sequence_identity_train"].iloc[0])
    lig_val = float(out_chain.loc[out_chain["ENTITY_TYPE"] == "ligand", "ligand_fingerprint_similarity_train"].iloc[0])

    assert prot_val == 100.0
    assert lig_val == 1.0
    assert float(out_system["protein_sequence_identity_train_max"].iloc[0]) == 100.0
    assert float(out_system["ligand_fingerprint_similarity_train_max"].iloc[0]) == 1.0

    assert len(bias_df) == 1
    assert float(bias_df["x_ligand_ecfp_similarity"].iloc[0]) == 1.0
    assert float(bias_df["y_protein_sequence_identity"].iloc[0]) == 100.0
