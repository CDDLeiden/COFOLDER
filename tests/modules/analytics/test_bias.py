"""Tests for bias metric analytics."""

from datetime import date

import pandas as pd
import pytest

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

    prot_val = float(out_chain.loc[out_chain["CHAIN_ID"] == "A", "bias_prot_sim_train"].iloc[0])
    lig_b = float(out_chain.loc[out_chain["CHAIN_ID"] == "B", "bias_lig_sim_train"].iloc[0])
    lig_c = float(out_chain.loc[out_chain["CHAIN_ID"] == "C", "bias_lig_sim_train"].iloc[0])

    assert prot_val == 100.0
    assert lig_b == 1.0
    assert lig_c == 1.0
    assert float(out_system["bias_prot_sim_train_max"].iloc[0]) == 100.0
    assert float(out_system["bias_lig_sim_train_max"].iloc[0]) == 1.0


def test_apply_bias_metrics_resolves_ccd_from_components_cif(temp_dir):
    protein_ref = temp_dir / "protein_training.csv"
    ligand_ref = temp_dir / "ligand_training.csv"
    components_cif = temp_dir / "components.cif"

    protein_ref.write_text(
        "pdb_id,release_date,sequence,sequence_similarity\n"
        "1ABC,2022-01-01,MAAA,99\n",
        encoding="utf-8",
    )
    ligand_ref.write_text(
        "pdb_id,release_date,ligand_id,smiles,ecfp_similarity\n"
        "1ABC,2022-01-01,NOTEDO,CCO,0.6\n",
        encoding="utf-8",
    )
    components_cif.write_text(
        "\n".join(
            [
                "data_EDO",
                "_chem_comp.id EDO",
                "loop_",
                "_pdbx_chem_comp_descriptor.comp_id",
                "_pdbx_chem_comp_descriptor.type",
                "_pdbx_chem_comp_descriptor.program",
                "_pdbx_chem_comp_descriptor.descriptor",
                "EDO SMILES_CANONICAL RDKit CCO",
                "",
            ]
        ),
        encoding="utf-8",
    )

    system_df = pd.DataFrame(
        [{"model_name": "system", "repeat": 1, "diffusion_sample": 0}]
    )
    chain_df = pd.DataFrame(
        [
            {"CHAIN_ID": "A", "ENTITY_TYPE": "protein"},
            {"CHAIN_ID": "B", "ENTITY_TYPE": "ligand", "ligand_molecule_id": "EDO"},
        ]
    )
    sys_obj = _MockSystem(
        sequences=[
            {"protein": {"id": "A", "sequence": "MAAA"}},
            {"ligand": {"id": "B", "ccd": "EDO"}},
        ]
    )

    out_system, out_chain = apply_bias_metrics(
        system_df=system_df,
        chain_df=chain_df,
        sys_obj=sys_obj,
        protein_training_data_path=protein_ref,
        ligand_training_data_path=ligand_ref,
        release_cutoff="2023-06-01",
        components_cif_path=components_cif,
        strict_query_resolution=True,
    )

    assert float(out_chain.loc[out_chain["CHAIN_ID"] == "B", "bias_lig_sim_train"].iloc[0]) == 1.0
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


def test_apply_bias_metrics_rejects_unresolved_ccd_when_strict(monkeypatch, temp_dir):
    home_dir = temp_dir / "home"
    home_dir.mkdir()
    monkeypatch.setenv("HOME", str(home_dir))

    protein_ref = temp_dir / "protein_training.csv"
    ligand_ref = temp_dir / "ligand_training.csv"

    protein_ref.write_text(
        "pdb_id,release_date,sequence\n"
        "1ABC,2022-01-01,MAAA\n",
        encoding="utf-8",
    )
    ligand_ref.write_text(
        "pdb_id,release_date,ligand_id,smiles\n"
        "1ABC,2022-01-01,NOTEDO,CCO\n",
        encoding="utf-8",
    )

    system_df = pd.DataFrame(
        [{"model_name": "system", "repeat": 1, "diffusion_sample": 0}]
    )
    chain_df = pd.DataFrame(
        [{"CHAIN_ID": "B", "ENTITY_TYPE": "ligand", "ligand_molecule_id": "EDO"}]
    )
    sys_obj = _MockSystem(
        sequences=[{"ligand": {"id": "B", "ccd": "EDO"}}]
    )

    with pytest.raises(ValueError, match="Could not resolve SMILES for CCD-backed ligand chain"):
        apply_bias_metrics(
            system_df=system_df,
            chain_df=chain_df,
            sys_obj=sys_obj,
            protein_training_data_path=protein_ref,
            ligand_training_data_path=ligand_ref,
            release_cutoff="2023-06-01",
            strict_query_resolution=True,
        )


def test_apply_bias_metrics_merges_public_and_custom_references_with_summary(temp_dir):
    protein_ref = temp_dir / "protein_training.csv"
    ligand_ref = temp_dir / "ligand_training.csv"
    custom_protein = temp_dir / "custom_protein.csv"
    custom_ligand = temp_dir / "custom_ligand.csv"

    protein_ref.write_text(
        "pdb_id,release_date,sequence,sequence_similarity\n"
        "1PUB,2022-01-01,MKRAAS,83.0\n",
        encoding="utf-8",
    )
    ligand_ref.write_text(
        "pdb_id,release_date,ligand_id,smiles,ecfp_similarity\n"
        "1PUB,2022-01-01,LIG,CCN,0.4\n",
        encoding="utf-8",
    )
    custom_protein.write_text(
        "sequence,dataset_name\n"
        "MKRAAT,private_proteins\n",
        encoding="utf-8",
    )
    custom_ligand.write_text(
        "smiles,dataset_name\n"
        "CCO,private_ligands\n",
        encoding="utf-8",
    )

    system_df = pd.DataFrame([{"model_name": "system", "repeat": 1, "diffusion_sample": 0}])
    chain_df = pd.DataFrame(
        [
            {"CHAIN_ID": "A", "ENTITY_TYPE": "protein"},
            {"CHAIN_ID": "B", "ENTITY_TYPE": "ligand", "ligand_molecule_id": "LIG"},
        ]
    )
    sys_obj = _MockSystem(
        sequences=[
            {"protein": {"id": "A", "sequence": "MKRAAT"}},
            {"ligand": {"id": "B", "smiles": "CCO"}},
        ]
    )

    apply_bias_metrics(
        system_df=system_df,
        chain_df=chain_df,
        sys_obj=sys_obj,
        protein_training_data_path=protein_ref,
        ligand_training_data_path=ligand_ref,
        custom_protein_reference_path=custom_protein,
        custom_ligand_reference_path=custom_ligand,
        release_cutoff="2023-06-01",
        output_dir=temp_dir / "results",
    )

    protein_view = pd.read_csv(temp_dir / "results" / "protein_training_data.csv")
    ligand_view = pd.read_csv(temp_dir / "results" / "ligand_training_data_B.csv")
    summary = pd.read_csv(temp_dir / "results" / "reference_landscape_summary.csv")

    assert set(protein_view["source"]) == {"public", "custom"}
    assert "dataset_name" in protein_view.columns
    assert set(ligand_view["source"]) == {"public", "custom"}
    assert set(summary["nearest_overall_source"]) == {"custom"}
    assert summary["custom_changes_nearest_reference"].tolist() == [True, True]


def test_apply_bias_metrics_rejects_missing_custom_provenance_paths(temp_dir):
    protein_ref = temp_dir / "protein_training.csv"
    ligand_ref = temp_dir / "ligand_training.csv"
    custom_protein = temp_dir / "custom_protein.csv"

    protein_ref.write_text(
        "pdb_id,release_date,sequence\n"
        "1ABC,2022-01-01,MAAA\n",
        encoding="utf-8",
    )
    ligand_ref.write_text(
        "pdb_id,release_date,ligand_id,smiles\n"
        "1ABC,2022-01-01,LIG,CCO\n",
        encoding="utf-8",
    )
    custom_protein.write_text(
        "sequence,source_structure_path\n"
        "MAAA,missing_reference.pdb\n",
        encoding="utf-8",
    )

    system_df = pd.DataFrame([{"model_name": "system", "repeat": 1, "diffusion_sample": 0}])
    chain_df = pd.DataFrame([{"CHAIN_ID": "A", "ENTITY_TYPE": "protein"}])
    sys_obj = _MockSystem(sequences=[{"protein": {"id": "A", "sequence": "MAAA"}}])

    with pytest.raises(ValueError, match="source_structure_path"):
        apply_bias_metrics(
            system_df=system_df,
            chain_df=chain_df,
            sys_obj=sys_obj,
            protein_training_data_path=protein_ref,
            ligand_training_data_path=ligand_ref,
            custom_protein_reference_path=custom_protein,
            release_cutoff="2023-06-01",
        )
