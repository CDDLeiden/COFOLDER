"""Tests for bias metric analytics."""

from datetime import date
import logging
from pathlib import Path
import subprocess
import sys
import tempfile

import pandas as pd
import pytest

from cofolder.modules.analytics import plots as plots_module
from cofolder.modules.analytics.bias import _references as bias_references
from cofolder.modules.analytics.bias import _similarity as bias_similarity
from cofolder.modules.analytics.bias import _enrichment as bias_enrichment
from cofolder.modules.analytics.bias import _artifacts as bias_artifacts
from cofolder.modules.analytics.bias import (
    BIAS_TRAINING_DATA_COLUMNS,
    apply_bias_metrics,
)
from cofolder.modules.analytics.bias._artifacts import _same_type_pair_dataset
from cofolder.modules.analytics.bias._datasets import _build_bias_training_dataset
from cofolder.modules.analytics.bias._references import _load_ligand_training
from cofolder.modules.analytics.bias._similarity import _build_protein_training_view
from cofolder.modules.analytics.plots import plot_bias_reference_overlap, plot_reference_overlap_scatter


class _MockSystem:
    def __init__(self, sequences):
        self._sequences = sequences

    def find_value(self, key=None, path=None):
        if key == "sequences":
            return self._sequences
        return None


def test_bias_package_imports_have_no_execution_side_effects(temp_dir):
    modules = (
        "cofolder.modules.analytics.bias",
        "cofolder.modules.analytics.bias._common",
        "cofolder.modules.analytics.bias._references",
        "cofolder.modules.analytics.bias._similarity",
        "cofolder.modules.analytics.bias._datasets",
        "cofolder.modules.analytics.bias._enrichment",
        "cofolder.modules.analytics.bias._artifacts",
        "cofolder.modules.analytics.bias._workflow",
    )
    code = f"""
import sys

def deny(event, args):
    if event in {{'subprocess.Popen', 'socket.connect'}}:
        raise RuntimeError(event)

sys.addaudithook(deny)
for module in {modules!r}:
    __import__(module)
"""
    before = tuple(temp_dir.iterdir())
    completed = subprocess.run(
        [sys.executable, "-c", code],
        cwd=temp_dir,
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    assert tuple(temp_dir.iterdir()) == before


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

    protein_row = out_chain.loc[out_chain["ENTITY_TYPE"] == "protein"].iloc[0]
    lig_val = float(out_chain.loc[out_chain["ENTITY_TYPE"] == "ligand", "bias_lig_sim_train"].iloc[0])

    assert pd.isna(protein_row["bias_prot_sim_train"])
    assert float(protein_row["bias_prot_sim_train_pairwise"]) == 100.0
    assert lig_val == 1.0
    assert float(out_system["bias_prot_sim_train_pairwise_max"].iloc[0]) == 100.0
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

    assert prot_val == 99.0
    assert lig_b == 1.0
    assert lig_c == 1.0
    assert float(out_system["bias_prot_sim_train_max"].iloc[0]) == 99.0
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


def test_smiles_from_components_cif_reuses_cached_document(monkeypatch, temp_dir):
    components_cif = temp_dir / "components.cif"
    components_cif.write_text(
        "\n".join(
            [
                "data_MG",
                "_chem_comp.id MG",
                "loop_",
                "_pdbx_chem_comp_descriptor.comp_id",
                "_pdbx_chem_comp_descriptor.type",
                "_pdbx_chem_comp_descriptor.program",
                "_pdbx_chem_comp_descriptor.descriptor",
                "MG SMILES_CANONICAL RDKit [Mg+2]",
                "",
                "data_K",
                "_chem_comp.id K",
                "loop_",
                "_pdbx_chem_comp_descriptor.comp_id",
                "_pdbx_chem_comp_descriptor.type",
                "_pdbx_chem_comp_descriptor.program",
                "_pdbx_chem_comp_descriptor.descriptor",
                "K SMILES_CANONICAL RDKit [K+]",
                "",
            ]
        ),
        encoding="utf-8",
    )

    read_calls = 0
    original_read_file = bias_references.gemmi.cif.read_file

    def _counting_read_file(path):
        nonlocal read_calls
        read_calls += 1
        return original_read_file(path)

    bias_references._components_cif_smiles_index.cache_clear()
    monkeypatch.setattr(bias_references.gemmi.cif, "read_file", _counting_read_file)

    assert bias_references._smiles_from_components_cif("MG", components_cif) == "[Mg+2]"
    assert bias_references._smiles_from_components_cif("K", components_cif) == "[K+]"
    assert read_calls == 1


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

    protein_row = out_chain.loc[out_chain["CHAIN_ID"] == "A"].iloc[0]
    lig_val = out_chain.loc[out_chain["CHAIN_ID"] == "B", "bias_lig_sim_train"].iloc[0]

    assert pd.isna(protein_row["bias_prot_sim_train"])
    assert float(protein_row["bias_prot_sim_train_pairwise"]) == 100.0
    assert pd.isna(lig_val)
    assert float(out_system["bias_prot_sim_train_pairwise_max"].iloc[0]) == 100.0


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
    structure_ref = temp_dir / "custom_reference.pdb"
    structure_ref.write_text("HEADER CUSTOM\n", encoding="utf-8")

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
        "sequence,dataset_name,source_structure_path\n"
        f"MKRAAT,private_proteins,{structure_ref.name}\n",
        encoding="utf-8",
    )
    custom_ligand.write_text(
        "smiles,dataset_name,source_reference_path\n"
        f"CCO,private_ligands,{structure_ref.name}\n",
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
    bias_training = pd.read_csv(temp_dir / "results" / "bias_training_data.csv")
    summary = pd.read_csv(temp_dir / "results" / "reference_landscape_summary.csv")

    assert set(protein_view["source"]) == {"public", "custom"}
    assert "dataset_name" in protein_view.columns
    assert set(ligand_view["source"]) == {"public", "custom"}
    assert bias_training.columns.tolist() == BIAS_TRAINING_DATA_COLUMNS
    assert bias_training["source"].tolist() == ["custom", "public"]
    assert bias_training["query_protein_chain_id"].tolist() == ["A", "A"]
    assert bias_training["query_ligand_chain_id"].tolist() == ["B", "B"]
    assert bias_training["pairing_status"].tolist() == ["paired", "paired"]
    assert pd.isna(bias_training["plot_sequence_similarity"].iloc[0])
    assert bias_training["plot_sequence_similarity"].iloc[1] == 0.83
    assert bias_training["plot_sequence_similarity_pairwise"].iloc[0] == 1.0
    assert pd.isna(bias_training["plot_sequence_similarity_pairwise"].iloc[1])
    assert bias_training["sequence_similarity_method"].tolist() == [
        "pairwise_aligner",
        "mmseqs_pident",
    ]
    assert bias_training["plot_ecfp_similarity"].tolist() == [1.0, 0.4]
    assert (temp_dir / "results" / "bias_reference_overlap_scatter.png").exists()
    assert (temp_dir / "results" / "bias_reference_overlap_scatter.pdf").exists()
    assert set(summary["nearest_overall_source"]) == {"custom"}
    summary_methods = summary.set_index("entity_type")[
        "nearest_overall_similarity_method"
    ].to_dict()
    assert summary_methods == {
        "protein": "pairwise_aligner",
        "ligand": "ecfp4_tanimoto",
    }
    assert summary["custom_changes_nearest_reference"].tolist() == [True, True]


def test_apply_bias_metrics_writes_pair_specific_outputs_for_unique_ligand_groups(temp_dir):
    protein_ref = temp_dir / "protein_training.csv"
    ligand_ref = temp_dir / "ligand_training.csv"

    protein_ref.write_text(
        "pdb_id,release_date,sequence,sequence_similarity\n"
        "1ABC,2022-01-01,MAAA,100.0\n",
        encoding="utf-8",
    )
    ligand_ref.write_text(
        "pdb_id,release_date,ligand_id,smiles,ecfp_similarity\n"
        "1ABC,2022-01-01,LIG,CCO,1.0\n"
        "1ABC,2022-01-01,MG,[Mg+2],1.0\n",
        encoding="utf-8",
    )

    system_df = pd.DataFrame([{"model_name": "system", "repeat": 1, "diffusion_sample": 0}])
    chain_df = pd.DataFrame(
        [
            {"CHAIN_ID": "A", "ENTITY_TYPE": "protein", "ligand_molecule_id": "protein_A"},
            {"CHAIN_ID": "B", "ENTITY_TYPE": "ligand", "ligand_molecule_id": "CCO"},
            {"CHAIN_ID": "C", "ENTITY_TYPE": "ligand", "ligand_molecule_id": "MG"},
            {"CHAIN_ID": "D", "ENTITY_TYPE": "ligand", "ligand_molecule_id": "MG"},
        ]
    )
    sys_obj = _MockSystem(
        sequences=[
            {"protein": {"id": "A", "sequence": "MAAA"}},
            {"ligand": {"id": "B", "smiles": "CCO"}},
            {"ligand": {"id": ["C", "D"], "ccd": "MG"}},
        ]
    )

    apply_bias_metrics(
        system_df=system_df,
        chain_df=chain_df,
        sys_obj=sys_obj,
        protein_training_data_path=protein_ref,
        ligand_training_data_path=ligand_ref,
        release_cutoff="2023-06-01",
        output_dir=temp_dir / "results",
    )

    output_dir = temp_dir / "results"
    assert (output_dir / "bias_training_data_A__B.csv").exists()
    assert (output_dir / "bias_training_data_A__MG.csv").exists()
    assert (output_dir / "bias_reference_overlap_scatter_A__B.png").exists()
    assert (output_dir / "bias_reference_overlap_scatter_A__B.pdf").exists()
    assert (output_dir / "bias_reference_overlap_scatter_A__MG.png").exists()
    assert (output_dir / "bias_reference_overlap_scatter_A__MG.pdf").exists()
    assert not (output_dir / "bias_training_data_A__C.csv").exists()
    assert not (output_dir / "bias_training_data_A__D.csv").exists()

    mg_pair = pd.read_csv(output_dir / "bias_training_data_A__MG.csv")
    assert mg_pair["query_pair_id"].tolist() == ["A__MG"]
    assert mg_pair["query_ligand_chain_id"].tolist() == ["MG"]

    bias_training = pd.read_csv(output_dir / "bias_training_data.csv")
    assert set(bias_training["query_pair_id"]) == {"A__B", "A__MG"}
    assert set(bias_training["query_ligand_chain_id"]) == {"B", "MG"}
    assert not set(bias_training["query_ligand_chain_id"]).intersection({"C", "D"})
    assert (output_dir / "ligand_training_data.csv").exists()


def test_apply_bias_metrics_collapses_duplicate_protein_queries_in_pair_specific_outputs(temp_dir):
    protein_ref = temp_dir / "protein_training.csv"
    ligand_ref = temp_dir / "ligand_training.csv"

    protein_ref.write_text(
        "pdb_id,release_date,sequence,sequence_similarity\n"
        "1ABC,2022-01-01,MAAA,100.0\n",
        encoding="utf-8",
    )
    ligand_ref.write_text(
        "pdb_id,release_date,ligand_id,smiles,ecfp_similarity\n"
        "1ABC,2022-01-01,LIG,CCO,1.0\n",
        encoding="utf-8",
    )

    system_df = pd.DataFrame([{"model_name": "system", "repeat": 1, "diffusion_sample": 0}])
    chain_df = pd.DataFrame(
        [
            {"CHAIN_ID": "A", "ENTITY_TYPE": "protein", "ligand_molecule_id": "protein_A"},
            {"CHAIN_ID": "C", "ENTITY_TYPE": "protein", "ligand_molecule_id": "protein_C"},
            {"CHAIN_ID": "B", "ENTITY_TYPE": "ligand", "ligand_molecule_id": "CCO"},
        ]
    )
    sys_obj = _MockSystem(
        sequences=[
            {"protein": {"id": ["A", "C"], "sequence": "MAAA"}},
            {"ligand": {"id": "B", "smiles": "CCO"}},
        ]
    )

    apply_bias_metrics(
        system_df=system_df,
        chain_df=chain_df,
        sys_obj=sys_obj,
        protein_training_data_path=protein_ref,
        ligand_training_data_path=ligand_ref,
        release_cutoff="2023-06-01",
        output_dir=temp_dir / "results",
    )

    output_dir = temp_dir / "results"
    assert (output_dir / "bias_training_data_A__B.csv").exists()
    assert (output_dir / "bias_reference_overlap_scatter_A__B.png").exists()
    assert not (output_dir / "bias_training_data_C__B.csv").exists()
    assert not (output_dir / "bias_reference_overlap_scatter_C__B.png").exists()

    collapsed_pair = pd.read_csv(output_dir / "bias_training_data_A__B.csv")
    assert set(collapsed_pair["query_protein_chain_id"]) == {"A"}


def test_build_bias_training_dataset_clears_ambiguous_public_detail_fields():
    protein_view = pd.DataFrame(
        [
            {
                "query_chain_id": "A",
                "pdb_id": "1ABC",
                "release_date": "2022-01-01",
                "source": "public",
                "dataset_name": "public",
                "sequence_similarity": 95.0,
                "sequence": "SEQ1",
            },
            {
                "query_chain_id": "A",
                "pdb_id": "1ABC",
                "release_date": "2022-01-01",
                "source": "public",
                "dataset_name": "public",
                "sequence_similarity": 90.0,
                "sequence": "SEQ2",
            },
        ]
    )
    ligand_view = pd.DataFrame(
        [
            {
                "query_chain_id": "B",
                "pdb_id": "1ABC",
                "release_date": "2022-01-01",
                "source": "public",
                "dataset_name": "public",
                "ligand_id": "L1",
                "ecfp_similarity": 0.8,
                "smiles": "CCO",
            },
            {
                "query_chain_id": "B",
                "pdb_id": "1ABC",
                "release_date": "2022-01-01",
                "source": "public",
                "dataset_name": "public",
                "ligand_id": "L2",
                "ecfp_similarity": 0.7,
                "smiles": "CCN",
            },
        ]
    )

    result = _build_bias_training_dataset(
        protein_view=protein_view,
        ligand_views={"B": ligand_view},
    )

    assert len(result) == 4
    assert result["reference_key"].tolist() == ["public:pdb:1ABC"] * 4
    assert result["sequence"].isna().all()
    assert result["ligand_id"].isna().all()
    assert result["smiles"].isna().all()
    assert result["pairing_status"].tolist() == ["paired"] * 4


def test_build_bias_training_dataset_leaves_disjoint_retained_rows_for_pdb_backfill():
    protein_view = pd.DataFrame(
        [
            {
                "query_chain_id": "A",
                "pdb_id": "5ZOD",
                "release_date": "2022-01-01",
                "source": "public",
                "dataset_name": "public",
                "sequence_similarity": 94.1,
                "sequence": "SEQ_A",
            },
            {
                "query_chain_id": "A",
                "pdb_id": "3Q8K",
                "release_date": "2022-01-01",
                "source": "public",
                "dataset_name": "public",
                "sequence_similarity": 94.0,
                "sequence": "SEQ_B",
            },
        ]
    )
    ligand_view = pd.DataFrame(
        [
            {
                "query_chain_id": "B",
                "pdb_id": "1XM1",
                "release_date": "2022-01-01",
                "source": "public",
                "dataset_name": "public",
                "ligand_id": "L1",
                "ecfp_similarity": 0.38,
                "smiles": "CCO",
            },
            {
                "query_chain_id": "B",
                "pdb_id": "6AGT",
                "release_date": "2022-01-01",
                "source": "public",
                "dataset_name": "public",
                "ligand_id": "L2",
                "ecfp_similarity": 0.36,
                "smiles": "CCN",
            },
        ]
    )
    protein_lookup_view = pd.DataFrame(
        [
            *protein_view.to_dict("records"),
            {
                "query_chain_id": "A",
                "pdb_id": "1XM1",
                "release_date": "2022-01-01",
                "source": "public",
                "dataset_name": "public",
                "sequence_similarity": 61.0,
                "sequence": "SEQ_C",
            },
            {
                "query_chain_id": "A",
                "pdb_id": "6AGT",
                "release_date": "2022-01-01",
                "source": "public",
                "dataset_name": "public",
                "sequence_similarity": 58.0,
                "sequence": "SEQ_D",
            },
        ]
    )
    ligand_lookup_view = pd.DataFrame(
        [
            *ligand_view.to_dict("records"),
            {
                "query_chain_id": "B",
                "pdb_id": "5ZOD",
                "release_date": "2022-01-01",
                "source": "public",
                "dataset_name": "public",
                "ligand_id": "L5",
                "ecfp_similarity": 0.12,
                "smiles": "CCC",
            },
            {
                "query_chain_id": "B",
                "pdb_id": "3Q8K",
                "release_date": "2022-01-01",
                "source": "public",
                "dataset_name": "public",
                "ligand_id": "L6",
                "ecfp_similarity": 0.09,
                "smiles": "CCCl",
            },
        ]
    )

    result = _build_bias_training_dataset(
        protein_view=protein_view,
        ligand_views={"B": ligand_view},
        protein_lookup_view=protein_lookup_view,
        ligand_lookup_views={"B": ligand_lookup_view},
    )

    assert len(result) == 4
    assert set(result["pairing_status"]) == {"protein_only", "ligand_only"}
    assert set(result["protein_pdb_id"].dropna()) == {"5ZOD", "3Q8K"}
    assert set(result["ligand_pdb_id"].dropna()) == {"1XM1", "6AGT"}
    assert set(result["reference_key"]) == {
        "public:pdb:5ZOD",
        "public:pdb:3Q8K",
        "public:pdb:1XM1",
        "public:pdb:6AGT",
    }


def test_mixed_pair_dataset_handles_direct_overlap_per_reference(monkeypatch):
    protein_view = pd.DataFrame(
        [
            {
                "query_chain_id": "A",
                "pdb_id": "1UL1",
                "release_date": "2022-01-01",
                "source": "public",
                "dataset_name": "public",
                "sequence_similarity": 99.0,
                "sequence": "SEQ_OVERLAP",
            },
            {
                "query_chain_id": "A",
                "pdb_id": "5ZOD",
                "release_date": "2022-01-01",
                "source": "public",
                "dataset_name": "public",
                "sequence_similarity": 94.1,
                "sequence": "SEQ_A",
            },
        ]
    )
    ligand_view = pd.DataFrame(
        [
            {
                "query_chain_id": "MG",
                "pdb_id": "1UL1",
                "release_date": "2022-01-01",
                "source": "public",
                "dataset_name": "public",
                "ligand_id": "MG",
                "ecfp_similarity": 1.0,
                "smiles": "[Mg+2]",
            },
            {
                "query_chain_id": "MG",
                "pdb_id": "1XM1",
                "release_date": "2022-01-01",
                "source": "public",
                "dataset_name": "public",
                "ligand_id": "MG",
                "ecfp_similarity": 0.91,
                "smiles": "[Mg+2]",
            },
        ]
    )

    monkeypatch.setattr(
        "cofolder.modules.analytics.bias._enrichment._pdb_ligand_similarity_rows",
        lambda **kwargs: [{"ligand_id": "MG", "smiles": "[Mg+2]", "ecfp_similarity": 0.05}]
        if kwargs["pdb_id"] == "5ZOD"
        else [],
    )
    monkeypatch.setattr(
        "cofolder.modules.analytics.bias._enrichment._pdb_protein_similarity_rows",
        lambda **kwargs: [{"sequence": "SEQ_X", "sequence_similarity": 12.5}]
        if kwargs["pdb_id"] == "1XM1"
        else [],
    )
    monkeypatch.setattr(
        "cofolder.modules.analytics.bias._enrichment._resolved_mmseqs_bin_token",
        lambda: "mmseqs",
    )

    result = bias_artifacts._mixed_pair_dataset(
        protein_label="A",
        ligand_label="MG",
        protein_query_sequence="SEQ_QUERY",
        ligand_query_smiles="[Mg+2]",
        protein_view=protein_view,
        ligand_view=ligand_view,
        protein_lookup_view=protein_view,
        ligand_lookup_view=ligand_view,
        boltz_cache_path=Path(tempfile.gettempdir()),
        components_cif_path=None,
    )

    assert len(result) == 3
    assert set(result["pairing_status"]) == {"paired"}
    assert set(result["pdb_id"]) == {"1UL1", "5ZOD", "1XM1"}
    overlap = result[result["pdb_id"] == "1UL1"].iloc[0]
    assert overlap["sequence_similarity"] == 99.0
    assert overlap["ecfp_similarity"] == 1.0
    protein_seed = result[result["pdb_id"] == "5ZOD"].iloc[0]
    assert protein_seed["sequence_similarity"] == 94.1
    assert protein_seed["ecfp_similarity"] == 0.05
    ligand_seed = result[result["pdb_id"] == "1XM1"].iloc[0]
    assert ligand_seed["sequence_similarity"] == 12.5
    assert ligand_seed["ecfp_similarity"] == 0.91


def test_same_type_pair_dataset_uses_other_ligands_from_same_pdb(monkeypatch):
    component_1_view = pd.DataFrame(
        [
            {
                "query_chain_id": "B",
                "pdb_id": "1XM1",
                "release_date": "2022-01-01",
                "source": "public",
                "dataset_name": "public",
                "ligand_id": "L1",
                "ecfp_similarity": 1.0,
                "smiles": "CCO",
            },
            {
                "query_chain_id": "B",
                "pdb_id": "6AGT",
                "release_date": "2022-01-01",
                "source": "public",
                "dataset_name": "public",
                "ligand_id": "L2",
                "ecfp_similarity": 0.9,
                "smiles": "CCO",
            },
        ]
    )
    component_2_view = pd.DataFrame(
        [
            {
                "query_chain_id": "MG",
                "pdb_id": "5ZOD",
                "release_date": "2022-01-01",
                "source": "public",
                "dataset_name": "public",
                "ligand_id": "MG",
                "ecfp_similarity": 1.0,
                "smiles": "[Mg+2]",
            },
            {
                "query_chain_id": "MG",
                "pdb_id": "3Q8K",
                "release_date": "2022-01-01",
                "source": "public",
                "dataset_name": "public",
                "ligand_id": "MG",
                "ecfp_similarity": 0.95,
                "smiles": "[Mg+2]",
            },
        ]
    )
    calls: list[tuple[str, tuple[str, ...]]] = []

    def _fake_pdb_ligand_rows(**kwargs):
        calls.append(
            (
                kwargs["pdb_id"],
                tuple(sorted(kwargs.get("exclude_ligand_ids") or set())),
            )
        )
        pdb_id = kwargs["pdb_id"]
        if pdb_id == "1XM1":
            return [
                {"ligand_id": "MG", "smiles": "[Mg+2]", "ecfp_similarity": 1.0},
                {"ligand_id": "SO4", "smiles": "[O-]S([O-])(=O)=O", "ecfp_similarity": 0.0},
            ]
        if pdb_id == "6AGT":
            return [{"ligand_id": "MG", "smiles": "[Mg+2]", "ecfp_similarity": 0.8}]
        if pdb_id == "5ZOD":
            return [{"ligand_id": "L5", "smiles": "CCO", "ecfp_similarity": 0.2}]
        if pdb_id == "3Q8K":
            return [{"ligand_id": "L6", "smiles": "CCO", "ecfp_similarity": 0.15}]
        return []

    monkeypatch.setattr(
        "cofolder.modules.analytics.bias._enrichment._pdb_ligand_similarity_rows",
        _fake_pdb_ligand_rows,
    )

    result = _same_type_pair_dataset(
        component_1_label="B",
        component_2_label="MG",
        component_type="ligand",
        component_1_view=component_1_view,
        component_2_view=component_2_view,
        component_1_lookup_view=component_1_view,
        component_2_lookup_view=component_2_view,
        ligand_queries={"B": "CCO", "MG": "[Mg+2]"},
        boltz_cache_path=Path(tempfile.gettempdir()),
        components_cif_path=None,
    )

    assert len(result) == 5
    assert set(result["pairing_status"]) == {"paired"}
    assert set(result["query_pair_id"]) == {"B__MG"}
    assert set(result["component_1_type"]) == {"ligand"}
    assert set(result["component_2_type"]) == {"ligand"}
    assert set(result["reference_pdb_id"]) == {"1XM1", "6AGT", "5ZOD", "3Q8K"}
    assert sorted(zip(result["reference_pdb_id"], result["component_1_similarity"], result["component_2_similarity"])) == [
        ("1XM1", 1.0, 0.0),
        ("1XM1", 1.0, 1.0),
        ("3Q8K", 0.15, 0.95),
        ("5ZOD", 0.2, 1.0),
        ("6AGT", 0.9, 0.8),
    ]
    assert sorted(calls) == [
        ("1XM1", ("L1",)),
        ("3Q8K", ("MG",)),
        ("5ZOD", ("MG",)),
        ("6AGT", ("L2",)),
    ]


def test_enrich_mixed_bias_dataset_with_pdb_backfill_fetches_missing_axes(monkeypatch):
    result = pd.DataFrame(
        [
            {
                "query_pair_id": "A__B",
                "query_protein_chain_id": "A",
                "query_ligand_chain_id": "B",
                "reference_key": "public:pdb:5ZOD",
                "reference_label": "5ZOD",
                "pairing_status": "protein_only",
                "source": "public",
                "dataset_name": "public",
                "pdb_id": "5ZOD",
                "protein_pdb_id": "5ZOD",
                "ligand_pdb_id": pd.NA,
                "protein_release_date": "2022-01-01",
                "ligand_release_date": pd.NA,
                "protein_source": "public",
                "protein_dataset_name": "public",
                "protein_source_structure_path": pd.NA,
                "protein_source_reference_path": pd.NA,
                "ligand_source": pd.NA,
                "ligand_dataset_name": pd.NA,
                "ligand_source_structure_path": pd.NA,
                "ligand_source_reference_path": pd.NA,
                "sequence_similarity": 94.1,
                "ecfp_similarity": 0.0,
                "plot_sequence_similarity": 0.941,
                "plot_ecfp_similarity": 0.0,
                "sequence": "SEQ_A",
                "ligand_id": pd.NA,
                "smiles": pd.NA,
            },
            {
                "query_pair_id": "A__B",
                "query_protein_chain_id": "A",
                "query_ligand_chain_id": "B",
                "reference_key": "public:pdb:1XM1",
                "reference_label": "1XM1",
                "pairing_status": "ligand_only",
                "source": "public",
                "dataset_name": "public",
                "pdb_id": "1XM1",
                "protein_pdb_id": pd.NA,
                "ligand_pdb_id": "1XM1",
                "protein_release_date": pd.NA,
                "ligand_release_date": "2022-01-01",
                "protein_source": pd.NA,
                "protein_dataset_name": pd.NA,
                "protein_source_structure_path": pd.NA,
                "protein_source_reference_path": pd.NA,
                "ligand_source": "public",
                "ligand_dataset_name": "public",
                "ligand_source_structure_path": pd.NA,
                "ligand_source_reference_path": pd.NA,
                "sequence_similarity": 0.0,
                "ecfp_similarity": 0.38,
                "plot_sequence_similarity": 0.0,
                "plot_ecfp_similarity": 0.38,
                "sequence": pd.NA,
                "ligand_id": "L1",
                "smiles": "CCO",
            },
        ]
    )

    monkeypatch.setattr(
        "cofolder.modules.analytics.bias._enrichment._pdb_ligand_similarity_rows",
        lambda **kwargs: [{"ligand_id": "MG", "smiles": "[Mg++]", "ecfp_similarity": 0.11}]
        if kwargs["pdb_id"] == "5ZOD"
        else [],
    )
    monkeypatch.setattr(
        "cofolder.modules.analytics.bias._enrichment._pdb_protein_similarity_rows",
        lambda **kwargs: [{"sequence": "SEQ_X", "sequence_similarity": 12.5}]
        if kwargs["pdb_id"] == "1XM1"
        else [],
    )
    monkeypatch.setattr(
        "cofolder.modules.analytics.bias._enrichment._resolved_mmseqs_bin_token",
        lambda: "mmseqs",
    )

    enriched = bias_enrichment._enrich_mixed_bias_dataset_with_pdb_backfill(
        result,
        protein_queries={"A": "SEQ_QUERY"},
        ligand_queries={"B": "CCO"},
        boltz_cache_path=Path(tempfile.gettempdir()),
        components_cif_path=None,
    )

    assert enriched["pairing_status"].tolist() == ["paired", "paired"]
    zod = enriched[enriched["pdb_id"] == "5ZOD"].iloc[0]
    assert zod["ligand_pdb_id"] == "5ZOD"
    assert zod["ecfp_similarity"] == 0.11
    xm1 = enriched[enriched["pdb_id"] == "1XM1"].iloc[0]
    assert xm1["protein_pdb_id"] == "1XM1"
    assert xm1["sequence_similarity"] == 12.5
    assert xm1["plot_sequence_similarity"] == 0.125
    assert pd.isna(xm1["sequence_similarity_pairwise"])
    assert xm1["sequence_similarity_method"] == "mmseqs_pident"


def test_enrich_mixed_bias_dataset_with_pdb_backfill_reuses_protein_lookup_rows(monkeypatch):
    rows = pd.DataFrame(
        [
            {
                "query_pair_id": "A__MG",
                "query_protein_chain_id": "A",
                "query_ligand_chain_id": "MG",
                "reference_key": "public:pdb:1XM1",
                "reference_label": "1XM1",
                "pairing_status": "ligand_only",
                "source": "public",
                "dataset_name": "public",
                "pdb_id": "1XM1",
                "protein_pdb_id": pd.NA,
                "ligand_pdb_id": "1XM1",
                "protein_release_date": pd.NA,
                "ligand_release_date": "2022-01-01",
                "protein_source": pd.NA,
                "protein_dataset_name": pd.NA,
                "protein_source_structure_path": pd.NA,
                "protein_source_reference_path": pd.NA,
                "ligand_source": "public",
                "ligand_dataset_name": "public",
                "ligand_source_structure_path": pd.NA,
                "ligand_source_reference_path": pd.NA,
                "sequence_similarity": 0.0,
                "ecfp_similarity": 0.91,
                "plot_sequence_similarity": 0.0,
                "plot_ecfp_similarity": 0.91,
                "sequence": pd.NA,
                "ligand_id": "MG",
                "smiles": "[Mg+2]",
            },
        ]
    )
    protein_lookup_view = pd.DataFrame(
        [
            {
                "query_chain_id": "A",
                "pdb_id": "1XM1",
                "release_date": "2021-12-31",
                "source": "custom",
                "dataset_name": "protein_lookup",
                "source_structure_path": "/tmp/1xm1.cif",
                "source_reference_path": "/tmp/protein.csv",
                "sequence_similarity": 12.5,
                "sequence": "SEQ_LOOKUP",
            },
        ]
    )

    monkeypatch.setattr(
        "cofolder.modules.analytics.bias._enrichment._pdb_protein_similarity_rows",
        lambda **kwargs: (_ for _ in ()).throw(AssertionError("MMseqs fallback should not run")),
    )

    enriched = bias_enrichment._enrich_mixed_bias_dataset_with_pdb_backfill(
        rows,
        protein_queries={"A": "SEQ_QUERY"},
        ligand_queries={"MG": "[Mg+2]"},
        protein_lookup_view=protein_lookup_view,
        boltz_cache_path=Path(tempfile.gettempdir()),
        components_cif_path=None,
    )

    assert enriched["pairing_status"].tolist() == ["paired"]
    row = enriched.iloc[0]
    assert row["protein_pdb_id"] == "1XM1"
    assert row["sequence_similarity"] == 12.5
    assert row["plot_sequence_similarity"] == 0.125
    assert row["sequence"] == "SEQ_LOOKUP"
    assert row["protein_source"] == "custom"
    assert row["protein_dataset_name"] == "protein_lookup"
    assert row["protein_source_structure_path"] == "/tmp/1xm1.cif"
    assert row["protein_source_reference_path"] == "/tmp/protein.csv"


def test_enrich_mixed_bias_dataset_expands_each_protein_over_all_recovered_ligands(
    monkeypatch,
):
    rows = pd.DataFrame(
        [
            {
                "query_pair_id": "A__B",
                "query_protein_chain_id": "A",
                "query_ligand_chain_id": "B",
                "reference_key": "public:pdb:1XM1",
                "reference_label": "1XM1",
                "pairing_status": "protein_only",
                "source": "public",
                "dataset_name": "public",
                "pdb_id": "1XM1",
                "protein_pdb_id": "1XM1",
                "ligand_pdb_id": pd.NA,
                "protein_release_date": "2022-01-01",
                "ligand_release_date": pd.NA,
                "protein_source": "public",
                "protein_dataset_name": "public",
                "ligand_source": pd.NA,
                "ligand_dataset_name": pd.NA,
                "sequence_similarity": similarity,
                "ecfp_similarity": 0.0,
                "plot_sequence_similarity": similarity / 100.0,
                "plot_ecfp_similarity": 0.0,
                "sequence": sequence,
                "ligand_id": pd.NA,
                "smiles": pd.NA,
            }
            for sequence, similarity in (("SEQ_CHAIN_1", 80.0), ("SEQ_CHAIN_2", 70.0))
        ]
    )

    monkeypatch.setattr(
        "cofolder.modules.analytics.bias._enrichment._pdb_ligand_similarity_rows",
        lambda **kwargs: [
            {"ligand_id": "LIG_1", "smiles": "CCO", "ecfp_similarity": 0.81},
            {"ligand_id": "LIG_2", "smiles": "CCN", "ecfp_similarity": 0.72},
        ],
    )

    enriched = bias_enrichment._enrich_mixed_bias_dataset_with_pdb_backfill(
        rows,
        protein_queries={"A": "QUERY_SEQUENCE"},
        ligand_queries={"B": "QUERY_SMILES"},
        boltz_cache_path=Path(tempfile.gettempdir()),
        components_cif_path=None,
    )

    assert len(enriched) == 4
    assert enriched["pairing_status"].tolist() == ["paired"] * 4
    assert set(zip(enriched["sequence"], enriched["ligand_id"])) == {
        ("SEQ_CHAIN_1", "LIG_1"),
        ("SEQ_CHAIN_1", "LIG_2"),
        ("SEQ_CHAIN_2", "LIG_1"),
        ("SEQ_CHAIN_2", "LIG_2"),
    }


def test_enrich_mixed_bias_dataset_expands_each_ligand_over_all_recovered_proteins(
    monkeypatch,
):
    rows = pd.DataFrame(
        [
            {
                "query_pair_id": "A__B",
                "query_protein_chain_id": "A",
                "query_ligand_chain_id": "B",
                "reference_key": "public:pdb:1XM1",
                "reference_label": "1XM1",
                "pairing_status": "ligand_only",
                "source": "public",
                "dataset_name": "public",
                "pdb_id": "1XM1",
                "protein_pdb_id": pd.NA,
                "ligand_pdb_id": "1XM1",
                "protein_release_date": pd.NA,
                "ligand_release_date": "2022-01-01",
                "protein_source": pd.NA,
                "protein_dataset_name": pd.NA,
                "ligand_source": "public",
                "ligand_dataset_name": "public",
                "sequence_similarity": pd.NA,
                "ecfp_similarity": similarity,
                "plot_sequence_similarity": pd.NA,
                "plot_ecfp_similarity": similarity,
                "sequence": pd.NA,
                "ligand_id": ligand_id,
                "smiles": smiles,
            }
            for ligand_id, smiles, similarity in (
                ("LIG_1", "CCO", 0.81),
                ("LIG_2", "CCN", 0.72),
            )
        ]
    )
    protein_lookup_view = pd.DataFrame(
        [
            {
                "query_chain_id": "A",
                "pdb_id": "1XM1",
                "release_date": "2022-01-01",
                "source": "public",
                "dataset_name": "public",
                "sequence_similarity": similarity,
                "sequence": sequence,
            }
            for sequence, similarity in (("SEQ_CHAIN_1", 80.0), ("SEQ_CHAIN_2", 70.0))
        ]
    )
    monkeypatch.setattr(
        "cofolder.modules.analytics.bias._enrichment._pdb_protein_similarity_rows",
        lambda **kwargs: (_ for _ in ()).throw(
            AssertionError("MMseqs fallback should not run")
        ),
    )

    enriched = bias_enrichment._enrich_mixed_bias_dataset_with_pdb_backfill(
        rows,
        protein_queries={"A": "QUERY_SEQUENCE"},
        ligand_queries={"B": "QUERY_SMILES"},
        protein_lookup_view=protein_lookup_view,
        boltz_cache_path=Path(tempfile.gettempdir()),
        components_cif_path=None,
    )

    assert len(enriched) == 4
    assert enriched["pairing_status"].tolist() == ["paired"] * 4
    assert set(zip(enriched["sequence"], enriched["ligand_id"])) == {
        ("SEQ_CHAIN_1", "LIG_1"),
        ("SEQ_CHAIN_1", "LIG_2"),
        ("SEQ_CHAIN_2", "LIG_1"),
        ("SEQ_CHAIN_2", "LIG_2"),
    }


def test_enrich_mixed_bias_dataset_with_pdb_backfill_warns_and_keeps_ligand_only_without_mmseqs(
    monkeypatch,
    caplog,
):
    rows = pd.DataFrame(
        [
            {
                "query_pair_id": "A__MG",
                "query_protein_chain_id": "A",
                "query_ligand_chain_id": "MG",
                "reference_key": "public:pdb:1XM1",
                "reference_label": "1XM1",
                "pairing_status": "ligand_only",
                "source": "public",
                "dataset_name": "public",
                "pdb_id": "1XM1",
                "protein_pdb_id": pd.NA,
                "ligand_pdb_id": "1XM1",
                "protein_release_date": pd.NA,
                "ligand_release_date": "2022-01-01",
                "protein_source": pd.NA,
                "protein_dataset_name": pd.NA,
                "protein_source_structure_path": pd.NA,
                "protein_source_reference_path": pd.NA,
                "ligand_source": "public",
                "ligand_dataset_name": "public",
                "ligand_source_structure_path": pd.NA,
                "ligand_source_reference_path": pd.NA,
                "sequence_similarity": 0.0,
                "ecfp_similarity": 0.91,
                "plot_sequence_similarity": 0.0,
                "plot_ecfp_similarity": 0.91,
                "sequence": pd.NA,
                "ligand_id": "MG",
                "smiles": "[Mg+2]",
            },
        ]
    )
    monkeypatch.setattr(
        "cofolder.modules.analytics.bias._enrichment._resolved_mmseqs_bin_token",
        lambda: "",
    )

    with caplog.at_level(logging.WARNING, logger="cofolder.modules.analytics.bias"):
        enriched = bias_enrichment._enrich_mixed_bias_dataset_with_pdb_backfill(
            rows,
            protein_queries={"A": "SEQ_QUERY"},
            ligand_queries={"MG": "[Mg+2]"},
            boltz_cache_path=Path(tempfile.gettempdir()),
            components_cif_path=None,
        )

    assert enriched["pairing_status"].tolist() == ["ligand_only"]
    assert "MMseqs protein fallback unavailable for mixed bias row A__MG" in caplog.text


def test_enrich_mixed_bias_dataset_with_pdb_backfill_warns_on_threshold_drift(
    monkeypatch,
    caplog,
):
    rows = pd.DataFrame(
        [
            {
                "query_pair_id": "A__MG",
                "query_protein_chain_id": "A",
                "query_ligand_chain_id": "MG",
                "reference_key": "public:pdb:1XM1",
                "reference_label": "1XM1",
                "pairing_status": "ligand_only",
                "source": "public",
                "dataset_name": "public",
                "pdb_id": "1XM1",
                "protein_pdb_id": pd.NA,
                "ligand_pdb_id": "1XM1",
                "protein_release_date": pd.NA,
                "ligand_release_date": "2022-01-01",
                "protein_source": pd.NA,
                "protein_dataset_name": pd.NA,
                "protein_source_structure_path": pd.NA,
                "protein_source_reference_path": pd.NA,
                "ligand_source": "public",
                "ligand_dataset_name": "public",
                "ligand_source_structure_path": pd.NA,
                "ligand_source_reference_path": pd.NA,
                "sequence_similarity": 0.0,
                "ecfp_similarity": 0.91,
                "plot_sequence_similarity": 0.0,
                "plot_ecfp_similarity": 0.91,
                "sequence": pd.NA,
                "ligand_id": "MG",
                "smiles": "[Mg+2]",
            },
        ]
    )
    monkeypatch.setattr(
        "cofolder.modules.analytics.bias._enrichment._resolved_mmseqs_bin_token",
        lambda: "mmseqs",
    )
    monkeypatch.setattr(
        "cofolder.modules.analytics.bias._enrichment._pdb_protein_similarity_rows",
        lambda **kwargs: [{"sequence": "SEQ_X", "sequence_similarity": 61.0}],
    )

    with caplog.at_level(logging.WARNING, logger="cofolder.modules.analytics.bias"):
        enriched = bias_enrichment._enrich_mixed_bias_dataset_with_pdb_backfill(
            rows,
            protein_queries={"A": "SEQ_QUERY"},
            ligand_queries={"MG": "[Mg+2]"},
            boltz_cache_path=Path(tempfile.gettempdir()),
            components_cif_path=None,
        )

    assert enriched["pairing_status"].tolist() == ["paired"]
    assert enriched["sequence_similarity"].tolist() == [61.0]
    assert "Protein MMseqs fallback exceeded reference threshold for mixed bias row A__MG" in caplog.text


def test_enrich_mixed_bias_dataset_with_pdb_backfill_writes_progress_checkpoints(
    monkeypatch,
    temp_dir,
    caplog,
):
    rows = []
    for index in range(51):
        pdb_id = f"{index:04d}"
        rows.append(
            {
                "query_pair_id": "A__B",
                "query_protein_chain_id": "A",
                "query_ligand_chain_id": "B",
                "reference_key": f"public:pdb:{pdb_id}",
                "reference_label": pdb_id,
                "pairing_status": "protein_only",
                "source": "public",
                "dataset_name": "public",
                "pdb_id": pdb_id,
                "protein_pdb_id": pdb_id,
                "ligand_pdb_id": pd.NA,
                "protein_release_date": "2022-01-01",
                "ligand_release_date": pd.NA,
                "protein_source": "public",
                "protein_dataset_name": "public",
                "protein_source_structure_path": pd.NA,
                "protein_source_reference_path": pd.NA,
                "ligand_source": pd.NA,
                "ligand_dataset_name": pd.NA,
                "ligand_source_structure_path": pd.NA,
                "ligand_source_reference_path": pd.NA,
                "sequence_similarity": 90.0,
                "ecfp_similarity": 0.0,
                "plot_sequence_similarity": 0.9,
                "plot_ecfp_similarity": 0.0,
                "sequence": "SEQ_A",
                "ligand_id": pd.NA,
                "smiles": pd.NA,
            }
        )

    monkeypatch.setattr(
        "cofolder.modules.analytics.bias._enrichment._pdb_ligand_similarity_rows",
        lambda **kwargs: [{"ligand_id": "L1", "smiles": "CCO", "ecfp_similarity": 0.12}],
    )

    output_path = temp_dir / "bias_training_data_A__B.csv"
    with caplog.at_level(logging.INFO, logger="cofolder.modules.analytics.bias"):
        enriched = bias_enrichment._enrich_mixed_bias_dataset_with_pdb_backfill(
            pd.DataFrame(rows),
            protein_queries={"A": "SEQ_QUERY"},
            ligand_queries={"B": "CCO"},
            boltz_cache_path=Path(tempfile.gettempdir()),
            components_cif_path=None,
            progress_output_path=output_path,
            progress_label="A__B",
        )

    assert output_path.exists()
    assert len(enriched) == 51
    assert set(enriched["pairing_status"]) == {"paired"}
    assert "Bias lookup progress for A__B: 0/51 lookups completed" in caplog.text
    assert "Bias lookup progress for A__B: 10/51 lookups completed" in caplog.text
    assert "Bias lookup progress for A__B: 50/51 lookups completed" in caplog.text
    assert "Bias lookup progress for A__B: 51/51 lookups completed" in caplog.text


def test_enrich_same_type_ligand_pair_dataset_writes_progress_checkpoints(
    monkeypatch,
    temp_dir,
    caplog,
):
    rows = []
    for index in range(51):
        pdb_id = f"{index:04d}"
        rows.append(
            {
                "query_pair_id": "B__MG",
                "component_1_id": "B",
                "component_1_type": "ligand",
                "component_2_id": "MG",
                "component_2_type": "ligand",
                "reference_key": f"public:pdb:{pdb_id}",
                "reference_label": pdb_id,
                "reference_pdb_id": pdb_id,
                "pairing_status": "component_1_only",
                "source": "public",
                "dataset_name": "public",
                "component_1_reference_key": f"public:pdb:{pdb_id}",
                "component_1_pdb_id": pdb_id,
                "component_1_release_date": "2022-01-01",
                "component_1_source": "public",
                "component_1_dataset_name": "public",
                "component_1_source_structure_path": pd.NA,
                "component_1_source_reference_path": pd.NA,
                "component_2_reference_key": pd.NA,
                "component_2_pdb_id": pd.NA,
                "component_2_release_date": pd.NA,
                "component_2_source": pd.NA,
                "component_2_dataset_name": pd.NA,
                "component_2_source_structure_path": pd.NA,
                "component_2_source_reference_path": pd.NA,
                "component_1_similarity": 1.0,
                "component_2_similarity": 0.0,
                "plot_component_1_similarity": 1.0,
                "plot_component_2_similarity": 0.0,
                "_component_1_ligand_id": f"L{index}",
                "_component_1_smiles": "CCO",
                "_component_2_ligand_id": pd.NA,
                "_component_2_smiles": pd.NA,
            }
        )

    monkeypatch.setattr(
        "cofolder.modules.analytics.bias._enrichment._pdb_ligand_similarity_rows",
        lambda **kwargs: [{"ligand_id": "MG", "smiles": "[Mg+2]", "ecfp_similarity": 0.2}],
    )

    output_path = temp_dir / "bias_ligand_pair_data_B__MG.csv"
    with caplog.at_level(logging.INFO, logger="cofolder.modules.analytics.bias"):
        enriched = bias_enrichment._enrich_same_type_ligand_pair_dataset_with_pdb_backfill(
            pd.DataFrame(rows),
            ligand_queries={"B": "CCO", "MG": "[Mg+2]"},
            boltz_cache_path=Path(tempfile.gettempdir()),
            components_cif_path=None,
            progress_output_path=output_path,
            progress_label="B__MG",
        )

    assert output_path.exists()
    assert len(enriched) == 51
    assert set(enriched["pairing_status"]) == {"paired"}
    assert "Bias lookup progress for B__MG: 0/51 lookups completed" in caplog.text
    assert "Bias lookup progress for B__MG: 10/51 lookups completed" in caplog.text
    assert "Bias lookup progress for B__MG: 50/51 lookups completed" in caplog.text
    assert "Bias lookup progress for B__MG: 51/51 lookups completed" in caplog.text


def test_build_bias_training_dataset_sets_apo_ligand_similarity_to_zero():
    protein_view = pd.DataFrame(
        [
            {
                "query_chain_id": "A",
                "pdb_id": "5ZOD",
                "release_date": "2022-01-01",
                "source": "public",
                "dataset_name": "public",
                "sequence_similarity": 94.1,
                "sequence": "SEQ_A",
            },
        ]
    )
    ligand_view = pd.DataFrame(
        [
            {
                "query_chain_id": "B",
                "pdb_id": pd.NA,
                "release_date": pd.NA,
                "source": "custom",
                "dataset_name": "private_ligands",
                "source_reference_path": "/tmp/custom_ligand.sdf",
                "ligand_id": "L1",
                "ecfp_similarity": 0.91,
                "smiles": "CCO",
            },
        ]
    )

    result = _build_bias_training_dataset(
        protein_view=protein_view,
        ligand_views={"B": ligand_view},
        protein_lookup_view=protein_view,
        ligand_lookup_views={"B": ligand_view},
    )

    apo_rows = result[result["reference_key"] == "public:pdb:5ZOD"].copy()
    ligand_only = result[result["pairing_status"] == "ligand_only"].copy()

    assert len(apo_rows) == 1
    assert apo_rows["pairing_status"].tolist() == ["protein_only"]
    assert apo_rows["protein_pdb_id"].tolist() == ["5ZOD"]
    assert apo_rows["ligand_pdb_id"].isna().all()
    assert apo_rows["ecfp_similarity"].tolist() == [0.0]
    assert apo_rows["plot_ecfp_similarity"].tolist() == [0.0]
    assert len(ligand_only) == 1


def test_build_bias_training_dataset_leaves_ligand_only_protein_similarity_unavailable():
    protein_view = pd.DataFrame(
        [
            {
                "query_chain_id": "A",
                "pdb_id": "5ZOD",
                "release_date": "2022-01-01",
                "source": "public",
                "dataset_name": "public",
                "sequence_similarity": 94.1,
                "sequence": "SEQ_A",
            },
        ]
    )
    ligand_view = pd.DataFrame(
        [
            {
                "query_chain_id": "B",
                "pdb_id": "1XM1",
                "release_date": "2022-01-01",
                "source": "public",
                "dataset_name": "public",
                "ligand_id": "L1",
                "ecfp_similarity": 0.38,
                "smiles": "CCO",
            },
        ]
    )

    result = _build_bias_training_dataset(
        protein_view=protein_view,
        ligand_views={"B": ligand_view},
        protein_lookup_view=protein_view,
        ligand_lookup_views={"B": ligand_view},
    )

    ligand_only = result[result["reference_key"] == "public:pdb:1XM1"].copy()

    assert len(ligand_only) == 1
    assert ligand_only["pairing_status"].tolist() == ["ligand_only"]
    assert ligand_only["sequence_similarity"].isna().all()
    assert ligand_only["sequence_similarity_pairwise"].isna().all()
    assert ligand_only["plot_sequence_similarity"].isna().all()
    assert ligand_only["plot_sequence_similarity_pairwise"].isna().all()
    assert ligand_only["sequence_similarity_method"].tolist() == ["unavailable"]
    assert ligand_only["protein_pdb_id"].isna().all()


def test_plot_bias_reference_overlap_keeps_one_sided_threshold_hits_and_removes_decorations(
    temp_dir,
    monkeypatch,
):
    captured: dict[str, object] = {}

    def _capture_save(fig, output_path):
        ax = fig.axes[0]
        offsets = []
        for collection in ax.collections:
            offsets.extend(tuple(round(float(value), 2) for value in point) for point in collection.get_offsets())
        captured["offsets"] = offsets
        captured["legend"] = ax.get_legend()
        captured["texts"] = [text.get_text() for text in ax.texts]
        path = temp_dir / Path(output_path).name
        path.write_text("plot", encoding="utf-8")
        return path

    monkeypatch.setattr(plots_module, "_save_figure", _capture_save)

    saved_paths = plot_bias_reference_overlap(
        pd.DataFrame(
            [
                {
                    "reference_label": "protein-only",
                    "pairing_status": "protein_only",
                    "source": "public",
                    "query_protein_chain_id": "A",
                    "query_ligand_chain_id": "B",
                    "plot_sequence_similarity": 0.94,
                    "plot_ecfp_similarity": pd.NA,
                },
                {
                    "reference_label": "ligand-only",
                    "pairing_status": "ligand_only",
                    "source": "custom",
                    "query_protein_chain_id": "A",
                    "query_ligand_chain_id": "B",
                    "plot_sequence_similarity": pd.NA,
                    "plot_ecfp_similarity": 0.38,
                },
                {
                    "reference_label": "below-threshold",
                    "pairing_status": "paired",
                    "source": "public",
                    "query_protein_chain_id": "A",
                    "query_ligand_chain_id": "B",
                    "plot_sequence_similarity": 0.24,
                    "plot_ecfp_similarity": 0.34,
                },
                {
                    "reference_label": "mmseqs-backfilled",
                    "pairing_status": "paired",
                    "source": "public",
                    "query_protein_chain_id": "A",
                    "query_ligand_chain_id": "B",
                    "plot_sequence_similarity": 0.125,
                    "plot_ecfp_similarity": 0.38,
                },
                {
                    "reference_label": "exact-thresholds",
                    "pairing_status": "paired",
                    "source": "public",
                    "query_protein_chain_id": "A",
                    "query_ligand_chain_id": "B",
                    "plot_sequence_similarity": 0.25,
                    "plot_ecfp_similarity": 0.35,
                },
                {
                    "reference_label": "protein-boundary-ligand-pass",
                    "pairing_status": "paired",
                    "source": "public",
                    "query_protein_chain_id": "A",
                    "query_ligand_chain_id": "B",
                    "plot_sequence_similarity": 0.25,
                    "plot_ecfp_similarity": 0.36,
                },
                {
                    "reference_label": "ligand-boundary-protein-pass",
                    "pairing_status": "paired",
                    "source": "public",
                    "query_protein_chain_id": "A",
                    "query_ligand_chain_id": "B",
                    "plot_sequence_similarity": 0.26,
                    "plot_ecfp_similarity": 0.35,
                },
            ]
        ),
        output_dir=temp_dir / "results",
    )

    assert len(saved_paths) == 2
    assert set(captured["offsets"]) == {
        (0.0, 0.94),
        (0.35, 0.26),
        (0.36, 0.25),
        (0.38, 0.0),
        (0.38, 0.12),
    }
    assert captured["legend"] is None
    assert captured["texts"] == []


def test_plot_reference_overlap_scatter_supports_same_type_axes(temp_dir, monkeypatch):
    captured: dict[str, object] = {}

    def _capture_save(fig, output_path):
        ax = fig.axes[0]
        offsets = []
        for collection in ax.collections:
            offsets.extend(tuple(round(float(value), 2) for value in point) for point in collection.get_offsets())
        captured["offsets"] = offsets
        captured["legend"] = ax.get_legend()
        captured["texts"] = [text.get_text() for text in ax.texts]
        path = temp_dir / Path(output_path).name
        path.write_text("plot", encoding="utf-8")
        return path

    monkeypatch.setattr(plots_module, "_save_figure", _capture_save)

    saved_paths = plot_reference_overlap_scatter(
        pd.DataFrame(
            [
                {
                    "query_pair_id": "A__C",
                    "component_1_id": "A",
                    "component_2_id": "C",
                    "source": "public",
                    "plot_component_1_similarity": 0.94,
                    "plot_component_2_similarity": pd.NA,
                },
                {
                    "query_pair_id": "A__C",
                    "component_1_id": "A",
                    "component_2_id": "C",
                    "source": "custom",
                    "plot_component_1_similarity": pd.NA,
                    "plot_component_2_similarity": 0.31,
                },
                {
                    "query_pair_id": "A__C",
                    "component_1_id": "A",
                    "component_2_id": "C",
                    "source": "public",
                    "plot_component_1_similarity": 0.24,
                    "plot_component_2_similarity": 0.2,
                },
            ]
        ),
        output_dir=temp_dir / "results",
        file_stem="bias_protein_pair_scatter_A__C",
        x_col="plot_component_1_similarity",
        y_col="plot_component_2_similarity",
        x_threshold=0.25,
        y_threshold=0.25,
        x_label="Protein reference overlap (component 1)",
        y_label="Protein reference overlap (component 2)",
        title="Bias protein-protein reference-overlap diagnostic",
        query_1_col="component_1_id",
        query_2_col="component_2_id",
    )

    assert len(saved_paths) == 2
    assert set(captured["offsets"]) == {(0.94, 0.0), (0.0, 0.31)}
    assert captured["legend"] is None
    assert captured["texts"] == []


def test_build_protein_training_view_keeps_distinct_custom_paths_for_same_sequence():
    proteins_df = pd.DataFrame(
        [
            {
                "pdb_id": pd.NA,
                "release_date": pd.NA,
                "source": "custom",
                "dataset_name": "private_set",
                "source_structure_path": "/tmp/reference_a.pdb",
                "source_reference_path": pd.NA,
                "sequence": "MKRAAT",
            },
            {
                "pdb_id": pd.NA,
                "release_date": pd.NA,
                "source": "custom",
                "dataset_name": "private_set",
                "source_structure_path": "/tmp/reference_b.pdb",
                "source_reference_path": pd.NA,
                "sequence": "MKRAAT",
            },
        ]
    )

    result = _build_protein_training_view(proteins_df, {"A": "MKRAAT"})

    assert len(result) == 2
    assert result["source_structure_path"].tolist() == [
        "/tmp/reference_a.pdb",
        "/tmp/reference_b.pdb",
    ]
    assert result["sequence_similarity"].isna().all()
    assert result["sequence_similarity_pairwise"].tolist() == [100.0, 100.0]
    assert result["sequence_similarity_method"].tolist() == [
        "pairwise_aligner",
        "pairwise_aligner",
    ]


def test_training_views_use_strict_manuscript_thresholds(monkeypatch):
    proteins_df = pd.DataFrame(
        [
            {
                "pdb_id": pdb_id,
                "release_date": "2022-01-01",
                "source": "public",
                "dataset_name": "public",
                "source_structure_path": pd.NA,
                "source_reference_path": pd.NA,
                "sequence": sequence,
                "sequence_similarity": similarity,
            }
            for pdb_id, sequence, similarity in (
                ("1AAA", "SEQ_BOUNDARY", 25.0),
                ("1AAB", "SEQ_ABOVE", 25.1),
            )
        ]
    )
    protein_view = _build_protein_training_view(proteins_df, {"A": "QUERY"})

    assert protein_view["pdb_id"].tolist() == ["1AAB"]

    chain_df = pd.DataFrame(
        [{"CHAIN_ID": "B", "ENTITY_TYPE": "ligand", "ligand_molecule_id": "QUERY"}]
    )
    ligands_df = pd.DataFrame(
        [
            {
                "pdb_id": "1AAA",
                "release_date": "2022-01-01",
                "source": "public",
                "dataset_name": "public",
                "ligand_id": "BOUNDARY",
                "smiles": "CCO",
            },
            {
                "pdb_id": "1AAB",
                "release_date": "2022-01-01",
                "source": "public",
                "dataset_name": "public",
                "ligand_id": "ABOVE",
                "smiles": "CCN",
            },
        ]
    )
    monkeypatch.setattr(
        bias_similarity,
        "_ligand_similarity_series",
        lambda *args, **kwargs: pd.Series([0.35, 0.36], index=ligands_df.index),
    )

    ligand_view = bias_similarity._build_ligand_training_views(
        chain_df=chain_df,
        ligands_df=ligands_df,
        ligand_queries={"B": "QUERY_SMILES"},
        boltz_cache_path=Path(tempfile.gettempdir()),
    )["B"]

    assert ligand_view["pdb_id"].tolist() == ["1AAB"]


def test_build_bias_training_dataset_marks_protein_only_rows_with_provenance():
    protein_view = pd.DataFrame(
        [
            {
                "query_chain_id": "A",
                "pdb_id": "1ABC",
                "release_date": "2022-01-01",
                "source": "public",
                "dataset_name": "public",
                "sequence_similarity": 95.0,
                "sequence": "SEQ1",
            },
        ]
    )

    result = _build_bias_training_dataset(
        protein_view=protein_view,
        ligand_views={},
    )

    assert result["pairing_status"].tolist() == ["protein_only"]
    assert result["source"].tolist() == ["public"]
    assert result["dataset_name"].tolist() == ["public"]
    assert result["protein_pdb_id"].tolist() == ["1ABC"]
    assert result["sequence_similarity"].tolist() == [95.0]


def test_protein_reference_rejects_mixed_mmseqs_and_pairwise_scores(temp_dir):
    reference_path = temp_dir / "mixed_protein_reference.csv"
    reference_path.write_text(
        "sequence,sequence_similarity,sequence_similarity_pairwise,"
        "sequence_similarity_method\n"
        "MAAA,24.0,75.0,mmseqs_pident\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="cannot contain both"):
        bias_references._load_custom_protein_training(reference_path)


def test_combined_bias_rows_keep_mmseqs_and_pairwise_in_separate_columns():
    protein_view = pd.DataFrame(
        [
            {
                "query_chain_id": "A",
                "pdb_id": "1ABC",
                "source": "public",
                "dataset_name": "public",
                "sequence_similarity": 24.0,
                "sequence_similarity_pairwise": pd.NA,
                "sequence_similarity_method": "mmseqs_pident",
                "sequence": "MAAA",
            },
            {
                "query_chain_id": "A",
                "pdb_id": pd.NA,
                "source": "custom",
                "dataset_name": "custom",
                "sequence_similarity": pd.NA,
                "sequence_similarity_pairwise": 75.0,
                "sequence_similarity_method": "pairwise_aligner",
                "sequence": "MATA",
            },
        ]
    )

    result = _build_bias_training_dataset(protein_view=protein_view, ligand_views={})

    assert not (
        result["sequence_similarity"].notna()
        & result["sequence_similarity_pairwise"].notna()
    ).any()
    assert set(
        result.loc[result["sequence_similarity"].notna(), "sequence_similarity_method"]
    ) == {"mmseqs_pident"}
    assert set(
        result.loc[
            result["sequence_similarity_pairwise"].notna(),
            "sequence_similarity_method",
        ]
    ) == {"pairwise_aligner"}


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


def test_apply_bias_metrics_clears_stale_plot_and_ligand_outputs_when_plot_is_skipped(temp_dir):
    protein_ref = temp_dir / "protein_training.csv"
    protein_ref.write_text(
        "pdb_id,release_date,sequence\n"
        "1ABC,2022-01-01,MAAA\n",
        encoding="utf-8",
    )

    output_dir = temp_dir / "results"
    output_dir.mkdir()
    (output_dir / "ligand_training_data_B.csv").write_text("stale\n", encoding="utf-8")
    (output_dir / "bias_reference_overlap_scatter.png").write_text("stale\n", encoding="utf-8")
    (output_dir / "bias_reference_overlap_scatter.pdf").write_text("stale\n", encoding="utf-8")

    system_df = pd.DataFrame([{"model_name": "system", "repeat": 1, "diffusion_sample": 0}])
    chain_df = pd.DataFrame([{"CHAIN_ID": "A", "ENTITY_TYPE": "protein"}])
    sys_obj = _MockSystem(sequences=[{"protein": {"id": "A", "sequence": "MAAA"}}])

    apply_bias_metrics(
        system_df=system_df,
        chain_df=chain_df,
        sys_obj=sys_obj,
        protein_training_data_path=protein_ref,
        ligand_training_data_path=None,
        release_cutoff="2023-06-01",
        output_dir=output_dir,
    )

    assert not (output_dir / "ligand_training_data_B.csv").exists()
    assert not (output_dir / "bias_reference_overlap_scatter.png").exists()
    assert not (output_dir / "bias_reference_overlap_scatter.pdf").exists()
    assert (output_dir / "bias_reference_overlap_scatter.skipped.txt").exists()
