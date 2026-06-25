"""Tests for bias metric analytics."""

from datetime import date
from pathlib import Path

import pandas as pd
import pytest

from cofolder.modules.analytics import plots as plots_module
from cofolder.modules.analytics.bias import (
    BIAS_TRAINING_DATA_COLUMNS,
    _build_bias_training_dataset,
    _build_protein_training_view,
    _load_ligand_training,
    _synthetic_reference_key,
    apply_bias_metrics,
)
from cofolder.modules.analytics.plots import plot_bias_reference_overlap


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
    assert bias_training["source"].tolist() == ["custom", "public", "mixed", "mixed"]
    assert bias_training["query_protein_chain_id"].tolist() == ["A", "A", "A", "A"]
    assert bias_training["query_ligand_chain_id"].tolist() == ["B", "B", "B", "B"]
    assert bias_training["pairing_status"].tolist() == ["paired", "paired", "synthetic_paired", "synthetic_paired"]
    assert bias_training["plot_sequence_similarity"].tolist()[:2] == [1.0, 0.83]
    assert bias_training["plot_ecfp_similarity"].tolist()[:2] == [1.0, 0.4]
    assert set(
        zip(
            bias_training["plot_sequence_similarity"].tolist()[2:],
            bias_training["plot_ecfp_similarity"].tolist()[2:],
        )
    ) == {(1.0, 0.4), (0.83, 1.0)}
    assert (temp_dir / "results" / "bias_reference_overlap_scatter.png").exists()
    assert (temp_dir / "results" / "bias_reference_overlap_scatter.pdf").exists()
    assert set(summary["nearest_overall_source"]) == {"custom"}
    assert summary["custom_changes_nearest_reference"].tolist() == [True, True]


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


def test_build_bias_training_dataset_adds_synthetic_rows_for_disjoint_public_hits():
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

    result = _build_bias_training_dataset(
        protein_view=protein_view,
        ligand_views={"B": ligand_view},
    )

    synthetic = result[result["pairing_status"] == "synthetic_paired"].copy()
    protein_only = result[result["pairing_status"] == "protein_only"].copy()
    ligand_only = result[result["pairing_status"] == "ligand_only"].copy()

    assert len(result) == 8
    assert len(synthetic) == 4
    assert len(protein_only) == 2
    assert len(ligand_only) == 2
    assert synthetic["reference_key"].str.startswith("synthetic:").all()
    assert synthetic["reference_label"].str.contains(" x ").all()
    assert synthetic["pdb_id"].tolist() == ["synthetic_pair"] * 4
    assert synthetic["source"].tolist() == ["public"] * 4
    assert synthetic["query_pair_id"].tolist() == ["A__B"] * 4
    assert set(zip(synthetic["protein_pdb_id"], synthetic["ligand_pdb_id"])) == {
        ("5ZOD", "1XM1"),
        ("5ZOD", "6AGT"),
        ("3Q8K", "1XM1"),
        ("3Q8K", "6AGT"),
    }


def test_build_bias_training_dataset_marks_mixed_provenance_synthetic_rows():
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
    )

    synthetic = result[result["pairing_status"] == "synthetic_paired"].copy()

    assert len(synthetic) == 1
    assert synthetic["source"].tolist() == ["mixed"]
    assert synthetic["protein_source"].tolist() == ["public"]
    assert synthetic["ligand_source"].tolist() == ["custom"]
    assert synthetic["protein_pdb_id"].tolist() == ["5ZOD"]
    assert synthetic["ligand_pdb_id"].isna().all()
    assert synthetic["pdb_id"].tolist() == ["synthetic_pair"]
    assert synthetic["reference_key"].str.startswith("synthetic:").all()


def test_synthetic_reference_key_avoids_separator_collisions():
    left_then_right = _synthetic_reference_key("public:a__b", "custom:c")
    split_other_way = _synthetic_reference_key("public:a", "custom:b__c")

    assert left_then_right != split_other_way
    assert left_then_right.startswith("synthetic:")
    assert split_other_way.startswith("synthetic:")


def test_build_bias_training_dataset_keeps_synthetic_rows_scoped_to_each_query_pair():
    protein_view = pd.DataFrame(
        [
            {
                "query_chain_id": "A",
                "pdb_id": "1AAA",
                "release_date": "2022-01-01",
                "source": "public",
                "dataset_name": "public",
                "sequence_similarity": 90.0,
                "sequence": "SEQ_A",
            },
            {
                "query_chain_id": "C",
                "pdb_id": "2CCC",
                "release_date": "2022-01-01",
                "source": "public",
                "dataset_name": "public",
                "sequence_similarity": 91.0,
                "sequence": "SEQ_C",
            },
        ]
    )
    ligand_views = {
        "B": pd.DataFrame(
            [
                {
                    "query_chain_id": "B",
                    "pdb_id": "1BBB",
                    "release_date": "2022-01-01",
                    "source": "public",
                    "dataset_name": "public",
                    "ligand_id": "LB",
                    "ecfp_similarity": 0.8,
                    "smiles": "CCO",
                }
            ]
        ),
        "D": pd.DataFrame(
            [
                {
                    "query_chain_id": "D",
                    "pdb_id": "2DDD",
                    "release_date": "2022-01-01",
                    "source": "public",
                    "dataset_name": "public",
                    "ligand_id": "LD",
                    "ecfp_similarity": 0.7,
                    "smiles": "CCN",
                }
            ]
        ),
    }

    result = _build_bias_training_dataset(
        protein_view=protein_view,
        ligand_views=ligand_views,
    )

    synthetic = result[result["pairing_status"] == "synthetic_paired"].copy()

    assert set(synthetic["query_pair_id"]) == {"A__B", "A__D", "C__B", "C__D"}
    for _, row in synthetic.iterrows():
        if row["query_pair_id"] == "A__B":
            assert row["protein_pdb_id"] == "1AAA"
            assert row["ligand_pdb_id"] == "1BBB"
        elif row["query_pair_id"] == "A__D":
            assert row["protein_pdb_id"] == "1AAA"
            assert row["ligand_pdb_id"] == "2DDD"
        elif row["query_pair_id"] == "C__B":
            assert row["protein_pdb_id"] == "2CCC"
            assert row["ligand_pdb_id"] == "1BBB"
        else:
            assert row["protein_pdb_id"] == "2CCC"
            assert row["ligand_pdb_id"] == "2DDD"


def test_plot_bias_reference_overlap_distinguishes_synthetic_rows(temp_dir, monkeypatch):
    captured: dict[str, object] = {}

    def _capture_save(fig, output_path):
        ax = fig.axes[0]
        _, labels = ax.get_legend_handles_labels()
        captured["labels"] = labels
        captured["annotations"] = [text.get_text() for text in ax.texts]
        path = Path(output_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("plot", encoding="utf-8")
        return path

    monkeypatch.setattr(plots_module, "_save_figure", _capture_save)

    saved_paths = plot_bias_reference_overlap(
        pd.DataFrame(
            [
                {
                    "reference_label": "1ABC",
                    "pairing_status": "paired",
                    "source": "public",
                    "plot_sequence_similarity": 0.95,
                    "plot_ecfp_similarity": 0.91,
                },
                {
                    "reference_label": "2XYZ x 3LMN",
                    "pairing_status": "synthetic_paired",
                    "source": "public",
                    "plot_sequence_similarity": 0.99,
                    "plot_ecfp_similarity": 0.97,
                },
            ]
        ),
        output_dir=temp_dir / "results",
    )

    assert len(saved_paths) == 2
    assert "Public references" in captured["labels"]
    assert "Public references synthetic combinations" in captured["labels"]
    assert plots_module._annotation_label(
        pd.Series(
            {
                "reference_label": "2XYZ x 3LMN",
                "pairing_status": "synthetic_paired",
            }
        )
    ) == "Synthetic: 2XYZ x 3LMN"


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
    assert result["sequence_similarity"].tolist() == [100.0, 100.0]


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
