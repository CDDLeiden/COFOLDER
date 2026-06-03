"""Tests for cofolder.recipes.bias module."""

from __future__ import annotations

import pandas as pd
import pytest
import yaml

from cofolder.modules.input.system import System
from cofolder.recipes.bias import Bias, build_bias_dataframes


class TestBuildBiasDataframes:
    def test_scaffolds_current_system_format(self):
        sys_obj = System(
            system={
                "sequences": [
                    {"protein": {"id": ["A", "C"], "sequence": "MKRAAT"}},
                    {"ligand": {"id": "B", "ccd": "EDO"}},
                ]
            }
        )

        system_df, chain_df = build_bias_dataframes(
            sys_obj,
            system_name="standalone_bias",
        )

        assert list(system_df.columns) == ["model_name", "repeat", "diffusion_sample"]
        assert system_df.iloc[0]["model_name"] == "standalone_bias"
        assert chain_df["CHAIN_ID"].tolist() == ["A", "C", "B"]
        assert chain_df["ENTITY_TYPE"].tolist() == ["protein", "protein", "ligand"]
        assert chain_df["ligand_molecule_id"].tolist() == ["protein_A", "protein_C", "EDO"]


class TestBiasRun:
    def test_run_supports_custom_only_references_and_writes_summary(self, temp_dir):
        system_path = temp_dir / "system.yaml"
        system_path.write_text(
            yaml.safe_dump(
                {
                    "sequences": [
                        {"protein": {"id": "A", "sequence": "MKRAAT"}},
                        {"ligand": {"id": "B", "smiles": "CCO"}},
                    ]
                }
            ),
            encoding="utf-8",
        )
        structure_ref = temp_dir / "custom_reference.pdb"
        structure_ref.write_text("HEADER CUSTOM\n", encoding="utf-8")

        custom_protein = temp_dir / "custom_protein.csv"
        custom_ligand = temp_dir / "custom_ligand.csv"
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

        bias = Bias(
            wrk_dir=str(temp_dir),
            system_path=str(system_path),
            custom_protein_reference_path=str(custom_protein),
            custom_ligand_reference_path=str(custom_ligand),
        )

        system_df, chain_df = bias.run()

        output_dir = temp_dir / "results" / "bias_train"
        protein_view = pd.read_csv(output_dir / "protein_training_data.csv")
        ligand_view = pd.read_csv(output_dir / "ligand_training_data_B.csv")
        bias_training = pd.read_csv(output_dir / "bias_training_data.csv")
        summary = pd.read_csv(output_dir / "reference_landscape_summary.csv")

        assert float(system_df["bias_prot_sim_train_max"].iloc[0]) == 100.0
        assert float(system_df["bias_lig_sim_train_max"].iloc[0]) == 1.0
        assert protein_view["source"].tolist() == ["custom"]
        assert protein_view["dataset_name"].tolist() == ["private_proteins"]
        assert ligand_view["source"].tolist() == ["custom"]
        assert ligand_view["dataset_name"].tolist() == ["private_ligands"]
        assert bias_training["source"].tolist() == ["custom"]
        assert bias_training["pairing_status"].tolist() == ["paired"]
        assert (output_dir / "bias_reference_overlap_scatter.png").exists()
        assert (output_dir / "bias_reference_overlap_scatter.pdf").exists()
        assert set(summary["nearest_overall_source"]) == {"custom"}
        assert summary["custom_changes_nearest_reference"].tolist() == [True, True]
        assert chain_df["CHAIN_ID"].tolist() == ["A", "B"]

    def test_run_writes_bias_outputs_for_ccd_backed_ligands(self, temp_dir):
        system_path = temp_dir / "system.yaml"
        system_path.write_text(
            yaml.safe_dump(
                {
                    "sequences": [
                        {"protein": {"id": "A", "sequence": "MKRAAT"}},
                        {"ligand": {"id": ["B", "C"], "ccd": "EDO"}},
                    ]
                }
            ),
            encoding="utf-8",
        )

        components_cif = temp_dir / "components.cif"
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

        protein_ref = temp_dir / "protein_training.csv"
        ligand_ref = temp_dir / "ligand_training.csv"
        protein_ref.write_text(
            "pdb_id,release_date,sequence,sequence_similarity\n"
            "1ABC,2022-01-01,MKRAAT,100.0\n",
            encoding="utf-8",
        )
        ligand_ref.write_text(
            "pdb_id,release_date,ligand_id,smiles,ecfp_similarity\n"
            "1ABC,2022-01-01,NOTEDO,CCO,0.6\n",
            encoding="utf-8",
        )

        bias = Bias(
            wrk_dir=str(temp_dir),
            system_path=str(system_path),
            protein_training_data_path=str(protein_ref),
            ligand_training_data_path=str(ligand_ref),
            bias_training_components_cif=str(components_cif),
        )

        system_df, chain_df = bias.run()

        output_dir = temp_dir / "results" / "bias_train"
        assert (output_dir / "system_metrics.csv").exists()
        assert (output_dir / "chain_metrics.csv").exists()
        assert (output_dir / "protein_training_data.csv").exists()
        assert (output_dir / "ligand_training_data_B.csv").exists()
        assert (output_dir / "bias_training_data.csv").exists()
        assert (output_dir / "bias_reference_overlap_scatter.png").exists()
        assert (output_dir / "bias_reference_overlap_scatter.pdf").exists()

        assert float(system_df["bias_prot_sim_train_max"].iloc[0]) == 100.0
        assert float(system_df["bias_lig_sim_train_max"].iloc[0]) == 1.0

        ligand_rows = chain_df[chain_df["ENTITY_TYPE"] == "ligand"].copy()
        assert ligand_rows["CHAIN_ID"].tolist() == ["B", "C"]
        assert ligand_rows["bias_lig_sim_train"].astype(float).tolist() == [1.0, 1.0]

    def test_run_build_mode_defaults_ligand_training_output(self, monkeypatch, temp_dir):
        system_path = temp_dir / "system.yaml"
        system_path.write_text(
            yaml.safe_dump(
                {
                    "sequences": [
                        {"protein": {"id": "A", "sequence": "MKRAAT"}},
                        {"ligand": {"id": "B", "smiles": "CCO", "ccd": "ETH"}},
                    ]
                }
            ),
            encoding="utf-8",
        )

        components_cif = temp_dir / "components.cif"
        components_cif.write_text("data_components\n", encoding="utf-8")
        protein_ref = temp_dir / "generated" / "protein_training.csv"

        calls: list[dict[str, object]] = []

        def _fake_build(**kwargs):
            calls.append(kwargs)
            kwargs["output_protein_csv"].write_text(
                "pdb_id,release_date,sequence,sequence_similarity\n"
                "1ABC,2022-01-01,MKRAAT,100.0\n",
                encoding="utf-8",
            )
            kwargs["output_ligand_csv"].write_text(
                "pdb_id,release_date,ligand_id,smiles,ecfp_similarity\n"
                "1ABC,2022-01-01,ETH,CCO,1.0\n",
                encoding="utf-8",
            )

        monkeypatch.setattr("cofolder.recipes.bias.run_build_bias_training_data", _fake_build)

        bias = Bias(
            wrk_dir=str(temp_dir),
            system_path=str(system_path),
            protein_training_data_path=str(protein_ref),
            build_bias_training_data=True,
            bias_training_components_cif=str(components_cif),
        )

        bias.run()

        assert len(calls) == 1
        assert calls[0]["output_protein_csv"] == protein_ref
        assert calls[0]["output_ligand_csv"] == temp_dir / "results" / "bias_train" / "ligand_training_data.csv"

        output_dir = temp_dir / "results" / "bias_train"
        system_out = pd.read_csv(output_dir / "system_metrics.csv")
        chain_out = pd.read_csv(output_dir / "chain_metrics.csv")

        assert float(system_out["bias_prot_sim_train_max"].iloc[0]) == 100.0
        assert float(system_out["bias_lig_sim_train_max"].iloc[0]) == 1.0
        assert {"CHAIN_ID", "ENTITY_TYPE", "ligand_molecule_id"}.issubset(chain_out.columns)
        assert (output_dir / "protein_training_data.csv").exists()
        assert (output_dir / "ligand_training_data_B.csv").exists()
        assert (output_dir / "bias_training_data.csv").exists()
        assert (output_dir / "reference_landscape_summary.csv").exists()
        assert (output_dir / "bias_reference_overlap_scatter.png").exists()
        assert (output_dir / "bias_reference_overlap_scatter.pdf").exists()

    def test_run_rejects_unresolved_ccd_ligands_without_components_or_matching_training_id(
        self,
        monkeypatch,
        temp_dir,
    ):
        home_dir = temp_dir / "home"
        home_dir.mkdir()
        monkeypatch.setenv("HOME", str(home_dir))

        system_path = temp_dir / "system.yaml"
        system_path.write_text(
            yaml.safe_dump(
                {
                    "sequences": [
                        {"protein": {"id": "A", "sequence": "MKRAAT"}},
                        {"ligand": {"id": "B", "ccd": "EDO"}},
                    ]
                }
            ),
            encoding="utf-8",
        )

        protein_ref = temp_dir / "protein_training.csv"
        ligand_ref = temp_dir / "ligand_training.csv"
        protein_ref.write_text(
            "pdb_id,release_date,sequence,sequence_similarity\n"
            "1ABC,2022-01-01,MKRAAT,100.0\n",
            encoding="utf-8",
        )
        ligand_ref.write_text(
            "pdb_id,release_date,ligand_id,smiles,ecfp_similarity\n"
            "1ABC,2022-01-01,NOTEDO,CCO,0.6\n",
            encoding="utf-8",
        )

        bias = Bias(
            wrk_dir=str(temp_dir),
            system_path=str(system_path),
            protein_training_data_path=str(protein_ref),
            ligand_training_data_path=str(ligand_ref),
        )

        with pytest.raises(ValueError, match="Could not resolve SMILES for CCD-backed ligand chain"):
            bias.run()

    def test_run_allows_ligand_only_bias_when_bias_chains_select_only_ligands(self, temp_dir):
        system_path = temp_dir / "system.yaml"
        system_path.write_text(
            yaml.safe_dump(
                {
                    "sequences": [
                        {"protein": {"id": "A", "sequence": "MKRAAT"}},
                        {"ligand": {"id": "B", "smiles": "CCO"}},
                    ]
                }
            ),
            encoding="utf-8",
        )
        custom_ligand = temp_dir / "custom_ligand.csv"
        custom_ligand.write_text("smiles\nCCO\n", encoding="utf-8")

        bias = Bias(
            wrk_dir=str(temp_dir),
            system_path=str(system_path),
            custom_ligand_reference_path=str(custom_ligand),
            bias_chains=["B"],
        )

        system_df, chain_df = bias.run()

        assert float(system_df["bias_lig_sim_train_max"].iloc[0]) == 1.0
        assert "bias_prot_sim_train_max" not in system_df.columns or pd.isna(
            system_df.get("bias_prot_sim_train_max").iloc[0]
        )
        ligand_rows = chain_df[chain_df["CHAIN_ID"] == "B"].copy()
        assert ligand_rows["bias_lig_sim_train"].astype(float).tolist() == [1.0]

    def test_run_merges_public_and_custom_sources_into_plotting_dataset(self, temp_dir):
        system_path = temp_dir / "system.yaml"
        system_path.write_text(
            yaml.safe_dump(
                {
                    "sequences": [
                        {"protein": {"id": "A", "sequence": "MKRAAT"}},
                        {"ligand": {"id": "B", "smiles": "CCO"}},
                    ]
                }
            ),
            encoding="utf-8",
        )

        structure_ref = temp_dir / "custom_reference.pdb"
        structure_ref.write_text("HEADER CUSTOM\n", encoding="utf-8")
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
            "sequence,dataset_name,source_structure_path\n"
            f"MKRAAT,private_proteins,{structure_ref.name}\n",
            encoding="utf-8",
        )
        custom_ligand.write_text(
            "smiles,dataset_name,source_reference_path\n"
            f"CCO,private_ligands,{structure_ref.name}\n",
            encoding="utf-8",
        )

        bias = Bias(
            wrk_dir=str(temp_dir),
            system_path=str(system_path),
            protein_training_data_path=str(protein_ref),
            ligand_training_data_path=str(ligand_ref),
            custom_protein_reference_path=str(custom_protein),
            custom_ligand_reference_path=str(custom_ligand),
        )

        bias.run()

        output_dir = temp_dir / "results" / "bias_train"
        bias_training = pd.read_csv(output_dir / "bias_training_data.csv")

        assert bias_training["source"].tolist() == ["custom", "public"]
        assert bias_training["query_protein_chain_id"].tolist() == ["A", "A"]
        assert bias_training["query_ligand_chain_id"].tolist() == ["B", "B"]
        assert (output_dir / "bias_reference_overlap_scatter.png").exists()
        assert (output_dir / "bias_reference_overlap_scatter.pdf").exists()
