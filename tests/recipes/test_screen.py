"""Tests for cofolder.recipes.screen module."""
from unittest.mock import patch

import pandas as pd
import pytest
import yaml

from cofolder.recipes.screen import Screen


class TestScreenInit:
    """Tests for Screen initialization."""

    def test_init_basic(self, sample_system_yaml, sample_options_yaml, sample_csv_file, temp_dir):
        screener = Screen(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            variable=["sequences,1,ligand,smiles"],
            variable_csv=str(sample_csv_file),
            col_variable=["smiles"],
            col_id="compound_id",
        )

        assert screener.variable_paths == [["sequences", 1, "ligand", "smiles"]]
        assert screener.col_variable == ["smiles"]

    def test_pair_count_mismatch_raises(self, sample_system_yaml, sample_options_yaml, sample_csv_file, temp_dir):
        with pytest.raises(ValueError, match="Number of --variable entries must match"):
            Screen(
                wrk_dir=str(temp_dir),
                system_path=str(sample_system_yaml),
                options_path=str(sample_options_yaml),
                variable=["sequences,1,ligand,smiles", "sequences,1,ligand,ccd"],
                variable_csv=str(sample_csv_file),
                col_variable=["smiles"],
                col_id="compound_id",
            )

    def test_missing_variable_csv_raises(self, sample_system_yaml, sample_options_yaml, temp_dir):
        with pytest.raises(ValueError, match="--variable_csv is required"):
            Screen(
                wrk_dir=str(temp_dir),
                system_path=str(sample_system_yaml),
                options_path=str(sample_options_yaml),
                variable=["sequences,1,ligand,smiles"],
                variable_csv=None,
                col_variable=["smiles"],
                col_id="compound_id",
            )

    def test_missing_col_id_raises(self, sample_system_yaml, sample_options_yaml, sample_csv_file, temp_dir):
        with pytest.raises(ValueError, match="--col_id is required"):
            Screen(
                wrk_dir=str(temp_dir),
                system_path=str(sample_system_yaml),
                options_path=str(sample_options_yaml),
                variable=["sequences,1,ligand,smiles"],
                variable_csv=str(sample_csv_file),
                col_variable=["smiles"],
                col_id=None,
            )


class TestScreenRun:
    @patch("cofolder.recipes.screen.Validate.run")
    def test_run_invalid_smiles_warns_and_continues(
        self,
        mock_validate_run,
        sample_system_yaml,
        sample_options_yaml,
        temp_dir,
    ):
        csv_path = temp_dir / "invalid_smiles.csv"
        csv_path.write_text(
            "compound_id,smiles\n"
            "CMPD_BAD,C1(\n"
            "CMPD_OK,CCO\n",
            encoding="utf-8",
        )

        screener = Screen(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            variable=["sequences,1,ligand,smiles"],
            variable_csv=str(csv_path),
            col_variable=["smiles"],
            col_id="compound_id",
        )
        screener.run()

        assert mock_validate_run.call_count == 1
        out_df = pd.read_csv(temp_dir / "screen_results.csv")
        assert len(out_df) == 2
        bad = out_df[out_df["compound_id"] == "CMPD_BAD"].iloc[0]
        good = out_df[out_df["compound_id"] == "CMPD_OK"].iloc[0]
        assert bad["status"] == "failed"
        assert "Invalid SMILES" in str(bad["error_message"])
        assert good["status"] == "success"

    @patch("cofolder.recipes.screen.Validate.run")
    def test_run_writes_merged_input_with_scores(
        self,
        mock_validate_run,
        sample_system_yaml,
        sample_options_yaml,
        sample_csv_file,
        temp_dir,
    ):
        input_df = pd.read_csv(sample_csv_file)
        for i, (_, in_row) in enumerate(input_df.iterrows(), 1):
            run_dir = temp_dir / f"{i}_{in_row['compound_id']}"
            results_dir = run_dir / "results"
            results_dir.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(
                [
                    {
                        "repeat": 1,
                        "diffusion_sample": 0,
                        "confidence_score": 0.9,
                    }
                ]
            ).to_csv(results_dir / "system_metrics.csv", index=False)
            pd.DataFrame(
                [
                    {
                        "CHAIN_ID": "A",
                        "ENTITY_TYPE": "protein",
                        "repeat": 1,
                        "diffusion_sample": 0,
                        "sasa_norm_heavy": 5.1,
                    },
                    {
                        "CHAIN_ID": "B",
                        "ENTITY_TYPE": "ligand",
                        "repeat": 1,
                        "diffusion_sample": 0,
                        "affinity_pred_value": 6.2,
                    },
                ]
            ).to_csv(results_dir / "chain_metrics.csv", index=False)

        screener = Screen(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            variable=["sequences,1,ligand,smiles"],
            variable_csv=str(sample_csv_file),
            col_variable=["smiles"],
            col_id="compound_id",
            merge_data="mw",
        )
        screener.run()

        merged_df = pd.read_csv(temp_dir / "screen_results_with_scores.csv")
        assert "system__confidence_score" in merged_df.columns
        assert "protein_A__sasa_norm_heavy" in merged_df.columns
        assert "ligand_B__affinity_pred_value" in merged_df.columns
        assert "mw" in merged_df.columns
        assert "compound_id" in merged_df.columns

    @patch("cofolder.recipes.screen.Validate.run")
    def test_run_calls_validate_per_row_and_writes_summary(
        self,
        mock_validate_run,
        sample_system_yaml,
        sample_options_yaml,
        sample_csv_file,
        temp_dir,
    ):
        screener = Screen(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            variable=["sequences,1,ligand,smiles"],
            variable_csv=str(sample_csv_file),
            col_variable=["smiles"],
            col_id="compound_id",
            merge_data="mw",
        )

        screener.run()

        assert mock_validate_run.call_count == 2

        out_csv = temp_dir / "screen_results.csv"
        assert out_csv.exists()
        out_df = pd.read_csv(out_csv)
        assert len(out_df) == 2
        assert set(["index", "compound_id", "status", "error_message", "run_dir", "smiles", "mw"]).issubset(
            out_df.columns
        )
        assert set(out_df["status"].tolist()) == {"success"}

        merged_csv = temp_dir / "screen_results_with_scores.csv"
        assert merged_csv.exists()
        merged_df = pd.read_csv(merged_csv)
        assert len(merged_df) == 2
        assert set(["compound_id", "smiles", "status", "error_message", "run_dir"]).issubset(merged_df.columns)

    @patch("cofolder.recipes.screen.Validate.run")
    def test_run_multi_variable_updates_system_yaml(
        self,
        mock_validate_run,
        sample_system_yaml,
        sample_options_yaml,
        temp_dir,
    ):
        csv_path = temp_dir / "multi.csv"
        csv_path.write_text(
            "compound_id,smiles,ccd\n"
            "CMPD001,CCO,EDO\n",
            encoding="utf-8",
        )

        screener = Screen(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            variable=["sequences,1,ligand,smiles", "sequences,1,ligand,ccd"],
            variable_csv=str(csv_path),
            col_variable=["smiles", "ccd"],
            col_id="compound_id",
        )

        screener.run()

        run_system_yaml = temp_dir / "1_CMPD001" / "screen_system.yaml"
        assert run_system_yaml.exists()
        data = yaml.safe_load(run_system_yaml.read_text(encoding="utf-8"))
        lig = data["sequences"][1]["ligand"]
        assert lig["smiles"] == "CCO"
        assert lig["ccd"] == "EDO"

    @patch("cofolder.recipes.screen.Validate.run")
    def test_run_continue_on_failure(
        self,
        mock_validate_run,
        sample_system_yaml,
        sample_options_yaml,
        sample_csv_file,
        temp_dir,
    ):
        def _run_side_effect():
            if mock_validate_run.call_count == 1:
                raise RuntimeError("boom")
            return None

        mock_validate_run.side_effect = _run_side_effect

        screener = Screen(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            variable=["sequences,1,ligand,smiles"],
            variable_csv=str(sample_csv_file),
            col_variable=["smiles"],
            col_id="compound_id",
        )

        screener.run()

        out_df = pd.read_csv(temp_dir / "screen_results.csv")
        assert len(out_df) == 2
        assert "failed" in set(out_df["status"].tolist())
        assert "success" in set(out_df["status"].tolist())

        merged_df = pd.read_csv(temp_dir / "screen_results_with_scores.csv")
        assert len(merged_df) == 2
        assert "failed" in set(merged_df["status"].tolist())
        assert "success" in set(merged_df["status"].tolist())
