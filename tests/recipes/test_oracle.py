"""Tests for cofolder.recipes.oracle module."""

from unittest.mock import MagicMock
from unittest.mock import patch

import pandas as pd
import pytest

from cofolder.recipes.oracle import Oracle
from cofolder.recipes.validate import DEFAULT_SCORING_FUNCTIONS


class TestOracleInit:
    def test_default_scoring_functions_remain_enabled(
        self,
        sample_system_yaml,
        sample_options_yaml,
        temp_dir,
    ):
        oracle = Oracle(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            input_smiles="CCO",
            output_metric="affinity_pred_value",
        )

        assert set(oracle.validate_kwargs["scoring_functions"]) == DEFAULT_SCORING_FUNCTIONS

    def test_requires_exactly_one_input(self, sample_system_yaml, sample_options_yaml, temp_dir):
        with pytest.raises(ValueError, match="exactly one input"):
            Oracle(
                wrk_dir=str(temp_dir),
                system_path=str(sample_system_yaml),
                options_path=str(sample_options_yaml),
                input_smiles=None,
                input_mol_file=None,
                output_metric="affinity_pred_value",
            )

    def test_metric_requires_enabled_scoring_functions(self, sample_system_yaml, sample_options_yaml, temp_dir):
        with pytest.raises(ValueError, match="requires one of scoring functions"):
            Oracle(
                wrk_dir=str(temp_dir),
                system_path=str(sample_system_yaml),
                options_path=str(sample_options_yaml),
                input_smiles="CCO",
                output_metric="ifp_distance",
                scoring_functions=["confidence_metrics"],
            )


class TestOracleRun:
    @patch("cofolder.recipes.oracle.Validate")
    def test_run_uses_boltz2_as_default_runner(
        self,
        mock_validate_cls,
        sample_system_yaml,
        sample_options_yaml,
        temp_dir,
    ):
        mock_validator = MagicMock()
        mock_validate_cls.return_value = mock_validator

        oracle = Oracle(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            input_smiles="CCO",
            output_metric="affinity_pred_value",
            scoring_functions=["affinity_metrics"],
        )

        run_dir = temp_dir / "oracle_run" / "results"
        run_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame([{"affinity_pred_value": 6.4}]).to_csv(run_dir / "chain_metrics.csv", index=False)

        oracle.run()

        assert mock_validate_cls.call_args.kwargs["runner"] == "boltz2"
        mock_validator.run.assert_called_once()

    @patch("cofolder.recipes.oracle.Validate")
    def test_run_passes_explicit_boltz_community_runner(
        self,
        mock_validate_cls,
        sample_system_yaml,
        sample_options_yaml,
        temp_dir,
    ):
        mock_validator = MagicMock()
        mock_validate_cls.return_value = mock_validator

        oracle = Oracle(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            runner="boltz-community",
            input_smiles="CCO",
            output_metric="affinity_pred_value",
            scoring_functions=["affinity_metrics"],
        )

        run_dir = temp_dir / "oracle_run" / "results"
        run_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame([{"affinity_pred_value": 6.4}]).to_csv(run_dir / "chain_metrics.csv", index=False)

        oracle.run()

        assert mock_validate_cls.call_args.kwargs["runner"] == "boltz-community"
        mock_validator.run.assert_called_once()

    @patch("cofolder.recipes.oracle.Validate.run")
    def test_run_returns_single_value_and_writes_output(
        self,
        mock_validate_run,
        sample_system_yaml,
        sample_options_yaml,
        temp_dir,
    ):
        oracle = Oracle(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            input_smiles="CCO",
            output_metric="affinity_pred_value",
            aggregate="first",
            scoring_functions=["affinity_metrics"],
        )

        run_dir = temp_dir / "oracle_run" / "results"
        run_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame([{"affinity_pred_value": 6.4}]).to_csv(run_dir / "chain_metrics.csv", index=False)
        pd.DataFrame([{"ptm": 0.5}]).to_csv(run_dir / "system_metrics.csv", index=False)

        value = oracle.run()

        assert mock_validate_run.call_count == 1
        assert value == pytest.approx(6.4)
        result_csv = temp_dir / "oracle_result.csv"
        assert result_csv.exists()
        out = pd.read_csv(result_csv)
        assert out.iloc[0]["output_metric"] == "affinity_pred_value"
        assert float(out.iloc[0]["value"]) == pytest.approx(6.4)

    @patch("cofolder.recipes.oracle.Validate.run")
    def test_run_aggregate_mean(
        self,
        mock_validate_run,
        sample_system_yaml,
        sample_options_yaml,
        temp_dir,
    ):
        oracle = Oracle(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            input_smiles="CCO",
            output_metric="sasa_norm_heavy",
            aggregate="mean",
            scoring_functions=["sasa_normalized"],
        )

        run_dir = temp_dir / "oracle_run" / "results"
        run_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(
            [
                {"sasa_norm_heavy": 5.0},
                {"sasa_norm_heavy": 7.0},
            ]
        ).to_csv(run_dir / "chain_metrics.csv", index=False)

        value = oracle.run()
        assert mock_validate_run.call_count == 1
        assert value == pytest.approx(6.0)
