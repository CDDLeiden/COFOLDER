"""Tests for cofolder.recipes.oracle module."""

import json
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from cofolder.recipes.oracle import (
    Oracle,
    OracleGate,
    OracleGatePolicy,
    OracleScoreContext,
)
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

        assert (
            set(oracle.validate_kwargs["scoring_functions"])
            == DEFAULT_SCORING_FUNCTIONS
        )

    def test_requires_exactly_one_input(
        self, sample_system_yaml, sample_options_yaml, temp_dir
    ):
        with pytest.raises(ValueError, match="exactly one input"):
            Oracle(
                wrk_dir=str(temp_dir),
                system_path=str(sample_system_yaml),
                options_path=str(sample_options_yaml),
                input_smiles=None,
                input_mol_file=None,
                output_metric="affinity_pred_value",
            )

    def test_metric_requires_enabled_scoring_functions(
        self, sample_system_yaml, sample_options_yaml, temp_dir
    ):
        with pytest.raises(ValueError, match="requires one of scoring functions"):
            Oracle(
                wrk_dir=str(temp_dir),
                system_path=str(sample_system_yaml),
                options_path=str(sample_options_yaml),
                input_smiles="CCO",
                output_metric="ifp_distance",
                scoring_functions=["confidence_metrics"],
            )

    @pytest.mark.parametrize(
        "kwargs",
        [
            {},
            {
                "output_metric": "confidence_score",
                "score_components": {"confidence_score": 1},
            },
            {
                "score_components": {"confidence_score": 1},
                "scoring_function": lambda context: 1,
            },
        ],
    )
    def test_requires_exactly_one_score_source(
        self, kwargs, sample_system_yaml, sample_options_yaml, temp_dir
    ):
        with pytest.raises(ValueError, match="exactly one score source"):
            Oracle(
                wrk_dir=str(temp_dir),
                system_path=str(sample_system_yaml),
                options_path=str(sample_options_yaml),
                input_smiles="CCO",
                **kwargs,
            )

    def test_gates_require_explicit_policy(
        self, sample_system_yaml, sample_options_yaml, temp_dir
    ):
        with pytest.raises(ValueError, match="gate_policy is required"):
            Oracle(
                wrk_dir=str(temp_dir),
                system_path=str(sample_system_yaml),
                options_path=str(sample_options_yaml),
                input_smiles="CCO",
                output_metric="confidence_score",
                score_gates=[OracleGate("confidence_score", "ge", 0.5)],
            )

    def test_custom_pocket_coverage_dependencies_are_checked(
        self, sample_system_yaml, sample_options_yaml, temp_dir
    ):
        with pytest.raises(ValueError, match="pocket_coverage_reference"):
            Oracle(
                wrk_dir=str(temp_dir),
                system_path=str(sample_system_yaml),
                options_path=str(sample_options_yaml),
                input_smiles="CCO",
                output_metric="ligand_B__pocket_coverage_custom",
                scoring_functions=["ifp_distance"],
            )

    def test_bias_metric_requires_bias_assessment(
        self, sample_system_yaml, sample_options_yaml, temp_dir
    ):
        with pytest.raises(ValueError, match="assess_bias=True"):
            Oracle(
                wrk_dir=str(temp_dir),
                system_path=str(sample_system_yaml),
                options_path=str(sample_options_yaml),
                input_smiles="CCO",
                output_metric="ligand_B__bias_lig_sim_train",
            )

    def test_reference_structural_metric_requires_reference(
        self, sample_system_yaml, sample_options_yaml, temp_dir
    ):
        with pytest.raises(ValueError, match="requires reference_path"):
            Oracle(
                wrk_dir=str(temp_dir),
                system_path=str(sample_system_yaml),
                options_path=str(sample_options_yaml),
                input_smiles="CCO",
                output_metric="ligand_B__ligand_rmsd_ref",
            )

    def test_vector_ifp_cannot_be_a_scalar_objective(
        self, sample_system_yaml, sample_options_yaml, temp_dir
    ):
        with pytest.raises(ValueError, match="vector-valued"):
            Oracle(
                wrk_dir=str(temp_dir),
                system_path=str(sample_system_yaml),
                options_path=str(sample_options_yaml),
                input_smiles="CCO",
                output_metric="ifp_distance",
                scoring_functions=["ifp_distance"],
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
        pd.DataFrame([{"affinity_pred_value": 6.4}]).to_csv(
            run_dir / "chain_metrics.csv", index=False
        )

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
        pd.DataFrame([{"affinity_pred_value": 6.4}]).to_csv(
            run_dir / "chain_metrics.csv", index=False
        )

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
        pd.DataFrame([{"affinity_pred_value": 6.4}]).to_csv(
            run_dir / "chain_metrics.csv", index=False
        )
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


def _write_current_metric_outputs(temp_dir, *, ligand=None, system=None, protein=None):
    results = temp_dir / "oracle_run" / "results"
    results.mkdir(parents=True, exist_ok=True)
    if system is not None:
        pd.DataFrame(system).to_csv(results / "system_metrics.csv", index=False)
    rows = []
    for values in protein or []:
        rows.append({"CHAIN_ID": "A", "ENTITY_TYPE": "protein", **values})
    for values in ligand or []:
        rows.append({"CHAIN_ID": "B", "ENTITY_TYPE": "ligand", **values})
    if rows:
        pd.DataFrame(rows).to_csv(results / "chain_metrics.csv", index=False)


class TestOracleScoring:
    @patch("cofolder.recipes.oracle.Validate.run")
    def test_bare_selector_prefers_query_ligand_and_qualified_selector_is_supported(
        self, mock_run, sample_system_yaml, sample_options_yaml, temp_dir
    ):
        _write_current_metric_outputs(
            temp_dir,
            protein=[{"sasa": 90.0}],
            ligand=[{"sasa": 12.0}],
            system=[{"confidence_score": 0.8}],
        )
        oracle = Oracle(
            str(temp_dir),
            str(sample_system_yaml),
            str(sample_options_yaml),
            input_smiles="CCO",
            output_metric="sasa",
            scoring_functions=["sasa"],
        )
        assert oracle.run() == pytest.approx(12.0)
        assert oracle._metric_value(
            "system__confidence_score",
            {"system__confidence_score": (0.8,)},
        ) == pytest.approx(0.8)

    @pytest.mark.parametrize(
        ("aggregate", "expected"),
        [("first", 4), ("mean", 5.5), ("max", 7), ("min", 4), ("median", 5.5)],
    )
    def test_all_aggregation_methods(
        self, aggregate, expected, sample_system_yaml, sample_options_yaml, temp_dir
    ):
        oracle = Oracle(
            str(temp_dir),
            str(sample_system_yaml),
            str(sample_options_yaml),
            input_smiles="CCO",
            output_metric="confidence_score",
            aggregate=aggregate,
            scoring_functions=["confidence_metrics"],
        )
        assert oracle._aggregate_values("system__confidence_score", [4, 7]) == expected

    def test_ambiguous_bare_selector_requires_qualification(
        self, sample_system_yaml, sample_options_yaml, temp_dir
    ):
        oracle = Oracle(
            str(temp_dir),
            str(sample_system_yaml),
            str(sample_options_yaml),
            input_smiles="CCO",
            output_metric="mystery",
        )
        with pytest.raises(ValueError, match="ambiguous"):
            oracle._resolve_selector(
                "mystery",
                {
                    "protein_A__mystery": (1.0,),
                    "protein_C__mystery": (2.0,),
                },
            )

    @patch("cofolder.recipes.oracle.Validate.run")
    def test_nonfinite_metric_is_rejected(
        self, mock_run, sample_system_yaml, sample_options_yaml, temp_dir
    ):
        _write_current_metric_outputs(
            temp_dir, system=[{"confidence_score": float("nan")}]
        )
        oracle = Oracle(
            str(temp_dir),
            str(sample_system_yaml),
            str(sample_options_yaml),
            input_smiles="CCO",
            output_metric="confidence_score",
            scoring_functions=["confidence_metrics"],
        )
        with pytest.raises(ValueError, match="finite numeric values"):
            oracle.run()

    @patch("cofolder.recipes.oracle.Validate.run")
    def test_weighted_composite_supports_negative_weights_and_writes_audit(
        self, mock_run, sample_system_yaml, sample_options_yaml, temp_dir
    ):
        _write_current_metric_outputs(
            temp_dir,
            ligand=[{"affinity_pred_value": 6.0}],
            system=[{"confidence_score": 0.8}],
        )
        oracle = Oracle(
            str(temp_dir),
            str(sample_system_yaml),
            str(sample_options_yaml),
            input_smiles="CCO",
            score_components={
                "ligand_B__affinity_pred_value": -1.0,
                "system__confidence_score": 2.0,
            },
            scoring_functions=["affinity_metrics", "confidence_metrics"],
        )
        assert oracle.run() == pytest.approx(-4.4)
        audit = pd.read_csv(temp_dir / "oracle_result.csv").iloc[0]
        assert audit["score_mode"] == "composite"
        assert audit["base_value"] == pytest.approx(-4.4)
        assert json.loads(audit["component_values"]) == {
            "ligand_B__affinity_pred_value": 6.0,
            "system__confidence_score": 0.8,
        }
        assert bool(audit["gate_pass"])
        assert audit["gate_action"] == "none"

    @patch("cofolder.recipes.oracle.Validate.run")
    def test_custom_function_receives_structured_context_once(
        self, mock_run, sample_system_yaml, sample_options_yaml, temp_dir
    ):
        _write_current_metric_outputs(
            temp_dir,
            ligand=[{"affinity_pred_value": 6.0}],
            system=[{"confidence_score": 0.75}],
        )
        scoring_function = MagicMock(
            side_effect=lambda context: context.aggregated_metrics[
                "system__confidence_score"
            ]
            * 10
        )
        oracle = Oracle(
            str(temp_dir),
            str(sample_system_yaml),
            str(sample_options_yaml),
            input_smiles="CCO",
            scoring_function=scoring_function,
            scoring_functions=["confidence_metrics"],
        )
        assert oracle.run() == pytest.approx(7.5)
        scoring_function.assert_called_once()
        context = scoring_function.call_args.args[0]
        assert isinstance(context, OracleScoreContext)
        assert context.query_smiles == "CCO"
        assert context.run_dir == temp_dir / "oracle_run"
        assert not context.system_metrics.empty
        assert not context.chain_metrics.empty

    @pytest.mark.parametrize(
        "invalid", [True, "not-a-score", float("nan"), float("inf")]
    )
    @patch("cofolder.recipes.oracle.Validate.run")
    def test_custom_function_must_return_finite_number(
        self, mock_run, invalid, sample_system_yaml, sample_options_yaml, temp_dir
    ):
        _write_current_metric_outputs(temp_dir, system=[{"confidence_score": 0.75}])
        oracle = Oracle(
            str(temp_dir),
            str(sample_system_yaml),
            str(sample_options_yaml),
            input_smiles="CCO",
            scoring_function=lambda context: invalid,
        )
        with pytest.raises((TypeError, ValueError), match="finite numeric scalar"):
            oracle.run()

    @pytest.mark.parametrize(
        ("comparison", "value", "passes"),
        [("ge", 0.5, True), ("gt", 0.5, False), ("le", 0.5, True), ("lt", 0.5, False)],
    )
    @patch("cofolder.recipes.oracle.Validate.run")
    def test_gate_comparison_boundaries(
        self,
        mock_run,
        comparison,
        value,
        passes,
        sample_system_yaml,
        sample_options_yaml,
        temp_dir,
    ):
        _write_current_metric_outputs(temp_dir, system=[{"confidence_score": value}])
        oracle = Oracle(
            str(temp_dir),
            str(sample_system_yaml),
            str(sample_options_yaml),
            input_smiles="CCO",
            output_metric="confidence_score",
            score_gates=[OracleGate("system__confidence_score", comparison, 0.5)],
            gate_policy=OracleGatePolicy("fixed_penalty", -2),
            scoring_functions=["confidence_metrics"],
        )
        expected = value if passes else -2
        assert oracle.run() == pytest.approx(expected)

    @pytest.mark.parametrize(
        ("mode", "policy_value", "expected"),
        [("downweight", 0.25, 1.5), ("fixed_penalty", -3, -3), ("non_binder", 0, 0)],
    )
    @patch("cofolder.recipes.oracle.Validate.run")
    def test_gate_failure_policies_are_applied_once(
        self,
        mock_run,
        mode,
        policy_value,
        expected,
        sample_system_yaml,
        sample_options_yaml,
        temp_dir,
    ):
        _write_current_metric_outputs(temp_dir, ligand=[{"affinity_pred_value": 6.0}])
        oracle = Oracle(
            str(temp_dir),
            str(sample_system_yaml),
            str(sample_options_yaml),
            input_smiles="CCO",
            output_metric="affinity_pred_value",
            score_gates=[
                OracleGate("missing_structural_metric", "ge", 0.5),
                OracleGate("another_missing_metric", "le", 2.0),
            ],
            gate_policy=OracleGatePolicy(mode, policy_value),
            scoring_functions=["affinity_metrics"],
        )
        assert oracle.run() == pytest.approx(expected)
        audit = pd.read_csv(temp_dir / "oracle_result.csv").iloc[0]
        assert not bool(audit["gate_pass"])
        assert json.loads(audit["failed_gates"]) == [
            "missing_structural_metric",
            "another_missing_metric",
        ]
        assert audit["gate_action"] == mode

    @patch("cofolder.recipes.oracle.Validate.run")
    def test_bias_similarity_scalar_regression(
        self, mock_run, sample_system_yaml, sample_options_yaml, temp_dir
    ):
        _write_current_metric_outputs(temp_dir, ligand=[{"bias_lig_sim_train": 0.42}])
        oracle = Oracle(
            str(temp_dir),
            str(sample_system_yaml),
            str(sample_options_yaml),
            input_smiles="CCO",
            output_metric="ligand_B__bias_lig_sim_train",
            assess_bias=True,
        )
        assert oracle.run() == pytest.approx(0.42)

    @patch("cofolder.recipes.oracle.Validate.run")
    def test_custom_pocket_coverage_scalar_regression(
        self, mock_run, sample_system_yaml, sample_options_yaml, temp_dir
    ):
        _write_current_metric_outputs(
            temp_dir, ligand=[{"pocket_coverage_custom": 0.75}]
        )
        oracle = Oracle(
            str(temp_dir),
            str(sample_system_yaml),
            str(sample_options_yaml),
            input_smiles="CCO",
            output_metric="ligand_B__pocket_coverage_custom",
            scoring_functions=["ifp_distance"],
            pocket_coverage_reference="A2 S8 T10",
            reproduction_metrics=["pocket_coverage"],
        )
        assert oracle.run() == pytest.approx(0.75)
