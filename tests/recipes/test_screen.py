"""Tests for cofolder.recipes.screen module."""

import json
import logging
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest
import yaml
from rdkit import Chem

from cofolder.modules.contracts import (
    ScreenExecutionCardinalityError,
    WorkflowExecutionError,
)
from cofolder.modules.input.system import System
from cofolder.modules.runners.boltz2_runner import Boltz2Runner
from cofolder.modules.runners.contracts import (
    RunnerExecutionResult,
    RunnerMetricOutcome,
    RunnerPreparationResult,
    RunnerRuntime,
)
from cofolder.modules.runners.msa import capture_generated_msas, inject_cached_msas
from cofolder.recipes.screen import Screen
from cofolder.recipes._screen_postprocess import ScreenPostprocessor
from tests.modules.analytics.test_reproduction import (
    _write_predicted_pdb,
    _write_reference_pdb,
)


def test_execution_key_validation_ignores_aggregate_metric_rows():
    processor = SimpleNamespace(
        execution_plan=SimpleNamespace(
            executions=(SimpleNamespace(repeat_id=1, sample_id=0),)
        )
    )
    system_df = pd.DataFrame(
        [
            {"repeat": 1, "diffusion_sample": 0, "confidence_score": 0.8},
            {"repeat": None, "diffusion_sample": None, "confidence_score": 0.8},
        ]
    )

    ScreenPostprocessor._validate_observed_execution_keys(processor, system_df)


def test_reference_complex_ifp_filter_annotates_primary_prediction(
    temp_dir, sample_system_yaml, sample_options_yaml, sample_csv_file
):
    reference_path = temp_dir / "reference.pdb"
    _write_reference_pdb(reference_path)
    run_dir = temp_dir / "compound"
    structures = run_dir / "results" / "structures"
    structures.mkdir(parents=True)
    predicted_name = "prediction.pdb"
    _write_predicted_pdb(
        structures / predicted_name,
        ((1.3, 1.2, 0.0), (2.7, 1.2, 0.0)),
    )
    pd.DataFrame(
        [
            {"CHAIN_ID": "A", "ENTITY_TYPE": "protein", "cif_file": predicted_name},
            {"CHAIN_ID": "Z", "ENTITY_TYPE": "ligand", "cif_file": predicted_name},
        ]
    ).to_csv(run_dir / "results" / "chain_metrics.csv", index=False)

    screen = Screen(
        wrk_dir=str(temp_dir / "screen"),
        system_path=str(sample_system_yaml),
        options_path=str(sample_options_yaml),
        library=str(sample_csv_file),
        ligand_chain="Z",
        smiles_column="smiles",
        col_id="compound_id",
        reference_path=str(reference_path),
        ifp_filter_threshold=1.0,
        ifp_filter_source="reference_complex",
        ifp_reference_ligand="L",
    )
    screen._validate_config()
    screen._prepare_reference_ifp()
    result = screen._evaluate_ifp_filter(run_dir)

    assert result["ifp_filter_status"] == "accepted"
    assert result["ifp_filter_similarity"] == 1.0
    assert result["ifp_filter_mapping_status"] == "mapped"
    assert json.loads(result["ifp_filter_missing_interactions"]) == []


class _ReusableScreenRunner:
    name = "boltz2"
    capabilities = {"confidence_metrics"}
    supports_msa_reuse = True

    def __init__(self, *, generate_msa=True, fail_after_msa_once=False):
        self.msa_missing_at_run: list[bool] = []
        self.msa_sequences_missing_at_run: list[tuple[str, ...]] = []
        self.generate_msa = generate_msa
        self.fail_after_msa_once = fail_after_msa_once

    def ensure_available(self):
        return None

    def validate_system(
        self,
        system_obj,
        options_obj,
        *,
        check_atom_names=True,
        source_path=None,
        requirements=None,
    ):
        return Boltz2Runner().validate_system(
            system_obj,
            options_obj,
            check_atom_names=check_atom_names,
            source_path=source_path,
            requirements=requirements,
        )

    def load_options(self, options_path):
        return {}

    def msa_reuse_settings(self, options_obj):
        return {}

    def prepare_system(
        self, system_obj, options_obj, wrk_dir, conformers, sdf_file, logger
    ):
        return RunnerPreparationResult(system_obj=system_obj, options_obj=options_obj)

    def inject_reusable_msas(self, system_obj, cache_dir, *, settings=None):
        return inject_cached_msas(system_obj, cache_dir, settings=settings)

    def capture_reusable_msas(
        self, system_obj, *, generated_dir, cache_dir, settings=None
    ):
        return capture_generated_msas(
            system_obj,
            generated_dir=generated_dir,
            cache_dir=cache_dir,
            runner_name=self.name,
            settings=settings,
        )

    def run(self, request):
        proteins = [
            entry["protein"]
            for entry in request.system_obj.system["sequences"]
            if "protein" in entry
        ]
        missing_proteins = [protein for protein in proteins if not protein.get("msa")]
        self.msa_missing_at_run.append(bool(missing_proteins))
        self.msa_sequences_missing_at_run.append(
            tuple(protein["sequence"] for protein in missing_proteins)
        )
        if missing_proteins and self.generate_msa:
            msa_dir = (
                request.repeat_dir / f"boltz_results_{request.system_name}" / "msa"
            )
            msa_dir.mkdir(parents=True, exist_ok=True)
            for index, protein in enumerate(missing_proteins):
                (msa_dir / f"generated_{index}.csv").write_text(
                    f"key,sequence\n-1,{protein['sequence']}\n",
                    encoding="utf-8",
                )
        if self.fail_after_msa_once:
            self.fail_after_msa_once = False
            raise RuntimeError("prediction failed after MSA generation")

        normalized_dir = request.repeat_dir / "normalized"
        structures_dir = normalized_dir / "structures"
        structures_dir.mkdir(parents=True, exist_ok=True)
        structure_name = f"{request.repeat}_{request.system_name}_model_0.cif"
        (structures_dir / structure_name).write_text("data", encoding="utf-8")
        system_metrics_path = normalized_dir / "system_metrics.csv"
        chain_metrics_path = normalized_dir / "chain_metrics.csv"
        manifest_path = normalized_dir / "manifest.json"
        pd.DataFrame(
            [
                {
                    "cif_file": structure_name,
                    "model_name": request.system_name,
                    "repeat": request.repeat,
                    "diffusion_sample": 0,
                    "ptm": 0.8,
                    "iptm": 0.7,
                    "confidence_score": 0.9,
                }
            ]
        ).to_csv(system_metrics_path, index=False)
        pd.DataFrame(
            [
                {
                    "conf_chain_id": 0,
                    "cif_file": structure_name,
                    "model_name": request.system_name,
                    "repeat": request.repeat,
                    "diffusion_sample": 0,
                    "chains_ptm": 0.8,
                    "pair_chains_iptm_1": 0.7,
                },
                {
                    "conf_chain_id": 1,
                    "cif_file": structure_name,
                    "model_name": request.system_name,
                    "repeat": request.repeat,
                    "diffusion_sample": 0,
                    "chains_ptm": 0.7,
                    "pair_chains_iptm_0": 0.7,
                },
            ]
        ).to_csv(chain_metrics_path, index=False)
        manifest_path.write_text(
            json.dumps({"runner": self.name, "capabilities": ["confidence_metrics"]}),
            encoding="utf-8",
        )
        return RunnerExecutionResult(
            runner_name=self.name,
            raw_output_dir=request.repeat_dir,
            normalized_dir=normalized_dir,
            structures_dir=structures_dir,
            system_metrics_path=system_metrics_path,
            chain_metrics_path=chain_metrics_path,
            manifest_path=manifest_path,
            capabilities=set(self.capabilities),
            runtime=RunnerRuntime(diffusion_samples=1),
            metric_outcomes={
                "confidence_metrics": RunnerMetricOutcome(state="computed")
            },
        )


class TestScreenInit:
    """Tests for Screen initialization."""

    def test_init_basic(
        self, sample_system_yaml, sample_options_yaml, sample_csv_file, temp_dir
    ):
        screener = Screen(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            ligand_chain="B",
            library=str(sample_csv_file),
            smiles_column="smiles",
            col_id="compound_id",
            bias_query_cache_path=str(temp_dir / "shared-bias-cache"),
        )

        assert screener.ligand_chain == "B"
        assert screener.smiles_column == "smiles"
        assert screener.smiles_column == "smiles"

    def test_unique_ligand_selector_is_inferred_during_preflight(
        self, sample_system_yaml, sample_options_yaml, sample_csv_file, temp_dir
    ):
        work_dir = temp_dir / "preflight-work"
        with patch(
            "cofolder.modules.runners.boltz2_runner.Boltz2Runner.check_availability",
            return_value=(True, None),
        ):
            screen = Screen(
                wrk_dir=str(temp_dir),
                system_path=str(sample_system_yaml),
                options_path=str(sample_options_yaml),
                library=str(sample_csv_file),
                smiles_column="smiles",
                col_id="compound_id",
            )
            screen.wrk_dir = work_dir
            report = screen.preflight()

        assert report.ready
        assert screen.ligand_chain == "B"
        assert not work_dir.exists()

    def test_mapped_preflight_validates_every_row_without_writing(
        self, sample_system_yaml, sample_options_yaml, temp_dir
    ):
        csv_path = temp_dir / "variants.csv"
        csv_path.write_text(
            "id,sequence\nvalid,MKRAAC\ninvalid,MKRAA*\n", encoding="utf-8"
        )
        work_dir = temp_dir / "mapped-preflight"
        with patch(
            "cofolder.modules.runners.boltz2_runner.Boltz2Runner.check_availability",
            return_value=(True, None),
        ):
            report = Screen(
                wrk_dir=str(work_dir),
                system_path=str(sample_system_yaml),
                options_path=str(sample_options_yaml),
                library=str(csv_path),
                col_id="id",
                mappings=["sequence=sequences.0.protein.sequence"],
            ).preflight()

        assert not report.ready
        assert "mapped_rows=2" in report.messages
        assert "mapped_rows_valid=1" in report.messages
        assert "mapped_rows_invalid=1" in report.messages
        assert any("row 2" in message and "Invalid" in message for message in report.messages)
        assert not work_dir.exists()

    def test_missing_library_raises(
        self, sample_system_yaml, sample_options_yaml, temp_dir
    ):
        with pytest.raises(WorkflowExecutionError, match="--library is required"):
            Screen(
                wrk_dir=str(temp_dir),
                system_path=str(sample_system_yaml),
                options_path=str(sample_options_yaml),
                ligand_chain="B",
                library=None,
                smiles_column="smiles",
                col_id="compound_id",
            ).run()

    def test_missing_col_id_raises(
        self, sample_system_yaml, sample_options_yaml, sample_csv_file, temp_dir
    ):
        with pytest.raises(WorkflowExecutionError, match="--col_id is required"):
            Screen(
                wrk_dir=str(temp_dir),
                system_path=str(sample_system_yaml),
                options_path=str(sample_options_yaml),
                ligand_chain="B",
                library=str(sample_csv_file),
                smiles_column="smiles",
                col_id=None,
            ).run()


class TestScreenRun:
    @patch("cofolder.recipes.screen.Validate.run")
    def test_mapped_screen_accepts_protein_only_system(
        self, mock_validate_run, sample_options_yaml, temp_dir
    ):
        system_path = temp_dir / "protein.yaml"
        system_path.write_text(
            yaml.safe_dump(
                {"sequences": [{"protein": {"id": "A", "sequence": "MKRAAT"}}]}
            ),
            encoding="utf-8",
        )
        csv_path = temp_dir / "proteins.csv"
        csv_path.write_text("id,sequence\nvariant,MKRAAC\n", encoding="utf-8")

        result = Screen(
            wrk_dir=str(temp_dir / "screen"),
            system_path=str(system_path),
            options_path=str(sample_options_yaml),
            library=str(csv_path),
            col_id="id",
            mappings=["sequence=sequences.0.protein.sequence"],
        ).run()

        assert result["status"].tolist() == ["success"]
        row_system = yaml.safe_load(
            (temp_dir / "screen" / "compound_000001" / "screen_system.yaml").read_text(
                encoding="utf-8"
            )
        )
        assert row_system["sequences"][0]["protein"]["sequence"] == "MKRAAC"

    @patch("cofolder.recipes.screen.Validate.run")
    def test_csv_and_sdf_publish_same_execution_schema(
        self,
        mock_validate_run,
        sample_system_yaml,
        sample_options_yaml,
        temp_dir,
    ):
        csv_path = temp_dir / "one.csv"
        csv_path.write_text("id,smiles\none,CCO\n", encoding="utf-8")
        molecule = Chem.MolFromSmiles("CCO")
        molecule.SetProp("_Name", "one")
        sdf_path = temp_dir / "one.sdf"
        sdf_path.write_text(
            Chem.MolToMolBlock(molecule) + "\n$$$$\n", encoding="utf-8"
        )

        csv_work = temp_dir / "csv_screen"
        sdf_work = temp_dir / "sdf_screen"
        Screen(
            wrk_dir=str(csv_work),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            ligand_chain="B",
            library=str(csv_path),
            smiles_column="smiles",
            col_id="id",
        ).run()
        Screen(
            wrk_dir=str(sdf_work),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            ligand_chain="B",
            library=str(sdf_path),
        ).run()

        csv_executions = pd.read_csv(csv_work / "results" / "executions.csv")
        sdf_executions = pd.read_csv(sdf_work / "results" / "executions.csv")
        assert csv_executions.columns.tolist() == sdf_executions.columns.tolist()
        assert csv_executions["status"].tolist() == ["success"]
        assert sdf_executions["status"].tolist() == ["success"]
        assert mock_validate_run.call_count == 2

    @patch("cofolder.recipes.screen.Validate.run", autospec=True)
    def test_duplicate_normalized_execution_raises_cardinality_error(
        self,
        mock_validate_run,
        sample_system_yaml,
        sample_options_yaml,
        sample_csv_file,
        temp_dir,
    ):
        def write_duplicates(validator):
            results_dir = validator.wrk_dir / "results"
            results_dir.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(
                [
                    {
                        "model_name": "boltz2",
                        "repeat": 1,
                        "diffusion_sample": 0,
                        "confidence_score": value,
                    }
                    for value in (0.8, 0.9)
                ]
            ).to_csv(results_dir / "system_metrics.csv", index=False)

        mock_validate_run.side_effect = write_duplicates
        with pytest.raises(ScreenExecutionCardinalityError, match="duplicate"):
            Screen(
                wrk_dir=str(temp_dir / "duplicate_screen"),
                system_path=str(sample_system_yaml),
                options_path=str(sample_options_yaml),
                ligand_chain="B",
                library=str(sample_csv_file),
                smiles_column="smiles",
                col_id="compound_id",
            ).run()

    @patch("cofolder.recipes.screen.Validate.run", autospec=True)
    def test_execution_matrix_expands_compounds_repeats_and_samples(
        self,
        mock_validate_run,
        sample_system_yaml,
        sample_csv_file,
        temp_dir,
    ):
        options_path = temp_dir / "matrix_options.yaml"
        options_path.write_text(
            yaml.safe_dump(
                {
                    "version": 1,
                    "runtime": {
                        "cache_path": "~/.boltz",
                        "diffusion_samples": 2,
                    },
                    "runner": {"recycling_steps": 3},
                }
            ),
            encoding="utf-8",
        )

        def write_matrix(validator):
            results_dir = validator.wrk_dir / "results"
            results_dir.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(
                [
                    {
                        "model_name": "boltz2",
                        "repeat": repeat,
                        "diffusion_sample": sample,
                        "confidence_score": 0.9,
                    }
                    for repeat in (1, 2)
                    for sample in (0, 1)
                ]
            ).to_csv(results_dir / "system_metrics.csv", index=False)

        mock_validate_run.side_effect = write_matrix
        work_dir = temp_dir / "matrix_screen"
        results = Screen(
            wrk_dir=str(work_dir),
            system_path=str(sample_system_yaml),
            options_path=str(options_path),
            ligand_chain="B",
            library=str(sample_csv_file),
            smiles_column="smiles",
            col_id="compound_id",
            repeats=2,
        ).run()

        assert len(results) == 8
        assert not results.duplicated(
            ["compound_id", "repeat_id", "model_id", "sample_id"]
        ).any()
        executions = pd.read_csv(work_dir / "results" / "executions.csv")
        assert len(executions) == 8
        assert set(executions["status"]) == {"success"}
        assert set(executions["schema_version"]) == {"1.0.0"}

    def test_screen_generates_msa_once_then_reuses_it_for_repeats_and_rows(
        self,
        monkeypatch,
        sample_options_yaml,
        sample_csv_file,
        temp_dir,
    ):
        system_path = temp_dir / "system_with_ids.yaml"
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
        runner = _ReusableScreenRunner()
        monkeypatch.setattr("cofolder.recipes.screen.get_runner", lambda name: runner)
        monkeypatch.setattr("cofolder.recipes.validate.get_runner", lambda name: runner)

        results = Screen(
            wrk_dir=str(temp_dir / "screen"),
            system_path=str(system_path),
            options_path=str(sample_options_yaml),
            ligand_chain="B",
            library=str(sample_csv_file),
            smiles_column="smiles",
            col_id="compound_id",
            repeats=2,
            scoring_functions=["confidence_metrics"],
            assess_robustness=False,
        ).run()

        assert runner.msa_missing_at_run == [True, False, False, False]
        assert results["status"].tolist() == ["success"] * 4
        assert results["ligand_B__pair_chains_iptm_A"].tolist() == [0.7] * 4
        assert results.groupby("compound_id")["repeat_id"].apply(list).tolist() == [
            [1, 2],
            [1, 2],
        ]
        for row_index, compound_id in enumerate(("CMPD001", "CMPD002"), 1):
            row_system = yaml.safe_load(
                (
                    temp_dir
                    / "screen"
                    / f"compound_{row_index:06d}"
                    / "screen_system.yaml"
                ).read_text(encoding="utf-8")
            )
            msa_path = row_system["sequences"][0]["protein"]["msa"]
            assert Path(msa_path).is_file()
            assert str(temp_dir / "screen" / "shared" / "msa") in msa_path

    def test_mapped_screen_generates_once_per_unique_sequence(
        self, monkeypatch, caplog, sample_options_yaml, temp_dir
    ):
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
        csv_path = temp_dir / "variants.csv"
        csv_path.write_text(
            "id,sequence\na1,MKRAAT\nb1,MKRAAC\na2,MKRAAT\nb2,MKRAAC\n",
            encoding="utf-8",
        )
        runner = _ReusableScreenRunner()
        caplog.set_level(logging.INFO, logger="cofolder.screen")
        monkeypatch.setattr("cofolder.recipes.screen.get_runner", lambda name: runner)
        monkeypatch.setattr("cofolder.recipes.validate.get_runner", lambda name: runner)

        results = Screen(
            wrk_dir=str(temp_dir / "screen"),
            system_path=str(system_path),
            options_path=str(sample_options_yaml),
            library=str(csv_path),
            col_id="id",
            mappings=["sequence=sequences.0.protein.sequence"],
            repeats=2,
            scoring_functions=["confidence_metrics"],
            assess_robustness=False,
        ).run()

        assert runner.msa_missing_at_run == [True, False, True, False, False, False, False, False]
        assert results["status"].tolist() == ["success"] * 8
        row_msa_paths = []
        for index in range(1, 5):
            row_system = yaml.safe_load(
                (
                    temp_dir
                    / "screen"
                    / f"compound_{index:06d}"
                    / "screen_system.yaml"
                ).read_text(encoding="utf-8")
            )
            row_msa_paths.append(row_system["sequences"][0]["protein"]["msa"])
        assert row_msa_paths[0] == row_msa_paths[2]
        assert row_msa_paths[1] == row_msa_paths[3]
        assert row_msa_paths[0] != row_msa_paths[1]
        manifest = json.loads(
            (temp_dir / "screen" / "shared" / "msa" / "boltz2" / "manifest.json").read_text(
                encoding="utf-8"
            )
        )
        assert len(manifest["proteins"]) == 2
        messages = [record.getMessage() for record in caplog.records]
        assert sum("first generation required" in message for message in messages) == 2
        assert sum("sequence-matched shared cache" in message for message in messages) == 2

    def test_mapped_screen_reuses_each_sequence_across_multiple_proteins(
        self, monkeypatch, sample_options_yaml, temp_dir
    ):
        system_path = temp_dir / "two_proteins.yaml"
        system_path.write_text(
            yaml.safe_dump(
                {
                    "sequences": [
                        {"protein": {"id": "A", "sequence": "MKRAAT"}},
                        {"protein": {"id": "C", "sequence": "GGGG"}},
                    ]
                }
            ),
            encoding="utf-8",
        )
        csv_path = temp_dir / "variants.csv"
        csv_path.write_text(
            "id,sequence_a,sequence_c\n"
            "one,MKRAAT,GGGG\n"
            "two,MKRAAC,GGGG\n"
            "three,MKRAAT,GGGH\n"
            "four,MKRAAC,GGGH\n",
            encoding="utf-8",
        )
        runner = _ReusableScreenRunner()
        monkeypatch.setattr("cofolder.recipes.screen.get_runner", lambda name: runner)
        monkeypatch.setattr("cofolder.recipes.validate.get_runner", lambda name: runner)

        Screen(
            wrk_dir=str(temp_dir / "screen"),
            system_path=str(system_path),
            options_path=str(sample_options_yaml),
            library=str(csv_path),
            col_id="id",
            mappings=[
                "sequence_a=sequences.0.protein.sequence",
                "sequence_c=sequences.1.protein.sequence",
            ],
            scoring_functions=["confidence_metrics"],
            assess_robustness=False,
        ).run()

        assert runner.msa_sequences_missing_at_run == [
            ("MKRAAT", "GGGG"),
            ("MKRAAC",),
            ("GGGH",),
            (),
        ]

    @patch("cofolder.recipes.screen.Validate.run", autospec=True)
    def test_changed_sequence_with_inherited_msa_stops_before_prediction(
        self, mock_validate_run, sample_options_yaml, temp_dir
    ):
        input_dir = temp_dir / "inputs"
        input_dir.mkdir()
        (input_dir / "protein.a3m").write_text(">query\nMKRAAT\n", encoding="utf-8")
        system_path = input_dir / "system.yaml"
        system_path.write_text(
            yaml.safe_dump(
                {
                    "sequences": [
                        {
                            "protein": {
                                "id": "A",
                                "sequence": "MKRAAT",
                                "msa": "protein.a3m",
                            }
                        },
                        {"ligand": {"id": "B", "smiles": "CCO"}},
                    ]
                }
            ),
            encoding="utf-8",
        )
        csv_path = input_dir / "variants.csv"
        csv_path.write_text(
            "id,sequence\nbase,MKRAAT\nvariant_2,MKRAAC\n", encoding="utf-8"
        )
        screen = Screen(
            wrk_dir=str(temp_dir / "screen"),
            system_path=str(system_path),
            options_path=str(sample_options_yaml),
            library=str(csv_path),
            col_id="id",
            mappings=["sequence=sequences.0.protein.sequence"],
        )

        with pytest.raises(WorkflowExecutionError, match="inherits the fixed MSA"):
            screen.run()
        mock_validate_run.assert_not_called()

        with patch(
            "cofolder.modules.runners.boltz2_runner.Boltz2Runner.check_availability",
            return_value=(True, None),
        ):
            report = screen.preflight()
        assert not report.ready
        assert any("variant_2" in message for message in report.messages)

    def test_mapped_matching_msas_avoid_generation(
        self, monkeypatch, sample_options_yaml, temp_dir
    ):
        input_dir = temp_dir / "inputs"
        input_dir.mkdir()
        (input_dir / "a.a3m").write_text(">query\nMKRAAT\n", encoding="utf-8")
        (input_dir / "b.a3m").write_text(">query\nMKRAAC\n", encoding="utf-8")
        system_path = input_dir / "system.yaml"
        system_path.write_text(
            yaml.safe_dump(
                {
                    "sequences": [
                        {
                            "protein": {
                                "id": "A",
                                "sequence": "MKRAAT",
                                "msa": "a.a3m",
                            }
                        },
                        {"ligand": {"id": "B", "smiles": "CCO"}},
                    ]
                }
            ),
            encoding="utf-8",
        )
        csv_path = input_dir / "variants.csv"
        csv_path.write_text(
            "id,sequence,msa\na,MKRAAT,a.a3m\nb,MKRAAC,b.a3m\n",
            encoding="utf-8",
        )
        runner = _ReusableScreenRunner()
        monkeypatch.setattr("cofolder.recipes.screen.get_runner", lambda name: runner)
        monkeypatch.setattr("cofolder.recipes.validate.get_runner", lambda name: runner)

        Screen(
            wrk_dir=str(temp_dir / "screen"),
            system_path=str(system_path),
            options_path=str(sample_options_yaml),
            library=str(csv_path),
            col_id="id",
            mappings=[
                "sequence=sequences.0.protein.sequence",
                "msa=sequences.0.protein.msa",
            ],
            scoring_functions=["confidence_metrics"],
            assess_robustness=False,
        ).run()

        assert runner.msa_missing_at_run == [False, False]

    @patch("cofolder.recipes.screen.Validate.run", autospec=True)
    def test_mapped_mismatched_msa_stops_before_prediction(
        self, mock_validate_run, sample_options_yaml, temp_dir
    ):
        input_dir = temp_dir / "inputs"
        input_dir.mkdir()
        (input_dir / "wrong.a3m").write_text(">query\nMKRAAT\n", encoding="utf-8")
        system_path = input_dir / "system.yaml"
        system_path.write_text(
            yaml.safe_dump(
                {
                    "sequences": [
                        {
                            "protein": {
                                "id": "A",
                                "sequence": "MKRAAT",
                                "msa": "wrong.a3m",
                            }
                        },
                        {"ligand": {"id": "B", "smiles": "CCO"}},
                    ]
                }
            ),
            encoding="utf-8",
        )
        csv_path = input_dir / "variants.csv"
        csv_path.write_text(
            "id,sequence,msa\nvariant,MKRAAC,wrong.a3m\n", encoding="utf-8"
        )

        with pytest.raises(WorkflowExecutionError, match="does not match protein chain") as caught:
            Screen(
                wrk_dir=str(temp_dir / "screen"),
                system_path=str(system_path),
                options_path=str(sample_options_yaml),
                library=str(csv_path),
                col_id="id",
                mappings=[
                    "sequence=sequences.0.protein.sequence",
                    "msa=sequences.0.protein.msa",
                ],
            ).run()
        mock_validate_run.assert_not_called()
        assert "sequences.0.protein.msa" in str(caught.value)

    def test_screen_preserves_precomputed_msa_without_generation(
        self,
        monkeypatch,
        sample_options_yaml,
        sample_csv_file,
        temp_dir,
    ):
        input_dir = temp_dir / "inputs"
        input_dir.mkdir()
        (input_dir / "protein.a3m").write_text(">query\nMKRAAT\n", encoding="utf-8")
        system_path = input_dir / "system.yaml"
        system_path.write_text(
            yaml.safe_dump(
                {
                    "sequences": [
                        {
                            "protein": {
                                "id": "A",
                                "sequence": "MKRAAT",
                                "msa": "protein.a3m",
                            }
                        },
                        {"ligand": {"id": "B", "smiles": "CCO"}},
                    ]
                }
            ),
            encoding="utf-8",
        )
        runner = _ReusableScreenRunner()
        monkeypatch.setattr("cofolder.recipes.screen.get_runner", lambda name: runner)
        monkeypatch.setattr("cofolder.recipes.validate.get_runner", lambda name: runner)

        Screen(
            wrk_dir=str(temp_dir / "screen"),
            system_path=str(system_path),
            options_path=str(sample_options_yaml),
            ligand_chain="B",
            library=str(sample_csv_file),
            smiles_column="smiles",
            col_id="compound_id",
            scoring_functions=["confidence_metrics"],
            assess_robustness=False,
        ).run()

        assert runner.msa_missing_at_run == [False, False]
        row_system = yaml.safe_load(
            (temp_dir / "screen" / "compound_000001" / "screen_system.yaml").read_text(
                encoding="utf-8"
            )
        )
        assert row_system["sequences"][0]["protein"]["msa"] == str(
            (input_dir / "protein.a3m").resolve()
        )

    def test_screen_does_not_retry_when_runner_omits_generated_msa(
        self,
        monkeypatch,
        sample_options_yaml,
        sample_csv_file,
        temp_dir,
    ):
        system_path = temp_dir / "system_with_ids.yaml"
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
        runner = _ReusableScreenRunner(generate_msa=False)
        monkeypatch.setattr("cofolder.recipes.screen.get_runner", lambda name: runner)
        monkeypatch.setattr("cofolder.recipes.validate.get_runner", lambda name: runner)

        with pytest.raises(WorkflowExecutionError) as caught:
            Screen(
                wrk_dir=str(temp_dir / "screen"),
                system_path=str(system_path),
                options_path=str(sample_options_yaml),
                ligand_chain="B",
                library=str(sample_csv_file),
                smiles_column="smiles",
                col_id="compound_id",
                scoring_functions=["confidence_metrics"],
                assess_robustness=False,
            ).run()

        assert runner.msa_missing_at_run == [True]
        assert len(caught.value.failures) == 2
        assert (
            "without producing valid reusable MSAs" in caught.value.failures[0].message
        )
        assert "already attempted" in caught.value.failures[1].message

    def test_generated_msa_is_reused_after_prediction_failure(
        self,
        monkeypatch,
        sample_options_yaml,
        sample_csv_file,
        temp_dir,
    ):
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
        runner = _ReusableScreenRunner(fail_after_msa_once=True)
        monkeypatch.setattr("cofolder.recipes.screen.get_runner", lambda name: runner)
        monkeypatch.setattr("cofolder.recipes.validate.get_runner", lambda name: runner)

        results = Screen(
            wrk_dir=str(temp_dir / "screen"),
            system_path=str(system_path),
            options_path=str(sample_options_yaml),
            ligand_chain="B",
            library=str(sample_csv_file),
            smiles_column="smiles",
            col_id="compound_id",
            scoring_functions=["confidence_metrics"],
            assess_robustness=False,
        ).run()

        assert runner.msa_missing_at_run == [True, False]
        assert results["status"].tolist() == ["failed", "success"]

    @patch("cofolder.recipes.screen.Validate.run")
    def test_run_preserves_constraints_unrelated_to_ligand_replacement(
        self, mock_validate_run, sample_options_yaml, temp_dir
    ):
        constraints = [
            {"contact": {"token1": ["A", 1], "token2": ["D", 1], "max_distance": 5.0}}
        ]
        system_path = temp_dir / "system.yaml"
        system_path.write_text(
            yaml.safe_dump(
                {
                    "sequences": [
                        {"protein": {"id": "A", "sequence": "AC"}},
                        {"dna": {"id": "D", "sequence": "AT"}},
                        {"ligand": {"id": "L", "smiles": "CCO"}},
                    ],
                    "constraints": constraints,
                }
            ),
            encoding="utf-8",
        )
        csv_path = temp_dir / "ligands.csv"
        csv_path.write_text("id,smiles\none,CCN\n", encoding="utf-8")

        Screen(
            wrk_dir=str(temp_dir / "screen"),
            system_path=str(system_path),
            options_path=str(sample_options_yaml),
            ligand_chain="L",
            library=str(csv_path),
            smiles_column="smiles",
            col_id="id",
        ).run()

        row_system = yaml.safe_load(
            (temp_dir / "screen" / "compound_000001" / "screen_system.yaml").read_text(
                encoding="utf-8"
            )
        )
        assert row_system["constraints"] == constraints

    def test_run_reports_ligand_replacement_that_invalidates_bond(
        self, monkeypatch, sample_options_yaml, temp_dir
    ):
        system_path = temp_dir / "system.yaml"
        system_path.write_text(
            yaml.safe_dump(
                {
                    "sequences": [
                        {"protein": {"id": "A", "sequence": "AC"}},
                        {"ligand": {"id": "L", "smiles": "CCO"}},
                    ],
                    "constraints": [
                        {"bond": {"atom1": ["A", 2, "SG"], "atom2": ["L", 1, "C1"]}}
                    ],
                }
            ),
            encoding="utf-8",
        )
        csv_path = temp_dir / "ligands.csv"
        csv_path.write_text("id,smiles\none,[Na+]\n", encoding="utf-8")

        class _PreflightOnlyValidate:
            def __init__(self, **kwargs):
                self.system_path = kwargs["system_path"]

            def run(self):
                Boltz2Runner().validate_system(
                    System(system_path=self.system_path), {}, check_atom_names=True
                )

        monkeypatch.setattr("cofolder.recipes.screen.Validate", _PreflightOnlyValidate)

        with pytest.raises(WorkflowExecutionError) as caught:
            Screen(
                wrk_dir=str(temp_dir / "screen"),
                system_path=str(system_path),
                options_path=str(sample_options_yaml),
                ligand_chain="L",
                library=str(csv_path),
                smiles_column="smiles",
                col_id="id",
            ).run()

        assert "atom 'C1' does not exist" in caught.value.failures[0].message
        assert (temp_dir / "screen" / "results" / "failures.csv").is_file()

    @patch("cofolder.recipes.screen.Validate")
    def test_run_uses_boltz2_as_default_runner(
        self,
        mock_validate_cls,
        sample_system_yaml,
        sample_options_yaml,
        sample_csv_file,
        temp_dir,
    ):
        mock_validator = MagicMock()
        mock_validate_cls.return_value = mock_validator

        screener = Screen(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            ligand_chain="B",
            library=str(sample_csv_file),
            smiles_column="smiles",
            col_id="compound_id",
            bias_query_cache_path=str(temp_dir / "shared-bias-cache"),
        )

        merged = screener.run()

        first_call_kwargs = mock_validate_cls.call_args_list[0].kwargs
        assert first_call_kwargs["runner"] == "boltz2"
        assert first_call_kwargs["scoring_functions"] is None
        assert {
            call.kwargs["bias_query_cache_path"]
            for call in mock_validate_cls.call_args_list
        } == {str(temp_dir / "shared-bias-cache")}
        assert mock_validator.run.call_count == 2

        assert merged["status"].tolist() == ["success", "success"]
        assert (temp_dir / "results" / "records.jsonl").is_file()

    @patch("cofolder.recipes.screen.Validate")
    def test_run_passes_explicit_boltz_community_runner(
        self,
        mock_validate_cls,
        sample_system_yaml,
        sample_options_yaml,
        sample_csv_file,
        temp_dir,
    ):
        mock_validator = MagicMock()
        mock_validate_cls.return_value = mock_validator

        with patch(
            "cofolder.modules.runners.boltz_community_runner.BoltzCommunityRunner.check_availability",
            return_value=(True, None),
        ):
            screener = Screen(
                wrk_dir=str(temp_dir),
                system_path=str(sample_system_yaml),
                options_path=str(sample_options_yaml),
                runner="boltz-community",
                ligand_chain="B",
                library=str(sample_csv_file),
                smiles_column="smiles",
                col_id="compound_id",
            )
            screener.run()

        first_call_kwargs = mock_validate_cls.call_args_list[0].kwargs
        assert first_call_kwargs["runner"] == "boltz-community"
        assert mock_validator.run.call_count == 2

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
            "compound_id,smiles\nCMPD_BAD,C1(\nCMPD_OK,CCO\n",
            encoding="utf-8",
        )

        screener = Screen(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            ligand_chain="B",
            library=str(csv_path),
            smiles_column="smiles",
            col_id="compound_id",
        )
        out_df = screener.run()

        assert mock_validate_run.call_count == 1
        assert len(out_df) == 2
        bad = out_df[out_df["compound_id"] == "CMPD_BAD"].iloc[0]
        good = out_df[out_df["compound_id"] == "CMPD_OK"].iloc[0]
        assert bad["status"] == "unavailable"
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
        input_df["pIC50"] = [5.4, 6.1]
        input_df.to_csv(sample_csv_file, index=False)
        for i, (_, in_row) in enumerate(input_df.iterrows(), 1):
            run_dir = temp_dir / f"compound_{i:06d}"
            results_dir = run_dir / "results"
            results_dir.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(
                [
                    {
                        "repeat": 1,
                        "diffusion_sample": 0,
                        "confidence_score": 0.9,
                        "bias_prot_sim_train_max": 42.0,
                        "bias_lig_sim_train_max": 0.35,
                        "pocket_coverage_custom_mean": 0.75,
                    }
                ]
            ).to_csv(results_dir / "system_metrics.csv", index=False)
            pd.DataFrame(
                [
                    {
                        "conf_chain_id": 0,
                        "CHAIN_ID": "A",
                        "ENTITY_TYPE": "protein",
                        "repeat": 1,
                        "diffusion_sample": 0,
                        "chains_ptm": 0.8,
                        "pair_chains_iptm_1": 0.72,
                        "bias_prot_sim_train": 42.0,
                        "sasa_norm_heavy": 5.1,
                    },
                    {
                        "conf_chain_id": 1,
                        "CHAIN_ID": "B",
                        "ENTITY_TYPE": "ligand",
                        "repeat": 1,
                        "diffusion_sample": 0,
                        "chains_ptm": 0.7,
                        "pair_chains_iptm_0": 0.72,
                        "affinity_pred_value": 6.2,
                        "affinity_probability_binary": 0.88,
                        "pIC50": 6.2,
                        "sasa": 12.5,
                        "sasa_norm_heavy": 1.25,
                        "ifp_distance": "[1, 0, 1]",
                        "bias_lig_sim_train": 0.35,
                        "pocket_coverage_custom": 0.75,
                    },
                ]
            ).to_csv(results_dir / "chain_metrics.csv", index=False)

        screener = Screen(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            ligand_chain="B",
            library=str(sample_csv_file),
            smiles_column="smiles",
            col_id="compound_id",
            merge_data="mw,pIC50",
        )
        returned_df = screener.run()

        merged_df = returned_df
        assert returned_df.columns.tolist() == merged_df.columns.tolist()
        assert returned_df["compound_id"].tolist() == merged_df["compound_id"].tolist()
        assert returned_df["ligand_B__affinity_pred_value"].tolist() == [6.2, 6.2]
        assert "system__confidence_score" in merged_df.columns
        assert "protein_A__sasa_norm_heavy" in merged_df.columns
        assert "ligand_B__affinity_pred_value" in merged_df.columns
        assert "ligand_B__affinity_probability_binary" in merged_df.columns
        assert "ligand_B__ifp_distance" in merged_df.columns
        assert "system__bias_prot_sim_train_max" in merged_df.columns
        assert "system__bias_lig_sim_train_max" in merged_df.columns
        assert "ligand_B__pocket_coverage_custom" in merged_df.columns
        assert merged_df["ligand_B__affinity_probability_binary"].tolist() == [
            0.88,
            0.88,
        ]
        assert merged_df["ligand_B__pair_chains_iptm_A"].tolist() == [0.72, 0.72]
        assert merged_df["ligand_B__sasa_norm_heavy"].tolist() == [1.25, 1.25]
        assert merged_df["ligand_B__ifp_distance"].tolist() == ["[1, 0, 1]"] * 2
        assert merged_df["system__bias_prot_sim_train_max"].tolist() == [42.0, 42.0]
        assert merged_df["system__bias_lig_sim_train_max"].tolist() == [0.35, 0.35]
        assert merged_df["ligand_B__pocket_coverage_custom"].tolist() == [0.75, 0.75]
        assert "mw" in merged_df.columns
        assert merged_df["pIC50"].tolist() == [5.4, 6.1]
        assert "compound_id" in merged_df.columns
        public_metrics = pd.read_csv(temp_dir / "results" / "metrics.csv")
        pic50 = public_metrics[public_metrics["metric_name"] == "pIC50"]
        assert not pic50["chain_id"].isna().any()
        assert set(pic50.loc[pic50["status"] == "computed", "chain_id"]) == {"B"}

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
            ligand_chain="B",
            library=str(sample_csv_file),
            smiles_column="smiles",
            col_id="compound_id",
            merge_data="mw",
        )

        out_df = screener.run()

        assert mock_validate_run.call_count == 2

        assert (temp_dir / "results" / "records.jsonl").exists()
        assert len(out_df) == 2
        assert set(
            [
                "index",
                "compound_id",
                "status",
                "error_message",
                "run_dir",
                "smiles",
                "mw",
            ]
        ).issubset(out_df.columns)
        assert set(out_df["status"].tolist()) == {"success"}

        merged_df = out_df
        assert len(merged_df) == 2
        assert set(
            ["compound_id", "smiles", "status", "error_message", "run_dir"]
        ).issubset(merged_df.columns)

    @patch("cofolder.recipes.screen.Validate.run")
    def test_run_mapped_csv_updates_multiple_system_fields(
        self,
        mock_validate_run,
        sample_system_yaml,
        sample_options_yaml,
        temp_dir,
    ):
        csv_path = temp_dir / "multi.csv"
        csv_path.write_text(
            "experiment,protein_sequence,ligand_smiles,ligand_id\n"
            "variant_1,MKRAAC,CCN,C\n",
            encoding="utf-8",
        )

        screener = Screen(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            library=str(csv_path),
            col_id="experiment",
            mappings=[
                "protein_sequence=sequences.0.protein.sequence",
                "ligand_smiles=sequences.1.ligand.smiles",
                "ligand_id=sequences.1.ligand.id",
            ],
        )

        result = screener.run()

        run_system_yaml = temp_dir / "compound_000001" / "screen_system.yaml"
        assert run_system_yaml.exists()
        data = yaml.safe_load(run_system_yaml.read_text(encoding="utf-8"))
        assert data["sequences"][0]["protein"]["sequence"] == "MKRAAC"
        lig = data["sequences"][1]["ligand"]
        assert lig["smiles"] == "CCN"
        assert lig["id"] == "C"
        assert "ccd" not in lig
        assert "ligand_C__chains_ptm" in result.columns
        mapping_artifact = temp_dir / "results" / "system_mappings.json"
        assert json.loads(mapping_artifact.read_text(encoding="utf-8")) == [
            {
                "column": "protein_sequence",
                "yaml_path": "sequences.0.protein.sequence",
            },
            {
                "column": "ligand_smiles",
                "yaml_path": "sequences.1.ligand.smiles",
            },
            {"column": "ligand_id", "yaml_path": "sequences.1.ligand.id"},
        ]

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

        mock_validate_run.side_effect = _run_side_effect

        screener = Screen(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            ligand_chain="B",
            library=str(sample_csv_file),
            smiles_column="smiles",
            col_id="compound_id",
        )

        out_df = screener.run()

        assert len(out_df) == 2
        assert "failed" in set(out_df["status"].tolist())
        assert "success" in set(out_df["status"].tolist())

        merged_df = out_df
        assert len(merged_df) == 2
        assert "failed" in set(merged_df["status"].tolist())
        assert "success" in set(merged_df["status"].tolist())

    @patch("cofolder.recipes.screen.Validate.run")
    def test_mixed_sdf_preserves_failures_and_source_mapping(
        self,
        mock_validate_run,
        sample_system_yaml,
        sample_options_yaml,
        temp_dir,
    ):
        sdf_path = temp_dir / "mixed.sdf"
        molecule = Chem.MolFromSmiles("CCO")
        molecule.SetProp("_Name", "valid")
        sdf_path.write_text(
            Chem.MolToMolBlock(molecule)
            + "\n$$$$\n"
            + "malformed record\n$$$$\n",
            encoding="utf-8",
        )

        results = Screen(
            wrk_dir=str(temp_dir / "screen"),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            ligand_chain="B",
            library=str(sdf_path),
        ).run()

        assert mock_validate_run.call_count == 1
        assert results["status"].tolist() == ["success", "unavailable"]
        assert results["execution_id"].tolist() == ["valid", "malformed record"]
        assert (temp_dir / "screen" / "compound_000001" / "source_ligand.sdf").is_file()
        members = pd.read_csv(temp_dir / "screen" / "results" / "compound_members.csv")
        assert members["source_record_id"].tolist() == [
            "record_000001",
            "record_000002",
        ]
        assert members["status"].tolist() == ["success", "failed"]
        failures = pd.read_csv(temp_dir / "screen" / "results" / "failures.csv")
        assert failures["error_code"].tolist() == ["molblock_record_parse_failed"]
        assert '"source_record_id": "record_000002"' in failures.loc[0, "details"]
