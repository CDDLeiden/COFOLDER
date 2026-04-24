"""Tests for cofolder.recipes.validate module."""

import logging
from pathlib import Path

import pandas as pd
import yaml

from cofolder.modules.runners.base import RunnerPreparation, RunnerResult
from cofolder.recipes.validate import Validate


def _make_runner_results(system_name: str = "system"):
    system_df = pd.DataFrame(
        [
            {
                "cif_file": f"1_{system_name}_model_0.cif",
                "model_name": system_name,
                "repeat": 1,
                "diffusion_sample": 0,
                "ptm": 0.8,
                "iptm": 0.7,
                "confidence_score": 0.9,
            }
        ]
    )
    chain_df = pd.DataFrame(
        [
            {
                "conf_chain_id": 0,
                "CHAIN_ID": "A",
                "ENTITY_TYPE": "protein",
                "ligand_molecule_id": "protein_A",
                "cif_file": f"1_{system_name}_model_0.cif",
                "model_name": system_name,
                "repeat": 1,
                "diffusion_sample": 0,
                "chains_ptm": 0.85,
            },
            {
                "conf_chain_id": 1,
                "CHAIN_ID": "B",
                "ENTITY_TYPE": "ligand",
                "ligand_molecule_id": "ETH",
                "cif_file": f"1_{system_name}_model_0.cif",
                "model_name": system_name,
                "repeat": 1,
                "diffusion_sample": 0,
                "chains_ptm": 0.65,
                "affinity_pred_value": 6.1,
                "affinity_probability_binary": 0.8,
                "pIC50": 5.0,
                "IC50_M": 1e-6,
                "pIC50_kcal_per_mol": 7.0,
            },
        ]
    )
    manifests = [
        {
            "runner": "boltz",
            "runtime_context": {"cache_path": "~/.boltz", "diffusion_samples": 1},
        }
    ]
    return system_df, chain_df, manifests


class _FakeRunner:
    name = "boltz"
    capabilities = {"confidence_metrics", "affinity_metrics", "affinity_metrics_ext"}

    def __init__(self, system_name: str = "system"):
        self.system_name = system_name

    def is_available(self):
        return True, None

    def load_options(self, options_path):
        return {"options_path": str(options_path)}

    def prepare_system(self, system_obj, options_obj, wrk_dir, conformers, sdf_file, logger):
        return RunnerPreparation(
            system_obj=system_obj,
            options_obj=options_obj,
            runtime_context={"cache_path": "~/.boltz", "diffusion_samples": 1},
        )

    def run(self, request):
        if request.timings is not None:
            request.timings.record(
                f"{request.label_prefix}.boltz.total",
                1.250,
                logger=logging.getLogger("cofolder.recipes.validate"),
            )
            request.timings.record(
                f"{request.label_prefix}.boltz.msa",
                0.400,
                logger=logging.getLogger("cofolder.recipes.validate"),
            )
            request.timings.record(
                f"{request.label_prefix}.boltz.affinity_prediction",
                0.300,
                logger=logging.getLogger("cofolder.recipes.validate"),
            )

        normalized_dir = request.repeat_dir / "normalized"
        structures_dir = normalized_dir / "structures"
        structures_dir.mkdir(parents=True, exist_ok=True)
        (structures_dir / f"{request.repeat}_{self.system_name}_model_0.cif").write_text("data", encoding="utf-8")

        return RunnerResult(
            runner_name="boltz",
            raw_output_dir=request.repeat_dir,
            normalized_dir=normalized_dir,
            structures_dir=structures_dir,
            system_metrics_path=normalized_dir / "system_metrics.csv",
            chain_metrics_path=normalized_dir / "chain_metrics.csv",
            manifest_path=normalized_dir / "manifest.json",
            diffusion_samples=1,
            capabilities=set(self.capabilities),
            runtime_context={"cache_path": "~/.boltz", "diffusion_samples": 1},
        )


class _NoMetricsRunner(_FakeRunner):
    capabilities = set()


def _patch_validate_pipeline(monkeypatch, system_name: str = "system"):
    monkeypatch.setattr(
        "cofolder.recipes.validate.helpers.get_seeds",
        lambda repeats, seed, logger: (123, [123]),
    )
    monkeypatch.setattr(
        "cofolder.recipes.validate.get_runner",
        lambda name: _FakeRunner(system_name=system_name),
    )
    monkeypatch.setattr(
        "cofolder.recipes.validate.gather.gather_structures",
        lambda base_dir, system_name, repeats, logger: None,
    )
    monkeypatch.setattr(
        "cofolder.recipes.validate.gather.merge_runner_results",
        lambda raw_dir, repeats, logger: _make_runner_results(system_name=system_name),
    )
    monkeypatch.setattr(
        "cofolder.recipes.validate.gather.add_chain_info",
        lambda chain_df, sys: chain_df,
    )
    monkeypatch.setattr(
        "cofolder.recipes.validate.scaffold_reproduction_metrics",
        lambda system_df, chain_df, reference_path, wrk_dir, pocket_coverage_reference, reproduction_metrics, logger: (
            system_df,
            chain_df,
        ),
    )

    class _FakeStructure:
        def __init__(self, wrk_dir, chain_df, cif_folder):
            self.chain_df = chain_df.copy()

        def add_ifp_distance(self):
            self.chain_df["ifp_distance"] = [None, "[]"]
            return self.chain_df

        def add_ifp_prolif(self):
            self.chain_df["ifp_prolif"] = [None, "ligand_B_ifp.pkl"]
            return self.chain_df

        def add_sasa(self, absolute=False, normalized=False):
            if absolute:
                self.chain_df["sasa"] = [10.0, 2.0]
            if normalized:
                self.chain_df["sasa_norm_heavy"] = [1.0, 0.2]
            return self.chain_df

    monkeypatch.setattr("cofolder.recipes.validate.Structure", _FakeStructure)


def _patch_validate_pipeline_no_metrics(monkeypatch, system_name: str = "system"):
    monkeypatch.setattr(
        "cofolder.recipes.validate.helpers.get_seeds",
        lambda repeats, seed, logger: (123, [123]),
    )
    monkeypatch.setattr(
        "cofolder.recipes.validate.get_runner",
        lambda name: _NoMetricsRunner(system_name=system_name),
    )
    monkeypatch.setattr(
        "cofolder.recipes.validate.gather.gather_structures",
        lambda base_dir, system_name, repeats, logger: None,
    )
    monkeypatch.setattr(
        "cofolder.recipes.validate.gather.merge_runner_results",
        lambda raw_dir, repeats, logger: (
            pd.DataFrame(
                [
                    {
                        "cif_file": f"1_{system_name}_model_0.cif",
                        "model_name": system_name,
                        "repeat": 1,
                        "diffusion_sample": 0,
                    }
                ]
            ),
            pd.DataFrame(
                [
                    {
                        "conf_chain_id": 0,
                        "CHAIN_ID": "A",
                        "ENTITY_TYPE": "protein",
                        "ligand_molecule_id": "protein_A",
                        "cif_file": f"1_{system_name}_model_0.cif",
                        "model_name": system_name,
                        "repeat": 1,
                        "diffusion_sample": 0,
                    }
                ]
            ),
            [{"runner": "boltz", "runtime_context": {"diffusion_samples": 1}}],
        ),
    )
    monkeypatch.setattr(
        "cofolder.recipes.validate.gather.add_chain_info",
        lambda chain_df, sys: chain_df,
    )
    monkeypatch.setattr(
        "cofolder.recipes.validate.scaffold_reproduction_metrics",
        lambda system_df, chain_df, reference_path, wrk_dir, pocket_coverage_reference, reproduction_metrics, logger: (
            system_df,
            chain_df,
        ),
    )


class TestValidateInit:
    """Tests for Validate initialization."""

    def test_init_basic(self, sample_system_yaml, sample_options_yaml, temp_dir):
        validator = Validate(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            scoring_functions=[],
        )

        assert str(validator.wrk_dir) == str(temp_dir)
        assert str(validator.system_path) == str(sample_system_yaml)
        assert str(validator.options_path) == str(sample_options_yaml)
        assert validator.base_system is not None
        assert validator.runner_name == "boltz"
        assert validator.scoring_functions == set()
        assert validator.reference_path is None
        assert validator.reproduction_metrics == {
            "protein_rmsd",
            "ligand_rmsd",
            "sucos",
            "pocket_coverage",
        }

    def test_init_with_reference_path(self, sample_system_yaml, sample_options_yaml, temp_dir):
        reference_path = temp_dir / "reference.pdb"
        reference_path.write_text("HEADER TEST\n")

        validator = Validate(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            scoring_functions=[],
            reference_path=str(reference_path),
            reproduction_metrics=["sucos"],
        )

        assert str(validator.reference_path) == str(reference_path)
        assert validator.reproduction_metrics == {"sucos"}

    def test_init_with_pocket_coverage_reference(
        self,
        sample_system_yaml,
        sample_options_yaml,
        temp_dir,
    ):
        validator = Validate(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            scoring_functions=[],
            pocket_coverage_reference="A2 S8 T10",
            reproduction_metrics=["pocket_coverage"],
        )

        assert validator.pocket_coverage_reference == "A2 S8 T10"
        assert validator.reproduction_metrics == {"pocket_coverage"}


class TestValidateRun:
    def test_debug_run_logs_timing_summary(
        self,
        monkeypatch,
        sample_system_yaml,
        sample_options_yaml,
        temp_dir,
        caplog,
    ):
        _patch_validate_pipeline(monkeypatch)

        validator = Validate(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            scoring_functions=["ifp_distance", "ifp_prolif", "sasa", "sasa_normalized"],
        )

        with caplog.at_level(logging.DEBUG, logger="cofolder.recipes.validate"):
            validator.run()

        assert "TIMER SUMMARY | validate.total" in caplog.text
        assert "TIMER SUMMARY | repeat_1.boltz.total" in caplog.text
        assert "TIMER SUMMARY | repeat_1.boltz.msa" in caplog.text
        assert "TIMER SUMMARY | repeat_1.boltz.affinity_prediction" in caplog.text
        assert "TIMER SUMMARY | scores.ifp_distance" in caplog.text
        assert "TIMER SUMMARY | scores.ifp_prolif" in caplog.text
        assert "TIMER SUMMARY | scores.sasa" in caplog.text

    def test_debug_run_logs_bias_similarity_timings(
        self,
        monkeypatch,
        sample_options_yaml,
        temp_dir,
        caplog,
    ):
        _patch_validate_pipeline(monkeypatch, system_name="bias_system")

        system_data = {
            "sequences": [
                {"protein": {"id": "A", "fasta": "MKRAAT"}},
                {"ligand": {"id": "B", "smiles": "CCO", "ccd": "ETH"}},
            ]
        }
        system_path = temp_dir / "bias_system.yaml"
        system_path.write_text(yaml.safe_dump(system_data), encoding="utf-8")

        protein_ref = temp_dir / "protein_training.csv"
        ligand_ref = temp_dir / "ligand_training.csv"
        protein_ref.write_text(
            "pdb_id,release_date,sequence,sequence_similarity\n"
            "1ABC,2022-01-01,MKRAAT,100.0\n",
            encoding="utf-8",
        )
        ligand_ref.write_text(
            "pdb_id,release_date,ligand_id,smiles,ecfp_similarity\n"
            "1ABC,2022-01-01,ETH,CCO,1.0\n",
            encoding="utf-8",
        )

        validator = Validate(
            wrk_dir=str(temp_dir),
            system_path=str(system_path),
            options_path=str(sample_options_yaml),
            scoring_functions=[],
            assess_bias=True,
            protein_training_data_path=str(protein_ref),
            ligand_training_data_path=str(ligand_ref),
        )

        with caplog.at_level(logging.DEBUG, logger="cofolder.recipes.validate"):
            validator.run()

        assert "TIMER SUMMARY | scores.bias_metrics.total" in caplog.text
        assert "TIMER SUMMARY | scores.bias_metrics.protein_similarity" in caplog.text
        assert "TIMER SUMMARY | scores.bias_metrics.ligand_similarity" in caplog.text

    def test_non_debug_run_does_not_log_timers(
        self,
        monkeypatch,
        sample_system_yaml,
        sample_options_yaml,
        temp_dir,
        caplog,
    ):
        _patch_validate_pipeline(monkeypatch)

        validator = Validate(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            scoring_functions=["sasa"],
        )

        with caplog.at_level(logging.INFO, logger="cofolder.recipes.validate"):
            validator.run()

        assert "TIMER |" not in caplog.text
        assert "TIMER SUMMARY |" not in caplog.text

    def test_unsupported_metric_groups_warn_and_write_empty_columns(
        self,
        monkeypatch,
        sample_system_yaml,
        sample_options_yaml,
        temp_dir,
        caplog,
    ):
        _patch_validate_pipeline_no_metrics(monkeypatch)

        validator = Validate(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            scoring_functions=["confidence_metrics", "affinity_metrics", "affinity_metrics_ext"],
        )

        with caplog.at_level(logging.WARNING, logger="cofolder.recipes.validate"):
            validator.run()

        assert "does not support scoring groups" in caplog.text

        system_df = pd.read_csv(Path(temp_dir) / "results" / "system_metrics.csv")
        chain_df = pd.read_csv(Path(temp_dir) / "results" / "chain_metrics.csv")

        assert {"ptm", "iptm", "confidence_score"}.issubset(system_df.columns)
        assert {"affinity_pred_value", "affinity_probability_binary", "pIC50", "IC50_M", "pIC50_kcal_per_mol"}.issubset(chain_df.columns)
