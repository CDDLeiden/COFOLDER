"""Tests for cofolder.recipes.validate module."""

import json
import logging
from pathlib import Path
from unittest.mock import Mock
from unittest.mock import patch

import pandas as pd
import pytest
import yaml

from cofolder.modules.runners.contracts import (
    RunnerExecutionResult,
    RunnerMetricOutcome,
    RunnerPreparationResult,
    RunnerRuntime,
)
from cofolder.modules.runners.boltz1_runner import Boltz1Runner
from cofolder.modules.runners.boltz_community_runner import BoltzCommunityRunner
from cofolder.modules.runners.boltz2_runner import Boltz2Runner
from cofolder.modules.runners.validators import RunnerBundleValidationError
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
            "runner": "boltz2",
            "runtime_context": {"cache_path": "~/.boltz", "diffusion_samples": 1},
        }
    ]
    return system_df, chain_df, manifests


def _write_normalized_bundle(
    repeat_dir: Path,
    *,
    system_name: str = "system",
    include_confidence: bool = True,
    include_affinity: bool = True,
    capabilities: set[str] | None = None,
    runtime_context: dict[str, object] | None = None,
) -> tuple[Path, Path, Path, Path, Path]:
    normalized_dir = repeat_dir / "normalized"
    structures_dir = normalized_dir / "structures"
    structures_dir.mkdir(parents=True, exist_ok=True)
    structure_name = f"1_{system_name}_model_0.cif"
    (structures_dir / structure_name).write_text("data", encoding="utf-8")

    system_row = {
        "cif_file": structure_name,
        "model_name": system_name,
        "repeat": 1,
        "diffusion_sample": 0,
    }
    if include_confidence:
        system_row.update(
            {
                "ptm": 0.8,
                "iptm": 0.7,
                "confidence_score": 0.9,
            }
        )

    chain_row = {
        "conf_chain_id": 0,
        "cif_file": structure_name,
        "model_name": system_name,
        "repeat": 1,
        "diffusion_sample": 0,
    }
    if include_confidence:
        chain_row["chains_ptm"] = 0.85
    if include_affinity:
        chain_row.update(
            {
                "affinity_pred_value": 6.1,
                "affinity_probability_binary": 0.8,
                "pIC50": 5.0,
                "IC50_M": 1e-6,
                "pIC50_kcal_per_mol": 7.0,
            }
        )

    system_metrics_path = normalized_dir / "system_metrics.csv"
    chain_metrics_path = normalized_dir / "chain_metrics.csv"
    manifest_path = normalized_dir / "manifest.json"
    pd.DataFrame([system_row]).to_csv(system_metrics_path, index=False)
    pd.DataFrame([chain_row]).to_csv(chain_metrics_path, index=False)
    manifest_path.write_text(
        json.dumps(
            {
                "runner": "boltz2",
                "capabilities": sorted(capabilities or []),
                "runtime_context": runtime_context or {"diffusion_samples": 1},
            }
        ),
        encoding="utf-8",
    )
    return normalized_dir, structures_dir, system_metrics_path, chain_metrics_path, manifest_path


class _FakeRunner:
    name = "boltz2"
    capabilities = {"confidence_metrics", "affinity_metrics", "affinity_metrics_ext"}

    def __init__(self, system_name: str = "system"):
        self.system_name = system_name

    def is_available(self):
        return True, None

    def check_availability(self):
        return True, None

    def ensure_available(self):
        return None

    def load_options(self, options_path):
        return {"options_path": str(options_path)}

    def prepare_system(self, system_obj, options_obj, wrk_dir, conformers, sdf_file, logger):
        return RunnerPreparationResult(
            system_obj=system_obj,
            options_obj=options_obj,
            runtime=RunnerRuntime(cache_path="~/.boltz", diffusion_samples=1, model_name="boltz2"),
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

        normalized_dir, structures_dir, system_metrics_path, chain_metrics_path, manifest_path = (
            _write_normalized_bundle(
                request.repeat_dir,
                system_name=self.system_name,
                capabilities=self.capabilities,
                runtime_context={"cache_path": "~/.boltz", "diffusion_samples": 1},
            )
        )

        return RunnerExecutionResult(
            runner_name="boltz2",
            raw_output_dir=request.repeat_dir,
            normalized_dir=normalized_dir,
            structures_dir=structures_dir,
            system_metrics_path=system_metrics_path,
            chain_metrics_path=chain_metrics_path,
            manifest_path=manifest_path,
            diffusion_samples=1,
            capabilities=set(self.capabilities),
            runtime=RunnerRuntime(cache_path="~/.boltz", diffusion_samples=1, model_name="boltz2"),
            metric_outcomes={
                group_name: RunnerMetricOutcome(
                    state="computed" if group_name in self.capabilities else "unsupported"
                )
                for group_name in (
                    "confidence_metrics",
                    "affinity_metrics",
                    "affinity_metrics_ext",
                )
            },
        )


class _NoMetricsRunner(_FakeRunner):
    capabilities = set()

    def run(self, request):
        normalized_dir, structures_dir, system_metrics_path, chain_metrics_path, manifest_path = (
            _write_normalized_bundle(
                request.repeat_dir,
                system_name=self.system_name,
                include_confidence=False,
                include_affinity=False,
                capabilities=self.capabilities,
                runtime_context={"cache_path": "~/.boltz", "diffusion_samples": 1},
            )
        )

        return RunnerExecutionResult(
            runner_name="boltz2",
            raw_output_dir=request.repeat_dir,
            normalized_dir=normalized_dir,
            structures_dir=structures_dir,
            system_metrics_path=system_metrics_path,
            chain_metrics_path=chain_metrics_path,
            manifest_path=manifest_path,
            diffusion_samples=1,
            capabilities=set(self.capabilities),
            runtime=RunnerRuntime(cache_path="~/.boltz", diffusion_samples=1, model_name="boltz2"),
            metric_outcomes={
                "confidence_metrics": RunnerMetricOutcome(state="unsupported"),
                "affinity_metrics": RunnerMetricOutcome(state="unsupported"),
                "affinity_metrics_ext": RunnerMetricOutcome(state="unsupported"),
            },
        )


class _ManifestOnlyRuntimeRunner(_FakeRunner):
    def prepare_system(self, system_obj, options_obj, wrk_dir, conformers, sdf_file, logger):
        return RunnerPreparationResult(
            system_obj=system_obj,
            options_obj=options_obj,
            runtime=RunnerRuntime(diffusion_samples=1, cache_path="/typed/cache", model_name="boltz2"),
        )

    def run(self, request):
        normalized_dir, structures_dir, system_metrics_path, chain_metrics_path, manifest_path = (
            _write_normalized_bundle(
                request.repeat_dir,
                system_name=self.system_name,
                include_confidence=False,
                include_affinity=False,
                capabilities=self.capabilities,
                runtime_context={"cache_path": "/legacy/cache", "diffusion_samples": 3},
            )
        )
        return RunnerExecutionResult(
            runner_name="boltz2",
            raw_output_dir=request.repeat_dir,
            normalized_dir=normalized_dir,
            structures_dir=structures_dir,
            system_metrics_path=system_metrics_path,
            chain_metrics_path=chain_metrics_path,
            manifest_path=manifest_path,
            capabilities=set(self.capabilities),
            runtime=RunnerRuntime(diffusion_samples=1, cache_path="/typed/cache", model_name="boltz2"),
        )


class _TypedDiffusionRuntimeRunner(_FakeRunner):
    def prepare_system(self, system_obj, options_obj, wrk_dir, conformers, sdf_file, logger):
        return RunnerPreparationResult(
            system_obj=system_obj,
            options_obj=options_obj,
        )

    def run(self, request):
        normalized_dir, structures_dir, system_metrics_path, chain_metrics_path, manifest_path = (
            _write_normalized_bundle(
                request.repeat_dir,
                system_name=self.system_name,
                include_confidence=False,
                include_affinity=False,
                capabilities=self.capabilities,
                runtime_context={"cache_path": "/legacy/cache", "diffusion_samples": 1},
            )
        )
        return RunnerExecutionResult(
            runner_name="boltz2",
            raw_output_dir=request.repeat_dir,
            normalized_dir=normalized_dir,
            structures_dir=structures_dir,
            system_metrics_path=system_metrics_path,
            chain_metrics_path=chain_metrics_path,
            manifest_path=manifest_path,
            capabilities=set(self.capabilities),
            runtime=RunnerRuntime(diffusion_samples=3, cache_path="/typed/cache", model_name="boltz2"),
        )


class _ImplicitOutcomeRunner(_FakeRunner):
    capabilities = {"confidence_metrics"}

    def run(self, request):
        normalized_dir, structures_dir, system_metrics_path, chain_metrics_path, manifest_path = (
            _write_normalized_bundle(
                request.repeat_dir,
                system_name=self.system_name,
                include_confidence=True,
                include_affinity=False,
                capabilities=self.capabilities,
                runtime_context={"cache_path": "~/.boltz", "diffusion_samples": 1},
            )
        )
        return RunnerExecutionResult(
            runner_name="boltz2",
            raw_output_dir=request.repeat_dir,
            normalized_dir=normalized_dir,
            structures_dir=structures_dir,
            system_metrics_path=system_metrics_path,
            chain_metrics_path=chain_metrics_path,
            manifest_path=manifest_path,
            diffusion_samples=1,
            capabilities=set(self.capabilities),
            runtime=RunnerRuntime(cache_path="~/.boltz", diffusion_samples=1, model_name="boltz2"),
        )


class _MissingConfidencePayloadRunner(_FakeRunner):
    capabilities = {"confidence_metrics"}

    def run(self, request):
        normalized_dir, structures_dir, system_metrics_path, chain_metrics_path, manifest_path = (
            _write_normalized_bundle(
                request.repeat_dir,
                system_name=self.system_name,
                include_confidence=False,
                include_affinity=False,
                capabilities=self.capabilities,
                runtime_context={"diffusion_samples": 1},
            )
        )
        return RunnerExecutionResult(
            runner_name="boltz2",
            raw_output_dir=request.repeat_dir,
            normalized_dir=normalized_dir,
            structures_dir=structures_dir,
            system_metrics_path=system_metrics_path,
            chain_metrics_path=chain_metrics_path,
            manifest_path=manifest_path,
            capabilities=set(self.capabilities),
            runtime=RunnerRuntime(diffusion_samples=1, cache_path="~/.boltz", model_name="boltz2"),
            metric_outcomes={
                "confidence_metrics": RunnerMetricOutcome(state="missing"),
            },
        )


class _MixedOutcomeRunner(_FakeRunner):
    capabilities = {"confidence_metrics", "affinity_metrics", "affinity_metrics_ext"}

    def run(self, request):
        normalized_dir, structures_dir, system_metrics_path, chain_metrics_path, manifest_path = (
            _write_normalized_bundle(
                request.repeat_dir,
                system_name=self.system_name,
                include_confidence=True,
                include_affinity=False,
                capabilities=self.capabilities,
                runtime_context={"cache_path": "~/.boltz", "diffusion_samples": 1},
            )
        )

        return RunnerExecutionResult(
            runner_name="boltz2",
            raw_output_dir=request.repeat_dir,
            normalized_dir=normalized_dir,
            structures_dir=structures_dir,
            system_metrics_path=system_metrics_path,
            chain_metrics_path=chain_metrics_path,
            manifest_path=manifest_path,
            diffusion_samples=1,
            capabilities=set(self.capabilities),
            runtime=RunnerRuntime(cache_path="~/.boltz", diffusion_samples=1, model_name="boltz2"),
            metric_outcomes={
                "confidence_metrics": RunnerMetricOutcome(state="computed"),
                "affinity_metrics": RunnerMetricOutcome(state="missing"),
                "affinity_metrics_ext": RunnerMetricOutcome(state="missing"),
            },
        )


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
            [{"runner": "boltz2", "runtime_context": {"diffusion_samples": 1}}],
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


def _patch_validate_pipeline_real_boltz1(monkeypatch, system_name: str = "system"):
    monkeypatch.setattr(
        "cofolder.recipes.validate.helpers.get_seeds",
        lambda repeats, seed, logger: (123, [123]),
    )
    runner = Boltz1Runner()
    runner.ensure_available = lambda: None
    monkeypatch.setattr(
        "cofolder.recipes.validate.get_runner",
        lambda name: runner,
    )
    monkeypatch.setattr(
        "cofolder.recipes.validate.scaffold_reproduction_metrics",
        lambda system_df, chain_df, reference_path, wrk_dir, pocket_coverage_reference, reproduction_metrics, logger: (
            system_df,
            chain_df,
        ),
    )

    def _fake_run_boltz(cmd, check=True, timings=None, label_prefix=None):
        out_dir = Path(cmd[cmd.index("--out_dir") + 1])
        system_path = Path(cmd[2])
        system_name = system_path.stem
        prediction_dir = out_dir / f"boltz_results_{system_name}" / "predictions" / system_name
        prediction_dir.mkdir(parents=True, exist_ok=True)
        (prediction_dir / f"{system_name}_model_0.cif").write_text("data_test", encoding="utf-8")
        (prediction_dir / f"confidence_{system_name}_model_0.json").write_text(
            json.dumps(
                {
                    "ptm": 0.8,
                    "iptm": 0.7,
                    "confidence_score": 0.9,
                    "chains_ptm": {"0": 0.85, "1": 0.65},
                    "pair_chains_iptm": {"0": {"1": 0.55}, "1": {"0": 0.55}},
                }
            ),
            encoding="utf-8",
        )

    monkeypatch.setattr("cofolder.modules.runners.boltz_runner.run_boltz", _fake_run_boltz)


def _patch_validate_pipeline_real_boltz2(monkeypatch, temp_dir: Path):
    monkeypatch.setattr(
        "cofolder.recipes.validate.helpers.get_seeds",
        lambda repeats, seed, logger: (123, [123]),
    )
    runner = Boltz2Runner()
    runner.ensure_available = lambda: None
    monkeypatch.setattr(
        "cofolder.recipes.validate.get_runner",
        lambda name: runner,
    )
    monkeypatch.setattr(
        "cofolder.recipes.validate.gather.gather_structures",
        lambda base_dir, system_name, repeats, logger: None,
    )
    monkeypatch.setattr(
        "cofolder.recipes.validate.gather.merge_runner_results",
        lambda raw_dir, repeats, logger: _make_runner_results(),
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

    def _fake_run_boltz(cmd, check=True, timings=None, label_prefix=None):
        out_dir = Path(cmd[cmd.index("--out_dir") + 1])
        system_path = Path(cmd[2])
        system_name = system_path.stem
        prediction_dir = out_dir / f"boltz_results_{system_name}" / "predictions" / system_name
        prediction_dir.mkdir(parents=True, exist_ok=True)
        (prediction_dir / f"{system_name}_model_0.cif").write_text("data_test", encoding="utf-8")
        (prediction_dir / f"confidence_{system_name}_model_0.json").write_text(
            json.dumps(
                {
                    "ptm": 0.8,
                    "iptm": 0.7,
                    "confidence_score": 0.9,
                    "chains_ptm": {"0": 0.85, "1": 0.65},
                    "pair_chains_iptm": {"0": {"1": 0.55}, "1": {"0": 0.55}},
                }
            ),
            encoding="utf-8",
        )
        (prediction_dir / f"affinity_{system_name}.json").write_text(
            json.dumps(
                {
                    "affinity_pred_value": 6.1,
                    "affinity_probability_binary": 0.8,
                }
            ),
            encoding="utf-8",
        )

    monkeypatch.setattr("cofolder.modules.runners.boltz_runner.run_boltz", _fake_run_boltz)


def _patch_validate_pipeline_real_boltz_community(
    monkeypatch,
    *,
    include_affinity: bool = True,
):
    monkeypatch.setattr(
        "cofolder.recipes.validate.helpers.get_seeds",
        lambda repeats, seed, logger: (123, [123]),
    )
    runner = BoltzCommunityRunner()
    runner.ensure_available = lambda: None
    monkeypatch.setattr(
        "cofolder.recipes.validate.get_runner",
        lambda name: runner,
    )
    monkeypatch.setattr(
        "cofolder.recipes.validate.scaffold_reproduction_metrics",
        lambda system_df, chain_df, reference_path, wrk_dir, pocket_coverage_reference, reproduction_metrics, logger: (
            system_df,
            chain_df,
        ),
    )

    def _fake_run_boltz(cmd, check=True, timings=None, label_prefix=None):
        out_dir = Path(cmd[cmd.index("--out_dir") + 1])
        system_path = Path(cmd[2])
        system_name = system_path.stem
        prediction_dir = out_dir / f"boltz_results_{system_name}" / "predictions" / system_name
        prediction_dir.mkdir(parents=True, exist_ok=True)
        (prediction_dir / f"{system_name}_model_0.cif").write_text("data_test", encoding="utf-8")
        (prediction_dir / f"confidence_{system_name}_model_0.json").write_text(
            json.dumps(
                {
                    "ptm": 0.8,
                    "iptm": 0.7,
                    "confidence_score": 0.9,
                    "chains_ptm": {"0": 0.85, "1": 0.65},
                    "pair_chains_iptm": {"0": {"1": 0.55}, "1": {"0": 0.55}},
                }
            ),
            encoding="utf-8",
        )
        if include_affinity:
            (prediction_dir / f"affinity_{system_name}.json").write_text(
                json.dumps(
                    {
                        "affinity_pred_value": 6.1,
                        "affinity_probability_binary": 0.8,
                    }
                ),
                encoding="utf-8",
            )

    monkeypatch.setattr("cofolder.modules.runners.boltz_runner.run_boltz", _fake_run_boltz)


class TestValidateInit:
    """Tests for Validate initialization."""

    def test_init_checks_runner_availability_before_loading_inputs(
        self,
        sample_system_yaml,
        sample_options_yaml,
        temp_dir,
    ):
        fake_runner = Mock()
        fake_runner.ensure_available.side_effect = RuntimeError("install boltz first")

        with patch("cofolder.recipes.validate.get_runner", return_value=fake_runner):
            with patch("cofolder.recipes.validate.read.read_yaml") as mock_read_yaml:
                with pytest.raises(RuntimeError, match="install boltz first"):
                    Validate(
                        wrk_dir=str(temp_dir),
                        system_path=str(sample_system_yaml),
                        options_path=str(sample_options_yaml),
                        scoring_functions=[],
                    )

        mock_read_yaml.assert_not_called()

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
        assert validator.runner_name == "boltz2"
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

    def test_reference_runner_contract_preserves_raw_and_results_layout(
        self,
        monkeypatch,
        sample_system_yaml,
        sample_options_yaml,
        temp_dir,
    ):
        _patch_validate_pipeline(monkeypatch)

        validator = Validate(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            scoring_functions=["confidence_metrics", "affinity_metrics", "affinity_metrics_ext"],
        )

        validator.run()

        raw_repeat_dir = Path(temp_dir) / "raw" / "repeat_1"
        normalized_dir = raw_repeat_dir / "normalized"
        assert raw_repeat_dir.exists()
        assert (normalized_dir / "system_metrics.csv").exists()
        assert (normalized_dir / "chain_metrics.csv").exists()
        assert (normalized_dir / "manifest.json").exists()
        assert (normalized_dir / "structures" / "1_system_model_0.cif").exists()
        assert (Path(temp_dir) / "results" / "system_metrics.csv").exists()
        assert (Path(temp_dir) / "results" / "chain_metrics.csv").exists()

    def test_reference_runner_validate_uses_real_boltz2_path(
        self,
        monkeypatch,
        sample_options_yaml,
        temp_dir,
    ):
        _patch_validate_pipeline_real_boltz2(monkeypatch, Path(temp_dir))
        system_path = Path(temp_dir) / "affinity_system.yaml"
        system_path.write_text(
            yaml.safe_dump(
                {
                    "sequences": [
                        {"protein": {"id": "A", "fasta": "MKRAAT"}},
                        {"ligand": {"id": "B", "smiles": "CCO", "ccd": "ETH"}},
                    ],
                    "properties": [{"affinity": {"binder": "B"}}],
                }
            ),
            encoding="utf-8",
        )

        validator = Validate(
            wrk_dir=str(temp_dir),
            system_path=str(system_path),
            options_path=str(sample_options_yaml),
            scoring_functions=["confidence_metrics", "affinity_metrics", "affinity_metrics_ext"],
        )

        validator.run()

        raw_repeat_dir = Path(temp_dir) / "raw" / "repeat_1"
        normalized_dir = raw_repeat_dir / "normalized"
        raw_system_df = pd.read_csv(normalized_dir / "system_metrics.csv")
        raw_chain_df = pd.read_csv(normalized_dir / "chain_metrics.csv")
        assert validator.runner.__class__ is Boltz2Runner
        assert (normalized_dir / "manifest.json").exists()
        assert {"ptm", "iptm", "confidence_score"}.issubset(raw_system_df.columns)
        assert {
            "affinity_pred_value",
            "affinity_probability_binary",
            "pIC50",
            "IC50_M",
            "pIC50_kcal_per_mol",
        }.issubset(raw_chain_df.columns)
        assert (Path(temp_dir) / "results" / "system_metrics.csv").exists()
        assert (Path(temp_dir) / "results" / "chain_metrics.csv").exists()

    def test_reference_runner_mixed_metric_outcomes_allow_unrequested_missing_groups(
        self,
        monkeypatch,
        sample_system_yaml,
        sample_options_yaml,
        temp_dir,
    ):
        monkeypatch.setattr(
            "cofolder.recipes.validate.helpers.get_seeds",
            lambda repeats, seed, logger: (123, [123]),
        )
        monkeypatch.setattr(
            "cofolder.recipes.validate.get_runner",
            lambda name: _MixedOutcomeRunner(),
        )
        monkeypatch.setattr(
            "cofolder.recipes.validate.gather.gather_structures",
            lambda base_dir, system_name, repeats, logger: None,
        )
        monkeypatch.setattr(
            "cofolder.recipes.validate.gather.merge_runner_results",
            lambda raw_dir, repeats, logger: _make_runner_results(),
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

        validator = Validate(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            scoring_functions=["confidence_metrics"],
        )

        validator.run()

        raw_chain_df = pd.read_csv(Path(temp_dir) / "raw" / "repeat_1" / "normalized" / "chain_metrics.csv")
        system_df = pd.read_csv(Path(temp_dir) / "results" / "system_metrics.csv")
        assert {"ptm", "iptm", "confidence_score"}.issubset(system_df.columns)
        assert "affinity_pred_value" not in raw_chain_df.columns
        assert (Path(temp_dir) / "results" / "chain_metrics.csv").exists()

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

        raw_system_df = pd.read_csv(Path(temp_dir) / "raw" / "repeat_1" / "normalized" / "system_metrics.csv")
        raw_chain_df = pd.read_csv(Path(temp_dir) / "raw" / "repeat_1" / "normalized" / "chain_metrics.csv")
        system_df = pd.read_csv(Path(temp_dir) / "results" / "system_metrics.csv")
        chain_df = pd.read_csv(Path(temp_dir) / "results" / "chain_metrics.csv")

        assert {"ptm", "iptm", "confidence_score"}.isdisjoint(raw_system_df.columns)
        assert {
            "affinity_pred_value",
            "affinity_probability_binary",
            "pIC50",
            "IC50_M",
            "pIC50_kcal_per_mol",
        }.isdisjoint(raw_chain_df.columns)
        assert {"ptm", "iptm", "confidence_score"}.issubset(system_df.columns)
        assert {"affinity_pred_value", "affinity_probability_binary", "pIC50", "IC50_M", "pIC50_kcal_per_mol"}.issubset(chain_df.columns)
        assert system_df[["ptm", "iptm", "confidence_score"]].isna().all().all()
        assert chain_df[
            [
                "affinity_pred_value",
                "affinity_probability_binary",
                "pIC50",
                "IC50_M",
                "pIC50_kcal_per_mol",
            ]
        ].isna().all().all()

    def test_boltz1_validate_warns_for_unsupported_affinity_groups_and_keeps_schemas_stable(
        self,
        monkeypatch,
        sample_options_yaml,
        temp_dir,
        caplog,
    ):
        _patch_validate_pipeline_real_boltz1(monkeypatch, system_name="affinity_system")
        system_path = Path(temp_dir) / "affinity_system.yaml"
        system_path.write_text(
            yaml.safe_dump(
                {
                    "sequences": [
                        {"protein": {"id": "A", "fasta": "MKRAAT"}},
                        {"ligand": {"id": "B", "smiles": "CCO", "ccd": "ETH"}},
                    ],
                    "properties": [{"affinity": {"binder": "B"}}],
                }
            ),
            encoding="utf-8",
        )

        validator = Validate(
            wrk_dir=str(temp_dir),
            system_path=str(system_path),
            options_path=str(sample_options_yaml),
            runner="boltz1",
            scoring_functions=["confidence_metrics", "affinity_metrics", "affinity_metrics_ext"],
        )

        with caplog.at_level(logging.WARNING, logger="cofolder.recipes.validate"):
            validator.run()

        raw_repeat_dir = Path(temp_dir) / "raw" / "repeat_1"
        normalized_dir = raw_repeat_dir / "normalized"
        raw_system_df = pd.read_csv(normalized_dir / "system_metrics.csv")
        raw_chain_df = pd.read_csv(normalized_dir / "chain_metrics.csv")
        system_df = pd.read_csv(Path(temp_dir) / "results" / "system_metrics.csv")
        chain_df = pd.read_csv(Path(temp_dir) / "results" / "chain_metrics.csv")

        assert validator.runner.__class__ is Boltz1Runner
        assert "does not support scoring groups" in caplog.text
        assert (normalized_dir / "manifest.json").exists()
        assert (Path(temp_dir) / "results" / "structures" / "1_affinity_system_model_0.cif").exists()
        assert {"ptm", "iptm", "confidence_score"}.issubset(raw_system_df.columns)
        assert {
            "affinity_pred_value",
            "affinity_probability_binary",
            "pIC50",
            "IC50_M",
            "pIC50_kcal_per_mol",
        }.isdisjoint(raw_chain_df.columns)
        assert {"CHAIN_ID", "ENTITY_TYPE", "ligand_molecule_id"}.issubset(chain_df.columns)
        assert {"ptm", "iptm", "confidence_score"}.issubset(system_df.columns)
        assert {"affinity_pred_value", "affinity_probability_binary", "pIC50", "IC50_M", "pIC50_kcal_per_mol"}.issubset(chain_df.columns)
        assert chain_df[
            [
                "affinity_pred_value",
                "affinity_probability_binary",
                "pIC50",
                "IC50_M",
                "pIC50_kcal_per_mol",
            ]
        ].isna().all().all()

    def test_boltz_community_validate_uses_real_path_and_preserves_affinity_outputs(
        self,
        monkeypatch,
        sample_options_yaml,
        temp_dir,
    ):
        _patch_validate_pipeline_real_boltz_community(monkeypatch, include_affinity=True)
        system_path = Path(temp_dir) / "affinity_system.yaml"
        system_path.write_text(
            yaml.safe_dump(
                {
                    "sequences": [
                        {"protein": {"id": "A", "fasta": "MKRAAT"}},
                        {"ligand": {"id": "B", "smiles": "CCO", "ccd": "ETH"}},
                    ],
                    "properties": [{"affinity": {"binder": "B"}}],
                }
            ),
            encoding="utf-8",
        )

        validator = Validate(
            wrk_dir=str(temp_dir),
            system_path=str(system_path),
            options_path=str(sample_options_yaml),
            runner="boltz-community",
            scoring_functions=["confidence_metrics", "affinity_metrics", "affinity_metrics_ext"],
        )

        validator.run()

        raw_repeat_dir = Path(temp_dir) / "raw" / "repeat_1"
        normalized_dir = raw_repeat_dir / "normalized"
        raw_system_df = pd.read_csv(normalized_dir / "system_metrics.csv")
        raw_chain_df = pd.read_csv(normalized_dir / "chain_metrics.csv")
        system_df = pd.read_csv(Path(temp_dir) / "results" / "system_metrics.csv")
        chain_df = pd.read_csv(Path(temp_dir) / "results" / "chain_metrics.csv")

        assert validator.runner.__class__ is BoltzCommunityRunner
        assert (normalized_dir / "manifest.json").exists()
        assert (Path(temp_dir) / "results" / "structures" / "1_affinity_system_model_0.cif").exists()
        assert {"ptm", "iptm", "confidence_score"}.issubset(raw_system_df.columns)
        assert {
            "affinity_pred_value",
            "affinity_probability_binary",
            "pIC50",
            "IC50_M",
            "pIC50_kcal_per_mol",
        }.issubset(raw_chain_df.columns)
        assert {"CHAIN_ID", "ENTITY_TYPE", "ligand_molecule_id"}.issubset(chain_df.columns)
        assert {"ptm", "iptm", "confidence_score"}.issubset(system_df.columns)
        assert {"affinity_pred_value", "affinity_probability_binary", "pIC50", "IC50_M", "pIC50_kcal_per_mol"}.issubset(chain_df.columns)

    def test_boltz_community_validate_fails_before_gather_on_missing_requested_affinity_payload(
        self,
        monkeypatch,
        sample_options_yaml,
        temp_dir,
    ):
        _patch_validate_pipeline_real_boltz_community(monkeypatch, include_affinity=False)
        system_path = Path(temp_dir) / "affinity_system.yaml"
        system_path.write_text(
            yaml.safe_dump(
                {
                    "sequences": [
                        {"protein": {"id": "A", "fasta": "MKRAAT"}},
                        {"ligand": {"id": "B", "smiles": "CCO", "ccd": "ETH"}},
                    ],
                    "properties": [{"affinity": {"binder": "B"}}],
                }
            ),
            encoding="utf-8",
        )

        called = {"gather_structures": False}

        def _fake_gather_structures(base_dir, system_name, repeats, logger):
            called["gather_structures"] = True

        monkeypatch.setattr(
            "cofolder.recipes.validate.gather.gather_structures",
            _fake_gather_structures,
        )

        validator = Validate(
            wrk_dir=str(temp_dir),
            system_path=str(system_path),
            options_path=str(sample_options_yaml),
            runner="boltz-community",
            scoring_functions=["confidence_metrics", "affinity_metrics", "affinity_metrics_ext"],
        )

        with pytest.raises(RunnerBundleValidationError, match="cannot satisfy requested metric groups"):
            validator.run()

        assert called["gather_structures"] is False

    def test_validate_ignores_manifest_runtime_context_for_robustness(
        self,
        monkeypatch,
        sample_system_yaml,
        sample_options_yaml,
        temp_dir,
    ):
        monkeypatch.setattr(
            "cofolder.recipes.validate.helpers.get_seeds",
            lambda repeats, seed, logger: (123, [123]),
        )
        monkeypatch.setattr(
            "cofolder.recipes.validate.get_runner",
            lambda name: _ManifestOnlyRuntimeRunner(),
        )
        monkeypatch.setattr(
            "cofolder.recipes.validate.gather.gather_structures",
            lambda base_dir, system_name, repeats, logger: None,
        )
        monkeypatch.setattr(
            "cofolder.recipes.validate.gather.merge_runner_results",
            lambda raw_dir, repeats, logger: (
                pd.DataFrame([{"cif_file": "1_system_model_0.cif", "model_name": "system", "repeat": 1}]),
                pd.DataFrame(),
                [{"runner": "boltz2", "runtime_context": {"diffusion_samples": 9}}],
            ),
        )
        monkeypatch.setattr(
            "cofolder.recipes.validate.scaffold_reproduction_metrics",
            lambda system_df, chain_df, reference_path, wrk_dir, pocket_coverage_reference, reproduction_metrics, logger: (
                system_df,
                chain_df,
            ),
        )

        called = {"robustness": False}

        def _fake_gather_robustness_results(system_df, chain_df, wrk_dir, reference_path):
            called["robustness"] = True
            return pd.DataFrame([{"repeat": 1, "diffusion_sample": 0}])

        monkeypatch.setattr(
            "cofolder.recipes.validate.gather.gather_robustness_results",
            _fake_gather_robustness_results,
        )

        validator = Validate(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            scoring_functions=[],
            assess_robustness=True,
        )

        validator.run()

        assert called["robustness"] is False

    def test_typed_runtime_controls_robustness_without_manifest_fallback(
        self,
        monkeypatch,
        sample_system_yaml,
        sample_options_yaml,
        temp_dir,
    ):
        monkeypatch.setattr(
            "cofolder.recipes.validate.helpers.get_seeds",
            lambda repeats, seed, logger: (123, [123]),
        )
        monkeypatch.setattr(
            "cofolder.recipes.validate.get_runner",
            lambda name: _TypedDiffusionRuntimeRunner(),
        )
        monkeypatch.setattr(
            "cofolder.recipes.validate.gather.gather_structures",
            lambda base_dir, system_name, repeats, logger: None,
        )
        monkeypatch.setattr(
            "cofolder.recipes.validate.gather.merge_runner_results",
            lambda raw_dir, repeats, logger: (
                pd.DataFrame([{"cif_file": "1_system_model_0.cif", "model_name": "system", "repeat": 1}]),
                pd.DataFrame(),
                [{"runner": "boltz2", "runtime_context": {"diffusion_samples": 1}}],
            ),
        )
        monkeypatch.setattr(
            "cofolder.recipes.validate.scaffold_reproduction_metrics",
            lambda system_df, chain_df, reference_path, wrk_dir, pocket_coverage_reference, reproduction_metrics, logger: (
                system_df,
                chain_df,
            ),
        )

        called = {"robustness": False}

        def _fake_gather_robustness_results(system_df, chain_df, wrk_dir, reference_path):
            called["robustness"] = True
            return pd.DataFrame([{"repeat": 1, "diffusion_sample": 0}])

        monkeypatch.setattr(
            "cofolder.recipes.validate.gather.gather_robustness_results",
            _fake_gather_robustness_results,
        )

        validator = Validate(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            scoring_functions=[],
            assess_robustness=True,
        )

        validator.run()

        assert called["robustness"] is True

    def test_validate_rejects_requested_metric_group_without_explicit_outcome(
        self,
        monkeypatch,
        sample_system_yaml,
        sample_options_yaml,
        temp_dir,
    ):
        monkeypatch.setattr(
            "cofolder.recipes.validate.helpers.get_seeds",
            lambda repeats, seed, logger: (123, [123]),
        )
        monkeypatch.setattr(
            "cofolder.recipes.validate.get_runner",
            lambda name: _ImplicitOutcomeRunner(),
        )

        called = {"gather_structures": False}
        monkeypatch.setattr(
            "cofolder.recipes.validate.gather.gather_structures",
            lambda base_dir, system_name, repeats, logger: called.__setitem__("gather_structures", True),
        )
        monkeypatch.setattr(
            "cofolder.recipes.validate.gather.merge_runner_results",
            lambda raw_dir, repeats, logger: (
                pd.DataFrame([{"cif_file": "1_system_model_0.cif", "model_name": "system", "repeat": 1}]),
                pd.DataFrame(),
                [{"runner": "boltz2", "runtime_context": {"diffusion_samples": 3, "cache_path": "/legacy/cache"}}],
            ),
        )

        validator = Validate(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            scoring_functions=["confidence_metrics"],
        )

        with pytest.raises(RunnerBundleValidationError, match="omits an explicit outcome"):
            validator.run()

        assert called["gather_structures"] is False

    def test_invalid_supported_metric_bundle_fails_before_gather(
        self,
        monkeypatch,
        sample_system_yaml,
        sample_options_yaml,
        temp_dir,
    ):
        monkeypatch.setattr(
            "cofolder.recipes.validate.helpers.get_seeds",
            lambda repeats, seed, logger: (123, [123]),
        )
        monkeypatch.setattr(
            "cofolder.recipes.validate.get_runner",
            lambda name: _MissingConfidencePayloadRunner(),
        )

        called = {"gather_structures": False}

        def _fake_gather_structures(base_dir, system_name, repeats, logger):
            called["gather_structures"] = True

        monkeypatch.setattr(
            "cofolder.recipes.validate.gather.gather_structures",
            _fake_gather_structures,
        )
        monkeypatch.setattr(
            "cofolder.recipes.validate.gather.merge_runner_results",
            lambda raw_dir, repeats, logger: _make_runner_results(),
        )
        monkeypatch.setattr(
            "cofolder.recipes.validate.scaffold_reproduction_metrics",
            lambda system_df, chain_df, reference_path, wrk_dir, pocket_coverage_reference, reproduction_metrics, logger: (
                system_df,
                chain_df,
            ),
        )

        validator = Validate(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            scoring_functions=["confidence_metrics"],
        )

        with pytest.raises(RunnerBundleValidationError, match="cannot satisfy requested metric groups"):
            validator.run()

        assert called["gather_structures"] is False
