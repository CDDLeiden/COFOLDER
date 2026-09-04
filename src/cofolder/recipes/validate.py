from __future__ import annotations

import copy
import logging
import os
from pathlib import Path

from cofolder.modules.analytics.reproduction import scaffold_reproduction_metrics
from cofolder.modules.analytics.structure import Structure
from cofolder.modules.input import system
from cofolder.modules.runners import (
    RunnerExecutionRequest,
    RunnerRuntime,
    get_runner,
    merge_runner_runtime,
)
from cofolder.modules.runners.msa import (
    msa_cache_key,
    resolve_declared_msa_paths,
    unresolved_protein_sequences,
)
from cofolder.modules.runners.validators import validate_runner_bundle
from cofolder.modules.utils import gather, helpers, read, write
from cofolder.modules.utils.timing import DebugTimingCollector
from cofolder.recipes.bias import BiasAssessmentWorkflow

logger = logging.getLogger(__name__)

DEFAULT_SCORING_FUNCTIONS = {
    "confidence_metrics",
    "affinity_metrics",
    "affinity_metrics_ext",
    "ifp_distance",
    "sasa",
    "sasa_normalized",
}
DEFAULT_REPRODUCTION_METRICS = {
    "protein_rmsd",
    "ligand_rmsd",
    "sucos",
    "pocket_coverage",
}
RUNNER_METRIC_GROUPS = {
    "confidence_metrics",
    "affinity_metrics",
    "affinity_metrics_ext",
}
AFFINITY_METRIC_GROUPS = {
    "affinity_metrics",
    "affinity_metrics_ext",
}


class Validate(object):
    """High-level orchestrator for validation workflow."""

    def __init__(
        self,
        wrk_dir: str,
        system_path: str,
        options_path: str,
        runner: str = "boltz2",
        repeats: int = 1,
        seed: int | None = None,
        scoring_functions: list[str] | None = None,
        assess_robustness: bool = True,
        assess_bias: bool = False,
        protein_training_data_path: str | None = None,
        ligand_training_data_path: str | None = None,
        bias_release_cutoff: str = "2023-06-01",
        bias_ligand_similarity_threshold: float = 0.35,
        bias_chains: list[str] | None = None,
        build_bias_training_data: bool = False,
        bias_training_components_cif: str | None = None,
        conformers: str | None = None,
        sdf_file: str | None = None,
        reference_path: str | None = None,
        pocket_coverage_reference: str | None = None,
        reproduction_metrics: list[str] | None = None,
        reusable_msa_dir: str | None = None,
    ):
        self.wrk_dir = Path(wrk_dir)
        self.system_path = Path(system_path)
        self.options_path = Path(options_path)
        self.runner_name = str(runner)
        self.repeats = repeats
        self.seed = seed
        self.assess_robustness = assess_robustness
        self.assess_bias = assess_bias
        self.protein_training_data_path = (
            Path(protein_training_data_path) if protein_training_data_path else None
        )
        self.ligand_training_data_path = (
            Path(ligand_training_data_path) if ligand_training_data_path else None
        )
        self.bias_release_cutoff = bias_release_cutoff
        self.bias_ligand_similarity_threshold = float(bias_ligand_similarity_threshold)
        self.build_bias_training_data = build_bias_training_data
        self.bias_training_components_cif = (
            Path(bias_training_components_cif) if bias_training_components_cif else None
        )
        normalized_bias_chains: set[str] = set()
        for value in bias_chains or []:
            for token in str(value).split(","):
                token = token.strip()
                if token:
                    normalized_bias_chains.add(token.upper())
        self.bias_chains = normalized_bias_chains or None
        self.conformers = conformers
        self.sdf_file = Path(sdf_file) if sdf_file else None
        self.reference_path = Path(reference_path) if reference_path else None
        self.pocket_coverage_reference = pocket_coverage_reference
        self.reproduction_metrics = set(
            reproduction_metrics or DEFAULT_REPRODUCTION_METRICS
        )
        self.reusable_msa_dir = Path(reusable_msa_dir) if reusable_msa_dir else None

        if self.reference_path is not None:
            if not self.reference_path.exists():
                raise ValueError(
                    f"reference_path does not exist: {self.reference_path}"
                )
            if not self.reference_path.is_file():
                raise ValueError(f"reference_path is not a file: {self.reference_path}")

        self.logger = logging.getLogger(__name__)
        self.explicit_scoring_functions = set(scoring_functions or ())
        self.scoring_functions = (
            set(scoring_functions)
            if scoring_functions is not None
            else set(DEFAULT_SCORING_FUNCTIONS)
        )

        self.runner = get_runner(self.runner_name)
        self.runner.ensure_available()

        self._system = read.read_yaml(path=self.system_path)
        self.base_system = system.System(system=self._system)
        resolve_declared_msa_paths(
            self.base_system,
            base_dir=self.system_path.parent,
        )

        self.logger.debug(
            "Initializing Validate with parameters: %s",
            {
                "wrk_dir": self.wrk_dir,
                "system_path": self.system_path,
                "options_path": self.options_path,
                "runner": self.runner_name,
                "repeats": self.repeats,
                "seed": self.seed,
                "conformers": self.conformers,
                "sdf_file": self.sdf_file,
                "assess_bias": self.assess_bias,
                "protein_training_data_path": self.protein_training_data_path,
                "ligand_training_data_path": self.ligand_training_data_path,
                "bias_release_cutoff": self.bias_release_cutoff,
                "bias_ligand_similarity_threshold": self.bias_ligand_similarity_threshold,
                "bias_chains": sorted(self.bias_chains) if self.bias_chains else None,
                "build_bias_training_data": self.build_bias_training_data,
                "bias_training_components_cif": self.bias_training_components_cif,
                "reference_path": self.reference_path,
                "pocket_coverage_reference": self.pocket_coverage_reference,
                "reproduction_metrics": sorted(self.reproduction_metrics),
                "scoring_functions": sorted(self.scoring_functions),
                "reusable_msa_dir": self.reusable_msa_dir,
            },
        )

        self.timings = DebugTimingCollector(logger=self.logger)

    def _log_timing_summary(self) -> None:
        self.timings.log_summary(logger=self.logger)

    def _debug_timer(self, label: str):
        return self.timings.measure(label, logger=self.logger)

    def run(self):
        with self._debug_timer("validate.total"):
            self.wrk_dir.mkdir(parents=True, exist_ok=True)
            self.raw_dir = self.wrk_dir / "raw"
            self.raw_dir.mkdir(parents=True, exist_ok=True)
            self.logger.debug("Using raw directory for outputs: %s", self.raw_dir)

            unsupported_metric_groups = sorted(
                (self.scoring_functions & RUNNER_METRIC_GROUPS)
                - set(self.runner.capabilities)
            )
            if unsupported_metric_groups:
                self.logger.warning(
                    "Runner '%s' does not support scoring groups %s. "
                    "The run will continue and the related output columns will be empty.",
                    self.runner_name,
                    unsupported_metric_groups,
                )

            with self._debug_timer("seeds.resolve"):
                self.seed, self.run_seeds = helpers.get_seeds(
                    repeats=self.repeats,
                    seed=self.seed,
                    logger=self.logger,
                )

            with self._debug_timer("runner.options.load"):
                runner_options = self.runner.load_options(self.options_path)
            resolve_msa_reuse_settings = getattr(
                self.runner, "msa_reuse_settings", None
            )
            msa_reuse_settings = (
                resolve_msa_reuse_settings(runner_options)
                if callable(resolve_msa_reuse_settings)
                else {}
            )

            self.sys = system.System(system=copy.deepcopy(self.base_system.system))
            if self.reusable_msa_dir is not None and getattr(
                self.runner, "supports_msa_reuse", False
            ):
                injected = self.runner.inject_reusable_msas(
                    self.sys,
                    self.reusable_msa_dir,
                    settings=msa_reuse_settings,
                )
                if injected:
                    self.logger.info(
                        "Injected %d reusable protein MSA(s) from %s.",
                        injected,
                        self.reusable_msa_dir,
                    )
            with self._debug_timer("runner.system.validate.preparation_input"):
                self.runner.validate_system(
                    self.sys,
                    runner_options,
                    check_atom_names=False,
                )
            with self._debug_timer("runner.prepare_system"):
                preparation = self.runner.prepare_system(
                    # keep backend-specific prep behind the runner boundary
                    system_obj=self.sys,
                    options_obj=runner_options,
                    wrk_dir=self.raw_dir,
                    conformers=self.conformers,
                    sdf_file=self.sdf_file,
                    logger=self.logger,
                )
            preparation_runtime = preparation.runtime
            self.sys = preparation.system_obj
            runner_options = preparation.options_obj
            with self._debug_timer("runner.system.validate.execution_input"):
                self.runner.validate_system(
                    self.sys,
                    runner_options,
                    check_atom_names=True,
                )
            for warning in preparation.warnings:
                self.logger.warning("%s", warning)

            yaml_path = os.path.join(self.raw_dir, str(self.system_path.name))
            with self._debug_timer("system_yaml.write"):
                write.write_yaml(self.sys, path=yaml_path)
            self.logger.info("Updated system saved to YAML: %s", yaml_path)

            runner_results = []
            for i, seed in enumerate(self.run_seeds, 1):
                self.logger.info(
                    "Running repeat %d/%d with seed %d using runner '%s'",
                    i,
                    self.repeats,
                    seed,
                    self.runner_name,
                )
                request = RunnerExecutionRequest(
                    runner_name=self.runner_name,
                    system_name=self.system_path.stem,
                    system_path=Path(yaml_path),
                    system_obj=self.sys,
                    options_path=self.options_path,
                    options_obj=runner_options,
                    repeat=i,
                    seed=seed,
                    repeat_dir=self.raw_dir / f"repeat_{i}",
                    raw_dir=self.raw_dir,
                    logger=self.logger,
                    timings=self.timings,
                    label_prefix=f"repeat_{i}",
                    runtime=preparation_runtime,
                )
                unresolved_before_run = []
                if self.reusable_msa_dir is not None and getattr(
                    self.runner, "supports_msa_reuse", False
                ):
                    unresolved_before_run = unresolved_protein_sequences(self.sys)
                    if unresolved_before_run:
                        self.reusable_msa_dir.mkdir(parents=True, exist_ok=True)
                        attempt_markers = [
                            self.reusable_msa_dir
                            / f".generation_attempted_{msa_cache_key(sequence, msa_reuse_settings)[:16]}"
                            for sequence in unresolved_before_run
                        ]
                        if any(marker.exists() for marker in attempt_markers):
                            raise RuntimeError(
                                "Protein MSA generation was already attempted for this screen, "
                                "but no valid reusable artifact resolved every protein. "
                                "Refusing to call the MSA server again silently."
                            )
                        for marker in attempt_markers:
                            marker.write_text(
                                "MSA generation was requested once by the screen workflow.\n",
                                encoding="utf-8",
                            )
                try:
                    result = self.runner.run(request)
                finally:
                    if self.reusable_msa_dir is not None and getattr(
                        self.runner, "supports_msa_reuse", False
                    ):
                        captured = self.runner.capture_reusable_msas(
                            self.sys,
                            generated_dir=request.repeat_dir,
                            cache_dir=self.reusable_msa_dir,
                            settings=msa_reuse_settings,
                        )
                        injected = self.runner.inject_reusable_msas(
                            self.sys,
                            self.reusable_msa_dir,
                            settings=msa_reuse_settings,
                        )
                        if captured or injected:
                            write.write_yaml(self.sys, path=yaml_path)
                            self.logger.info(
                                "Staged %d and injected %d reusable protein MSA(s); "
                                "later repeats will not call the MSA server.",
                                captured,
                                injected,
                            )
                        if unresolved_before_run and unresolved_protein_sequences(
                            self.sys
                        ):
                            raise RuntimeError(
                                "The runner completed without producing valid reusable MSAs "
                                "for every unresolved protein. Refusing to recalculate them "
                                "on a later repeat or ligand."
                            )
                runner_results.append(result)

            selected_runner_metric_groups = (
                self.scoring_functions & RUNNER_METRIC_GROUPS
            )
            requested_runner_metric_groups = self._required_runner_metric_groups()
            optional_runner_metric_groups = (
                selected_runner_metric_groups - requested_runner_metric_groups
            )
            unavailable_metric_groups = set(unsupported_metric_groups)
            with self._debug_timer("runner.bundle_validation"):
                for repeat, result in enumerate(runner_results, 1):
                    bundle = validate_runner_bundle(
                        result,
                        requested_metric_groups=requested_runner_metric_groups,
                    )
                    for group_name in sorted(optional_runner_metric_groups):
                        outcome = bundle.metric_outcomes.get(group_name)
                        if outcome is None or outcome.state not in {
                            "missing",
                            "failed",
                        }:
                            continue
                        unavailable_metric_groups.add(group_name)
                        reason = (
                            outcome.message
                            or "runner did not produce the normalized payload"
                        )
                        self.logger.warning(
                            "Runner '%s' did not produce optional scoring group '%s' for "
                            "repeat %d: %s The run will continue and the related output "
                            "columns will be empty.",
                            self.runner_name,
                            group_name,
                            repeat,
                            reason,
                        )

            with self._debug_timer("structures.gather"):
                gather.gather_structures(
                    base_dir=self.wrk_dir,
                    system_name=self.system_path.stem,
                    repeats=self.repeats,
                    logger=self.logger,
                )

            with self._debug_timer("results.merge_runner_outputs"):
                system_df, chain_df, _manifests = gather.merge_runner_results(
                    raw_dir=self.raw_dir,
                    repeats=self.repeats,
                    logger=self.logger,
                )

            if not chain_df.empty:
                with self._debug_timer("results.add_chain_info"):
                    chain_df = gather.add_chain_info(chain_df, self.sys)

            system_df, chain_df = self._apply_unavailable_metric_columns(
                system_df=system_df,
                chain_df=chain_df,
                unavailable_metric_groups=unavailable_metric_groups,
            )

            structure_metrics = {
                "ifp_distance",
                "ifp_prolif",
                "sasa",
                "sasa_normalized",
            }
            if self.scoring_functions & structure_metrics and not chain_df.empty:
                with self._debug_timer("structure.init"):
                    structure = Structure(
                        wrk_dir=self.wrk_dir,
                        chain_df=chain_df,
                        cif_folder=Path(self.wrk_dir / "results" / "structures"),
                    )

                if "ifp_distance" in self.scoring_functions:
                    with self._debug_timer("scores.ifp_distance"):
                        chain_df = structure.add_ifp_distance()

                if "ifp_prolif" in self.scoring_functions:
                    with self._debug_timer("scores.ifp_prolif"):
                        chain_df = structure.add_ifp_prolif()

                if {"sasa", "sasa_normalized"} & self.scoring_functions:
                    with self._debug_timer("scores.sasa"):
                        chain_df = structure.add_sasa(
                            absolute="sasa" in self.scoring_functions,
                            normalized="sasa_normalized" in self.scoring_functions,
                        )

            with self._debug_timer("scores.reproduction_metrics"):
                system_df, chain_df = scaffold_reproduction_metrics(
                    system_df=system_df,
                    chain_df=chain_df,
                    reference_path=self.reference_path,
                    wrk_dir=self.wrk_dir,
                    pocket_coverage_reference=self.pocket_coverage_reference,
                    reproduction_metrics=sorted(self.reproduction_metrics),
                    logger=self.logger,
                )

            # Story 2.4 retires recipe-side runtime reconstruction from manifest/runtime_context
            # fallbacks. Shared workflow behavior now depends only on typed runner metadata.
            runtime = merge_runner_runtime(
                preparation_runtime,
                *(result.runtime for result in runner_results),
            )

            if self.assess_bias:
                system_df, chain_df = self._apply_bias_metrics(
                    system_df=system_df,
                    chain_df=chain_df,
                    runtime=runtime,
                )

            with self._debug_timer("results.write.system_chain"):
                write.write_csv(
                    system_df,
                    output_path=self.wrk_dir / "results" / "system_metrics.csv",
                )
                write.write_csv(
                    chain_df, output_path=self.wrk_dir / "results" / "chain_metrics.csv"
                )

            diffusion_samples = self._resolve_diffusion_samples(system_df, runtime)
            if self.assess_robustness and (self.repeats > 1 or diffusion_samples > 1):
                with self._debug_timer("scores.robustness_metrics"):
                    results_df = gather.gather_robustness_results(
                        system_df=system_df,
                        chain_df=chain_df,
                        wrk_dir=self.wrk_dir,
                        reference_path=self.reference_path,
                    )
                    write.write_csv(
                        results_df,
                        output_path=self.wrk_dir / "results" / "robustness_metrics.csv",
                    )

        self._log_timing_summary()

    def _required_runner_metric_groups(self) -> set[str]:
        selected = self.scoring_functions & RUNNER_METRIC_GROUPS
        required = selected - AFFINITY_METRIC_GROUPS
        affinity_is_activated = self._system_declares_affinity()
        for group_name in selected & AFFINITY_METRIC_GROUPS:
            if affinity_is_activated or group_name in self.explicit_scoring_functions:
                required.add(group_name)
        return required

    def _system_declares_affinity(self) -> bool:
        properties = self.sys.find_value(key="properties") or []
        if not isinstance(properties, list):
            return False
        return any(
            isinstance(property_entry, dict) and "affinity" in property_entry
            for property_entry in properties
        )

    def _apply_unavailable_metric_columns(
        self,
        system_df,
        chain_df,
        unavailable_metric_groups: set[str],
    ):
        if "confidence_metrics" in unavailable_metric_groups:
            for column in ("ptm", "iptm", "confidence_score"):
                if column not in system_df.columns:
                    system_df[column] = None
            if "chains_ptm" not in chain_df.columns:
                chain_df["chains_ptm"] = None

        if "affinity_metrics" in unavailable_metric_groups:
            for column in ("affinity_pred_value", "affinity_probability_binary"):
                if column not in chain_df.columns:
                    chain_df[column] = None

        if "affinity_metrics_ext" in unavailable_metric_groups:
            for column in ("pIC50", "IC50_M", "pIC50_kcal_per_mol"):
                if column not in chain_df.columns:
                    chain_df[column] = None

        return system_df, chain_df

    def _apply_bias_metrics(self, system_df, chain_df, runtime: RunnerRuntime):
        workflow = BiasAssessmentWorkflow(
            wrk_dir=self.wrk_dir,
            system_path=self.system_path,
            sys_obj=self.sys,
            protein_training_data_path=self.protein_training_data_path,
            ligand_training_data_path=self.ligand_training_data_path,
            bias_release_cutoff=self.bias_release_cutoff,
            bias_ligand_similarity_threshold=self.bias_ligand_similarity_threshold,
            bias_chains=self.bias_chains,
            build_bias_training_data=self.build_bias_training_data,
            bias_training_components_cif=self.bias_training_components_cif,
            logger=self.logger,
            timings=self.timings,
            strict_training_data=False,
        )
        return workflow.apply(
            system_df=system_df,
            chain_df=chain_df,
            boltz_cache_path=runtime.cache_path or "~/.boltz",
        )

    @staticmethod
    def _resolve_diffusion_samples(system_df, runtime: RunnerRuntime) -> int:
        if not system_df.empty and "diffusion_sample" in system_df.columns:
            try:
                return int(system_df["diffusion_sample"].max()) + 1
            except Exception:
                pass
        try:
            return int(runtime.diffusion_samples or 1)
        except Exception:
            return 1
