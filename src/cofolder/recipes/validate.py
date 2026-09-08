from __future__ import annotations

import copy
import logging
import os
from contextlib import contextmanager
from dataclasses import replace
from pathlib import Path
from uuid import uuid4

from cofolder.modules.analytics.reproduction import scaffold_reproduction_metrics
from cofolder.modules.analytics.structure import Structure
from cofolder.modules.contracts import (
    PUBLIC_SCHEMA_VERSION,
    ArtifactReference,
    EvidenceSource,
    FailureStage,
    OutputIdentity,
    PublicManifest,
    PublicOutputBundle,
    PublicSerializationError,
    WorkflowExecutionError,
    WorkflowKind,
    bundle_from_frames,
    failure_from_exception,
    write_public_bundle,
)
from cofolder.modules.input import load_yaml_document, system
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
from cofolder.modules.utils import gather, helpers, write
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


class Validate:
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
        self.run_id = str(uuid4())
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
        self.output_identity = OutputIdentity(
            workflow=WorkflowKind.VALIDATE,
            run_id=self.run_id,
            system_id=self.system_path.stem,
            runner_id=self.runner_name,
        )

        self._system = None
        self.base_system: system.System | None = None

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
        self._failure_stage = FailureStage.INPUT_VALIDATION

    def _log_timing_summary(self) -> None:
        self.timings.log_summary(logger=self.logger)

    @contextmanager
    def _debug_timer(self, label: str):
        previous = self._failure_stage
        self._failure_stage = self._stage_for_label(label)
        try:
            with self.timings.measure(label, logger=self.logger):
                yield
        except Exception:
            raise
        else:
            self._failure_stage = previous

    @staticmethod
    def _stage_for_label(label: str) -> FailureStage:
        if label.startswith("runner.prepare"):
            return FailureStage.PREPARATION
        if label.startswith("runner.bundle"):
            return FailureStage.OUTPUT_VALIDATION
        if label.startswith(
            ("structures.gather", "results.merge", "results.add_chain")
        ):
            return FailureStage.GATHER
        if label.startswith(("scores.", "structure.")):
            return FailureStage.ANALYTICS
        if label.startswith("results.write"):
            return FailureStage.SERIALIZATION
        return FailureStage.INPUT_VALIDATION

    def run(self):
        try:
            return self._run_impl()
        except WorkflowExecutionError:
            raise
        except PublicSerializationError:
            raise
        except Exception as exc:
            failure = failure_from_exception(
                exc,
                identity=self.output_identity,
                stage=self._failure_stage,
                error_code=getattr(
                    exc,
                    "error_code",
                    f"validate_{self._failure_stage.value}_failed",
                ),
            )
            output_dir = self.wrk_dir / "results"
            write_public_bundle(
                PublicOutputBundle(
                    manifest=PublicManifest(
                        schema_version=PUBLIC_SCHEMA_VERSION,
                        identity=self.output_identity,
                        status="failed",
                    ),
                    records=(failure,),
                ),
                output_dir,
            )
            raise WorkflowExecutionError(
                str(exc), failures=(failure,), output_dir=output_dir
            ) from exc
        finally:
            self._log_timing_summary()

    def _run_impl(self):
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
            if (
                self.conformers
                and self.conformers
                not in self.runner.ligand_preparation_capabilities.conformer_modes
            ):
                from cofolder.modules.input import LigandValidationError

                raise LigandValidationError(
                    f"Runner '{self.runner_name}' does not support "
                    f"{self.conformers} ligand conformers.",
                    source_path=self.system_path,
                )
            with self._debug_timer("system_yaml.load"):
                document = load_yaml_document(self.system_path)
                self._system = document.value
                if not isinstance(self._system, dict):
                    from cofolder.modules.input import SystemInputValidationError

                    raise SystemInputValidationError(
                        "System YAML root must be a mapping.",
                        source_path=self.system_path,
                    )
                self.base_system = system.System(system=self._system)
                resolve_declared_msa_paths(
                    self.base_system,
                    base_dir=self.system_path.parent,
                )
            resolve_msa_reuse_settings = getattr(
                self.runner, "msa_reuse_settings", None
            )
            msa_reuse_settings = (
                resolve_msa_reuse_settings(runner_options)
                if callable(resolve_msa_reuse_settings)
                else {}
            )

            assert self.base_system is not None
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
                validated_system = self.runner.validate_system(
                    self.sys,
                    runner_options,
                    check_atom_names=False,
                    source_path=self.system_path,
                )
            with self._debug_timer("runner.prepare.availability"):
                self.runner.ensure_available()
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
            chain_identities = validated_system.chain_identities
            with self._debug_timer("runner.system.validate.execution_input"):
                self.runner.validate_system(
                    self.sys,
                    runner_options,
                    check_atom_names=True,
                    source_path=self.system_path,
                )
            for warning in preparation.warnings:
                self.logger.warning("%s", warning)

            yaml_path = os.path.join(self.raw_dir, str(self.system_path.name))
            with self._debug_timer("system_yaml.write"):
                write.write_yaml(self.sys, path=yaml_path)
            self.logger.info("Updated system saved to YAML: %s", yaml_path)

            runner_results = []
            workflow_failures = []
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
                    identity=replace(self.output_identity, repeat_id=i),
                    chain_identities=chain_identities,
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
                result = None
                try:
                    result = self.runner.run(request)
                except Exception as exc:
                    workflow_failures.append(
                        failure_from_exception(
                            exc,
                            identity=replace(self.output_identity, repeat_id=i),
                            stage=FailureStage.BACKEND_EXECUTION,
                            error_code="runner_backend_execution_failed",
                            details={"seed": seed},
                        )
                    )
                    self.logger.exception("Runner repeat %d failed: %s", i, exc)
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
                if result is not None:
                    runner_results.append((i, result))

            selected_runner_metric_groups = (
                self.scoring_functions & RUNNER_METRIC_GROUPS
            )
            requested_runner_metric_groups = self._required_runner_metric_groups()
            optional_runner_metric_groups = (
                selected_runner_metric_groups - requested_runner_metric_groups
            )
            unavailable_metric_groups = set(unsupported_metric_groups)
            valid_runner_results = []
            valid_repeat_ids = []
            with self._debug_timer("runner.bundle_validation"):
                for repeat, result in runner_results:
                    try:
                        bundle = validate_runner_bundle(
                            result,
                            requested_metric_groups=requested_runner_metric_groups,
                        )
                    except Exception as exc:
                        workflow_failures.append(
                            failure_from_exception(
                                exc,
                                identity=replace(
                                    self.output_identity, repeat_id=repeat
                                ),
                                stage=FailureStage.OUTPUT_VALIDATION,
                                error_code="runner_output_validation_failed",
                            )
                        )
                        self.logger.exception(
                            "Runner output validation failed for repeat %d: %s",
                            repeat,
                            exc,
                        )
                        continue
                    valid_runner_results.append(result)
                    valid_repeat_ids.append(repeat)
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

            if not valid_runner_results:
                output_dir = self.wrk_dir / "results"
                write_public_bundle(
                    PublicOutputBundle(
                        manifest=PublicManifest(
                            schema_version=PUBLIC_SCHEMA_VERSION,
                            identity=self.output_identity,
                            status="failed",
                        ),
                        records=tuple(workflow_failures),
                    ),
                    output_dir,
                )
                raise WorkflowExecutionError(
                    "No runner repeat produced a usable normalized result.",
                    failures=tuple(workflow_failures),
                    output_dir=output_dir,
                )

            with self._debug_timer("structures.gather"):
                gather_kwargs = {
                    "base_dir": self.wrk_dir,
                    "system_name": self.system_path.stem,
                    "repeats": self.repeats,
                    "logger": self.logger,
                }
                if valid_repeat_ids != list(range(1, self.repeats + 1)):
                    gather_kwargs["repeat_ids"] = valid_repeat_ids
                gather.gather_structures(**gather_kwargs)

            with self._debug_timer("results.merge_runner_outputs"):
                merge_kwargs = {
                    "raw_dir": self.raw_dir,
                    "repeats": self.repeats,
                    "logger": self.logger,
                }
                if valid_repeat_ids != list(range(1, self.repeats + 1)):
                    merge_kwargs["repeat_ids"] = valid_repeat_ids
                system_df, chain_df, _manifests = gather.merge_runner_results(
                    **merge_kwargs
                )

            if not chain_df.empty:
                with self._debug_timer("results.add_chain_info"):
                    chain_df = gather.add_chain_info(chain_df, chain_identities)

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
                try:
                    with self._debug_timer("structure.init"):
                        structure = Structure(
                            wrk_dir=self.wrk_dir,
                            chain_df=chain_df,
                            cif_folder=Path(self.wrk_dir / "results" / "structures"),
                        )
                except Exception as exc:
                    workflow_failures.append(
                        self._analytics_failure(exc, "structure_initialization_failed")
                    )
                    self.logger.exception("Structure analytics initialization failed")
                    structure = None

                if structure is not None and "ifp_distance" in self.scoring_functions:
                    try:
                        with self._debug_timer("scores.ifp_distance"):
                            chain_df = structure.add_ifp_distance()
                    except Exception as exc:
                        workflow_failures.append(
                            self._analytics_failure(exc, "ifp_distance_failed")
                        )
                        self.logger.exception("Distance IFP analytics failed")

                if structure is not None and "ifp_prolif" in self.scoring_functions:
                    try:
                        with self._debug_timer("scores.ifp_prolif"):
                            chain_df = structure.add_ifp_prolif()
                    except Exception as exc:
                        workflow_failures.append(
                            self._analytics_failure(exc, "ifp_prolif_failed")
                        )
                        self.logger.exception("ProLIF analytics failed")

                if (
                    structure is not None
                    and {
                        "sasa",
                        "sasa_normalized",
                    }
                    & self.scoring_functions
                ):
                    try:
                        with self._debug_timer("scores.sasa"):
                            chain_df = structure.add_sasa(
                                absolute="sasa" in self.scoring_functions,
                                normalized="sasa_normalized" in self.scoring_functions,
                            )
                    except Exception as exc:
                        workflow_failures.append(
                            self._analytics_failure(exc, "sasa_failed")
                        )
                        self.logger.exception("SASA analytics failed")

            try:
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
            except Exception as exc:
                workflow_failures.append(
                    self._analytics_failure(exc, "reproduction_metrics_failed")
                )
                self.logger.exception("Reproduction analytics failed")

            # Story 2.4 retires recipe-side runtime reconstruction from manifest/runtime_context
            # fallbacks. Shared workflow behavior now depends only on typed runner metadata.
            runtime = merge_runner_runtime(
                preparation_runtime,
                *(result.runtime for result in valid_runner_results),
            )

            if self.assess_bias:
                try:
                    system_df, chain_df = self._apply_bias_metrics(
                        system_df=system_df,
                        chain_df=chain_df,
                        runtime=runtime,
                    )
                except Exception as exc:
                    workflow_failures.append(
                        self._analytics_failure(exc, "bias_metrics_failed")
                    )
                    self.logger.exception("Bias analytics failed")

            diffusion_samples = self._resolve_diffusion_samples(system_df, runtime)
            robustness_df = None
            if self.assess_robustness and (self.repeats > 1 or diffusion_samples > 1):
                try:
                    with self._debug_timer("scores.robustness_metrics"):
                        robustness_df = gather.gather_robustness_results(
                            system_df=system_df,
                            chain_df=chain_df,
                            wrk_dir=self.wrk_dir,
                            reference_path=self.reference_path,
                        )
                except Exception as exc:
                    workflow_failures.append(
                        self._analytics_failure(
                            exc,
                            "robustness_metrics_failed",
                            stage=FailureStage.AGGREGATION,
                        )
                    )
                    self.logger.exception("Robustness aggregation failed")

            with self._debug_timer("results.write.public_contract"):
                self._write_public_results(
                    system_df=system_df,
                    chain_df=chain_df,
                    unavailable_metric_groups=unavailable_metric_groups,
                    failures=workflow_failures,
                    robustness_df=robustness_df,
                )

    def _analytics_failure(
        self,
        exc: BaseException,
        error_code: str,
        *,
        stage: FailureStage = FailureStage.ANALYTICS,
    ):
        return failure_from_exception(
            exc,
            identity=self.output_identity,
            stage=stage,
            error_code=error_code,
        )

    def _write_public_results(
        self,
        *,
        system_df,
        chain_df,
        unavailable_metric_groups: set[str],
        failures=(),
        robustness_df=None,
    ) -> None:
        evidence: list[EvidenceSource] = []
        if self.reference_path is not None:
            evidence.append(
                EvidenceSource(
                    kind="reference_structure",
                    identifier=self.reference_path.name,
                    path=str(self.reference_path),
                )
            )
        if self.pocket_coverage_reference:
            evidence.append(
                EvidenceSource(
                    kind="custom_pocket",
                    identifier="configured_custom_pocket",
                )
            )
        requested = set(self.scoring_functions)
        if self.reproduction_metrics:
            requested.add("reproduction_metrics")
        if self.assess_bias:
            requested.add("bias_metrics")
        artifact_refs: list[ArtifactReference] = []
        for label, relative_path in (
            ("predicted_structures", "structures"),
            ("aligned_structures", "structures_aligned"),
            ("robustness_matrices", "matrices"),
            ("interaction_fingerprints", "ifp"),
            ("bias_supporting_outputs", "bias_train"),
        ):
            if (self.wrk_dir / "results" / relative_path).exists():
                artifact_refs.append(
                    ArtifactReference(
                        label=label,
                        relative_path=relative_path,
                        kind="directory",
                    )
                )
        bundle = bundle_from_frames(
            system_df,
            chain_df,
            identity=self.output_identity,
            evidence=evidence,
            requested_metrics=requested,
            unavailable_groups=unavailable_metric_groups,
            failures=failures,
            artifacts=artifact_refs,
            robustness_df=robustness_df,
        )
        write_public_bundle(bundle, self.wrk_dir / "results")

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
