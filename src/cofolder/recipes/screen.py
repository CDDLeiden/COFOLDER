"""Screening workflow over canonical CSV, SDF, and MOL library members."""

from __future__ import annotations

import logging
import math
import re
from dataclasses import replace
from pathlib import Path
from typing import Any
from uuid import uuid4

import pandas as pd

from cofolder.modules.analytics.ifp_filtering import (
    IFPFilterMode,
    ReferenceIFPFilterPolicy,
)
from cofolder.modules.analytics.reference_ifp import (
    IFPSimilarityMetric,
    IFPTaxonomy,
    InteractionFingerprint,
    InteractionKey,
    LigandSelector,
)
from cofolder.modules.analytics.reproduction import (
    _parse_custom_pocket_reference,
)
from cofolder.modules.contracts import (
    METRIC_CATALOG,
    PUBLIC_SCHEMA_VERSION,
    ArtifactReference,
    EvidenceRegime,
    ExecutionRecord,
    ExecutionStatus,
    FailureStage,
    MetricRecord,
    OutputIdentity,
    PublicManifest,
    PublicOutputBundle,
    PublicSerializationError,
    RecordKind,
    ScreenOutputNormalizationError,
    StructuredExecutionError,
    SuccessRecord,
    WorkflowExecutionError,
    WorkflowFailureRecord,
    WorkflowKind,
    make_envelope,
    metric_records_from_frames,
    read_public_records,
    write_public_bundle,
)
from cofolder.modules.input import (
    SystemInputValidationError,
    WorkflowInputRequirements,
    load_yaml_document,
    system,
)
from cofolder.modules.input.compound_library import (
    CompoundLibraryFormat,
    CompoundMember,
    CompoundMemberFailure,
    DuplicateIdPolicy,
    load_compound_library,
)
from cofolder.modules.input.ligand import (
    replace_ligand_smiles,
    resolve_ligand_target,
)
from cofolder.modules.input.system import iter_system_chains
from cofolder.modules.runners import (
    PlannedExecution,
    RunnerExecutionPlan,
    get_runner,
)
from cofolder.modules.runners.msa import resolve_declared_msa_paths
from cofolder.modules.utils import write
from cofolder.recipes._metrics import read_metric_frames
from cofolder.recipes._completion import report_completion
from cofolder.recipes._diagnostics import workflow_failure
from cofolder.recipes._execution import build_execution_plan
from cofolder.recipes._results import standard_evidence, write_failure_bundle
from cofolder.recipes._screen_postprocess import (
    ScreenPostprocessConfig,
    ScreenPostprocessor,
)
from cofolder.recipes.validate import Validate

logger = logging.getLogger(__name__)


class Screen:
    """Run Validate for each canonical library member and consolidate analyses.

    ``cluster_ifps`` applies deterministic average-linkage clustering to compatible
    binary distance IFPs after prediction. ``ifp_cluster_similarity_threshold`` is
    the inclusive Jaccard-similarity cut and defaults to ``0.5``. Both consolidated
    Returned rows always contain cluster ID/status columns; an enabled run additionally
    publishes ``results/ifp_cluster_summary.csv``.
    """

    def __init__(
        self,
        wrk_dir: str,
        system_path: str,
        options_path: str,
        runner: str = "boltz2",
        library: str | None = None,
        library_format: str | None = None,
        ligand_chain: str | None = None,
        smiles_column: str | None = None,
        col_id: str | None = None,
        id_property: str = "_Name",
        duplicate_id_policy: str = "reject",
        merge_data: str | None = None,
        repeats: int = 1,
        seed: int | None = None,
        scoring_functions: list[str] | None = None,
        assess_robustness: bool = True,
        assess_bias: bool = False,
        protein_training_data_path: str | None = None,
        ligand_training_data_path: str | None = None,
        bias_training_data_protein_path: str | None = None,
        bias_training_data_ligand_path: str | None = None,
        bias_query_cache_path: str | None = None,
        custom_bias_reference_path: str | None = None,
        custom_protein_reference_path: str | None = None,
        custom_ligand_reference_path: str | None = None,
        bias_release_cutoff: str = "2023-06-01",
        bias_protein_similarity_threshold: float = 0.25,
        bias_ligand_similarity_threshold: float = 0.35,
        bias_chains: list[str] | None = None,
        build_bias_training_data: bool = False,
        bias_training_components_cif: str | None = None,
        conformers: str | None = None,
        sdf_file: str | None = None,
        reference_path: str | None = None,
        pocket_coverage_reference: str | None = None,
        reproduction_metrics: list[str] | None = None,
        ifp_filter_threshold: float | None = None,
        ifp_filter_source: str = "auto",
        ifp_taxonomy: str = "distance",
        ifp_similarity_metric: str | None = None,
        ifp_filter_policy: str = "similarity",
        ifp_required_interactions: list[str] | None = None,
        ifp_reference_ligand: str | None = None,
        ifp_reference_receptor_chains: list[str] | None = None,
        cluster_ifps: bool = False,
        ifp_cluster_similarity_threshold: float = 0.5,
    ):
        self.wrk_dir = Path(wrk_dir)
        self.system_path = Path(system_path)
        self.options_path = Path(options_path)
        self.runner = str(runner)
        self.run_id = str(uuid4())

        self.ligand_chain = str(ligand_chain).strip() if ligand_chain else ""
        self._smiles_column_configured = smiles_column is not None
        self._col_id_configured = col_id is not None
        self.smiles_column = str(smiles_column).strip() if smiles_column else "smiles"
        self.col_id = str(col_id).strip() if col_id else "execution_id"
        self.library = Path(library) if library else None
        self.library_format = (
            CompoundLibraryFormat(library_format) if library_format else None
        )
        self.id_property = str(id_property)
        self.duplicate_id_policy = DuplicateIdPolicy(duplicate_id_policy)
        self.merge_data = self._parse_list(merge_data)
        self.ifp_filter_threshold = ifp_filter_threshold
        self.ifp_filter_source = str(ifp_filter_source)
        self.ifp_taxonomy = IFPTaxonomy(ifp_taxonomy)
        self.ifp_similarity_metric = (
            IFPSimilarityMetric(ifp_similarity_metric)
            if ifp_similarity_metric is not None else None
        )
        self.ifp_filter_mode = IFPFilterMode(ifp_filter_policy)
        self.ifp_required_interaction_values = tuple(ifp_required_interactions or ())
        self.ifp_reference_ligand_value = ifp_reference_ligand
        self.ifp_reference_receptor_chains = tuple(ifp_reference_receptor_chains or ())
        self.pocket_coverage_reference = pocket_coverage_reference
        self.cluster_ifps = bool(cluster_ifps)
        self.ifp_cluster_similarity_threshold = ifp_cluster_similarity_threshold
        self._ifp_filter_reference_spec: dict[str, Any] | None = None
        self._reference_ifp: InteractionFingerprint | None = None
        self._reference_ifp_failure: str | None = None
        self._resolved_ifp_filter_source: str | None = None
        self._reference_ligand_selector: LigandSelector | None = None
        self._ifp_filter_policy_config: ReferenceIFPFilterPolicy | None = None
        self._prediction_ifp_cache: dict[
            tuple[Path, int | None, int | None],
            tuple[InteractionFingerprint | None, str],
        ] = {}
        self._postprocessor: ScreenPostprocessor | None = None

        self.validate_kwargs: dict[str, Any] = {
            "repeats": repeats,
            "seed": seed,
            "scoring_functions": scoring_functions,
            "assess_robustness": assess_robustness,
            "assess_bias": assess_bias,
            "protein_training_data_path": protein_training_data_path,
            "ligand_training_data_path": ligand_training_data_path,
            "bias_training_data_protein_path": bias_training_data_protein_path,
            "bias_training_data_ligand_path": bias_training_data_ligand_path,
            "bias_query_cache_path": bias_query_cache_path,
            "custom_bias_reference_path": custom_bias_reference_path,
            "custom_protein_reference_path": custom_protein_reference_path,
            "custom_ligand_reference_path": custom_ligand_reference_path,
            "bias_release_cutoff": bias_release_cutoff,
            "bias_protein_similarity_threshold": bias_protein_similarity_threshold,
            "bias_ligand_similarity_threshold": bias_ligand_similarity_threshold,
            "bias_chains": bias_chains,
            "build_bias_training_data": build_bias_training_data,
            "bias_training_components_cif": bias_training_components_cif,
            "conformers": conformers,
            "sdf_file": sdf_file,
            "reference_path": reference_path,
            "pocket_coverage_reference": pocket_coverage_reference,
            "reproduction_metrics": reproduction_metrics,
        }

        self.base_system: dict[str, Any] = {}
        self.base_system_obj: system.System | None = None
        self.ligand_target = None
        self.runner_impl = get_runner(self.runner)
        self.execution_plan: RunnerExecutionPlan | None = None
        self._failure_stage = FailureStage.INPUT_VALIDATION
        self.reusable_msa_dir = (
            self.wrk_dir / "shared" / "msa" / self.runner
            if getattr(self.runner_impl, "supports_msa_reuse", False)
            else None
        )
        self.msa_reuse_settings = {}
        self.logger = logging.getLogger("cofolder.screen")

    def preflight(self):
        """Validate and describe the screen without searches or inference."""
        from cofolder.recipes.preflight import prediction_preflight

        try:
            self._validate_config()
        except Exception as exc:
            from cofolder.recipes.preflight import PreflightReport

            return PreflightReport(
                workflow="screen",
                ready=False,
                runner=self.runner,
                output_dir=str(self.wrk_dir / "results"),
                messages=(str(exc),),
            )
        report, system_obj, _options = prediction_preflight(
            workflow="screen",
            system_path=self.system_path,
            options_path=self.options_path,
            runner_name=self.runner,
            repeats=int(self.validate_kwargs["repeats"]),
            output_dir=self.wrk_dir,
            require_ligand=True,
            assess_bias=bool(self.validate_kwargs["assess_bias"]),
            use_bias_databases=bool(
                self.validate_kwargs["bias_training_data_protein_path"]
                or self.validate_kwargs["bias_training_data_ligand_path"]
                or (
                    self.validate_kwargs["protein_training_data_path"] is None
                    and self.validate_kwargs["ligand_training_data_path"] is None
                    and self.validate_kwargs["custom_bias_reference_path"] is None
                    and self.validate_kwargs["custom_protein_reference_path"] is None
                    and self.validate_kwargs["custom_ligand_reference_path"] is None
                    and not self.validate_kwargs["build_bias_training_data"]
                )
            ),
            protein_database_path=self.validate_kwargs["bias_training_data_protein_path"],
            ligand_database_path=self.validate_kwargs["bias_training_data_ligand_path"],
            release_cutoff=self.validate_kwargs["bias_release_cutoff"],
            bias_chains=self.validate_kwargs["bias_chains"],
            bias_query_cache_path=self.validate_kwargs["bias_query_cache_path"],
            custom_bias_reference_path=self.validate_kwargs["custom_bias_reference_path"],
            protein_similarity_threshold=self.validate_kwargs[
                "bias_protein_similarity_threshold"
            ],
            ligand_similarity_threshold=self.validate_kwargs[
                "bias_ligand_similarity_threshold"
            ],
        )
        if not report.ready:
            return report
        try:
            target = resolve_ligand_target(system_obj, self.ligand_chain or None)
        except Exception as exc:
            return replace(report, ready=False, messages=report.messages + (str(exc),))
        self.ligand_chain = target.chain_ids[0]
        return replace(
            report,
            messages=report.messages + (f"selected_ligand_chain={self.ligand_chain}",),
        )

    def _build_execution_plan(self, options_obj: Any) -> RunnerExecutionPlan:
        """Resolve immutable runner, seed, repeat, model, and sample axes once."""
        return build_execution_plan(
            self.runner_impl,
            self.runner,
            options_obj,
            repeats=int(self.validate_kwargs["repeats"]),
            seed=self.validate_kwargs.get("seed"),
            logger=self.logger,
        )


    def _validate_config(self) -> None:
        if self.library is None:
            raise ValueError("--library is required.")
        if not self.library.exists():
            raise ValueError(f"--library does not exist: {self.library}")
        if not self.library.is_file():
            raise ValueError(f"--library is not a file: {self.library}")
        inferred = self.library_format
        if inferred is None:
            inferred = {
                ".csv": CompoundLibraryFormat.CSV,
                ".sdf": CompoundLibraryFormat.SDF,
                ".sd": CompoundLibraryFormat.SDF,
                ".mol": CompoundLibraryFormat.MOL,
            }.get(self.library.suffix.lower())
        if inferred is CompoundLibraryFormat.CSV:
            if not self._col_id_configured or not self.col_id:
                raise ValueError("--col_id is required for CSV libraries.")
            if not self._smiles_column_configured or not self.smiles_column:
                raise ValueError("--smiles_column is required for CSV libraries.")
        elif inferred in {CompoundLibraryFormat.SDF, CompoundLibraryFormat.MOL}:
            if self._col_id_configured or self._smiles_column_configured:
                raise ValueError(
                    "--col_id and --smiles_column are only valid for CSV libraries."
                )
            if self.validate_kwargs.get("sdf_file") is not None:
                raise ValueError(
                    "--sdf_file cannot be combined with an SDF/MOL screening library."
                )

        try:
            self.ifp_cluster_similarity_threshold = float(
                self.ifp_cluster_similarity_threshold
            )
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "--ifp_cluster_similarity_threshold must be a number in [0, 1]."
            ) from exc
        if (
            not math.isfinite(self.ifp_cluster_similarity_threshold)
            or not 0 <= self.ifp_cluster_similarity_threshold <= 1
        ):
            raise ValueError("--ifp_cluster_similarity_threshold must be in [0, 1].")

        filter_enabled = (
            self.ifp_filter_threshold is not None
            or self.ifp_filter_mode is IFPFilterMode.REQUIRED
            or bool(self.ifp_required_interaction_values)
        )
        if not filter_enabled:
            return
        if self.ifp_filter_mode is IFPFilterMode.SIMILARITY:
            try:
                self.ifp_filter_threshold = float(self.ifp_filter_threshold)
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    "Similarity IFP filtering requires --ifp_filter_threshold in [0, 1]."
                ) from exc
            if not math.isfinite(self.ifp_filter_threshold) or not 0 <= self.ifp_filter_threshold <= 1:
                raise ValueError("--ifp_filter_threshold must be in [0, 1].")
        elif self.ifp_filter_threshold is not None:
            raise ValueError("Required-interaction filtering does not accept --ifp_filter_threshold.")
        self._validate_filter_reference()

    def _validate_filter_reference(self) -> None:
        """Parse the filter reference before any prediction work starts."""
        if self.ifp_filter_source not in {"auto", "reference_complex", "custom_pocket"}:
            raise ValueError("--ifp_filter_source must be auto, reference_complex, or custom_pocket.")
        reference_path = self.validate_kwargs.get("reference_path")
        self._resolved_ifp_filter_source = (
            "custom_pocket" if self.ifp_filter_source == "auto" and self.pocket_coverage_reference
            else "reference_complex" if self.ifp_filter_source == "auto" and reference_path
            else self.ifp_filter_source
        )
        if self._resolved_ifp_filter_source == "auto":
            raise ValueError("IFP filtering requires a reference complex or custom pocket.")

        metric = self.ifp_similarity_metric or (
            IFPSimilarityMetric.REFERENCE_COVERAGE
            if self._resolved_ifp_filter_source == "custom_pocket"
            else IFPSimilarityMetric.JACCARD
        )
        required = frozenset(
            InteractionKey.parse(item) for item in self.ifp_required_interaction_values
        ) or None
        self._ifp_filter_policy_config = ReferenceIFPFilterPolicy(
            mode=self.ifp_filter_mode,
            similarity_metric=metric,
            threshold=self.ifp_filter_threshold,
            required_interactions=required,
        )

        if self._resolved_ifp_filter_source == "reference_complex":
            if reference_path is None:
                raise ValueError("Reference-complex IFP filtering requires --reference_path.")
            self._reference_ligand_selector = self._parse_ligand_selector(
                self.ifp_reference_ligand_value
            )
            return
        if self.ifp_taxonomy is IFPTaxonomy.PROLIF:
            raise ValueError("ProLIF filtering requires --ifp_filter_source reference_complex.")
        if self.ifp_filter_mode is IFPFilterMode.REQUIRED:
            raise ValueError("Required-interaction filtering requires a reference complex.")
        value = self.pocket_coverage_reference
        if value is None or not str(value).strip():
            raise ValueError(
                "--ifp_filter_threshold requires --pocket_coverage_reference."
            )

        raw = str(value).strip()
        path = Path(raw).expanduser()
        looks_like_path = bool(path.suffix) or "/" in raw or "\\" in raw
        if looks_like_path and not path.exists():
            raise ValueError(f"--pocket_coverage_reference file does not exist: {path}")
        if path.exists() and not path.is_file():
            raise ValueError(f"--pocket_coverage_reference is not a file: {path}")

        try:
            reference_value = str(path) if path.exists() else raw
            self._ifp_filter_reference_spec = _parse_custom_pocket_reference(
                reference_value
            )
        except (OSError, ValueError) as exc:
            raise ValueError(f"Invalid --pocket_coverage_reference: {exc}") from exc
        if self._ifp_filter_reference_spec is None:
            raise ValueError("--pocket_coverage_reference must not be empty.")

    @staticmethod
    def _parse_ligand_selector(value: str | None) -> LigandSelector | None:
        if value is None or not str(value).strip():
            return None
        tokens = str(value).strip().split(":")
        if len(tokens) == 1:
            return LigandSelector(chain_id=tokens[0])
        if len(tokens) != 2 or not tokens[0] or not tokens[1]:
            raise ValueError("--ifp_reference_ligand expects CHAIN or CHAIN:RESNUM[ICODE].")
        residue = tokens[1]
        index = 1 if residue.startswith("-") else 0
        while index < len(residue) and residue[index].isdigit():
            index += 1
        if not residue[:index] or len(residue[index:]) > 1:
            raise ValueError("--ifp_reference_ligand expects CHAIN or CHAIN:RESNUM[ICODE].")
        return LigandSelector(tokens[0], int(residue[:index]), residue[index:])

    @report_completion(WorkflowKind.SCREEN)
    def run(self) -> pd.DataFrame:
        try:
            return self._run_impl()
        except (
            WorkflowExecutionError,
            PublicSerializationError,
            ScreenOutputNormalizationError,
        ):
            raise
        except Exception as exc:
            identity = OutputIdentity(
                workflow=WorkflowKind.SCREEN,
                run_id=self.run_id,
                system_id=self.system_path.stem,
                runner_id=self.runner,
            )
            output_dir = self.wrk_dir / "results"
            failure = write_failure_bundle(
                output_dir=output_dir,
                exc=exc,
                identity=identity,
                stage=self._failure_stage,
                error_code=getattr(exc, "error_code", "screen_input_validation_failed"),
            )
            raise WorkflowExecutionError(
                str(exc), failures=(failure,), output_dir=output_dir
            ) from exc

    def _screen_postprocessor(self) -> ScreenPostprocessor:
        if self._postprocessor is None:
            self._postprocessor = ScreenPostprocessor(
                ScreenPostprocessConfig.from_screen(self),
                logger=self.logger,
            )
        return self._postprocessor

    def _ifp_filter_enabled(self) -> bool:
        return self._ifp_filter_policy_config is not None

    def _prepare_reference_ifp(self) -> None:
        self._screen_postprocessor()._prepare_reference_ifp()

    def _run_impl(self) -> pd.DataFrame:
        """Run the screen, publish contract records, and return the merged results."""

        self.wrk_dir.mkdir(parents=True, exist_ok=True)
        self._validate_config()
        document = load_yaml_document(self.system_path)
        if not isinstance(document.value, dict):
            raise SystemInputValidationError(
                "System YAML root must be a mapping.", source_path=self.system_path
            )
        self.base_system_obj = system.System(system=document.value)
        resolve_declared_msa_paths(
            self.base_system_obj, base_dir=self.system_path.parent
        )
        options_obj = self.runner_impl.load_options(self.options_path)
        self.execution_plan = self._build_execution_plan(options_obj)
        validated = self.runner_impl.validate_system(
            self.base_system_obj,
            options_obj,
            check_atom_names=False,
            source_path=self.system_path,
            requirements=WorkflowInputRequirements(
                require_protein=True, require_ligand=True
            ),
        )
        self.base_system_obj = validated.system
        self.base_system = self.base_system_obj.system
        self.ligand_target = resolve_ligand_target(
            self.base_system_obj, self.ligand_chain or None
        )
        self.ligand_chain = self.ligand_target.chain_ids[0]
        self._postprocessor = None
        self._prepare_reference_ifp()
        if self.reusable_msa_dir is not None:
            self.msa_reuse_settings = self.runner_impl.msa_reuse_settings(options_obj)
        assert self.ligand_target is not None
        assert self.library is not None
        compound_library = load_compound_library(
            self.library,
            target=self.ligand_target,
            source_format=self.library_format,
            smiles_column=self.smiles_column,
            id_column=self.col_id if self.col_id != "execution_id" else None,
            id_property=self.id_property,
            metadata_fields=self.merge_data,
            duplicate_policy=self.duplicate_id_policy,
        )
        if self.validate_kwargs.get("conformers") == "sdf":
            if compound_library.source_format is CompoundLibraryFormat.CSV:
                raise ValueError("--conformers sdf requires an SDF/MOL screening library.")
            capabilities = self.runner_impl.ligand_preparation_capabilities
            if "sdf" not in capabilities.conformer_modes:
                raise ValueError(
                    f"Runner '{self.runner}' does not support SDF ligand conformers."
                )

        records: list[dict[str, Any]] = []
        records_with_scores: list[dict[str, Any]] = []
        public_failures = []
        member_manifest: list[dict[str, Any]] = []
        total = len(compound_library.outcomes)
        for outcome in compound_library.outcomes:
            source = outcome.source
            i = source.source_record_index
            compound_id = outcome.execution_id
            run_dir = self.wrk_dir / outcome.execution_directory
            run_dir.mkdir(parents=True, exist_ok=True)
            row_system_path = run_dir / "screen_system.yaml"

            summary = {
                "index": i,
                self.col_id: compound_id,
                "source_format": source.source_format.value,
                "source_path": str(source.source_path),
                "source_record_index": i,
                "source_record_id": source.source_record_id,
                "original_id": source.original_id,
                "execution_id": compound_id,
                "execution_directory": outcome.execution_directory,
                "status": "success",
                "error_message": "",
                "run_dir": str(run_dir),
            }
            filter_result = self._default_filter_result()
            summary.update(filter_result)
            detailed: dict[str, Any] = dict(source.metadata)
            detailed.update(summary)
            self._ensure_screen_metric_schema(detailed)
            detailed.update(self._default_cluster_result())
            summary.update(self._default_cluster_result())
            summary[self.smiles_column] = (
                outcome.ligand.source_smiles
                if isinstance(outcome, CompoundMember)
                else source.metadata.get(self.smiles_column)
            )
            for col in self.merge_data:
                summary[col] = source.metadata.get(col)

            self.logger.info("(%d/%d) screening %s", i, total, compound_id)

            sys_obj: system.System | None = None
            row_stage = FailureStage.INPUT_VALIDATION
            coordinate_mode = (
                outcome.coordinate_mode
                if isinstance(outcome, CompoundMember)
                else None
            )
            row_exception: BaseException | None = None
            try:
                assert self.ligand_target is not None
                if isinstance(outcome, CompoundMemberFailure):
                    raise outcome.exception
                normalized_ligand = outcome.ligand
                sys_obj = replace_ligand_smiles(
                    self.base_system_obj,
                    target=self.ligand_target,
                    ligand=normalized_ligand,
                )

                if self.reusable_msa_dir is not None:
                    injected = self.runner_impl.inject_reusable_msas(
                        sys_obj,
                        self.reusable_msa_dir,
                        settings=self.msa_reuse_settings,
                    )
                    if injected:
                        self.logger.info(
                            "Reusing %d shared protein MSA(s) for %s.",
                            injected,
                            compound_id,
                        )

                write.write_yaml(sys_obj, path=row_system_path)

                row_validate_kwargs = dict(self.validate_kwargs)
                if outcome.conformer_molblock is not None:
                    source_sdf = run_dir / "source_ligand.sdf"
                    source_sdf.write_text(
                        outcome.conformer_molblock.rstrip() + "\n$$$$\n",
                        encoding="utf-8",
                    )
                    requested_mode = row_validate_kwargs.get("conformers")
                    supports_source = (
                        "sdf"
                        in self.runner_impl.ligand_preparation_capabilities.conformer_modes
                    )
                    if requested_mode in {"2D", "3D"}:
                        coordinate_mode = "generated"
                        row_validate_kwargs["sdf_file"] = None
                    elif supports_source:
                        coordinate_mode = "source"
                        row_validate_kwargs["conformers"] = "sdf"
                        row_validate_kwargs["sdf_file"] = str(source_sdf)
                    else:
                        coordinate_mode = "native_smiles"
                        row_validate_kwargs["conformers"] = None
                        row_validate_kwargs["sdf_file"] = None
                if row_validate_kwargs.get("assess_bias") and row_validate_kwargs.get(
                    "build_bias_training_data"
                ):
                    bias_train_dir = run_dir / "results" / "bias_train"
                    bias_train_dir.mkdir(parents=True, exist_ok=True)
                    # Build ligand references per-screened system after row system YAML exists.
                    row_validate_kwargs["ligand_training_data_path"] = str(
                        bias_train_dir / "ligand_training_data.csv"
                    )

                row_stage = FailureStage.PREPARATION
                validator = Validate(
                    wrk_dir=str(run_dir),
                    system_path=str(row_system_path),
                    options_path=str(self.options_path),
                    runner=self.runner,
                    execution_plan=self.execution_plan,
                    reusable_msa_dir=(
                        str(self.reusable_msa_dir)
                        if self.reusable_msa_dir is not None
                        else None
                    ),
                    **row_validate_kwargs,
                )
                row_stage = FailureStage.BACKEND_EXECUTION
                validator.run()
                if self.reusable_msa_dir is not None:
                    self.runner_impl.inject_reusable_msas(
                        sys_obj,
                        self.reusable_msa_dir,
                        settings=self.msa_reuse_settings,
                    )
                    write.write_yaml(sys_obj, path=row_system_path)
                row_stage = FailureStage.ANALYTICS
            except Exception as exc:
                row_exception = exc
                if self.reusable_msa_dir is not None and sys_obj is not None:
                    try:
                        injected = self.runner_impl.inject_reusable_msas(
                            sys_obj,
                            self.reusable_msa_dir,
                            settings=self.msa_reuse_settings,
                        )
                        if injected and row_system_path.exists():
                            write.write_yaml(sys_obj, path=row_system_path)
                    except Exception as cache_exc:
                        self.logger.warning(
                            "Could not update failed row YAML from reusable MSA cache: %s",
                            cache_exc,
                        )
                summary["status"] = "failed"
                summary["error_message"] = str(exc)
                detailed["status"] = "failed"
                detailed["error_message"] = str(exc)
                if self._ifp_filter_enabled():
                    filter_result = self._not_evaluable_filter_result("row_failed")
                    summary.update(filter_result)
                    detailed.update(filter_result)
                self.logger.exception("Screen row failed (%s): %s", compound_id, exc)
                source_exc = exc if isinstance(outcome, CompoundMemberFailure) else exc.__cause__ or exc
                if isinstance(exc, WorkflowExecutionError) and exc.failures:
                    row_stage = exc.failures[0].stage
                failure_details = {
                    "row_index": i,
                    "source_format": source.source_format.value,
                    "source_path": str(source.source_path),
                    "source_record_index": i,
                    "source_record_id": source.source_record_id,
                    "original_id": source.original_id,
                    "execution_id": compound_id,
                    "execution_directory": outcome.execution_directory,
                }
                child_failures = (
                    exc.failures
                    if isinstance(exc, WorkflowExecutionError) and exc.failures
                    else ()
                )
                if child_failures:
                    for child_failure in child_failures:
                        child_identity = child_failure.envelope.identity
                        rebased_identity = replace(
                            child_identity,
                            workflow=WorkflowKind.SCREEN,
                            run_id=self.run_id,
                            system_id=self.system_path.stem,
                            compound_id=compound_id,
                            execution_directory=outcome.execution_directory,
                            runner_id=self.execution_plan.backend.runner_name,
                            runner_version=self.execution_plan.backend.version,
                            backend_name=self.execution_plan.backend.backend_name,
                            backend_version_status=(
                                self.execution_plan.backend.version_status
                            ),
                        )
                        public_failures.append(
                            replace(
                                child_failure,
                                envelope=make_envelope(
                                    RecordKind.FAILURE,
                                    rebased_identity,
                                    child_failure.stage.value,
                                    child_failure.error_code,
                                ),
                                details={**child_failure.details, **failure_details},
                            )
                        )
                else:
                    public_failures.append(
                        workflow_failure(
                            source_exc,
                            identity=OutputIdentity(
                                workflow=WorkflowKind.SCREEN,
                                run_id=self.run_id,
                                system_id=self.system_path.stem,
                                compound_id=compound_id,
                                execution_directory=outcome.execution_directory,
                                runner_id=self.runner,
                                runner_version=self.execution_plan.backend.version,
                                backend_name=self.execution_plan.backend.backend_name,
                                backend_version_status=(
                                    self.execution_plan.backend.version_status
                                ),
                            ),
                            stage=row_stage,
                            error_code=getattr(
                                source_exc, "error_code", "screen_compound_failed"
                            ),
                            details=failure_details,
                        )
                    )

            detailed["coordinate_mode"] = coordinate_mode
            summary["coordinate_mode"] = coordinate_mode
            assert self.execution_plan is not None
            system_frame, chain_frame = read_metric_frames(run_dir)
            has_normalized_output = not system_frame.empty or not chain_frame.empty
            child_bundle_path = run_dir / "results" / "records.jsonl"
            child_bundle_failures = (
                tuple(
                    record
                    for record in read_public_records(child_bundle_path)
                    if isinstance(record, WorkflowFailureRecord)
                )
                if child_bundle_path.is_file()
                else ()
            )
            if row_exception is None and has_normalized_output:
                self._validate_observed_execution_keys(system_frame)
            known_child_failures = (
                row_exception.failures
                if isinstance(row_exception, WorkflowExecutionError)
                else child_bundle_failures
            )
            failure_repeat_ids = {
                failure.envelope.identity.repeat_id
                for failure in known_child_failures
                if failure.envelope.identity.repeat_id is not None
            }
            for slot in self.execution_plan.executions:
                matching_child_failures = [
                    failure
                    for failure in known_child_failures
                    if failure.envelope.identity.repeat_id in {None, slot.repeat_id}
                    and (
                        failure.envelope.identity.model_id in {None, slot.model_id}
                    )
                    and (
                        failure.envelope.identity.sample_id in {None, slot.sample_id}
                    )
                ]
                slot_failure = (
                    matching_child_failures[0]
                    if matching_child_failures
                    else None
                )
                execution_row = dict(detailed)
                execution_row.update(
                    {
                        "repeat": slot.repeat_id,
                        "repeat_id": slot.repeat_id,
                        "model_name": slot.model_id,
                        "model_id": slot.model_id,
                        "diffusion_sample": slot.sample_id,
                        "sample_id": slot.sample_id,
                        "execution_key": (
                            f"{compound_id}|repeat={slot.repeat_id}|"
                            f"model={slot.model_id}|sample={slot.sample_id}"
                        ),
                        "effective_seed": slot.effective_seed,
                        "runner_id": self.execution_plan.backend.runner_name,
                        "backend_name": self.execution_plan.backend.backend_name,
                        "runner_version": self.execution_plan.backend.version,
                        "backend_version_status": (
                            self.execution_plan.backend.version_status.value
                        ),
                    }
                )
                if row_exception is None and slot_failure is None:
                    slot_values, observed = self._collect_execution_score_columns(
                        system_frame,
                        chain_frame,
                        slot,
                    )
                    execution_row.update(slot_values)
                    execution_row.update(
                        self._evaluate_ifp_filter(
                            run_dir,
                            repeat_id=slot.repeat_id,
                            sample_id=slot.sample_id,
                        )
                    )
                    if has_normalized_output and not observed:
                        execution_row["status"] = "unavailable"
                        execution_row["error_stage"] = FailureStage.OUTPUT_VALIDATION.value
                        execution_row["exception_type"] = "MissingExecutionOutput"
                        execution_row["error_code"] = "screen_execution_output_missing"
                        execution_row["error_message"] = (
                            "The runner did not produce this planned model/sample output."
                        )
                        execution_row["error_details"] = {
                            "repeat_id": slot.repeat_id,
                            "model_id": slot.model_id,
                            "sample_id": slot.sample_id,
                        }
                    else:
                        execution_row["status"] = "success"
                        execution_row["error_stage"] = ""
                        execution_row["exception_type"] = ""
                        execution_row["error_code"] = ""
                        execution_row["error_message"] = ""
                        execution_row["error_details"] = {}
                else:
                    attempted = (
                        not isinstance(outcome, CompoundMemberFailure)
                        and (slot_failure.stage if slot_failure else row_stage)
                        in {
                            FailureStage.BACKEND_EXECUTION,
                            FailureStage.OUTPUT_VALIDATION,
                            FailureStage.ANALYTICS,
                            FailureStage.AGGREGATION,
                        }
                        and (
                            not failure_repeat_ids
                            or slot.repeat_id in failure_repeat_ids
                        )
                    )
                    execution_row["status"] = "failed" if attempted else "unavailable"
                    execution_row["error_stage"] = (
                        slot_failure.stage.value if slot_failure else row_stage.value
                    )
                    execution_row["exception_type"] = (
                        slot_failure.exception_type
                        if slot_failure
                        else type(row_exception).__name__
                    )
                    execution_row["error_code"] = (
                        slot_failure.error_code
                        if slot_failure
                        else getattr(
                            row_exception,
                            "error_code",
                            "screen_compound_failed"
                            if attempted
                            else "screen_execution_unavailable",
                        )
                    )
                    execution_row["error_message"] = (
                        slot_failure.message if slot_failure else str(row_exception)
                    )
                    execution_row["error_details"] = (
                        dict(slot_failure.details) if slot_failure else {}
                    )
                    if self._ifp_filter_enabled():
                        execution_row.update(
                            self._not_evaluable_filter_result("row_failed")
                        )
                records.append(dict(execution_row))
                records_with_scores.append(execution_row)
            member_manifest.append(
                {
                    "source_format": source.source_format.value,
                    "source_path": str(source.source_path),
                    "source_record_index": i,
                    "source_record_id": source.source_record_id,
                    "original_id": source.original_id,
                    "execution_id": compound_id,
                    "execution_directory": outcome.execution_directory,
                    "coordinate_mode": coordinate_mode,
                    "status": summary["status"],
                    "error_code": (
                        public_failures[-1].error_code
                        if summary["status"] == "failed"
                        else ""
                    ),
                    "error_message": summary["error_message"],
                }
            )

        summary_df = pd.DataFrame(records)
        results_df = pd.DataFrame(records_with_scores)
        self._apply_ifp_clustering(summary_df, results_df)

        self._set_filter_dtypes(summary_df)
        self._set_filter_dtypes(results_df)
        results_dir = self.wrk_dir / "results"
        results_dir.mkdir(parents=True, exist_ok=True)
        members_path = results_dir / "compound_members.csv"
        temporary_members_path = results_dir / f".{members_path.name}.{uuid4().hex}.tmp"
        try:
            pd.DataFrame(member_manifest).to_csv(temporary_members_path, index=False)
            temporary_members_path.replace(members_path)
        finally:
            temporary_members_path.unlink(missing_ok=True)
        has_success = self._write_public_results(results_df, public_failures)

        if not has_success:
            raise WorkflowExecutionError(
                "No screened compound produced a usable result.",
                failures=tuple(public_failures),
                output_dir=self.wrk_dir / "results",
            )
        return results_df

    def _write_public_results(self, results_df: pd.DataFrame, failures) -> bool:
        public_records = []
        reference_path = self.validate_kwargs.get("reference_path")
        evidence = standard_evidence(
            reference_path=Path(reference_path) if reference_path else None,
            pocket_coverage_reference=self.pocket_coverage_reference,
        )
        requested = set(self.validate_kwargs.get("scoring_functions") or ())
        if self.validate_kwargs.get("reproduction_metrics") is not None:
            requested.add("reproduction_metrics")
        if self.validate_kwargs.get("assess_bias"):
            requested.add("bias_metrics")
        if self._ifp_filter_enabled():
            requested.add("screen_metrics")
        if self.cluster_ifps:
            requested.add("screen_metrics")
        chain_specs = [
            (chain.entity_type, chain.chain_id, chain.sequence_index)
            for chain in iter_system_chains(self.base_system_obj)
        ]
        child_record_cache: dict[Path, tuple[Any, ...]] = {}
        for _, row in results_df.iterrows():
            compound_id = str(row[self.col_id])
            identity = OutputIdentity(
                workflow=WorkflowKind.SCREEN,
                run_id=self.run_id,
                system_id=self.system_path.stem,
                compound_id=compound_id,
                execution_directory=str(row["execution_directory"]),
                runner_id=str(row.get("runner_id") or self.runner),
                runner_version=(
                    str(row["runner_version"])
                    if pd.notna(row.get("runner_version"))
                    else None
                ),
                backend_name=str(row.get("backend_name") or self.runner),
                backend_version_status=(
                    self.execution_plan.backend.version_status
                    if self.execution_plan is not None
                    else None
                ),
                effective_seed=int(row["effective_seed"]),
                repeat_id=int(row["repeat_id"]),
                model_id=str(row["model_id"]),
                sample_id=(
                    int(row["sample_id"])
                    if pd.notna(row.get("sample_id"))
                    else None
                ),
            )
            execution_status = ExecutionStatus(str(row.get("status")))
            execution_error = None
            if execution_status is not ExecutionStatus.SUCCESS:
                execution_error = StructuredExecutionError(
                    stage=FailureStage(str(row.get("error_stage"))),
                    exception_type=str(row.get("exception_type") or "ExecutionUnavailable"),
                    error_code=str(row.get("error_code") or "screen_execution_unavailable"),
                    message=str(row.get("error_message") or "Execution unavailable."),
                    details={
                        "source_format": str(row.get("source_format")),
                        "source_record_index": int(row.get("source_record_index")),
                        "source_record_id": str(row.get("source_record_id")),
                        **(
                            row.get("error_details")
                            if isinstance(row.get("error_details"), dict)
                            else {}
                        ),
                    },
                )
            public_records.append(
                ExecutionRecord(
                    envelope=make_envelope(
                        RecordKind.EXECUTION, identity, "execution"
                    ),
                    status=execution_status,
                    execution_directory=str(row["execution_directory"]),
                    error=execution_error,
                )
            )
            child_records_path = Path(row["run_dir"]) / "results" / "records.jsonl"
            child_records: tuple[Any, ...] = ()
            if child_records_path.is_file():
                if child_records_path not in child_record_cache:
                    child_record_cache[child_records_path] = read_public_records(
                        child_records_path
                    )
                child_records = child_record_cache[child_records_path]
                for child_failure in child_records:
                    if not isinstance(child_failure, WorkflowFailureRecord):
                        continue
                    child_identity = child_failure.envelope.identity
                    if child_identity.repeat_id not in {None, identity.repeat_id}:
                        continue
                    if child_identity.model_id not in {None, identity.model_id}:
                        continue
                    if child_identity.sample_id not in {None, identity.sample_id}:
                        continue
                    rebased_identity = replace(
                        child_identity,
                        workflow=WorkflowKind.SCREEN,
                        run_id=self.run_id,
                        system_id=self.system_path.stem,
                        compound_id=compound_id,
                        execution_directory=str(row["execution_directory"]),
                        runner_id=identity.runner_id,
                        runner_version=identity.runner_version,
                        backend_name=identity.backend_name,
                        backend_version_status=identity.backend_version_status,
                        effective_seed=identity.effective_seed,
                        repeat_id=identity.repeat_id,
                        model_id=identity.model_id,
                        sample_id=identity.sample_id,
                    )
                    public_records.append(
                        replace(
                            child_failure,
                            envelope=make_envelope(
                                RecordKind.FAILURE,
                                rebased_identity,
                                child_failure.stage.value,
                                child_failure.error_code,
                                created_at=child_failure.envelope.created_at,
                            ),
                        )
                    )
            if execution_status is not ExecutionStatus.SUCCESS:
                continue
            preserved_metric_ids: set[str] = set()
            if child_records:
                for child_record in child_records:
                    if not isinstance(child_record, MetricRecord):
                        continue
                    child_identity = child_record.envelope.identity
                    if child_identity.repeat_id != identity.repeat_id:
                        continue
                    if child_identity.sample_id != identity.sample_id:
                        continue
                    rebased_identity = replace(
                        child_identity,
                        workflow=WorkflowKind.SCREEN,
                        run_id=self.run_id,
                        system_id=self.system_path.stem,
                        compound_id=compound_id,
                        execution_directory=str(row["execution_directory"]),
                        runner_id=identity.runner_id,
                        runner_version=identity.runner_version,
                        backend_name=identity.backend_name,
                        backend_version_status=identity.backend_version_status,
                        effective_seed=identity.effective_seed,
                        repeat_id=identity.repeat_id,
                        model_id=identity.model_id,
                        sample_id=identity.sample_id,
                    )
                    rebased_record = replace(
                        child_record,
                        envelope=make_envelope(
                            RecordKind.METRIC,
                            rebased_identity,
                            child_record.metric_name,
                            child_record.statistic,
                            created_at=child_record.envelope.created_at,
                        ),
                    )
                    public_records.append(rebased_record)
                    preserved_metric_ids.add(rebased_record.envelope.record_id)
            system_values = {
                "model_name": row["model_id"],
                "repeat": row["repeat_id"],
                "diffusion_sample": row.get("sample_id"),
                "runner_id": row.get("runner_id"),
                "runner_version": row.get("runner_version"),
                "backend_name": row.get("backend_name"),
                "backend_version_status": row.get("backend_version_status"),
                "effective_seed": row.get("effective_seed"),
            }
            chain_values: list[dict[str, Any]] = []
            for entity_type, chain_id, entity_position in chain_specs:
                prefix = f"{entity_type}_{chain_id}__"
                values = {
                    "CHAIN_ID": chain_id,
                    "ENTITY_ID": f"entity:{entity_position}",
                    "ENTITY_TYPE": entity_type,
                    "model_name": row["model_id"],
                    "repeat": row["repeat_id"],
                    "diffusion_sample": row.get("sample_id"),
                    "runner_id": row.get("runner_id"),
                    "runner_version": row.get("runner_version"),
                    "backend_name": row.get("backend_name"),
                    "backend_version_status": row.get("backend_version_status"),
                    "effective_seed": row.get("effective_seed"),
                }
                for key, value in row.items():
                    if str(key).startswith(prefix):
                        values[str(key).removeprefix(prefix)] = value
                chain_values.append(values)
            for key, value in row.items():
                name = str(key).removeprefix("system__")
                if str(key).startswith("system__") and name in METRIC_CATALOG:
                    system_values[name] = value
                elif str(key) in METRIC_CATALOG:
                    system_values[str(key)] = value
            generated_records = metric_records_from_frames(
                pd.DataFrame([system_values]),
                pd.DataFrame(chain_values),
                base_identity=identity,
                evidence=evidence,
                requested_metrics=requested or None,
                evidence_regime_overrides=(
                    {
                        name: (
                            EvidenceRegime.REFERENCE_STRUCTURE
                            if self._resolved_ifp_filter_source
                            == "reference_complex"
                            else EvidenceRegime.CUSTOM_POCKET
                        )
                        for name in METRIC_CATALOG
                        if name.startswith("ifp_filter_")
                    }
                    if self._ifp_filter_enabled()
                    else None
                ),
            )
            public_records.extend(
                record
                for record in generated_records
                if not isinstance(record, SuccessRecord)
                and record.envelope.record_id not in preserved_metric_ids
            )
        existing_record_ids = {
            record.envelope.record_id for record in public_records
        }
        public_records.extend(
            failure
            for failure in failures
            if failure.envelope.record_id not in existing_record_ids
        )
        has_success = any(
            isinstance(record, ExecutionRecord)
            and record.status is ExecutionStatus.SUCCESS
            for record in public_records
        )
        has_unsuccessful_execution = any(
            isinstance(record, ExecutionRecord)
            and record.status is not ExecutionStatus.SUCCESS
            for record in public_records
        )
        status = (
            "partial"
            if has_success and (failures or has_unsuccessful_execution)
            else "success"
            if has_success
            else "failed"
        )
        backend = self.execution_plan.backend if self.execution_plan else None
        manifest_identity = OutputIdentity(
            workflow=WorkflowKind.SCREEN,
            run_id=self.run_id,
            system_id=self.system_path.stem,
            runner_id=self.runner,
            runner_version=backend.version if backend else None,
            backend_name=backend.backend_name if backend else None,
            backend_version_status=backend.version_status if backend else None,
        )
        bundle = PublicOutputBundle(
            manifest=PublicManifest(
                schema_version=PUBLIC_SCHEMA_VERSION,
                identity=manifest_identity,
                status=status,
                evidence=tuple(evidence),
                requested_metrics=tuple(sorted(requested)),
                artifacts=(
                    ArtifactReference(
                        "executions",
                        "executions.csv",
                        "table",
                        "One terminal record per compound/repeat/model/sample.",
                    ),
                    ArtifactReference(
                        "compound_members",
                        "compound_members.csv",
                        "table",
                        "Stable source-record to execution-member mapping.",
                    ),
                    *(
                        (
                            ArtifactReference(
                                "ifp_cluster_summary",
                                "ifp_cluster_summary.csv",
                                "table",
                            ),
                        )
                        if self.cluster_ifps
                        else ()
                    ),
                ),
                backend=(
                    backend
                ),
                seed_plan=(
                    self.execution_plan.seed_plan
                    if self.execution_plan is not None
                    else None
                ),
            ),
            records=tuple(public_records),
        )
        write_public_bundle(bundle, self.wrk_dir / "results")
        for legacy_name in ("screen_results.csv", "screen_results_with_scores.csv"):
            (self.wrk_dir / legacy_name).unlink(missing_ok=True)
        return has_success

    def _default_cluster_result(self) -> dict[str, Any]:
        return self._screen_postprocessor()._default_cluster_result()

    def _apply_ifp_clustering(
        self,
        summary_df: pd.DataFrame,
        results_df: pd.DataFrame,
    ) -> None:
        self._screen_postprocessor()._apply_ifp_clustering(summary_df, results_df)

    def _default_filter_result(self) -> dict[str, Any]:
        return self._screen_postprocessor()._default_filter_result()

    def _not_evaluable_filter_result(self, reason: str) -> dict[str, Any]:
        return self._screen_postprocessor()._not_evaluable_filter_result(reason)

    def _evaluate_ifp_filter(
        self,
        run_dir: Path,
        *,
        repeat_id: int | None = None,
        sample_id: int | None = None,
    ) -> dict[str, Any]:
        return self._screen_postprocessor()._evaluate_ifp_filter(
            run_dir,
            repeat_id=repeat_id,
            sample_id=sample_id,
        )

    @staticmethod
    def _set_filter_dtypes(df: pd.DataFrame) -> None:
        ScreenPostprocessor._set_filter_dtypes(df)

    def _collect_execution_score_columns(
        self,
        system_df: pd.DataFrame,
        chain_df: pd.DataFrame,
        slot: PlannedExecution,
    ) -> tuple[dict[str, Any], bool]:
        return self._screen_postprocessor()._collect_execution_score_columns(
            system_df,
            chain_df,
            slot,
        )

    def _validate_observed_execution_keys(self, system_df: pd.DataFrame) -> None:
        self._screen_postprocessor()._validate_observed_execution_keys(system_df)

    def _ensure_screen_metric_schema(self, output: dict[str, Any]) -> None:
        self._screen_postprocessor()._ensure_screen_metric_schema(output)
    @staticmethod
    def _parse_list(input_str: str | None) -> list[str]:
        if not input_str:
            return []
        return [v.strip() for v in str(input_str).split(",") if v.strip()]

    @staticmethod
    def _safe_name(value: str) -> str:
        token = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())
        return token or "item"
