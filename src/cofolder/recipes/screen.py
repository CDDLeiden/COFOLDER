"""Screening workflow over canonical CSV, SDF, and MOL library members."""

from __future__ import annotations

import json
import logging
import math
import re
from dataclasses import replace
from pathlib import Path
from typing import Any
from uuid import uuid4

import pandas as pd

from cofolder.modules.analytics.ifp_clustering import (
    SUMMARY_COLUMNS as IFP_CLUSTER_SUMMARY_COLUMNS,
)
from cofolder.modules.analytics.ifp_clustering import cluster_interaction_fingerprints
from cofolder.modules.analytics.ifp_filtering import (
    IFPFilterMode,
    ReferenceIFPFilterPolicy,
    evaluate_reference_ifp_filter,
)
from cofolder.modules.analytics.reference_ifp import (
    IFPExtractionConfig,
    IFPSimilarityMetric,
    IFPTaxonomy,
    InteractionFingerprint,
    InteractionKey,
    LigandIdentity,
    LigandSelector,
    ReferenceEntitySelectionError,
    ReferenceIFPError,
    ResidueIdentity,
    compare_interaction_fingerprints,
    extract_interaction_fingerprint,
    map_reference_identities,
)
from cofolder.modules.analytics.reproduction import (
    _build_reference_ifp_from_custom,
    _parse_custom_pocket_reference,
)
from cofolder.modules.contracts import (
    METRIC_CATALOG,
    PUBLIC_SCHEMA_VERSION,
    SCREEN_METRIC_PROFILES,
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
    ScreenExecutionCardinalityError,
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
from cofolder.recipes._metrics import primary_metric_values, read_metric_frames
from cofolder.recipes._completion import report_completion
from cofolder.recipes._diagnostics import workflow_failure
from cofolder.recipes._execution import build_execution_plan
from cofolder.recipes._results import standard_evidence, write_failure_bundle
from cofolder.recipes.validate import Validate

logger = logging.getLogger(__name__)


SYSTEM_SCREEN_METRICS = SCREEN_METRIC_PROFILES["system"]
PROTEIN_SCREEN_METRICS = SCREEN_METRIC_PROFILES["protein"]
LIGAND_SCREEN_METRICS = SCREEN_METRIC_PROFILES["ligand"]


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

    def _ifp_filter_enabled(self) -> bool:
        return self._ifp_filter_policy_config is not None

    def _prepare_reference_ifp(self) -> None:
        if not self._ifp_filter_enabled() or self._resolved_ifp_filter_source != "reference_complex":
            return
        reference_path = Path(self.validate_kwargs["reference_path"])
        try:
            self._reference_ifp = extract_interaction_fingerprint(
                reference_path,
                ligand=self._reference_ligand_selector,
                receptor_chains=self.ifp_reference_receptor_chains or None,
                config=IFPExtractionConfig(taxonomy=self.ifp_taxonomy),
            )
        except ReferenceEntitySelectionError as exc:
            self._reference_ifp_failure = str(exc)
            return
        except ReferenceIFPError as exc:
            if self.ifp_taxonomy is not IFPTaxonomy.PROLIF:
                raise
            message = str(exc)
            self._reference_ifp_failure = (
                "prolif_worker_timeout" if "timeout" in message else "prolif_worker_crashed"
            )
            return
        try:
            policy = self._ifp_filter_policy_config
            if policy is not None and policy.required_interactions:
                unknown = policy.required_interactions - self._reference_ifp.interactions
                if unknown:
                    raise ValueError(
                        "Required interactions are absent from the reference fingerprint: "
                        + ", ".join(str(item) for item in sorted(unknown))
                    )
        except AttributeError as exc:
            raise ValueError("Reference IFP was not prepared.") from exc

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
        return {
            "ifp_cluster_id": None,
            "ifp_cluster_status": (
                "not_evaluable" if self.cluster_ifps else "not_applied"
            ),
        }

    def _apply_ifp_clustering(
        self,
        summary_df: pd.DataFrame,
        results_df: pd.DataFrame,
    ) -> None:
        """Annotate both outputs and write a deterministic cluster summary."""

        if not self.cluster_ifps:
            (self.wrk_dir / "results" / "ifp_cluster_summary.csv").unlink(
                missing_ok=True
            )
            return

        parsed_rows: list[tuple[int, InteractionFingerprint, str]] = []
        legacy_width: int | None = None
        for position, row in summary_df.iterrows():
            if row.get("status") != "success":
                continue
            fingerprint, _ = self._load_selected_interaction_fingerprint(
                Path(row["run_dir"]),
                repeat_id=int(row["repeat_id"]),
                sample_id=(
                    int(row["sample_id"])
                    if pd.notna(row.get("sample_id"))
                    else None
                ),
            )
            if fingerprint is None and self.ifp_taxonomy is IFPTaxonomy.DISTANCE:
                vector, _ = self._load_selected_ligand_ifp(
                    Path(row["run_dir"]),
                    repeat_id=int(row["repeat_id"]),
                    sample_id=(
                        int(row["sample_id"])
                        if pd.notna(row.get("sample_id"))
                        else None
                    ),
                )
                if vector is not None:
                    if legacy_width is None:
                        legacy_width = len(vector)
                    if len(vector) != legacy_width:
                        continue
                    fingerprint = InteractionFingerprint(
                        taxonomy=IFPTaxonomy.DISTANCE,
                        ligand=LigandIdentity(self.ligand_chain, 0),
                        receptor_chains=("_legacy",),
                        interactions=frozenset(
                            InteractionKey(
                                ResidueIdentity("_legacy", index + 1),
                                "distance_contact",
                            )
                            for index, active in enumerate(vector)
                            if active
                        ),
                        source_path=Path(row["run_dir"]),
                    )
            if fingerprint is None:
                continue
            parsed_rows.append((position, fingerprint, str(row["execution_key"])))

        if parsed_rows:
            clustered = cluster_interaction_fingerprints(
                [row[1] for row in parsed_rows],
                [row[2] for row in parsed_rows],
                similarity_threshold=self.ifp_cluster_similarity_threshold,
            )
            for (position, _, _), cluster_id in zip(parsed_rows, clustered.cluster_ids):
                for frame in (summary_df, results_df):
                    frame.at[position, "ifp_cluster_id"] = cluster_id
                    frame.at[position, "ifp_cluster_status"] = "clustered"
            cluster_summary = clustered.summary.copy()
            cluster_summary["ifp_taxonomy"] = self.ifp_taxonomy.value
            if self._ifp_filter_enabled():
                annotations = []
                for cluster_id in cluster_summary["ifp_cluster_id"]:
                    members = summary_df[summary_df["ifp_cluster_id"] == cluster_id]
                    similarities = pd.to_numeric(
                        members.get("ifp_filter_similarity", pd.Series(dtype=float)),
                        errors="coerce",
                    ).dropna()
                    accepted = members[
                        members.get("ifp_filter_status", pd.Series(index=members.index, dtype=object))
                        == "accepted"
                    ]
                    annotations.append({
                        "ifp_cluster_id": cluster_id,
                        "reference_evaluable_count": int(similarities.size),
                        "reference_accepted_count": len(accepted),
                        "reference_accepted_member_ids": json.dumps(
                            accepted[self.col_id].astype(str).tolist()
                        ),
                        "mean_reference_similarity": (
                            float(similarities.mean()) if not similarities.empty else None
                        ),
                        "max_reference_similarity": (
                            float(similarities.max()) if not similarities.empty else None
                        ),
                    })
                cluster_summary = cluster_summary.merge(
                    pd.DataFrame(annotations), on="ifp_cluster_id", how="left"
                )
        else:
            cluster_summary = pd.DataFrame(columns=IFP_CLUSTER_SUMMARY_COLUMNS)

        public_results_dir = self.wrk_dir / "results"
        public_results_dir.mkdir(parents=True, exist_ok=True)
        cluster_summary.to_csv(
            public_results_dir / "ifp_cluster_summary.csv", index=False
        )

    def _load_selected_interaction_fingerprint(
        self,
        run_dir: Path,
        *,
        repeat_id: int | None = None,
        sample_id: int | None = None,
    ) -> tuple[InteractionFingerprint | None, str]:
        run_dir = Path(run_dir)
        cache_key = (run_dir, repeat_id, sample_id)
        if cache_key in self._prediction_ifp_cache:
            return self._prediction_ifp_cache[cache_key]
        _, chain_df = read_metric_frames(run_dir)
        if chain_df.empty or not {"CHAIN_ID", "ENTITY_TYPE", "cif_file"}.issubset(chain_df.columns):
            result = (None, "missing_chain_metrics")
            self._prediction_ifp_cache[cache_key] = result
            return result
        selected = self._select_chain_metrics_rows(
            chain_df, repeat_id=repeat_id, sample_id=sample_id
        )
        ligand_rows = selected[
            (selected["ENTITY_TYPE"].astype(str).str.lower() == "ligand")
            & (selected["CHAIN_ID"].astype(str) == self.ligand_chain)
        ]
        if ligand_rows.empty:
            result = (None, "ligand_chain_not_found")
            self._prediction_ifp_cache[cache_key] = result
            return result
        cif_name = ligand_rows.iloc[0].get("cif_file")
        path = run_dir / "results" / "structures" / str(cif_name)
        receptor_chains = tuple(
            selected.loc[
                selected["ENTITY_TYPE"].astype(str).str.lower() == "protein", "CHAIN_ID"
            ].dropna().astype(str).unique()
        )
        try:
            result = (extract_interaction_fingerprint(
                path,
                ligand=LigandSelector(chain_id=self.ligand_chain),
                receptor_chains=receptor_chains or None,
                config=IFPExtractionConfig(taxonomy=self.ifp_taxonomy),
            ), "")
        except ReferenceIFPError as exc:
            result = (None, str(exc))
        self._prediction_ifp_cache[cache_key] = result
        return result

    def _load_selected_ligand_ifp(
        self,
        run_dir: Path,
        *,
        repeat_id: int | None = None,
        sample_id: int | None = None,
    ) -> tuple[list[int] | None, str]:
        _, chain_df = read_metric_frames(run_dir)
        if chain_df.empty:
            return None, "missing_chain_metrics"
        required = {"CHAIN_ID", "ENTITY_TYPE", "ifp_distance"}
        if chain_df.empty or not required.issubset(chain_df.columns):
            return None, "missing_ifp"
        chain_df = self._select_chain_metrics_rows(
            chain_df, repeat_id=repeat_id, sample_id=sample_id
        )
        ligand_rows = chain_df[
            chain_df["ENTITY_TYPE"].astype(str).str.lower() == "ligand"
        ]
        if self.ligand_chain:
            ligand_rows = ligand_rows[
                ligand_rows["CHAIN_ID"].astype(str) == self.ligand_chain
            ]
        elif ligand_rows["CHAIN_ID"].astype(str).nunique() != 1:
            return None, "ambiguous_ligand_chain"
        if ligand_rows.empty:
            return None, "ligand_chain_not_found"
        return self._parse_binary_ifp(ligand_rows.iloc[0].get("ifp_distance"))

    def _default_filter_result(self) -> dict[str, Any]:
        if self._ifp_filter_enabled():
            return self._not_evaluable_filter_result("missing_ifp")
        return {
            "ifp_filter_pass": pd.NA,
            "ifp_filter_status": "not_applied",
            "ifp_filter_reason": "filtering_disabled",
            "ifp_filter_overlap": None,
            "ifp_filter_similarity": None,
            "ifp_filter_threshold": None,
            "ifp_filter_reference": None,
            "ifp_filter_similarity_metric": None,
            "ifp_filter_policy": None,
            "ifp_filter_taxonomy": None,
            "ifp_filter_required_interactions": None,
            "ifp_filter_missing_interactions": None,
            "ifp_filter_mapping_status": None,
            "ifp_filter_mapping_failures": None,
        }

    def _not_evaluable_filter_result(self, reason: str) -> dict[str, Any]:
        return {
            "ifp_filter_pass": pd.NA,
            "ifp_filter_status": "not_evaluable",
            "ifp_filter_reason": reason,
            "ifp_filter_overlap": None,
            "ifp_filter_similarity": None,
            "ifp_filter_threshold": self.ifp_filter_threshold,
            "ifp_filter_reference": (
                str(self.validate_kwargs.get("reference_path"))
                if self._resolved_ifp_filter_source == "reference_complex"
                else self.pocket_coverage_reference
            ),
            "ifp_filter_similarity_metric": (
                self._ifp_filter_policy_config.similarity_metric.value
                if self._ifp_filter_policy_config else None
            ),
            "ifp_filter_policy": (
                self._ifp_filter_policy_config.mode.value
                if self._ifp_filter_policy_config else None
            ),
            "ifp_filter_taxonomy": self.ifp_taxonomy.value,
            "ifp_filter_required_interactions": json.dumps(
                [str(item) for item in sorted(
                    self._ifp_filter_policy_config.required_interactions or ()
                )]
            ) if self._ifp_filter_policy_config else None,
            "ifp_filter_missing_interactions": json.dumps([]),
            "ifp_filter_mapping_status": "unmappable",
            "ifp_filter_mapping_failures": json.dumps([reason]),
        }

    def _evaluate_ifp_filter(
        self,
        run_dir: Path,
        *,
        repeat_id: int | None = None,
        sample_id: int | None = None,
    ) -> dict[str, Any]:
        """Evaluate strict reference overlap for one completed execution."""
        if not self._ifp_filter_enabled():
            return self._default_filter_result()
        if self._resolved_ifp_filter_source == "reference_complex":
            return self._evaluate_reference_complex_filter(
                run_dir, repeat_id=repeat_id, sample_id=sample_id
            )

        _, chain_df = read_metric_frames(run_dir)
        if chain_df.empty:
            return self._not_evaluable_filter_result("missing_chain_metrics")
        required = {"CHAIN_ID", "ENTITY_TYPE", "ifp_distance"}
        if chain_df.empty or not required.issubset(chain_df.columns):
            return self._not_evaluable_filter_result("missing_ifp")

        chain_df = self._select_chain_metrics_rows(
            chain_df, repeat_id=repeat_id, sample_id=sample_id
        )
        ligand_rows = chain_df[
            chain_df["ENTITY_TYPE"].astype(str).str.lower() == "ligand"
        ]
        if self.ligand_chain:
            ligand_rows = ligand_rows[
                ligand_rows["CHAIN_ID"].astype(str) == self.ligand_chain
            ]
        elif ligand_rows["CHAIN_ID"].astype(str).nunique() != 1:
            return self._not_evaluable_filter_result("ambiguous_ligand_chain")
        if ligand_rows.empty:
            return self._not_evaluable_filter_result("ligand_chain_not_found")

        ligand_row = ligand_rows.iloc[0]
        pred_ifp, parse_reason = self._parse_binary_ifp(ligand_row.get("ifp_distance"))
        if pred_ifp is None:
            return self._not_evaluable_filter_result(parse_reason)

        ref_ifp = self._resolve_filter_reference(
            run_dir=run_dir,
            chain_df=chain_df,
            ligand_row=ligand_row,
        )
        if ref_ifp is None:
            return self._not_evaluable_filter_result("reference_resolution_failed")
        if len(pred_ifp) != len(ref_ifp):
            return self._not_evaluable_filter_result("incompatible_vector_lengths")
        if not any(ref_ifp):
            return self._not_evaluable_filter_result("empty_reference")

        intersection = sum(p and r for p, r in zip(pred_ifp, ref_ifp))
        overlap = intersection / sum(ref_ifp)
        union = sum(bool(p or r) for p, r in zip(pred_ifp, ref_ifp))
        jaccard = intersection / union if union else 0.0
        assert self._ifp_filter_policy_config is not None
        score = (
            jaccard
            if self._ifp_filter_policy_config.similarity_metric is IFPSimilarityMetric.JACCARD
            else overlap
        )
        passed = score >= float(self.ifp_filter_threshold)
        return {
            "ifp_filter_pass": bool(passed),
            "ifp_filter_status": "accepted" if passed else "rejected",
            "ifp_filter_reason": "threshold_met" if passed else "below_threshold",
            "ifp_filter_overlap": float(overlap),
            "ifp_filter_similarity": float(score),
            "ifp_filter_threshold": self.ifp_filter_threshold,
            "ifp_filter_reference": self.pocket_coverage_reference,
            "ifp_filter_similarity_metric": self._ifp_filter_policy_config.similarity_metric.value,
            "ifp_filter_policy": self._ifp_filter_policy_config.mode.value,
            "ifp_filter_taxonomy": self.ifp_taxonomy.value,
            "ifp_filter_required_interactions": json.dumps([]),
            "ifp_filter_missing_interactions": json.dumps([]),
            "ifp_filter_mapping_status": "not_applicable",
            "ifp_filter_mapping_failures": json.dumps([]),
        }

    def _evaluate_reference_complex_filter(
        self,
        run_dir: Path,
        *,
        repeat_id: int | None = None,
        sample_id: int | None = None,
    ) -> dict[str, Any]:
        if self._reference_ifp_failure:
            return self._not_evaluable_filter_result(self._reference_ifp_failure)
        if self._reference_ifp is None or self._ifp_filter_policy_config is None:
            return self._not_evaluable_filter_result("reference_ifp_extraction_failed")
        try:
            prediction, extraction_reason = self._load_selected_interaction_fingerprint(
                run_dir, repeat_id=repeat_id, sample_id=sample_id
            )
            if prediction is None:
                reason = (
                    "prolif_worker_timeout" if "timeout" in extraction_reason
                    else "prolif_worker_crashed" if "crashed" in extraction_reason
                    else extraction_reason or "prediction_ifp_extraction_failed"
                )
                return self._not_evaluable_filter_result(reason)
            predicted_path = prediction.source_path
            mapping = map_reference_identities(
                self._reference_ifp,
                prediction,
                reference_structure_path=Path(self.validate_kwargs["reference_path"]),
                predicted_structure_path=predicted_path,
            )
            comparison = compare_interaction_fingerprints(
                self._reference_ifp, prediction, mapping
            )
            outcome = evaluate_reference_ifp_filter(
                self._reference_ifp, comparison, self._ifp_filter_policy_config
            )
        except ReferenceIFPError as exc:
            message = str(exc)
            reason = (
                "prolif_worker_timeout" if "prolif_worker_timeout" in message
                else "prolif_worker_crashed" if "prolif_worker_crashed" in message
                else "prediction_ifp_extraction_failed"
            )
            return self._not_evaluable_filter_result(reason)

        selected_score = comparison.similarities.get(
            self._ifp_filter_policy_config.similarity_metric
        )
        return {
            "ifp_filter_pass": outcome.passed if outcome.passed is not None else pd.NA,
            "ifp_filter_status": outcome.status,
            "ifp_filter_reason": outcome.reason,
            "ifp_filter_overlap": comparison.similarities.get(
                IFPSimilarityMetric.REFERENCE_COVERAGE
            ),
            "ifp_filter_similarity": selected_score,
            "ifp_filter_threshold": self._ifp_filter_policy_config.threshold,
            "ifp_filter_reference": str(self.validate_kwargs["reference_path"]),
            "ifp_filter_similarity_metric": self._ifp_filter_policy_config.similarity_metric.value,
            "ifp_filter_policy": self._ifp_filter_policy_config.mode.value,
            "ifp_filter_taxonomy": self.ifp_taxonomy.value,
            "ifp_filter_required_interactions": json.dumps([
                str(item) for item in sorted(
                    self._ifp_filter_policy_config.required_interactions
                    or self._reference_ifp.interactions
                    if self._ifp_filter_policy_config.mode is IFPFilterMode.REQUIRED
                    else ()
                )
            ]),
            "ifp_filter_missing_interactions": json.dumps([
                str(item) for item in comparison.missing_interactions
            ]),
            "ifp_filter_mapping_status": comparison.mapping.status.value,
            "ifp_filter_mapping_failures": json.dumps(list(comparison.mapping.failures)),
        }

    def _resolve_filter_reference(
        self,
        run_dir: Path,
        chain_df: pd.DataFrame,
        ligand_row: pd.Series,
    ) -> list[int] | None:
        spec = self._ifp_filter_reference_spec
        if spec is None:
            return None
        if spec["kind"] == "bits":
            return list(spec["bits"])

        cif_name = ligand_row.get("cif_file")
        if pd.isna(cif_name) or "cif_file" not in chain_df.columns:
            return None
        model_rows = chain_df[chain_df["cif_file"] == cif_name]
        receptor_rows = model_rows[
            model_rows["ENTITY_TYPE"].astype(str).str.lower() == "protein"
        ]
        receptor_ids = receptor_rows["CHAIN_ID"].dropna().astype(str).unique()
        if len(receptor_ids) != 1:
            return None

        structure_path = run_dir / "results" / "structures" / str(cif_name)
        if not structure_path.exists():
            return None
        try:
            import gemmi

            model = gemmi.read_structure(str(structure_path))[0]
            receptor_chain = model[receptor_ids[0]]
            residues = list(receptor_chain)
            residue_order = [(" ", int(res.seqid.num), " ") for res in residues]
            residue_names = {
                residue_id: str(res.name)
                for residue_id, res in zip(residue_order, residues)
            }
        except Exception as exc:
            self.logger.warning(
                "Unable to resolve filter reference from %s: %s", structure_path, exc
            )
            return None

        return _build_reference_ifp_from_custom(
            spec,
            residue_order=residue_order,
            residue_names=residue_names,
            logger=self.logger,
            warn_key=f"cif={cif_name},chain={ligand_row.get('CHAIN_ID')}",
        )

    @staticmethod
    def _parse_binary_ifp(value: Any) -> tuple[list[int] | None, str]:
        if value is None or (not isinstance(value, (list, tuple)) and pd.isna(value)):
            return None, "missing_ifp"
        parsed = value
        if isinstance(value, str):
            if not value.strip():
                return None, "missing_ifp"
            try:
                parsed = json.loads(value)
            except (TypeError, json.JSONDecodeError):
                return None, "malformed_ifp"
        if not isinstance(parsed, (list, tuple)):
            return None, "malformed_ifp"
        if not parsed:
            return None, "missing_ifp"
        if any(item not in (0, 1, False, True) for item in parsed):
            return None, "malformed_ifp"
        return [int(item) for item in parsed], ""

    @staticmethod
    def _set_filter_dtypes(df: pd.DataFrame) -> None:
        if "ifp_filter_pass" in df.columns:
            df["ifp_filter_pass"] = df["ifp_filter_pass"].astype("boolean")

    def _collect_score_columns(self, run_dir: Path) -> dict[str, Any]:
        """Collect computed scores from per-row validate outputs."""
        out: dict[str, Any] = {}
        try:
            system_df, chain_df = read_metric_frames(run_dir)
            out.update(primary_metric_values(system_df, chain_df))
        except (OSError, UnicodeError, pd.errors.ParserError) as exc:
            self.logger.warning("Failed reading metrics from %s: %s", run_dir, exc)

        self._ensure_screen_metric_schema(out)
        return out

    def _collect_execution_score_columns(
        self,
        system_df: pd.DataFrame,
        chain_df: pd.DataFrame,
        slot: PlannedExecution,
    ) -> tuple[dict[str, Any], bool]:
        """Collect metrics for exactly one planned repeat/model/sample slot."""

        def select(frame: pd.DataFrame) -> pd.DataFrame:
            if frame.empty:
                return frame
            mask = pd.Series(True, index=frame.index)
            if "repeat" in frame:
                mask &= pd.to_numeric(frame["repeat"], errors="coerce").eq(
                    slot.repeat_id
                )
            if "diffusion_sample" in frame and slot.sample_id is not None:
                mask &= pd.to_numeric(
                    frame["diffusion_sample"], errors="coerce"
                ).eq(slot.sample_id)
            if "model_name" in frame:
                observed_models = set(frame["model_name"].dropna().astype(str))
                if slot.model_id in observed_models:
                    mask &= frame["model_name"].astype(str).eq(slot.model_id)
            return frame.loc[mask]

        selected_system = select(system_df)
        selected_chain = select(chain_df)
        if len(selected_system) > 1:
            raise ScreenExecutionCardinalityError(
                "Validate produced duplicate system rows for execution "
                f"repeat={slot.repeat_id}, model={slot.model_id!r}, "
                f"sample={slot.sample_id!r}."
            )
        observed = not selected_system.empty or not selected_chain.empty
        values = (
            primary_metric_values(selected_system, selected_chain)
            if observed
            else {}
        )
        self._ensure_screen_metric_schema(values)
        return values, observed

    def _validate_observed_execution_keys(self, system_df: pd.DataFrame) -> None:
        """Reject normalized system outputs outside the planned execution matrix."""
        if system_df.empty or self.execution_plan is None:
            return
        required = {"repeat", "diffusion_sample"}
        if not required.issubset(system_df.columns):
            raise ScreenExecutionCardinalityError(
                "Validate system output lacks repeat/diffusion_sample identity."
            )
        planned = {
            (slot.repeat_id, slot.sample_id) for slot in self.execution_plan.executions
        }
        observed = {
            (int(row["repeat"]), int(row["diffusion_sample"]))
            for _, row in system_df.iterrows()
        }
        extras = sorted(observed - planned)
        if extras:
            raise ScreenExecutionCardinalityError(
                f"Validate produced executions outside the planned matrix: {extras}."
            )

    def _ensure_screen_metric_schema(self, output: dict[str, Any]) -> None:
        """Populate stable manuscript-facing score columns, using nulls when unavailable."""

        for column in SYSTEM_SCREEN_METRICS:
            output.setdefault(f"system__{column}", None)

        chains = list(iter_system_chains(self.base_system_obj))
        protein_ids = [
            chain.chain_id for chain in chains if chain.entity_type == "protein"
        ]
        for chain in chains:
            prefix = f"{chain.entity_type}_{chain.chain_id}"
            metrics: tuple[str, ...] = ()
            if chain.entity_type == "protein":
                metrics = PROTEIN_SCREEN_METRICS
            elif chain.entity_type == "ligand":
                metrics = LIGAND_SCREEN_METRICS
            for column in metrics:
                output.setdefault(f"{prefix}__{column}", None)
            if chain.entity_type == "ligand":
                for protein_id in protein_ids:
                    output.setdefault(
                        f"{prefix}__pair_chains_iptm_{protein_id}",
                        None,
                    )

    @staticmethod
    def _select_metrics_row(df: pd.DataFrame) -> pd.Series:
        """Prefer repeat=1/sample=0 row if available; otherwise first row."""
        if {"repeat", "diffusion_sample"}.issubset(df.columns):
            sub = df[(df["repeat"] == 1) & (df["diffusion_sample"] == 0)]
            if not sub.empty:
                return sub.iloc[0]
        return df.iloc[0]

    @staticmethod
    def _select_chain_metrics_rows(
        df: pd.DataFrame,
        *,
        repeat_id: int | None = None,
        sample_id: int | None = None,
    ) -> pd.DataFrame:
        """Select one execution's chain rows, retaining legacy primary selection."""
        if repeat_id is not None and "repeat" in df.columns:
            selected = df[pd.to_numeric(df["repeat"], errors="coerce").eq(repeat_id)]
            if sample_id is not None and "diffusion_sample" in selected.columns:
                selected = selected[
                    pd.to_numeric(
                        selected["diffusion_sample"], errors="coerce"
                    ).eq(sample_id)
                ]
            return selected
        if {"repeat", "diffusion_sample"}.issubset(df.columns):
            sub = df[(df["repeat"] == 1) & (df["diffusion_sample"] == 0)]
            if not sub.empty:
                return sub
        if "CHAIN_ID" in df.columns:
            return df.drop_duplicates(subset=["CHAIN_ID"], keep="first")
        return df

    @staticmethod
    def _parse_list(input_str: str | None) -> list[str]:
        if not input_str:
            return []
        return [v.strip() for v in str(input_str).split(",") if v.strip()]

    @staticmethod
    def _safe_name(value: str) -> str:
        token = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())
        return token or "item"
