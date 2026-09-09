"""Screening workflow over canonical CSV, SDF, and MOL library members."""

from __future__ import annotations

import json
import logging
import math
import re
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
    EvidenceSource,
    FailureStage,
    OutputIdentity,
    PublicManifest,
    PublicOutputBundle,
    PublicSerializationError,
    SuccessRecord,
    WorkflowExecutionError,
    WorkflowKind,
    failure_from_exception,
    metric_records_from_frames,
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
from cofolder.modules.runners import get_runner
from cofolder.modules.runners.msa import resolve_declared_msa_paths
from cofolder.modules.utils import write
from cofolder.recipes._metrics import primary_metric_values, read_metric_frames
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
            Path, tuple[InteractionFingerprint | None, str]
        ] = {}

        self.validate_kwargs: dict[str, Any] = {
            "repeats": repeats,
            "seed": seed,
            "scoring_functions": scoring_functions,
            "assess_robustness": assess_robustness,
            "assess_bias": assess_bias,
            "protein_training_data_path": protein_training_data_path,
            "ligand_training_data_path": ligand_training_data_path,
            "bias_release_cutoff": bias_release_cutoff,
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
        self._failure_stage = FailureStage.INPUT_VALIDATION
        self.reusable_msa_dir = (
            self.wrk_dir / "shared" / "msa" / self.runner
            if getattr(self.runner_impl, "supports_msa_reuse", False)
            else None
        )
        self.msa_reuse_settings = {}
        self.logger = logging.getLogger("cofolder.screen")


    def _validate_config(self) -> None:
        if self.library is None:
            raise ValueError("--library is required.")
        if not self.library.exists():
            raise ValueError(f"--library does not exist: {self.library}")
        if not self.library.is_file():
            raise ValueError(f"--library is not a file: {self.library}")
        if not self.ligand_chain:
            raise ValueError("--ligand_chain is required.")

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

    def run(self) -> pd.DataFrame:
        try:
            return self._run_impl()
        except (WorkflowExecutionError, PublicSerializationError):
            raise
        except Exception as exc:
            identity = OutputIdentity(
                workflow=WorkflowKind.SCREEN,
                run_id=self.run_id,
                system_id=self.system_path.stem,
                runner_id=self.runner,
            )
            failure = failure_from_exception(
                exc,
                identity=identity,
                stage=self._failure_stage,
                error_code=getattr(exc, "error_code", "screen_input_validation_failed"),
            )
            output_dir = self.wrk_dir / "results"
            write_public_bundle(
                PublicOutputBundle(
                    manifest=PublicManifest(
                        schema_version=PUBLIC_SCHEMA_VERSION,
                        identity=identity,
                        status="failed",
                    ),
                    records=(failure,),
                ),
                output_dir,
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
            self.base_system_obj, self.ligand_chain
        )
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
        valid_members = [
            outcome
            for outcome in compound_library.outcomes
            if isinstance(outcome, CompoundMember)
        ]
        if self.validate_kwargs.get("conformers") == "sdf":
            if compound_library.source_format is CompoundLibraryFormat.CSV:
                raise ValueError("--conformers sdf requires an SDF/MOL screening library.")
            capabilities = self.runner_impl.ligand_preparation_capabilities
            if "sdf" not in capabilities.conformer_modes:
                raise ValueError(
                    f"Runner '{self.runner}' does not support SDF ligand conformers."
                )

        if valid_members:
            self._failure_stage = FailureStage.PREPARATION
            self.runner_impl.ensure_available()
            self._failure_stage = FailureStage.INPUT_VALIDATION

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
                detailed.update(self._collect_score_columns(run_dir=run_dir))
                filter_result = self._evaluate_ifp_filter(run_dir=run_dir)
                summary.update(filter_result)
                detailed.update(filter_result)
            except Exception as exc:
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
                public_failures.append(
                    failure_from_exception(
                        source_exc,
                        identity=OutputIdentity(
                            workflow=WorkflowKind.SCREEN,
                            run_id=self.run_id,
                            system_id=self.system_path.stem,
                            compound_id=compound_id,
                            runner_id=self.runner,
                        ),
                        stage=row_stage,
                        error_code=getattr(
                            source_exc, "error_code", "screen_compound_failed"
                        ),
                        details={
                            "row_index": i,
                            "source_format": source.source_format.value,
                            "source_path": str(source.source_path),
                            "source_record_index": i,
                            "source_record_id": source.source_record_id,
                            "original_id": source.original_id,
                            "execution_id": compound_id,
                            "execution_directory": outcome.execution_directory,
                        },
                    )
                )

            detailed["coordinate_mode"] = coordinate_mode
            summary["coordinate_mode"] = coordinate_mode
            records.append(summary)
            records_with_scores.append(detailed)
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

        failures = sum(1 for r in records if r["status"] == "failed")
        successes = len(records) - failures
        self.logger.info(
            "Screen complete: total=%d success=%d failed=%d public_results=%s",
            len(records),
            successes,
            failures,
            self.wrk_dir / "results",
        )
        if not has_success:
            raise WorkflowExecutionError(
                "No screened compound produced a usable result.",
                failures=tuple(public_failures),
                output_dir=self.wrk_dir / "results",
            )
        return results_df

    def _write_public_results(self, results_df: pd.DataFrame, failures) -> bool:
        public_records = []
        evidence = []
        reference_path = self.validate_kwargs.get("reference_path")
        if reference_path:
            reference = Path(reference_path)
            evidence.append(
                EvidenceSource(
                    kind="reference_structure",
                    identifier=reference.name,
                    path=str(reference),
                )
            )
        if self.pocket_coverage_reference:
            evidence.append(
                EvidenceSource(
                    kind="custom_pocket",
                    identifier="configured_custom_pocket",
                )
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
        for _, row in results_df.iterrows():
            compound_id = str(row[self.col_id])
            identity = OutputIdentity(
                workflow=WorkflowKind.SCREEN,
                run_id=self.run_id,
                system_id=self.system_path.stem,
                compound_id=compound_id,
                runner_id=self.runner,
            )
            if str(row.get("status")) == "failed":
                continue
            system_values = {
                "model_name": None,
                "repeat": None,
                "diffusion_sample": None,
            }
            chain_values: list[dict[str, Any]] = []
            for entity_type, chain_id, entity_position in chain_specs:
                prefix = f"{entity_type}_{chain_id}__"
                values = {
                    "CHAIN_ID": chain_id,
                    "ENTITY_ID": f"entity:{entity_position}",
                    "ENTITY_TYPE": entity_type,
                    "model_name": None,
                    "repeat": None,
                    "diffusion_sample": None,
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
            public_records.extend(
                metric_records_from_frames(
                    pd.DataFrame([system_values]),
                    pd.DataFrame(chain_values),
                    base_identity=identity,
                    evidence=evidence,
                    requested_metrics=requested or None,
                    evidence_regime_overrides={
                        name: (
                            EvidenceRegime.REFERENCE_STRUCTURE
                            if self._resolved_ifp_filter_source == "reference_complex"
                            else EvidenceRegime.CUSTOM_POCKET
                        )
                        for name in METRIC_CATALOG
                        if name.startswith("ifp_filter_")
                    } if self._ifp_filter_enabled() else None,
                )
            )
        public_records.extend(failures)
        has_success = any(
            isinstance(record, SuccessRecord) for record in public_records
        )
        status = (
            "partial"
            if failures and has_success
            else "success"
            if has_success
            else "failed"
        )
        manifest_identity = OutputIdentity(
            workflow=WorkflowKind.SCREEN,
            run_id=self.run_id,
            system_id=self.system_path.stem,
            runner_id=self.runner,
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
                Path(row["run_dir"])
            )
            if fingerprint is None and self.ifp_taxonomy is IFPTaxonomy.DISTANCE:
                vector, _ = self._load_selected_ligand_ifp(Path(row["run_dir"]))
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
            parsed_rows.append((position, fingerprint, str(row[self.col_id])))

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
        self, run_dir: Path
    ) -> tuple[InteractionFingerprint | None, str]:
        run_dir = Path(run_dir)
        if run_dir in self._prediction_ifp_cache:
            return self._prediction_ifp_cache[run_dir]
        _, chain_df = read_metric_frames(run_dir)
        if chain_df.empty or not {"CHAIN_ID", "ENTITY_TYPE", "cif_file"}.issubset(chain_df.columns):
            result = (None, "missing_chain_metrics")
            self._prediction_ifp_cache[run_dir] = result
            return result
        selected = self._select_chain_metrics_rows(chain_df)
        ligand_rows = selected[
            (selected["ENTITY_TYPE"].astype(str).str.lower() == "ligand")
            & (selected["CHAIN_ID"].astype(str) == self.ligand_chain)
        ]
        if ligand_rows.empty:
            result = (None, "ligand_chain_not_found")
            self._prediction_ifp_cache[run_dir] = result
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
        self._prediction_ifp_cache[run_dir] = result
        return result

    def _load_selected_ligand_ifp(self, run_dir: Path) -> tuple[list[int] | None, str]:
        _, chain_df = read_metric_frames(run_dir)
        if chain_df.empty:
            return None, "missing_chain_metrics"
        required = {"CHAIN_ID", "ENTITY_TYPE", "ifp_distance"}
        if chain_df.empty or not required.issubset(chain_df.columns):
            return None, "missing_ifp"
        chain_df = self._select_chain_metrics_rows(chain_df)
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

    def _evaluate_ifp_filter(self, run_dir: Path) -> dict[str, Any]:
        """Evaluate strict reference overlap for one completed screen row."""
        if not self._ifp_filter_enabled():
            return self._default_filter_result()
        if self._resolved_ifp_filter_source == "reference_complex":
            return self._evaluate_reference_complex_filter(run_dir)

        _, chain_df = read_metric_frames(run_dir)
        if chain_df.empty:
            return self._not_evaluable_filter_result("missing_chain_metrics")
        required = {"CHAIN_ID", "ENTITY_TYPE", "ifp_distance"}
        if chain_df.empty or not required.issubset(chain_df.columns):
            return self._not_evaluable_filter_result("missing_ifp")

        chain_df = self._select_chain_metrics_rows(chain_df)
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

    def _evaluate_reference_complex_filter(self, run_dir: Path) -> dict[str, Any]:
        if self._reference_ifp_failure:
            return self._not_evaluable_filter_result(self._reference_ifp_failure)
        if self._reference_ifp is None or self._ifp_filter_policy_config is None:
            return self._not_evaluable_filter_result("reference_ifp_extraction_failed")
        try:
            prediction, extraction_reason = self._load_selected_interaction_fingerprint(
                run_dir
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
    def _select_chain_metrics_rows(df: pd.DataFrame) -> pd.DataFrame:
        """Prefer repeat=1/sample=0 chain rows when available; otherwise de-duplicate by chain."""
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
