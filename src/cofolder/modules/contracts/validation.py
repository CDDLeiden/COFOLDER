from __future__ import annotations

import math
import re
from collections import defaultdict
from collections.abc import Collection, Mapping
from dataclasses import replace
from datetime import UTC, datetime
from statistics import mean, stdev
from uuid import NAMESPACE_URL, uuid5

from packaging.version import InvalidVersion, Version

from .metrics import get_metric_definition, validate_metric_record
from .models import (
    PUBLIC_SCHEMA_VERSION,
    AmbiguousIdentityError,
    ArtifactReference,
    BackendVersionStatus,
    EvidenceRegime,
    EvidenceSource,
    ExecutionRecord,
    ExecutionStatus,
    FailureStage,
    MetricClass,
    MetricRecord,
    OptimizationDirection,
    OutputIdentity,
    PublicOutputBundle,
    PublicRecord,
    PublicRecordEnvelope,
    PublicSchemaValidationError,
    RecordKind,
    RecordStatus,
    RunnerBackendIdentity,
    RunnerProvenanceError,
    SeedAdjustment,
    SeedOrigin,
    SeedPlan,
    StructuredExecutionError,
    SuccessRecord,
    WorkflowExecutionError,
    WorkflowFailureRecord,
    WorkflowKind,
)


def utc_now() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


def make_record_id(
    record_kind: RecordKind,
    identity: OutputIdentity,
    *discriminators: object,
) -> str:
    parts = [record_kind.value]
    parts.extend(
        str(getattr(identity, field)) for field in identity.__dataclass_fields__
    )
    parts.extend(str(value) for value in discriminators)
    return str(uuid5(NAMESPACE_URL, "cofolder:" + "|".join(parts)))


def make_envelope(
    record_kind: RecordKind,
    identity: OutputIdentity,
    *discriminators: object,
    created_at: str | None = None,
) -> PublicRecordEnvelope:
    return PublicRecordEnvelope(
        schema_version=PUBLIC_SCHEMA_VERSION,
        record_id=make_record_id(record_kind, identity, *discriminators),
        record_kind=record_kind,
        identity=identity,
        created_at=created_at or utc_now(),
    )


def validate_identity(
    identity: OutputIdentity,
    *,
    record_kind: RecordKind,
    enforce_workflow_scope: bool = True,
) -> None:
    if not isinstance(identity.workflow, WorkflowKind):
        raise PublicSchemaValidationError("identity.workflow must be a WorkflowKind.")
    if not str(identity.run_id).strip():
        raise PublicSchemaValidationError(
            "Public record identity requires a non-empty run_id."
        )
    if not str(identity.system_id).strip():
        raise PublicSchemaValidationError(
            "Public record identity requires a non-empty system_id."
        )
    if identity.workflow == WorkflowKind.BIAS:
        if any(
            value is not None
            for value in (
                identity.runner_id,
                identity.runner_version,
                identity.backend_name,
                identity.backend_version_status,
                identity.effective_seed,
            )
        ):
            raise PublicSchemaValidationError(
                "Bias record identity must not declare a runner."
            )
    elif identity.runner_id is None or not str(identity.runner_id).strip():
        raise PublicSchemaValidationError(
            f"{identity.workflow.value} record identity requires runner_id."
        )
    if identity.backend_version_status is not None and not isinstance(
        identity.backend_version_status, BackendVersionStatus
    ):
        raise PublicSchemaValidationError(
            "identity.backend_version_status must be a BackendVersionStatus."
        )
    if identity.backend_version_status == BackendVersionStatus.DETECTED:
        if not identity.runner_version:
            raise PublicSchemaValidationError(
                "Detected backend identity requires runner_version."
            )
    elif (
        identity.backend_version_status
        in {
            BackendVersionStatus.UNAVAILABLE,
            BackendVersionStatus.UNPARSEABLE,
        }
        and identity.runner_version is not None
    ):
        raise PublicSchemaValidationError(
            "Unavailable or unparseable backend version must serialize as null."
        )
    if identity.effective_seed is not None:
        if (
            isinstance(identity.effective_seed, bool)
            or not isinstance(identity.effective_seed, int)
            or not (0 <= identity.effective_seed <= 2**32 - 1)
        ):
            raise PublicSchemaValidationError(
                "identity.effective_seed must be an integer between 0 and 4294967295."
            )
        if identity.repeat_id is None:
            raise PublicSchemaValidationError(
                "identity.effective_seed requires a repeat_id."
            )
    if (
        enforce_workflow_scope
        and identity.workflow
        in {
            WorkflowKind.SCREEN,
            WorkflowKind.ORACLE,
        }
        and not identity.compound_id
        and record_kind != RecordKind.FAILURE
    ):
        raise PublicSchemaValidationError(
            f"{identity.workflow.value.title()} record identity requires compound_id."
        )
    entity_values = (identity.entity_id, identity.entity_type, identity.chain_id)
    if any(value is not None for value in entity_values) and not all(
        value is not None and str(value).strip() for value in entity_values
    ):
        raise PublicSchemaValidationError(
            "Entity-scoped identity requires entity_id, entity_type, and chain_id together."
        )
    if identity.related_entity_id is None and identity.related_chain_id is not None:
        raise PublicSchemaValidationError(
            "related_chain_id requires related_entity_id."
        )
    if identity.related_entity_id is not None and identity.related_chain_id is None:
        raise PublicSchemaValidationError(
            "related_entity_id requires related_chain_id."
        )
    if identity.related_entity_id is not None and identity.entity_id is None:
        raise PublicSchemaValidationError(
            "Pairwise identity requires a primary entity identity."
        )
    if identity.repeat_id is not None and identity.repeat_id < 1:
        raise PublicSchemaValidationError("repeat_id must be one-based when present.")
    if identity.sample_id is not None and identity.sample_id < 0:
        raise PublicSchemaValidationError(
            "sample_id must be non-negative when present."
        )


def _validate_json_value(value: object, *, path: str) -> None:
    if value is None or isinstance(value, (bool, str)):
        return
    if isinstance(value, int):
        return
    if isinstance(value, float):
        if math.isfinite(value):
            return
        raise PublicSchemaValidationError(f"{path} contains a non-finite number.")
    if isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            _validate_json_value(item, path=f"{path}[{index}]")
        return
    if isinstance(value, Mapping):
        for key, item in value.items():
            if not isinstance(key, str):
                raise PublicSchemaValidationError(f"{path} contains a non-string key.")
            _validate_json_value(item, path=f"{path}.{key}")
        return
    raise PublicSchemaValidationError(
        f"{path} contains non-JSON value of type {type(value).__name__}."
    )


def _validate_evidence(sources: Collection[EvidenceSource]) -> None:
    for source in sources:
        if not isinstance(source, EvidenceSource):
            raise PublicSchemaValidationError(
                "Evidence entries must be EvidenceSource instances."
            )
        if not source.kind.strip() or not source.identifier.strip():
            raise PublicSchemaValidationError(
                "Evidence entries require non-empty kind and identifier."
            )
        if (
            source.sha256 is not None
            and re.fullmatch(r"[0-9a-fA-F]{64}", source.sha256) is None
        ):
            raise PublicSchemaValidationError(
                "Evidence sha256 must contain exactly 64 hexadecimal characters."
            )


def _validate_artifacts(artifacts: Collection[ArtifactReference]) -> None:
    for artifact in artifacts:
        if not isinstance(artifact, ArtifactReference):
            raise PublicSchemaValidationError(
                "Artifact entries must be ArtifactReference instances."
            )
        if not artifact.label.strip() or not artifact.relative_path.strip():
            raise PublicSchemaValidationError(
                "Artifact entries require non-empty label and relative_path."
            )
        path = artifact.relative_path.replace("\\", "/")
        if path.startswith("/") or ".." in path.split("/"):
            raise PublicSchemaValidationError(
                f"Artifact path must remain relative to results/: {path!r}."
            )


def validate_public_record(record: PublicRecord) -> None:
    envelope = record.envelope
    if not isinstance(envelope.record_kind, RecordKind):
        raise PublicSchemaValidationError("envelope.record_kind must be a RecordKind.")
    if envelope.schema_version != PUBLIC_SCHEMA_VERSION:
        raise PublicSchemaValidationError(
            f"Unsupported public schema version {envelope.schema_version!r}."
        )
    try:
        datetime.fromisoformat(envelope.created_at)
    except (AttributeError, ValueError) as exc:
        raise PublicSchemaValidationError(
            "Record created_at must be an ISO-8601 timestamp."
        ) from exc
    validate_identity(envelope.identity, record_kind=envelope.record_kind)
    expected_id: str
    if isinstance(record, ExecutionRecord):
        if envelope.record_kind != RecordKind.EXECUTION:
            raise PublicSchemaValidationError(
                "ExecutionRecord requires record_kind='execution'."
            )
        identity = envelope.identity
        if identity.workflow != WorkflowKind.SCREEN:
            raise PublicSchemaValidationError(
                "ExecutionRecord is currently defined only for Screen outputs."
            )
        if not all(
            (
                identity.compound_id,
                identity.execution_directory,
                identity.repeat_id is not None,
                identity.model_id,
            )
        ):
            raise PublicSchemaValidationError(
                "Screen execution identity requires compound_id, execution_directory, "
                "repeat_id, and model_id."
            )
        if record.execution_directory != identity.execution_directory:
            raise PublicSchemaValidationError(
                "Execution record directory must match its identity."
            )
        if not isinstance(record.status, ExecutionStatus):
            raise PublicSchemaValidationError(
                "Execution record status must be an ExecutionStatus."
            )
        if record.status == ExecutionStatus.SUCCESS and record.error is not None:
            raise PublicSchemaValidationError(
                "Successful execution records must not contain an error."
            )
        if record.status != ExecutionStatus.SUCCESS and record.error is None:
            raise PublicSchemaValidationError(
                "Failed or unavailable execution records require an error."
            )
        if record.error is not None:
            if not isinstance(record.error, StructuredExecutionError):
                raise PublicSchemaValidationError(
                    "Execution error must be StructuredExecutionError."
                )
            if not isinstance(record.error.stage, FailureStage):
                raise PublicSchemaValidationError(
                    "Execution error stage must be a FailureStage."
                )
            if not record.error.error_code.strip() or not record.error.message.strip():
                raise PublicSchemaValidationError(
                    "Execution errors require error_code and message."
                )
            _validate_json_value(record.error.details, path="execution.error.details")
        _validate_artifacts(record.artifacts)
        expected_id = make_record_id(RecordKind.EXECUTION, identity, "execution")
    elif isinstance(record, MetricRecord):
        if envelope.record_kind != RecordKind.METRIC:
            raise PublicSchemaValidationError(
                "MetricRecord requires record_kind='metric'."
            )
        expected_id = make_record_id(
            RecordKind.METRIC, envelope.identity, record.metric_name, record.statistic
        )
        validate_metric_record(record)
        _validate_json_value(record.value, path="metric.value")
        _validate_evidence(record.evidence)
        if not isinstance(record.metric_class, MetricClass):
            raise PublicSchemaValidationError("metric_class must be a MetricClass.")
        if not isinstance(record.direction, OptimizationDirection):
            raise PublicSchemaValidationError(
                "direction must be an OptimizationDirection."
            )
        if not isinstance(record.evidence_regime, EvidenceRegime):
            raise PublicSchemaValidationError(
                "evidence_regime must be an EvidenceRegime."
            )
        if record.envelope.identity.repeat_id is not None and (
            record.envelope.identity.model_id is None
            and record.envelope.identity.sample_id is None
        ):
            raise PublicSchemaValidationError(
                "Repeat-scoped metric identity requires model_id or sample_id."
            )
    elif isinstance(record, WorkflowFailureRecord):
        if envelope.record_kind != RecordKind.FAILURE:
            raise PublicSchemaValidationError(
                "WorkflowFailureRecord requires record_kind='failure'."
            )
        if not isinstance(record.stage, FailureStage):
            raise PublicSchemaValidationError("failure.stage must be a FailureStage.")
        expected_id = make_record_id(
            RecordKind.FAILURE, envelope.identity, record.stage.value, record.error_code
        )
        if not record.error_code.strip() or not record.message.strip():
            raise PublicSchemaValidationError(
                "Failure records require error_code and message."
            )
        if (
            not isinstance(record.status, RecordStatus)
            or record.status != RecordStatus.FAILED
        ):
            raise PublicSchemaValidationError(
                "Failure records require status='failed'."
            )
        _validate_json_value(record.details, path="failure.details")
    elif isinstance(record, SuccessRecord):
        if envelope.record_kind != RecordKind.SUCCESS:
            raise PublicSchemaValidationError(
                "SuccessRecord requires record_kind='success'."
            )
        expected_id = make_record_id(RecordKind.SUCCESS, envelope.identity, "success")
        if (
            not isinstance(record.status, RecordStatus)
            or record.status != RecordStatus.SUCCESS
        ):
            raise PublicSchemaValidationError(
                "Success records require status='success'."
            )
        if record.envelope.identity.repeat_id is not None and (
            record.envelope.identity.model_id is None
            and record.envelope.identity.sample_id is None
        ):
            raise PublicSchemaValidationError(
                "Repeat-scoped success identity requires model_id or sample_id."
            )
        _validate_artifacts(record.artifacts)
    else:
        raise PublicSchemaValidationError(
            f"Unsupported public record type {type(record).__name__}."
        )
    if envelope.record_id != expected_id:
        raise PublicSchemaValidationError(
            f"Record {envelope.record_id!r} does not match its deterministic identity."
        )


def validate_public_bundle(bundle: PublicOutputBundle) -> PublicOutputBundle:
    if bundle.manifest.schema_version != PUBLIC_SCHEMA_VERSION:
        raise PublicSchemaValidationError(
            "Manifest schema_version does not match the public contract."
        )
    validate_identity(
        bundle.manifest.identity,
        record_kind=RecordKind.SUCCESS,
        enforce_workflow_scope=False,
    )
    _validate_evidence(bundle.manifest.evidence)
    _validate_artifacts(bundle.manifest.artifacts)
    _validate_execution_provenance(bundle)
    seen: set[str] = set()
    counts = defaultdict(int)
    for record in bundle.records:
        validate_public_record(record)
        record_id = record.envelope.record_id
        if record_id in seen:
            raise AmbiguousIdentityError(f"Duplicate public record_id {record_id!r}.")
        seen.add(record_id)
        counts[record.envelope.record_kind.value] += 1
    declared = {key: int(value) for key, value in bundle.manifest.record_counts.items()}
    actual = {kind.value: counts[kind.value] for kind in RecordKind}
    if declared and declared != actual:
        raise PublicSchemaValidationError(
            f"Manifest record_counts {declared!r} do not match records {actual!r}."
        )
    execution_success = any(
        isinstance(record, ExecutionRecord)
        and record.status == ExecutionStatus.SUCCESS
        for record in bundle.records
    )
    execution_failure = any(
        isinstance(record, ExecutionRecord)
        and record.status != ExecutionStatus.SUCCESS
        for record in bundle.records
    )
    has_success = bool(actual["success"]) or execution_success
    has_failure = bool(actual["failure"]) or execution_failure
    expected_status = (
        "partial"
        if has_success and has_failure
        else "success"
        if has_success
        else "failed"
    )
    if bundle.manifest.status != expected_status:
        raise PublicSchemaValidationError(
            f"Manifest status {bundle.manifest.status!r} does not match record outcomes "
            f"({expected_status!r})."
        )
    return bundle


def _validate_backend_identity(backend: RunnerBackendIdentity) -> None:
    if not backend.runner_name.strip() or not backend.backend_name.strip():
        raise RunnerProvenanceError(
            "Backend provenance requires non-empty runner and backend names."
        )
    if not isinstance(backend.version_status, BackendVersionStatus):
        raise RunnerProvenanceError(
            "Backend provenance version_status must be a BackendVersionStatus."
        )
    if backend.version_status == BackendVersionStatus.DETECTED:
        if not backend.version:
            raise RunnerProvenanceError(
                "Detected backend provenance requires a normalized version."
            )
        try:
            normalized_version = str(Version(backend.version))
        except InvalidVersion as exc:
            raise RunnerProvenanceError(
                "Detected backend provenance requires a valid PEP 440 version."
            ) from exc
        if normalized_version != backend.version:
            raise RunnerProvenanceError(
                "Detected backend provenance version must be normalized."
            )
    elif backend.version is not None:
        raise RunnerProvenanceError(
            "Unavailable or unparseable backend provenance must have a null version."
        )
    if (
        backend.version_status == BackendVersionStatus.UNAVAILABLE
        and backend.raw_version is not None
    ):
        raise RunnerProvenanceError(
            "Unavailable backend provenance must not contain a raw version."
        )
    if (
        backend.version_status == BackendVersionStatus.UNPARSEABLE
        and backend.raw_version is None
    ):
        raise RunnerProvenanceError(
            "Unparseable backend provenance must retain the raw version."
        )


def _validate_seed_plan(seed_plan: SeedPlan) -> dict[int, int]:
    def valid_seed(value: object) -> bool:
        return (
            not isinstance(value, bool)
            and isinstance(value, int)
            and 0 <= value <= 2**32 - 1
        )

    if not isinstance(seed_plan.origin, SeedOrigin):
        raise RunnerProvenanceError("Seed plan origin must be a SeedOrigin.")
    if not valid_seed(seed_plan.resolved_base_seed):
        raise RunnerProvenanceError(
            "Resolved base seed must be an integer between 0 and 4294967295."
        )
    if seed_plan.requested_base_seed is not None and not valid_seed(
        seed_plan.requested_base_seed
    ):
        raise RunnerProvenanceError(
            "Requested base seed must be null or an integer between 0 and 4294967295."
        )
    if seed_plan.requested_base_seed is not None and (
        seed_plan.requested_base_seed != seed_plan.resolved_base_seed
    ):
        raise RunnerProvenanceError(
            "A requested base seed must equal the resolved base seed."
        )
    if (
        seed_plan.origin == SeedOrigin.USER_SPECIFIED
        and seed_plan.requested_base_seed is None
    ):
        raise RunnerProvenanceError(
            "User-specified seed provenance requires a requested base seed."
        )
    if (
        seed_plan.origin == SeedOrigin.GENERATED
        and seed_plan.requested_base_seed is not None
    ):
        raise RunnerProvenanceError(
            "Generated seed provenance must have a null requested base seed."
        )
    if not seed_plan.repeats:
        raise RunnerProvenanceError(
            "Runner-backed seed plans require at least one repeat."
        )
    seeds: dict[int, int] = {}
    derived_seeds: set[int] = set()
    effective_seeds: set[int] = set()
    for item in seed_plan.repeats:
        if isinstance(item.repeat_id, bool) or not isinstance(item.repeat_id, int):
            raise RunnerProvenanceError("Seed-plan repeat IDs must be integers.")
        if not isinstance(item.adjustment, SeedAdjustment):
            raise RunnerProvenanceError(
                "Repeat seed adjustment must be a SeedAdjustment."
            )
        for label, value in (
            ("resolved_base_seed", item.resolved_base_seed),
            ("derived_seed", item.derived_seed),
            ("effective_seed", item.effective_seed),
        ):
            if not valid_seed(value):
                raise RunnerProvenanceError(
                    f"Seed provenance {label} must be between 0 and 4294967295."
                )
        if item.repeat_id in seeds:
            raise RunnerProvenanceError(
                f"Seed plan contains duplicate repeat {item.repeat_id}."
            )
        if item.resolved_base_seed != seed_plan.resolved_base_seed:
            raise RunnerProvenanceError(
                "Repeat seed provenance does not match the seed-plan base seed."
            )
        if (
            item.requested_base_seed != seed_plan.requested_base_seed
            or item.origin != seed_plan.origin
        ):
            raise RunnerProvenanceError(
                "Repeat seed provenance does not match the seed-plan origin."
            )
        if item.adjustment == SeedAdjustment.UNCHANGED and (
            item.effective_seed != item.derived_seed
            or item.adjustment_reason is not None
        ):
            raise RunnerProvenanceError(
                "Unchanged seed provenance must retain the derived seed and omit "
                "a reason."
            )
        if item.adjustment == SeedAdjustment.BACKEND_ADJUSTED and (
            item.effective_seed == item.derived_seed or not item.adjustment_reason
        ):
            raise RunnerProvenanceError(
                "Backend-adjusted seed provenance requires a changed seed and reason."
            )
        seeds[item.repeat_id] = item.effective_seed
        if (
            item.derived_seed in derived_seeds
            or item.effective_seed in effective_seeds
        ):
            raise RunnerProvenanceError(
                "Repeat seed provenance must contain distinct derived and effective "
                "seeds."
            )
        derived_seeds.add(item.derived_seed)
        effective_seeds.add(item.effective_seed)
    if sorted(seeds) != list(range(1, len(seeds) + 1)):
        raise RunnerProvenanceError(
            "Seed plan repeat IDs must be contiguous and start at one."
        )
    return seeds


def _validate_execution_provenance(bundle: PublicOutputBundle) -> None:
    backend = bundle.manifest.backend
    seed_plan = bundle.manifest.seed_plan
    declares_execution = any(
        value is not None
        for value in (
            bundle.manifest.identity.runner_version,
            bundle.manifest.identity.backend_name,
            bundle.manifest.identity.backend_version_status,
        )
    ) or any(
        record.envelope.identity.effective_seed is not None
        for record in bundle.records
    )
    if declares_execution and backend is None and seed_plan is None:
        raise RunnerProvenanceError(
            "Runner-backed execution provenance requires backend and seed_plan."
        )
    if (backend is None) != (seed_plan is None):
        raise RunnerProvenanceError(
            "Public manifest backend and seed_plan must be declared together."
        )
    if backend is None or seed_plan is None:
        return
    if not isinstance(backend, RunnerBackendIdentity) or not isinstance(
        seed_plan, SeedPlan
    ):
        raise RunnerProvenanceError(
            "Public manifest execution provenance must use typed contract values."
        )
    _validate_backend_identity(backend)
    seeds = _validate_seed_plan(seed_plan)
    manifest_identity = bundle.manifest.identity
    if manifest_identity.runner_id != backend.runner_name:
        raise RunnerProvenanceError(
            "Public manifest runner identity does not match backend provenance."
        )
    if (
        manifest_identity.runner_version != backend.version
        or manifest_identity.backend_name != backend.backend_name
        or manifest_identity.backend_version_status != backend.version_status
    ):
        raise RunnerProvenanceError(
            "Public manifest identity does not match backend provenance."
        )
    for record in bundle.records:
        identity = record.envelope.identity
        if (
            identity.runner_id != backend.runner_name
            or identity.runner_version != backend.version
            or identity.backend_name != backend.backend_name
            or identity.backend_version_status != backend.version_status
        ):
            raise RunnerProvenanceError(
                "Public record identity does not match backend provenance."
            )
        if identity.repeat_id is not None and identity.effective_seed != seeds.get(
            identity.repeat_id
        ):
            raise RunnerProvenanceError(
                f"Public record seed does not match repeat {identity.repeat_id}."
            )


def aggregate_metric_records(
    records: Collection[MetricRecord],
    *,
    statistic: str,
) -> tuple[MetricRecord, ...]:
    if statistic not in {"mean", "std"}:
        raise ValueError("statistic must be 'mean' or 'std'.")
    groups: dict[tuple[object, ...], list[MetricRecord]] = defaultdict(list)
    for record in records:
        identity = record.envelope.identity
        key = (
            identity.workflow,
            identity.run_id,
            identity.system_id,
            identity.compound_id,
            identity.runner_id,
            identity.runner_version,
            identity.entity_id,
            identity.entity_type,
            identity.chain_id,
            identity.related_entity_id,
            identity.related_chain_id,
            record.metric_name,
        )
        groups[key].append(record)
    output: list[MetricRecord] = []
    for grouped in groups.values():
        first = grouped[0]
        definition = get_metric_definition(first.metric_name)
        observations: set[tuple[object, ...]] = set()
        for item in grouped:
            if item.statistic != "value":
                raise AmbiguousIdentityError(
                    "Only value records can be aggregated into robustness statistics."
                )
            item_identity = item.envelope.identity
            observation = (
                item_identity.repeat_id,
                item_identity.model_id,
                item_identity.sample_id,
            )
            if observation in observations:
                raise AmbiguousIdentityError(
                    f"Metric {first.metric_name!r} has duplicate observations for "
                    f"repeat/model/sample identity {observation!r}."
                )
            observations.add(observation)
        if any(
            (
                item.metric_group,
                item.metric_class,
                item.unit,
                item.direction,
                item.evidence_regime,
                item.evidence,
            )
            != (
                first.metric_group,
                first.metric_class,
                first.unit,
                first.direction,
                first.evidence_regime,
                first.evidence,
            )
            for item in grouped[1:]
        ):
            raise AmbiguousIdentityError(
                f"Metric {first.metric_name!r} has conflicting metadata within one aggregate identity."
            )
        values = [
            float(item.value)
            for item in grouped
            if item.status == RecordStatus.COMPUTED
            and isinstance(item.value, (int, float))
            and not isinstance(item.value, bool)
            and math.isfinite(float(item.value))
        ]
        value = (
            mean(values)
            if statistic == "mean" and values
            else stdev(values)
            if statistic == "std" and len(values) > 1
            else None
        )
        status = RecordStatus.COMPUTED if value is not None else RecordStatus.MISSING
        identity = replace(
            first.envelope.identity, repeat_id=None, model_id=None, sample_id=None
        )
        output.append(
            MetricRecord(
                envelope=make_envelope(
                    RecordKind.METRIC, identity, first.metric_name, statistic
                ),
                metric_name=first.metric_name,
                metric_group="robustness_metrics",
                metric_class=MetricClass.ROBUSTNESS,
                status=status,
                value=value,
                unit=definition.unit,
                direction=definition.direction,
                evidence_regime=EvidenceRegime.REFERENCE_FREE,
                evidence=first.evidence,
                statistic=statistic,  # type: ignore[arg-type]
                reason_code=None if value is not None else "insufficient_observations",
                message=None
                if value is not None
                else "No finite observations were available for aggregation.",
            )
        )
    return tuple(output)


def failure_from_exception(
    exc: BaseException,
    *,
    identity: OutputIdentity,
    stage: FailureStage,
    error_code: str | None = None,
    retryable: bool = False,
    details: Mapping[str, object] | None = None,
) -> WorkflowFailureRecord:
    # Input-validation exceptions own their stable public code. Callers may still
    # add workflow-local details (for example a screening row number), but cannot
    # accidentally replace the validation contract with a generic code.
    from cofolder.modules.input.config import InputValidationError

    exception_details = getattr(exc, "details", {})
    combined_details = {**exception_details, **(details or {})}
    safe_details = {str(key): value for key, value in combined_details.items()}
    resolved_error_code = (
        exc.error_code
        if isinstance(exc, InputValidationError)
        else error_code or getattr(exc, "error_code", "workflow_execution_failed")
    )
    return WorkflowFailureRecord(
        envelope=make_envelope(
            RecordKind.FAILURE, identity, stage.value, resolved_error_code
        ),
        stage=stage,
        exception_type=type(exc).__name__,
        error_code=resolved_error_code,
        message=str(exc) or type(exc).__name__,
        retryable=retryable,
        details=safe_details,  # type: ignore[arg-type]
    )


__all__ = [
    "WorkflowExecutionError",
    "aggregate_metric_records",
    "failure_from_exception",
    "make_envelope",
    "make_record_id",
    "utc_now",
    "validate_identity",
    "validate_public_bundle",
    "validate_public_record",
]
