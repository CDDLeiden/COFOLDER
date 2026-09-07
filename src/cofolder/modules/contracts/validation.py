from __future__ import annotations

import math
import re
from collections import defaultdict
from collections.abc import Collection, Mapping
from dataclasses import replace
from datetime import UTC, datetime
from statistics import mean, stdev
from uuid import NAMESPACE_URL, uuid5

from .metrics import get_metric_definition, validate_metric_record
from .models import (
    PUBLIC_SCHEMA_VERSION,
    AmbiguousIdentityError,
    ArtifactReference,
    EvidenceRegime,
    EvidenceSource,
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
        if identity.runner_id is not None or identity.runner_version is not None:
            raise PublicSchemaValidationError(
                "Bias record identity must not declare a runner."
            )
    elif identity.runner_id is None or not str(identity.runner_id).strip():
        raise PublicSchemaValidationError(
            f"{identity.workflow.value} record identity requires runner_id."
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
    if isinstance(record, MetricRecord):
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
    expected_status = (
        "partial"
        if actual["success"] and actual["failure"]
        else "success"
        if actual["success"]
        else "failed"
    )
    if bundle.manifest.status != expected_status:
        raise PublicSchemaValidationError(
            f"Manifest status {bundle.manifest.status!r} does not match record outcomes "
            f"({expected_status!r})."
        )
    return bundle


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
    error_code: str,
    retryable: bool = False,
    details: Mapping[str, object] | None = None,
) -> WorkflowFailureRecord:
    safe_details = {str(key): value for key, value in (details or {}).items()}
    return WorkflowFailureRecord(
        envelope=make_envelope(RecordKind.FAILURE, identity, stage.value, error_code),
        stage=stage,
        exception_type=type(exc).__name__,
        error_code=error_code,
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
