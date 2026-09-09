from __future__ import annotations

import csv
import json
import os
import tempfile
from collections.abc import Iterable
from dataclasses import asdict, replace
from enum import Enum
from pathlib import Path
from typing import Any

from .models import (
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
    PublicSerializationError,
    RecordKind,
    RecordStatus,
    SerializedPublicOutput,
    StructuredExecutionError,
    SuccessRecord,
    WorkflowFailureRecord,
    WorkflowKind,
)
from .validation import validate_public_bundle, validate_public_record


def _jsonable(value: Any) -> Any:
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def record_to_dict(
    record: ExecutionRecord | SuccessRecord | MetricRecord | WorkflowFailureRecord,
) -> dict[str, Any]:
    raw = _jsonable(asdict(record))
    envelope = raw.pop("envelope")
    identity = envelope.pop("identity")
    return {**envelope, **identity, **raw}


def record_from_dict(raw: dict[str, Any]) -> PublicRecord:
    """Deserialize and validate one flattened public JSONL record."""
    identity_fields = OutputIdentity.__dataclass_fields__
    identity_values = {name: raw.get(name) for name in identity_fields}
    identity_values["workflow"] = WorkflowKind(identity_values["workflow"])
    if identity_values.get("backend_version_status") is not None:
        identity_values["backend_version_status"] = BackendVersionStatus(
            identity_values["backend_version_status"]
        )
    identity = OutputIdentity(**identity_values)
    kind = RecordKind(raw["record_kind"])
    envelope = PublicRecordEnvelope(
        schema_version=str(raw["schema_version"]),
        record_id=str(raw["record_id"]),
        record_kind=kind,
        identity=identity,
        created_at=str(raw["created_at"]),
    )
    artifacts = tuple(
        ArtifactReference(**item) for item in raw.get("artifacts") or ()
    )
    if kind is RecordKind.EXECUTION:
        error_value = raw.get("error")
        error = None
        if error_value is not None:
            error = StructuredExecutionError(
                stage=FailureStage(error_value["stage"]),
                exception_type=str(error_value["exception_type"]),
                error_code=str(error_value["error_code"]),
                message=str(error_value["message"]),
                retryable=bool(error_value.get("retryable", False)),
                details=error_value.get("details") or {},
            )
        record: PublicRecord = ExecutionRecord(
            envelope=envelope,
            status=ExecutionStatus(raw["status"]),
            execution_directory=str(raw["execution_directory"]),
            error=error,
            artifacts=artifacts,
        )
    elif kind is RecordKind.SUCCESS:
        record = SuccessRecord(envelope=envelope, artifacts=artifacts)
    elif kind is RecordKind.METRIC:
        record = MetricRecord(
            envelope=envelope,
            metric_name=str(raw["metric_name"]),
            metric_group=str(raw["metric_group"]),
            metric_class=MetricClass(raw["metric_class"]),
            status=RecordStatus(raw["status"]),
            value=raw.get("value"),
            unit=raw.get("unit"),
            direction=OptimizationDirection(raw["direction"]),
            evidence_regime=EvidenceRegime(raw["evidence_regime"]),
            evidence=tuple(
                EvidenceSource(**item) for item in raw.get("evidence") or ()
            ),
            statistic=raw.get("statistic", "value"),
            reason_code=raw.get("reason_code"),
            message=raw.get("message"),
        )
    elif kind is RecordKind.FAILURE:
        record = WorkflowFailureRecord(
            envelope=envelope,
            stage=FailureStage(raw["stage"]),
            exception_type=str(raw["exception_type"]),
            error_code=str(raw["error_code"]),
            message=str(raw["message"]),
            retryable=bool(raw.get("retryable", False)),
            details=raw.get("details") or {},
        )
    else:  # pragma: no cover - exhaustive over RecordKind
        raise ValueError(f"Unsupported public record kind: {kind!r}.")
    validate_public_record(record)
    return record


def read_public_records(path: str | Path) -> tuple[PublicRecord, ...]:
    """Read typed public records from a records.jsonl file."""
    records: list[PublicRecord] = []
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                records.append(record_from_dict(json.loads(line)))
    return tuple(records)


def _write_json(path: Path, value: object) -> None:
    path.write_text(
        json.dumps(_jsonable(value), indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def _write_jsonl(path: Path, records: Iterable[object]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(
                json.dumps(record_to_dict(record), sort_keys=True, allow_nan=False)
                + "\n"
            )


def _write_csv(
    path: Path, records: list[object], *, default_fields: tuple[str, ...]
) -> None:
    rows = [record_to_dict(record) for record in records]
    fields = list(default_fields)
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: json.dumps(value, allow_nan=False)
                    if isinstance(value, (dict, list))
                    else value
                    for key, value in row.items()
                }
            )


def write_public_bundle(
    bundle: PublicOutputBundle, output_dir: Path
) -> SerializedPublicOutput:
    output_dir = Path(output_dir)
    counts = {kind.value: 0 for kind in RecordKind}
    for record in bundle.records:
        counts[record.envelope.record_kind.value] += 1
    bundle = replace(bundle, manifest=replace(bundle.manifest, record_counts=counts))
    validate_public_bundle(bundle)
    output_dir.mkdir(parents=True, exist_ok=True)
    names = (
        "records.jsonl",
        "executions.csv",
        "successes.csv",
        "metrics.csv",
        "failures.csv",
        "manifest.json",
    )
    temp_paths: dict[str, Path] = {}
    try:
        for name in names:
            fd, raw_path = tempfile.mkstemp(prefix=f".{name}.", dir=output_dir)
            os.close(fd)
            temp_paths[name] = Path(raw_path)
        successes = [
            record for record in bundle.records if isinstance(record, SuccessRecord)
        ]
        executions = [
            record for record in bundle.records if isinstance(record, ExecutionRecord)
        ]
        metrics = [
            record for record in bundle.records if isinstance(record, MetricRecord)
        ]
        failures = [
            record
            for record in bundle.records
            if isinstance(record, WorkflowFailureRecord)
        ]
        _write_jsonl(temp_paths["records.jsonl"], bundle.records)
        common = (
            "schema_version",
            "record_id",
            "record_kind",
            "workflow",
            "run_id",
            "system_id",
            "compound_id",
            "execution_directory",
            "runner_id",
            "runner_version",
            "backend_name",
            "backend_version_status",
            "effective_seed",
            "repeat_id",
            "model_id",
            "sample_id",
            "entity_id",
            "entity_type",
            "chain_id",
            "related_entity_id",
            "related_chain_id",
            "created_at",
            "status",
        )
        _write_csv(
            temp_paths["executions.csv"],
            executions,
            default_fields=common + ("error", "artifacts"),
        )
        _write_csv(temp_paths["successes.csv"], successes, default_fields=common)
        _write_csv(
            temp_paths["metrics.csv"],
            metrics,
            default_fields=common
            + (
                "metric_name",
                "metric_group",
                "metric_class",
                "value",
                "unit",
                "direction",
                "evidence_regime",
                "evidence",
                "statistic",
                "reason_code",
                "message",
            ),
        )
        _write_csv(
            temp_paths["failures.csv"],
            failures,
            default_fields=common
            + (
                "stage",
                "exception_type",
                "error_code",
                "message",
                "retryable",
                "details",
            ),
        )
        _write_json(temp_paths["manifest.json"], asdict(bundle.manifest))
        for name in names:
            os.replace(temp_paths[name], output_dir / name)
        for legacy_name in (
            "system_metrics.csv",
            "chain_metrics.csv",
            "robustness_metrics.csv",
            "oracle_result.csv",
            "screen_results.csv",
            "screen_results_with_scores.csv",
        ):
            (output_dir / legacy_name).unlink(missing_ok=True)
    except Exception as exc:
        for path in temp_paths.values():
            path.unlink(missing_ok=True)
        if isinstance(exc, PublicSerializationError):
            raise
        raise PublicSerializationError(
            f"Could not write public output bundle to {output_dir}: {exc}"
        ) from exc
    return SerializedPublicOutput(
        output_dir=output_dir,
        records_path=output_dir / "records.jsonl",
        executions_path=output_dir / "executions.csv",
        successes_path=output_dir / "successes.csv",
        metrics_path=output_dir / "metrics.csv",
        failures_path=output_dir / "failures.csv",
        manifest_path=output_dir / "manifest.json",
    )


__all__ = [
    "read_public_records",
    "record_from_dict",
    "record_to_dict",
    "write_public_bundle",
]
