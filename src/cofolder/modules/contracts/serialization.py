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
    MetricRecord,
    PublicOutputBundle,
    PublicSerializationError,
    RecordKind,
    SerializedPublicOutput,
    SuccessRecord,
    WorkflowFailureRecord,
)
from .validation import validate_public_bundle


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
    record: SuccessRecord | MetricRecord | WorkflowFailureRecord,
) -> dict[str, Any]:
    raw = _jsonable(asdict(record))
    envelope = raw.pop("envelope")
    identity = envelope.pop("identity")
    return {**envelope, **identity, **raw}


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
            "runner_id",
            "runner_version",
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
        successes_path=output_dir / "successes.csv",
        metrics_path=output_dir / "metrics.csv",
        failures_path=output_dir / "failures.csv",
        manifest_path=output_dir / "manifest.json",
    )


__all__ = ["record_to_dict", "write_public_bundle"]
