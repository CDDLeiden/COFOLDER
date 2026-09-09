from __future__ import annotations

import json
import math
import re
from collections.abc import Collection, Iterable, Mapping
from dataclasses import replace
from pathlib import Path
from typing import Any

import pandas as pd

from .metrics import get_metric_definition, resolve_evidence_regime
from .models import (
    ArtifactReference,
    BackendVersionStatus,
    EvidenceRegime,
    EvidenceSource,
    MetricClass,
    MetricRecord,
    OutputIdentity,
    PublicManifest,
    PublicOutputBundle,
    RecordKind,
    RecordStatus,
    RunnerBackendIdentity,
    SeedPlan,
    SuccessRecord,
    WorkflowFailureRecord,
)
from .validation import make_envelope

EXECUTION_METADATA_COLUMNS = frozenset(
    {
        "runner_id",
        "backend_name",
        "runner_version",
        "backend_version_status",
        "effective_seed",
    }
)
SYSTEM_METADATA_COLUMNS = frozenset(
    {"idx", "cif_file", "model_name", "repeat", "diffusion_sample"}
    | EXECUTION_METADATA_COLUMNS
)
CHAIN_METADATA_COLUMNS = frozenset(
    SYSTEM_METADATA_COLUMNS
    | {"conf_chain_id", "CHAIN_ID", "ENTITY_TYPE", "ENTITY_ID", "ligand_molecule_id"}
)


def _present(value: object) -> bool:
    if value is None:
        return False
    try:
        result = pd.isna(value)
    except (TypeError, ValueError):
        return True
    return bool(not result) if isinstance(result, bool) else True


def _json_value(value: Any, *, value_type: str) -> Any:
    if not _present(value):
        return None
    if hasattr(value, "item"):
        value = value.item()
    if value_type == "boolean":
        if isinstance(value, str):
            return value.strip().lower() in {"1", "true", "yes"}
        return bool(value)
    if value_type == "integer":
        return int(value)
    if value_type == "float":
        number = float(value)
        return number if math.isfinite(number) else None
    if value_type == "json" and isinstance(value, str):
        try:
            return json.loads(value)
        except json.JSONDecodeError:
            return value
    return value


def _entity_ids(chain_df: pd.DataFrame) -> dict[str, str]:
    output: dict[str, str] = {}
    if chain_df.empty or "CHAIN_ID" not in chain_df.columns:
        return output
    for _, row in chain_df.iterrows():
        if not _present(row.get("CHAIN_ID")):
            continue
        chain_id = str(row["CHAIN_ID"])
        if chain_id in output:
            continue
        output[chain_id] = (
            str(row["ENTITY_ID"])
            if _present(row.get("ENTITY_ID"))
            else f"entity:{len(set(output.values()))}"
        )
    return output


def _metric_name_and_related_chain(
    column: str,
    conf_to_chain: dict[str, str],
) -> tuple[str, str | None]:
    match = re.fullmatch(r"pair_chains_iptm_(.+)", column)
    if match:
        token = match.group(1)
        return "pair_chains_iptm", conf_to_chain.get(token, token)
    if column.startswith("chain_pair_iptm_"):
        token = column.rsplit("_", 1)[-1]
        return "chain_pair_iptm", conf_to_chain.get(token, token)
    if column.startswith("bespoke_iptm_"):
        token = column.rsplit("_", 1)[-1]
        return "bespoke_iptm", conf_to_chain.get(token, token)
    return column, None


def _metric_status(
    *,
    value: object,
    metric_name: str,
    requested_metrics: Collection[str] | None,
    unavailable_groups: Collection[str],
    evidence: tuple[EvidenceSource, ...],
) -> tuple[RecordStatus, str | None, str | None]:
    definition = get_metric_definition(metric_name)
    if (
        requested_metrics is not None
        and metric_name not in requested_metrics
        and definition.group not in requested_metrics
    ):
        return RecordStatus.NOT_REQUESTED, None, None
    if definition.group in unavailable_groups:
        return (
            RecordStatus.UNSUPPORTED,
            "runner_capability_unavailable",
            "The selected runner does not support this metric group.",
        )
    evidence_kinds = {source.kind for source in evidence}
    if (
        definition.allowed_evidence_regimes
        == frozenset({EvidenceRegime.REFERENCE_STRUCTURE})
        and "reference_structure" not in evidence_kinds
    ):
        return (
            RecordStatus.UNSUPPORTED,
            "reference_structure_unavailable",
            "A reference structure is required for this metric.",
        )
    if (
        definition.allowed_evidence_regimes == frozenset({EvidenceRegime.CUSTOM_POCKET})
        and "custom_pocket" not in evidence_kinds
    ):
        return (
            RecordStatus.UNSUPPORTED,
            "custom_pocket_unavailable",
            "A custom pocket reference is required for this metric.",
        )
    if (
        EvidenceRegime.REFERENCE_FREE not in definition.allowed_evidence_regimes
        and not evidence_kinds.intersection(
            regime.value for regime in definition.allowed_evidence_regimes
        )
    ):
        return (
            RecordStatus.UNSUPPORTED,
            "required_evidence_unavailable",
            "The metric requires reference evidence that was not supplied.",
        )
    if not _present(value):
        return (
            RecordStatus.MISSING,
            "metric_value_missing",
            "The metric was applicable but no value was produced.",
        )
    return RecordStatus.COMPUTED, None, None


def metric_records_from_frames(
    system_df: pd.DataFrame,
    chain_df: pd.DataFrame,
    *,
    base_identity: OutputIdentity,
    evidence: Iterable[EvidenceSource] = (),
    requested_metrics: Collection[str] | None = None,
    unavailable_groups: Collection[str] = (),
    evidence_regime_overrides: Mapping[str, EvidenceRegime] | None = None,
) -> tuple[SuccessRecord | MetricRecord, ...]:
    evidence_tuple = tuple(evidence)
    entity_ids = _entity_ids(chain_df)
    conf_to_chain: dict[str, str] = {}
    if not chain_df.empty and {"conf_chain_id", "CHAIN_ID"}.issubset(chain_df.columns):
        for _, row in chain_df.iterrows():
            if _present(row.get("conf_chain_id")) and _present(row.get("CHAIN_ID")):
                conf_to_chain.setdefault(
                    str(int(row["conf_chain_id"])), str(row["CHAIN_ID"])
                )

    records: list[SuccessRecord | MetricRecord] = []
    for _, row in system_df.iterrows():
        identity = replace(
            base_identity,
            runner_id=str(row["runner_id"])
            if _present(row.get("runner_id"))
            else base_identity.runner_id,
            runner_version=str(row["runner_version"])
            if _present(row.get("runner_version"))
            else base_identity.runner_version,
            backend_name=str(row["backend_name"])
            if _present(row.get("backend_name"))
            else base_identity.backend_name,
            backend_version_status=BackendVersionStatus(
                str(row["backend_version_status"])
            )
            if _present(row.get("backend_version_status"))
            else base_identity.backend_version_status,
            effective_seed=int(row["effective_seed"])
            if _present(row.get("repeat")) and _present(row.get("effective_seed"))
            else base_identity.effective_seed,
            repeat_id=int(row["repeat"]) if _present(row.get("repeat")) else None,
            model_id=str(row["model_name"])
            if _present(row.get("model_name"))
            else None,
            sample_id=int(row["diffusion_sample"])
            if _present(row.get("diffusion_sample"))
            else None,
        )
        artifacts: tuple[ArtifactReference, ...] = ()
        if _present(row.get("cif_file")):
            name = str(row["cif_file"])
            artifacts = (
                ArtifactReference(
                    "predicted_structure", f"structures/{name}", "structure"
                ),
            )
        records.append(
            SuccessRecord(
                envelope=make_envelope(RecordKind.SUCCESS, identity, "success"),
                artifacts=artifacts,
            )
        )
        records.extend(
            _row_metric_records(
                row,
                SYSTEM_METADATA_COLUMNS,
                identity,
                evidence_tuple,
                requested_metrics,
                unavailable_groups,
                conf_to_chain,
                entity_ids,
                evidence_regime_overrides or {},
            )
        )

    for _, row in chain_df.iterrows():
        chain_id = str(row["CHAIN_ID"]) if _present(row.get("CHAIN_ID")) else None
        entity_type = (
            str(row["ENTITY_TYPE"]) if _present(row.get("ENTITY_TYPE")) else None
        )
        identity = replace(
            base_identity,
            runner_id=str(row["runner_id"])
            if _present(row.get("runner_id"))
            else base_identity.runner_id,
            runner_version=str(row["runner_version"])
            if _present(row.get("runner_version"))
            else base_identity.runner_version,
            backend_name=str(row["backend_name"])
            if _present(row.get("backend_name"))
            else base_identity.backend_name,
            backend_version_status=BackendVersionStatus(
                str(row["backend_version_status"])
            )
            if _present(row.get("backend_version_status"))
            else base_identity.backend_version_status,
            effective_seed=int(row["effective_seed"])
            if _present(row.get("repeat")) and _present(row.get("effective_seed"))
            else base_identity.effective_seed,
            repeat_id=int(row["repeat"]) if _present(row.get("repeat")) else None,
            model_id=str(row["model_name"])
            if _present(row.get("model_name"))
            else None,
            sample_id=int(row["diffusion_sample"])
            if _present(row.get("diffusion_sample"))
            else None,
            entity_id=(
                str(row["ENTITY_ID"])
                if _present(row.get("ENTITY_ID"))
                else entity_ids.get(chain_id)
                if chain_id
                else None
            ),
            entity_type=entity_type,
            chain_id=chain_id,
        )
        records.extend(
            _row_metric_records(
                row,
                CHAIN_METADATA_COLUMNS,
                identity,
                evidence_tuple,
                requested_metrics,
                unavailable_groups,
                conf_to_chain,
                entity_ids,
                evidence_regime_overrides or {},
            )
        )
    return tuple(records)


def _row_metric_records(
    row: pd.Series,
    metadata_columns: Collection[str],
    identity: OutputIdentity,
    evidence: tuple[EvidenceSource, ...],
    requested_metrics: Collection[str] | None,
    unavailable_groups: Collection[str],
    conf_to_chain: dict[str, str],
    entity_ids: dict[str, str],
    evidence_regime_overrides: Mapping[str, EvidenceRegime],
) -> list[MetricRecord]:
    records: list[MetricRecord] = []
    for column, raw_value in row.items():
        if column in metadata_columns:
            continue
        metric_name, related_chain_id = _metric_name_and_related_chain(
            str(column), conf_to_chain
        )
        definition = get_metric_definition(metric_name)
        record_identity = identity
        if related_chain_id is not None:
            record_identity = replace(
                identity,
                related_chain_id=related_chain_id,
                related_entity_id=entity_ids.get(
                    related_chain_id, f"entity:related:{related_chain_id}"
                ),
            )
        status, reason, message = _metric_status(
            value=raw_value,
            metric_name=metric_name,
            requested_metrics=requested_metrics,
            unavailable_groups=unavailable_groups,
            evidence=evidence,
        )
        value = (
            _json_value(raw_value, value_type=definition.value_type)
            if status == RecordStatus.COMPUTED
            else None
        )
        if value is None and status == RecordStatus.COMPUTED:
            status, reason, message = (
                RecordStatus.MISSING,
                "metric_value_non_finite",
                "The metric value was not finite.",
            )
        records.append(
            MetricRecord(
                envelope=make_envelope(
                    RecordKind.METRIC, record_identity, metric_name, "value"
                ),
                metric_name=metric_name,
                metric_group=definition.group,
                metric_class=definition.metric_class,
                status=status,
                value=value,
                unit=definition.unit,
                direction=definition.direction,
                evidence_regime=evidence_regime_overrides.get(
                    metric_name, resolve_evidence_regime(metric_name, evidence)
                ),
                evidence=evidence,
                reason_code=reason,
                message=message,
            )
        )
    return records


def bundle_from_frames(
    system_df: pd.DataFrame,
    chain_df: pd.DataFrame,
    *,
    identity: OutputIdentity,
    evidence: Iterable[EvidenceSource] = (),
    requested_metrics: Collection[str] = (),
    unavailable_groups: Collection[str] = (),
    failures: Iterable[WorkflowFailureRecord] = (),
    artifacts: Iterable[ArtifactReference] = (),
    robustness_df: pd.DataFrame | None = None,
    backend: RunnerBackendIdentity | None = None,
    seed_plan: SeedPlan | None = None,
) -> PublicOutputBundle:
    evidence_tuple = tuple(evidence)
    records = list(
        metric_records_from_frames(
            system_df,
            chain_df,
            base_identity=identity,
            evidence=evidence_tuple,
            requested_metrics=requested_metrics or None,
            unavailable_groups=unavailable_groups,
        )
    )
    records.extend(failures)
    if robustness_df is not None and not robustness_df.empty:
        records.extend(
            robustness_records_from_frame(
                robustness_df,
                base_identity=identity,
            )
        )
    has_success = any(isinstance(record, SuccessRecord) for record in records)
    has_failure = any(isinstance(record, WorkflowFailureRecord) for record in records)
    status = (
        "partial"
        if has_success and has_failure
        else "success"
        if has_success
        else "failed"
    )
    return PublicOutputBundle(
        manifest=PublicManifest(
            schema_version="1.0.0",
            identity=identity,
            status=status,
            evidence=evidence_tuple,
            requested_metrics=tuple(sorted(requested_metrics)),
            artifacts=tuple(artifacts),
            backend=backend,
            seed_plan=seed_plan,
        ),
        records=tuple(records),
    )


def robustness_records_from_frame(
    robustness_df: pd.DataFrame,
    *,
    base_identity: OutputIdentity,
) -> tuple[MetricRecord, ...]:
    records: list[MetricRecord] = []
    for _, row in robustness_df.iterrows():
        entity_type = str(row.get("ENTITY_TYPE", "system"))
        entity_id = (
            str(row.get("ENTITY_ID")) if _present(row.get("ENTITY_ID")) else None
        )
        identity = base_identity
        if entity_type != "system" and entity_id is not None:
            identity = replace(
                base_identity,
                entity_id=entity_id
                if entity_id.startswith("entity:")
                else f"entity:{entity_id}",
                entity_type=entity_type,
                chain_id=entity_id,
            )
        for column, raw_value in row.items():
            if column in {"ENTITY_TYPE", "ENTITY_ID"}:
                continue
            match = re.fullmatch(r"(.+)_(mean|std)", str(column))
            if match is None:
                continue
            raw_metric_name, statistic = match.groups()
            metric_name, related_chain_id = _metric_name_and_related_chain(
                raw_metric_name, {}
            )
            record_identity = identity
            if related_chain_id is not None:
                record_identity = replace(
                    identity,
                    related_entity_id=f"entity:related:{related_chain_id}",
                    related_chain_id=related_chain_id,
                )
            definition = get_metric_definition(metric_name)
            value = _json_value(raw_value, value_type="float")
            status = (
                RecordStatus.COMPUTED if value is not None else RecordStatus.MISSING
            )
            records.append(
                MetricRecord(
                    envelope=make_envelope(
                        RecordKind.METRIC, record_identity, metric_name, statistic
                    ),
                    metric_name=metric_name,
                    metric_group="robustness_metrics",
                    metric_class=MetricClass.ROBUSTNESS,
                    status=status,
                    value=value,
                    unit=definition.unit,
                    direction=definition.direction,
                    evidence_regime=EvidenceRegime.REFERENCE_FREE,
                    statistic=statistic,  # type: ignore[arg-type]
                    reason_code=None
                    if value is not None
                    else "insufficient_observations",
                    message=None
                    if value is not None
                    else "No finite observations were available for robustness aggregation.",
                )
            )
    return tuple(records)


def read_public_metric_frames(run_dir: str | Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    records_path = Path(run_dir) / "results" / "records.jsonl"
    if not records_path.exists():
        return pd.DataFrame(), pd.DataFrame()
    system_rows: dict[tuple[Any, ...], dict[str, Any]] = {}
    chain_rows: dict[tuple[Any, ...], dict[str, Any]] = {}
    with records_path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            record = json.loads(line)
            if record.get("record_kind") != "metric":
                continue
            key = (
                record.get("repeat_id"),
                record.get("model_id"),
                record.get("sample_id"),
                record.get("chain_id"),
            )
            target = chain_rows if record.get("chain_id") is not None else system_rows
            row = target.setdefault(
                key,
                {
                    "model_name": record.get("model_id"),
                    "repeat": record.get("repeat_id"),
                    "diffusion_sample": record.get("sample_id"),
                    "effective_seed": record.get("effective_seed"),
                },
            )
            if target is chain_rows:
                row.update(
                    {
                        "CHAIN_ID": record.get("chain_id"),
                        "ENTITY_TYPE": record.get("entity_type"),
                    }
                )
            name = record["metric_name"]
            if record.get("related_chain_id"):
                name = f"{name}_{record['related_chain_id']}"
            row[name] = (
                record.get("value") if record.get("status") == "computed" else None
            )
    return pd.DataFrame(system_rows.values()), pd.DataFrame(chain_rows.values())


__all__ = [
    "bundle_from_frames",
    "metric_records_from_frames",
    "read_public_metric_frames",
    "robustness_records_from_frame",
]
