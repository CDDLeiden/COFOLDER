from __future__ import annotations

import json

import pandas as pd
import pytest

from cofolder.modules.contracts import (
    METRIC_CATALOG,
    PUBLIC_SCHEMA_VERSION,
    AmbiguousIdentityError,
    EvidenceRegime,
    EvidenceSource,
    FailureStage,
    MetricClass,
    MetricRecord,
    OptimizationDirection,
    OutputIdentity,
    PublicManifest,
    PublicOutputBundle,
    PublicSchemaValidationError,
    RecordKind,
    RecordStatus,
    UnregisteredMetricError,
    WorkflowKind,
    aggregate_metric_records,
    bundle_from_frames,
    convert_metric_value,
    failure_from_exception,
    make_envelope,
    metric_records_from_frames,
    validate_public_bundle,
    write_public_bundle,
)


def _identity(**changes):
    values = {
        "workflow": WorkflowKind.VALIDATE,
        "run_id": "run-1",
        "system_id": "system-1",
        "runner_id": "boltz2",
    }
    values.update(changes)
    return OutputIdentity(**values)


def test_dataframe_adapter_emits_complete_identity_and_evidence():
    system_df = pd.DataFrame(
        [
            {
                "model_name": "model-0",
                "repeat": 1,
                "diffusion_sample": 0,
                "confidence_score": 0.8,
            }
        ]
    )
    chain_df = pd.DataFrame(
        [
            {
                "model_name": "model-0",
                "repeat": 1,
                "diffusion_sample": 0,
                "conf_chain_id": 0,
                "CHAIN_ID": "A",
                "ENTITY_ID": "entity:0",
                "ENTITY_TYPE": "protein",
                "chains_ptm": 0.7,
            }
        ]
    )

    records = metric_records_from_frames(
        system_df,
        chain_df,
        base_identity=_identity(),
        evidence=(EvidenceSource("predicted_structure", "model-0"),),
    )

    metrics = [
        record for record in records if record.envelope.record_kind == RecordKind.METRIC
    ]
    assert {record.metric_name for record in metrics} == {
        "confidence_score",
        "chains_ptm",
    }
    assert all(record.envelope.identity.repeat_id == 1 for record in metrics)
    chain_metric = next(
        record for record in metrics if record.metric_name == "chains_ptm"
    )
    assert chain_metric.envelope.identity.entity_id == "entity:0"
    assert chain_metric.envelope.identity.chain_id == "A"
    assert chain_metric.evidence_regime == EvidenceRegime.REFERENCE_FREE


def test_incompatible_reference_metric_is_explicitly_unsupported():
    records = metric_records_from_frames(
        pd.DataFrame(
            [
                {
                    "model_name": "model-0",
                    "repeat": 1,
                    "diffusion_sample": 0,
                    "ligand_rmsd_ref_mean": None,
                }
            ]
        ),
        pd.DataFrame(),
        base_identity=_identity(),
        requested_metrics={"reproduction_metrics"},
    )
    metric = next(
        record for record in records if record.envelope.record_kind == RecordKind.METRIC
    )
    assert metric.status == RecordStatus.UNSUPPORTED
    assert metric.reason_code == "reference_structure_unavailable"


def test_unregistered_metric_cannot_enter_public_output():
    with pytest.raises(UnregisteredMetricError, match="not registered"):
        metric_records_from_frames(
            pd.DataFrame([{"model_name": "m", "repeat": 1, "mystery": 2.0}]),
            pd.DataFrame(),
            base_identity=_identity(),
        )


def test_writer_creates_canonical_files_and_round_trippable_jsonl(tmp_path):
    bundle = bundle_from_frames(
        pd.DataFrame(
            [
                {
                    "model_name": "model-0",
                    "repeat": 1,
                    "diffusion_sample": 0,
                    "ptm": 0.5,
                }
            ]
        ),
        pd.DataFrame(),
        identity=_identity(),
        requested_metrics={"confidence_metrics"},
    )
    output = write_public_bundle(bundle, tmp_path)

    assert output.records_path.exists()
    assert output.successes_path.exists()
    assert output.metrics_path.exists()
    assert output.failures_path.exists()
    assert output.manifest_path.exists()
    assert not (tmp_path / "system_metrics.csv").exists()
    records = [
        json.loads(line) for line in output.records_path.read_text().splitlines()
    ]
    assert {record["record_kind"] for record in records} == {"success", "metric"}
    assert json.loads(output.manifest_path.read_text())["schema_version"] == "1.0.0"


def test_failure_factory_preserves_stage_and_actionable_context(tmp_path):
    identity = _identity(repeat_id=2)
    failure = failure_from_exception(
        RuntimeError("backend exploded"),
        identity=identity,
        stage=FailureStage.BACKEND_EXECUTION,
        error_code="runner_backend_execution_failed",
        details={"seed": 42},
    )
    assert failure.exception_type == "RuntimeError"
    assert failure.envelope.identity.repeat_id == 2
    assert failure.details["seed"] == 42


def test_record_ids_are_deterministic():
    identity = _identity()
    first = make_envelope(RecordKind.METRIC, identity, "ptm", "value")
    second = make_envelope(RecordKind.METRIC, identity, "ptm", "value")
    assert first.record_id == second.record_id


def test_pairwise_metric_uses_related_identity_instead_of_dynamic_name():
    records = metric_records_from_frames(
        pd.DataFrame(),
        pd.DataFrame(
            [
                {
                    "model_name": "model-0",
                    "repeat": 1,
                    "diffusion_sample": 0,
                    "conf_chain_id": 0,
                    "CHAIN_ID": "A",
                    "ENTITY_ID": "entity:0",
                    "ENTITY_TYPE": "dna",
                    "pair_chains_iptm_1": 0.75,
                },
                {
                    "model_name": "model-0",
                    "repeat": 1,
                    "diffusion_sample": 0,
                    "conf_chain_id": 1,
                    "CHAIN_ID": "B",
                    "ENTITY_ID": "entity:0",
                    "ENTITY_TYPE": "dna",
                },
            ]
        ),
        base_identity=_identity(),
    )
    metric = next(record for record in records if isinstance(record, MetricRecord))
    assert metric.metric_name == "pair_chains_iptm"
    assert metric.envelope.identity.chain_id == "A"
    assert metric.envelope.identity.related_chain_id == "B"
    assert metric.envelope.identity.related_entity_id == "entity:0"


def test_aggregation_preserves_identity_and_emits_robustness_metadata():
    records = metric_records_from_frames(
        pd.DataFrame(
            [
                {"model_name": "m1", "repeat": 1, "ptm": 0.4},
                {"model_name": "m2", "repeat": 2, "ptm": 0.8},
            ]
        ),
        pd.DataFrame(),
        base_identity=_identity(),
    )
    metrics = [record for record in records if isinstance(record, MetricRecord)]
    aggregate = aggregate_metric_records(metrics, statistic="mean")[0]
    assert aggregate.value == pytest.approx(0.6)
    assert aggregate.metric_class == MetricClass.ROBUSTNESS
    assert aggregate.envelope.identity.repeat_id is None
    validate_public_bundle(
        PublicOutputBundle(
            manifest=PublicManifest(
                schema_version=PUBLIC_SCHEMA_VERSION,
                identity=_identity(),
                status="failed",
            ),
            records=(aggregate,),
        )
    )


def test_aggregation_rejects_conflicting_metric_metadata():
    records = metric_records_from_frames(
        pd.DataFrame(
            [
                {"model_name": "m1", "repeat": 1, "ptm": 0.4},
                {"model_name": "m2", "repeat": 2, "ptm": 0.8},
            ]
        ),
        pd.DataFrame(),
        base_identity=_identity(),
    )
    metrics = [record for record in records if isinstance(record, MetricRecord)]
    conflicting = MetricRecord(
        envelope=metrics[1].envelope,
        metric_name="ptm",
        metric_group="confidence_metrics",
        metric_class=MetricClass.CONFIDENCE,
        status=RecordStatus.COMPUTED,
        value=0.8,
        unit="wrong-unit",
        direction=OptimizationDirection.MAXIMIZE,
        evidence_regime=EvidenceRegime.REFERENCE_FREE,
    )
    with pytest.raises(AmbiguousIdentityError, match="conflicting metadata"):
        aggregate_metric_records((metrics[0], conflicting), statistic="mean")


def test_aggregation_rejects_duplicate_observation_identity():
    records = metric_records_from_frames(
        pd.DataFrame([{"model_name": "m1", "repeat": 1, "ptm": 0.4}]),
        pd.DataFrame(),
        base_identity=_identity(),
    )
    metric = next(record for record in records if isinstance(record, MetricRecord))
    duplicate = MetricRecord(
        envelope=metric.envelope,
        metric_name=metric.metric_name,
        metric_group=metric.metric_group,
        metric_class=metric.metric_class,
        status=metric.status,
        value=0.8,
        unit=metric.unit,
        direction=metric.direction,
        evidence_regime=metric.evidence_regime,
    )
    with pytest.raises(AmbiguousIdentityError, match="duplicate observations"):
        aggregate_metric_records((metric, duplicate), statistic="mean")


def test_unrequested_metric_uses_null_not_requested_state():
    records = metric_records_from_frames(
        pd.DataFrame([{"model_name": "m1", "repeat": 1, "ptm": 0.4}]),
        pd.DataFrame(),
        base_identity=_identity(),
        requested_metrics={"sasa"},
    )
    metric = next(record for record in records if isinstance(record, MetricRecord))
    assert metric.status == RecordStatus.NOT_REQUESTED
    assert metric.value is None


def test_failure_details_reject_non_finite_json_before_writing(tmp_path):
    failure = failure_from_exception(
        RuntimeError("bad details"),
        identity=_identity(),
        stage=FailureStage.ANALYTICS,
        error_code="bad_details",
        details={"score": float("nan")},
    )
    bundle = PublicOutputBundle(
        manifest=PublicManifest(
            schema_version=PUBLIC_SCHEMA_VERSION,
            identity=_identity(),
            status="failed",
        ),
        records=(failure,),
    )
    with pytest.raises(PublicSchemaValidationError, match="non-finite"):
        write_public_bundle(bundle, tmp_path)
    assert not (tmp_path / "records.jsonl").exists()


def test_writer_removes_stale_legacy_public_files(tmp_path):
    legacy = tmp_path / "system_metrics.csv"
    legacy.write_text("old,data\n", encoding="utf-8")
    bundle = bundle_from_frames(
        pd.DataFrame([{"model_name": "m", "repeat": 1, "ptm": 0.5}]),
        pd.DataFrame(),
        identity=_identity(),
        requested_metrics={"confidence_metrics"},
    )
    write_public_bundle(bundle, tmp_path)
    assert not legacy.exists()


def test_metric_catalog_has_complete_metadata_and_validated_conversions():
    assert METRIC_CATALOG
    for name, definition in METRIC_CATALOG.items():
        assert definition.name == name
        assert definition.group
        assert definition.description
        assert definition.scopes
        assert definition.allowed_evidence_regimes

    assert convert_metric_value("affinity_pred_value", "pIC50", 1.0) == 5.0
    assert convert_metric_value("affinity_pred_value", "IC50_uM", 2.0) == 100.0
    with pytest.raises(ValueError, match="cannot change class"):
        convert_metric_value("affinity_pred_value", "confidence_score", 1.0)
