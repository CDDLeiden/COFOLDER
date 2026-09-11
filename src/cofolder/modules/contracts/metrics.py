from __future__ import annotations

import math
from collections.abc import Collection
from dataclasses import dataclass
from types import MappingProxyType
from typing import Literal

from .models import (
    EvidenceCompatibilityError,
    EvidenceRegime,
    EvidenceSource,
    MetricClass,
    MetricRecord,
    MetricRecordValidationError,
    OptimizationDirection,
    RecordStatus,
    UnregisteredMetricError,
)


@dataclass(frozen=True, slots=True)
class NumericRange:
    minimum: float | None = None
    maximum: float | None = None
    minimum_inclusive: bool = True
    maximum_inclusive: bool = True


@dataclass(frozen=True, slots=True)
class MetricDefinition:
    name: str
    group: str
    metric_class: MetricClass
    value_type: Literal["float", "integer", "boolean", "string", "json"]
    unit: str | None
    direction: OptimizationDirection
    valid_range: NumericRange | None
    description: str
    scopes: frozenset[str]
    allowed_evidence_regimes: frozenset[EvidenceRegime]


RF = frozenset({EvidenceRegime.REFERENCE_FREE})
REF = frozenset({EvidenceRegime.REFERENCE_STRUCTURE})
POCKET = frozenset({EvidenceRegime.CUSTOM_POCKET})
ANY_EVIDENCE = frozenset(EvidenceRegime)


def _metric(
    name: str,
    group: str,
    metric_class: MetricClass,
    *,
    unit: str | None = None,
    direction: OptimizationDirection = OptimizationDirection.NEUTRAL,
    valid_range: NumericRange | None = None,
    value_type: Literal["float", "integer", "boolean", "string", "json"] = "float",
    scopes: tuple[str, ...] = ("system", "chain"),
    evidence: frozenset[EvidenceRegime] = RF,
    description: str = "COFOLDER public workflow metric.",
) -> MetricDefinition:
    return MetricDefinition(
        name=name,
        group=group,
        metric_class=metric_class,
        value_type=value_type,
        unit=unit,
        direction=direction,
        valid_range=valid_range,
        description=description,
        scopes=frozenset(scopes),
        allowed_evidence_regimes=evidence,
    )


UNIT_INTERVAL = NumericRange(0.0, 1.0)
NONNEGATIVE = NumericRange(0.0, None)

_DEFINITIONS = [
    _metric(
        name,
        "confidence_metrics",
        MetricClass.CONFIDENCE,
        direction=OptimizationDirection.MAXIMIZE,
        valid_range=UNIT_INTERVAL,
    )
    for name in (
        "confidence_score",
        "ptm",
        "iptm",
        "chains_ptm",
        "chain_ptm",
        "sample_ranking_score",
        "avg_plddt",
    )
] + [
    *[
        _metric(
            name,
            "confidence_metrics",
            MetricClass.CONFIDENCE,
            direction=OptimizationDirection.MAXIMIZE,
            valid_range=UNIT_INTERVAL,
            scopes=("system",),
        )
        for name in (
            "ligand_iptm",
            "protein_iptm",
            "complex_plddt",
            "complex_iplddt",
        )
    ],
    *[
        _metric(
            name,
            "confidence_metrics",
            MetricClass.CONFIDENCE,
            direction=OptimizationDirection.MINIMIZE,
            valid_range=NONNEGATIVE,
            scopes=("system",),
        )
        for name in ("complex_pde", "complex_ipde", "complex_pae", "complex_ipae")
    ],
    _metric(
        "gpde",
        "confidence_metrics",
        MetricClass.CONFIDENCE,
        direction=OptimizationDirection.MINIMIZE,
        valid_range=NONNEGATIVE,
    ),
    _metric(
        "disorder",
        "confidence_metrics",
        MetricClass.CONFIDENCE,
        direction=OptimizationDirection.MINIMIZE,
        valid_range=UNIT_INTERVAL,
    ),
    _metric(
        "has_clash",
        "confidence_metrics",
        MetricClass.CONFIDENCE,
        value_type="boolean",
        direction=OptimizationDirection.MINIMIZE,
    ),
    _metric(
        "pair_chains_iptm",
        "confidence_metrics",
        MetricClass.CONFIDENCE,
        direction=OptimizationDirection.MAXIMIZE,
        valid_range=UNIT_INTERVAL,
        scopes=("chain_pair",),
    ),
    _metric(
        "chain_pair_iptm",
        "confidence_metrics",
        MetricClass.CONFIDENCE,
        direction=OptimizationDirection.MAXIMIZE,
        valid_range=UNIT_INTERVAL,
        scopes=("chain_pair",),
    ),
    _metric(
        "bespoke_iptm",
        "confidence_metrics",
        MetricClass.CONFIDENCE,
        direction=OptimizationDirection.MAXIMIZE,
        valid_range=UNIT_INTERVAL,
        scopes=("chain_pair",),
    ),
    _metric(
        "affinity_pred_value",
        "affinity_metrics",
        MetricClass.AFFINITY,
        unit="log10(µM)",
        direction=OptimizationDirection.MINIMIZE,
        scopes=("chain",),
    ),
    _metric(
        "affinity_probability_binary",
        "affinity_metrics",
        MetricClass.BINDING,
        direction=OptimizationDirection.MAXIMIZE,
        valid_range=UNIT_INTERVAL,
        scopes=("chain",),
    ),
    _metric(
        "pIC50",
        "affinity_metrics_ext",
        MetricClass.AFFINITY,
        unit="pIC50",
        direction=OptimizationDirection.MAXIMIZE,
        scopes=("chain",),
    ),
    _metric(
        "IC50_M",
        "affinity_metrics_ext",
        MetricClass.AFFINITY,
        unit="M",
        direction=OptimizationDirection.MINIMIZE,
        valid_range=NONNEGATIVE,
        scopes=("chain",),
    ),
    _metric(
        "IC50_uM",
        "affinity_metrics_ext",
        MetricClass.AFFINITY,
        unit="µM",
        direction=OptimizationDirection.MINIMIZE,
        valid_range=NONNEGATIVE,
        scopes=("chain",),
    ),
    _metric(
        "pIC50_kcal_per_mol",
        "affinity_metrics_ext",
        MetricClass.AFFINITY,
        unit="kcal/mol",
        direction=OptimizationDirection.MAXIMIZE,
        scopes=("chain",),
    ),
    _metric(
        "sasa",
        "structure_metrics",
        MetricClass.STRUCTURAL,
        unit="angstrom^2",
        direction=OptimizationDirection.NEUTRAL,
        valid_range=NONNEGATIVE,
        scopes=("chain",),
    ),
    _metric(
        "sasa_norm_heavy",
        "structure_metrics",
        MetricClass.STRUCTURAL,
        unit="angstrom^2/heavy_atom",
        direction=OptimizationDirection.NEUTRAL,
        valid_range=NONNEGATIVE,
        scopes=("chain",),
    ),
    _metric(
        "ifp_distance",
        "structure_metrics",
        MetricClass.STRUCTURAL,
        value_type="json",
        scopes=("chain",),
    ),
    _metric(
        "ifp_distance_features",
        "structure_metrics",
        MetricClass.STRUCTURAL,
        value_type="json",
        scopes=("chain",),
    ),
    _metric(
        "ifp_prolif",
        "structure_metrics",
        MetricClass.STRUCTURAL,
        value_type="json",
        scopes=("chain",),
    ),
    _metric(
        "ifp_prolif_features",
        "structure_metrics",
        MetricClass.STRUCTURAL,
        value_type="json",
        scopes=("chain",),
    ),
]

for name in (
    "ligand_rmsd_ref",
    "protein_rmsd_ref",
    "ligand_rmsd_ref_mean",
    "protein_rmsd_ref_mean",
):
    _DEFINITIONS.append(
        _metric(
            name,
            "reproduction_metrics",
            MetricClass.REPRODUCTION,
            unit="angstrom",
            direction=OptimizationDirection.MINIMIZE,
            valid_range=NONNEGATIVE,
            evidence=REF,
        )
    )
for name in (
    "sucos_ref",
    "sucos_shape_ref",
    "sucos_feature_ref",
    "sucos_ref_mean",
    "ligand_pose_overlap_ref",
):
    _DEFINITIONS.append(
        _metric(
            name,
            "reproduction_metrics",
            MetricClass.REPRODUCTION,
            direction=OptimizationDirection.MAXIMIZE,
            valid_range=UNIT_INTERVAL,
            evidence=REF,
        )
    )
for name in ("pocket_coverage_ref", "pocket_coverage_ref_mean"):
    _DEFINITIONS.append(
        _metric(
            name,
            "reproduction_metrics",
            MetricClass.REPRODUCTION,
            direction=OptimizationDirection.MAXIMIZE,
            valid_range=UNIT_INTERVAL,
            evidence=REF,
        )
    )
for name in ("pocket_coverage_custom", "pocket_coverage_custom_mean"):
    _DEFINITIONS.append(
        _metric(
            name,
            "reproduction_metrics",
            MetricClass.REPRODUCTION,
            direction=OptimizationDirection.MAXIMIZE,
            valid_range=UNIT_INTERVAL,
            evidence=POCKET,
        )
    )
for name in (
    "bias_prot_sim_train",
    "bias_prot_sim_train_max",
    "bias_prot_sim_train_pairwise",
    "bias_prot_sim_train_pairwise_max",
):
    _DEFINITIONS.append(
        _metric(
            name,
            "bias_metrics",
            MetricClass.TRAINING_PROXIMITY,
            unit="percent_identity",
            direction=OptimizationDirection.MINIMIZE,
            valid_range=NumericRange(0.0, 100.0),
            evidence=RF,
        )
    )
for name in ("bias_lig_sim_train", "bias_lig_sim_train_max"):
    _DEFINITIONS.append(
        _metric(
            name,
            "bias_metrics",
            MetricClass.TRAINING_PROXIMITY,
            unit="tanimoto",
            direction=OptimizationDirection.MINIMIZE,
            valid_range=UNIT_INTERVAL,
            evidence=RF,
        )
    )
IFP_FILTER_EVIDENCE = frozenset(
    {EvidenceRegime.REFERENCE_STRUCTURE, EvidenceRegime.CUSTOM_POCKET}
)
for name, value_type in (
    ("ifp_filter_pass", "boolean"),
    ("ifp_filter_overlap", "float"),
    ("ifp_filter_similarity", "float"),
    ("ifp_filter_threshold", "float"),
):
    _DEFINITIONS.append(
        _metric(
            name,
            "screen_metrics",
            MetricClass.FILTER,
            value_type=value_type,
            direction=(
                OptimizationDirection.MAXIMIZE
                if name in {"ifp_filter_overlap", "ifp_filter_similarity"}
                else OptimizationDirection.NEUTRAL
            ),
            valid_range=UNIT_INTERVAL
            if "overlap" in name or "threshold" in name
            else None,
            scopes=("compound",),
            evidence=IFP_FILTER_EVIDENCE,
        )
    )
for name in (
    "ifp_filter_status",
    "ifp_filter_reason",
    "ifp_filter_reference",
    "ifp_filter_similarity_metric",
    "ifp_filter_policy",
    "ifp_filter_taxonomy",
    "ifp_filter_mapping_status",
):
    _DEFINITIONS.append(
        _metric(
            name,
            "screen_metrics",
            MetricClass.FILTER,
            value_type="string",
            scopes=("compound",),
            evidence=IFP_FILTER_EVIDENCE,
        )
    )
for name in (
    "ifp_filter_required_interactions",
    "ifp_filter_missing_interactions",
    "ifp_filter_mapping_failures",
):
    _DEFINITIONS.append(
        _metric(
            name,
            "screen_metrics",
            MetricClass.FILTER,
            value_type="json",
            scopes=("compound",),
            evidence=IFP_FILTER_EVIDENCE,
        )
    )
for name in ("ifp_cluster_id", "ifp_cluster_status"):
    _DEFINITIONS.append(
        _metric(
            name,
            "screen_metrics",
            MetricClass.FILTER,
            value_type="string",
            scopes=("compound",),
            evidence=RF,
        )
    )
for name in ("oracle_score", "oracle_raw_score", "oracle_gate_adjusted_score"):
    _DEFINITIONS.append(
        _metric(
            name,
            "oracle_metrics",
            MetricClass.ORACLE,
            scopes=("compound",),
            evidence=ANY_EVIDENCE,
        )
    )
_DEFINITIONS.append(
    _metric(
        "struct_rmsd",
        "robustness_metrics",
        MetricClass.ROBUSTNESS,
        unit="angstrom",
        direction=OptimizationDirection.MINIMIZE,
        valid_range=NONNEGATIVE,
        evidence=RF,
    )
)

METRIC_CATALOG = MappingProxyType(
    {definition.name: definition for definition in _DEFINITIONS}
)

SCREEN_METRIC_PROFILES = MappingProxyType(
    {
        "system": (
            "confidence_score",
            "ptm",
            "iptm",
            "bias_prot_sim_train_max",
            "bias_prot_sim_train_pairwise_max",
            "bias_lig_sim_train_max",
            "pocket_coverage_ref",
            "pocket_coverage_ref_mean",
            "pocket_coverage_custom",
            "pocket_coverage_custom_mean",
        ),
        "protein": (
            "chains_ptm",
            "bias_prot_sim_train",
            "bias_prot_sim_train_pairwise",
        ),
        "ligand": (
            "chains_ptm",
            "affinity_pred_value",
            "affinity_probability_binary",
            "pIC50",
            "IC50_M",
            "pIC50_kcal_per_mol",
            "sasa",
            "sasa_norm_heavy",
            "ifp_distance",
            "ifp_distance_features",
            "ifp_prolif",
            "ifp_prolif_features",
            "bias_lig_sim_train",
            "pocket_coverage_ref",
            "pocket_coverage_custom",
            "ligand_pose_overlap_ref",
        ),
    }
)


def get_metric_definition(name: str) -> MetricDefinition:
    try:
        return METRIC_CATALOG[name]
    except KeyError as exc:
        raise UnregisteredMetricError(
            f"Metric {name!r} is not registered in the public catalog."
        ) from exc


def validate_metric_record(record: MetricRecord) -> None:
    definition = get_metric_definition(record.metric_name)
    identity = record.envelope.identity
    if definition.scopes == frozenset({"chain"}) and identity.entity_id is None:
        raise MetricRecordValidationError(
            f"Metric {record.metric_name!r} requires chain-scoped identity."
        )
    if definition.scopes == frozenset({"chain_pair"}) and (
        identity.entity_id is None or identity.related_entity_id is None
    ):
        raise MetricRecordValidationError(
            f"Metric {record.metric_name!r} requires primary and related chain identity."
        )
    if definition.scopes == frozenset({"compound"}) and identity.compound_id is None:
        raise MetricRecordValidationError(
            f"Metric {record.metric_name!r} requires compound identity."
        )
    if not isinstance(record.status, RecordStatus):
        raise MetricRecordValidationError(
            f"Metric {record.metric_name!r} status must be a RecordStatus."
        )
    if not isinstance(record.metric_class, MetricClass):
        raise MetricRecordValidationError(
            f"Metric {record.metric_name!r} class must be a MetricClass."
        )
    if not isinstance(record.direction, OptimizationDirection):
        raise MetricRecordValidationError(
            f"Metric {record.metric_name!r} direction must be an OptimizationDirection."
        )
    if not isinstance(record.evidence_regime, EvidenceRegime):
        raise MetricRecordValidationError(
            f"Metric {record.metric_name!r} evidence_regime must be an EvidenceRegime."
        )
    if record.status not in {
        RecordStatus.COMPUTED,
        RecordStatus.MISSING,
        RecordStatus.UNSUPPORTED,
        RecordStatus.FAILED,
        RecordStatus.NOT_REQUESTED,
    }:
        raise MetricRecordValidationError(
            f"Metric {record.metric_name!r} uses invalid status {record.status!r}."
        )
    expected_group = (
        "robustness_metrics" if record.statistic != "value" else definition.group
    )
    expected_class = (
        MetricClass.ROBUSTNESS
        if record.statistic != "value"
        else definition.metric_class
    )
    if (record.metric_group, record.metric_class, record.unit, record.direction) != (
        expected_group,
        expected_class,
        definition.unit,
        definition.direction,
    ):
        raise MetricRecordValidationError(
            f"Metric {record.metric_name!r} metadata does not match the catalog."
        )
    if record.status == RecordStatus.COMPUTED:
        if record.value is None:
            raise MetricRecordValidationError(
                f"Computed metric {record.metric_name!r} requires a value."
            )
        _validate_value(record.metric_name, record.value, definition)
        allowed_evidence = (
            RF if record.statistic != "value" else definition.allowed_evidence_regimes
        )
        if record.evidence_regime not in allowed_evidence:
            raise EvidenceCompatibilityError(
                f"Metric {record.metric_name!r} is incompatible with evidence regime {record.evidence_regime.value!r}."
            )
        if record.evidence_regime != EvidenceRegime.REFERENCE_FREE and not any(
            source.kind == record.evidence_regime.value for source in record.evidence
        ):
            raise EvidenceCompatibilityError(
                f"Metric {record.metric_name!r} does not include provenance for "
                f"evidence regime {record.evidence_regime.value!r}."
            )
    elif record.value is not None:
        raise MetricRecordValidationError(
            f"Non-computed metric {record.metric_name!r} must use a null value."
        )
    if (
        record.status
        in {RecordStatus.MISSING, RecordStatus.UNSUPPORTED, RecordStatus.FAILED}
        and not record.reason_code
    ):
        raise MetricRecordValidationError(
            f"Metric {record.metric_name!r} with status {record.status.value!r} requires reason_code."
        )


def _validate_value(name: str, value: object, definition: MetricDefinition) -> None:
    expected = definition.value_type
    valid = {
        "float": isinstance(value, (int, float)) and not isinstance(value, bool),
        "integer": isinstance(value, int) and not isinstance(value, bool),
        "boolean": isinstance(value, bool),
        "string": isinstance(value, str),
        "json": isinstance(value, (dict, list, str, int, float, bool)),
    }[expected]
    if not valid:
        raise MetricRecordValidationError(
            f"Metric {name!r} requires value type {expected!r}."
        )
    if isinstance(value, float) and not math.isfinite(value):
        raise MetricRecordValidationError(f"Metric {name!r} must be finite.")
    if (
        definition.valid_range is not None
        and isinstance(value, (int, float))
        and not isinstance(value, bool)
    ):
        bounds = definition.valid_range
        if bounds.minimum is not None and (
            value < bounds.minimum
            or (value == bounds.minimum and not bounds.minimum_inclusive)
        ):
            raise MetricRecordValidationError(
                f"Metric {name!r} is below its valid range."
            )
        if bounds.maximum is not None and (
            value > bounds.maximum
            or (value == bounds.maximum and not bounds.maximum_inclusive)
        ):
            raise MetricRecordValidationError(
                f"Metric {name!r} is above its valid range."
            )


def resolve_evidence_regime(
    metric_name: str, available_evidence: Collection[EvidenceSource]
) -> EvidenceRegime:
    definition = get_metric_definition(metric_name)
    kinds = {source.kind for source in available_evidence}
    preferred = (
        EvidenceRegime.REFERENCE_STRUCTURE
        if "reference_structure" in kinds
        else EvidenceRegime.CUSTOM_POCKET
        if "custom_pocket" in kinds
        else EvidenceRegime.REFERENCE_FREE
    )
    if preferred in definition.allowed_evidence_regimes:
        return preferred
    if EvidenceRegime.REFERENCE_FREE in definition.allowed_evidence_regimes:
        return EvidenceRegime.REFERENCE_FREE
    return next(iter(sorted(definition.allowed_evidence_regimes, key=str)))


def convert_metric_value(source_metric: str, target_metric: str, value: float) -> float:
    source = get_metric_definition(source_metric)
    target = get_metric_definition(target_metric)
    if source.metric_class != target.metric_class:
        raise ValueError(
            f"Metric conversion cannot change class from {source.metric_class.value!r} "
            f"to {target.metric_class.value!r}."
        )
    number = float(value)
    if source_metric == "affinity_pred_value" and target_metric == "IC50_uM":
        return 10**number
    if source_metric == "affinity_pred_value" and target_metric == "pIC50":
        return 6 - number
    if source_metric == "affinity_pred_value" and target_metric == "IC50_M":
        return 10 ** -(6 - number)
    if source_metric == "affinity_pred_value" and target_metric == "pIC50_kcal_per_mol":
        return (6 - number) * 1.364
    if source_metric == "pIC50" and target_metric == "IC50_M":
        return 10 ** (-number)
    raise ValueError(
        f"No registered conversion from {source_metric!r} to {target_metric!r}."
    )
