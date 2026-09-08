from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path
from typing import Literal, TypeAlias

PUBLIC_SCHEMA_VERSION = "1.0.0"
JSONValue: TypeAlias = (
    None | bool | int | float | str | list["JSONValue"] | dict[str, "JSONValue"]
)


class WorkflowKind(StrEnum):
    VALIDATE = "validate"
    SCREEN = "screen"
    ORACLE = "oracle"
    BIAS = "bias"


class RecordKind(StrEnum):
    SUCCESS = "success"
    METRIC = "metric"
    FAILURE = "failure"


class RecordStatus(StrEnum):
    SUCCESS = "success"
    COMPUTED = "computed"
    MISSING = "missing"
    UNSUPPORTED = "unsupported"
    FAILED = "failed"
    NOT_REQUESTED = "not_requested"


class EvidenceRegime(StrEnum):
    REFERENCE_STRUCTURE = "reference_structure"
    CUSTOM_POCKET = "custom_pocket"
    REFERENCE_FREE = "reference_free"


class FailureStage(StrEnum):
    INPUT_VALIDATION = "input_validation"
    PREPARATION = "preparation"
    BACKEND_EXECUTION = "backend_execution"
    OUTPUT_VALIDATION = "output_validation"
    GATHER = "gather"
    ANALYTICS = "analytics"
    AGGREGATION = "aggregation"
    SERIALIZATION = "serialization"


class MetricClass(StrEnum):
    CONFIDENCE = "confidence"
    BINDING = "binding"
    AFFINITY = "affinity"
    STRUCTURAL = "structural"
    REPRODUCTION = "reproduction"
    ROBUSTNESS = "robustness"
    TRAINING_PROXIMITY = "training_proximity"
    FILTER = "filter"
    ORACLE = "oracle"


class OptimizationDirection(StrEnum):
    MAXIMIZE = "maximize"
    MINIMIZE = "minimize"
    NEUTRAL = "neutral"


class BackendVersionStatus(StrEnum):
    DETECTED = "detected"
    UNAVAILABLE = "unavailable"
    UNPARSEABLE = "unparseable"


class SeedOrigin(StrEnum):
    USER_SPECIFIED = "user_specified"
    GENERATED = "generated"


class SeedAdjustment(StrEnum):
    UNCHANGED = "unchanged"
    BACKEND_ADJUSTED = "backend_adjusted"


@dataclass(frozen=True, slots=True)
class RunnerBackendIdentity:
    runner_name: str
    backend_name: str
    version: str | None
    version_status: BackendVersionStatus
    raw_version: str | None = None
    detail: str | None = None


@dataclass(frozen=True, slots=True)
class RepeatSeedProvenance:
    repeat_id: int
    requested_base_seed: int | None
    resolved_base_seed: int
    derived_seed: int
    effective_seed: int
    origin: SeedOrigin
    adjustment: SeedAdjustment = SeedAdjustment.UNCHANGED
    adjustment_reason: str | None = None


@dataclass(frozen=True, slots=True)
class SeedPlan:
    requested_base_seed: int | None
    resolved_base_seed: int
    origin: SeedOrigin
    repeats: tuple[RepeatSeedProvenance, ...]


@dataclass(frozen=True, slots=True)
class OutputIdentity:
    workflow: WorkflowKind
    run_id: str
    system_id: str
    compound_id: str | None = None
    runner_id: str | None = None
    runner_version: str | None = None
    backend_name: str | None = None
    backend_version_status: BackendVersionStatus | None = None
    effective_seed: int | None = None
    repeat_id: int | None = None
    model_id: str | None = None
    sample_id: int | None = None
    entity_id: str | None = None
    entity_type: str | None = None
    chain_id: str | None = None
    related_entity_id: str | None = None
    related_chain_id: str | None = None


@dataclass(frozen=True, slots=True)
class EvidenceSource:
    kind: str
    identifier: str
    path: str | None = None
    sha256: str | None = None


@dataclass(frozen=True, slots=True)
class ArtifactReference:
    label: str
    relative_path: str
    kind: str | None = None
    description: str | None = None


@dataclass(frozen=True, slots=True)
class PublicRecordEnvelope:
    schema_version: str
    record_id: str
    record_kind: RecordKind
    identity: OutputIdentity
    created_at: str


@dataclass(frozen=True, slots=True)
class SuccessRecord:
    envelope: PublicRecordEnvelope
    status: Literal[RecordStatus.SUCCESS] = RecordStatus.SUCCESS
    artifacts: tuple[ArtifactReference, ...] = ()


@dataclass(frozen=True, slots=True)
class MetricRecord:
    envelope: PublicRecordEnvelope
    metric_name: str
    metric_group: str
    metric_class: MetricClass
    status: RecordStatus
    value: JSONValue
    unit: str | None
    direction: OptimizationDirection
    evidence_regime: EvidenceRegime
    evidence: tuple[EvidenceSource, ...] = ()
    statistic: Literal["value", "mean", "std"] = "value"
    reason_code: str | None = None
    message: str | None = None


@dataclass(frozen=True, slots=True)
class WorkflowFailureRecord:
    envelope: PublicRecordEnvelope
    stage: FailureStage
    exception_type: str
    error_code: str
    message: str
    status: Literal[RecordStatus.FAILED] = RecordStatus.FAILED
    retryable: bool = False
    details: Mapping[str, JSONValue] = field(default_factory=dict)


PublicRecord: TypeAlias = SuccessRecord | MetricRecord | WorkflowFailureRecord


@dataclass(frozen=True, slots=True)
class PublicManifest:
    schema_version: str
    identity: OutputIdentity
    status: Literal["success", "partial", "failed"]
    evidence: tuple[EvidenceSource, ...] = ()
    requested_metrics: tuple[str, ...] = ()
    artifacts: tuple[ArtifactReference, ...] = ()
    record_counts: Mapping[str, int] = field(default_factory=dict)
    backend: RunnerBackendIdentity | None = None
    seed_plan: SeedPlan | None = None


@dataclass(frozen=True, slots=True)
class PublicOutputBundle:
    manifest: PublicManifest
    records: tuple[PublicRecord, ...]


@dataclass(frozen=True, slots=True)
class SerializedPublicOutput:
    output_dir: Path
    records_path: Path
    successes_path: Path
    metrics_path: Path
    failures_path: Path
    manifest_path: Path


class PublicContractError(ValueError):
    """Base class for invalid public-contract data."""


class PublicSchemaValidationError(PublicContractError):
    """Raised when a public record or bundle violates the schema."""


class AmbiguousIdentityError(PublicSchemaValidationError):
    """Raised when records cannot be associated with one unambiguous identity."""


class UnregisteredMetricError(PublicSchemaValidationError):
    """Raised when a metric is absent from the public catalog."""


class MetricRecordValidationError(PublicSchemaValidationError):
    """Raised when a metric value or state contradicts its catalog definition."""


class EvidenceCompatibilityError(PublicSchemaValidationError):
    """Raised when a computed metric uses an incompatible evidence regime."""


class PublicSerializationError(RuntimeError):
    """Raised when a validated public bundle cannot be serialized atomically."""


class SeedResolutionError(ValueError):
    """Raised when seed resolution or backend adjustment is invalid."""


class RunnerProvenanceError(ValueError):
    """Raised when execution provenance is contradictory or malformed."""


class WorkflowExecutionError(RuntimeError):
    """Raised after a workflow has persisted a terminal failure bundle."""

    def __init__(
        self,
        message: str,
        *,
        failures: tuple[WorkflowFailureRecord, ...],
        output_dir: Path,
    ) -> None:
        super().__init__(message)
        self.failures = failures
        self.output_dir = output_dir
