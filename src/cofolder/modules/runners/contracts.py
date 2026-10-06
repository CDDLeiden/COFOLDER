from __future__ import annotations

import logging
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Literal

from cofolder.modules.contracts.models import (
    OutputIdentity,
    PublicRecord,
    RepeatSeedProvenance,
    RunnerBackendIdentity,
    RunnerProvenanceError,
    SeedPlan,
    SeedResolutionError,
)

RUNNER_PROVENANCE_COLUMNS = (
    "runner_id",
    "backend_name",
    "runner_version",
    "backend_version_status",
    "effective_seed",
)


@dataclass(frozen=True, slots=True)
class RunnerModelSlot:
    model_id: str
    sample_id: int | None = None


@dataclass(frozen=True, slots=True)
class PlannedExecution:
    repeat_id: int
    model_id: str
    sample_id: int | None
    effective_seed: int


@dataclass(frozen=True, slots=True)
class RunnerExecutionPlan:
    backend: RunnerBackendIdentity
    seed_plan: SeedPlan
    executions: tuple[PlannedExecution, ...]


def runner_provenance_values(
    request: RunnerExecutionRequest,
) -> dict[str, str | int | None]:
    backend = request.backend_identity
    return {
        "runner_id": backend.runner_name,
        "backend_name": backend.backend_name,
        "runner_version": backend.version,
        "backend_version_status": backend.version_status.value,
        "effective_seed": request.seed_provenance.effective_seed,
    }


def attach_runner_provenance(frame: Any, request: RunnerExecutionRequest) -> Any:
    """Attach canonical execution provenance to every normalized frame row."""
    frame = frame.copy()
    for column, value in runner_provenance_values(request).items():
        frame[column] = value
    return frame


def attach_sample_provenance(
    records: Sequence[Mapping[str, Any]], request: RunnerExecutionRequest
) -> list[dict[str, Any]]:
    values = runner_provenance_values(request)
    return [{**dict(record), **values} for record in records]


def backend_manifest_value(identity: RunnerBackendIdentity) -> dict[str, Any]:
    value = asdict(identity)
    value["version_status"] = identity.version_status.value
    return value


def seed_manifest_value(seed: RepeatSeedProvenance) -> dict[str, Any]:
    value = asdict(seed)
    value["origin"] = seed.origin.value
    value["adjustment"] = seed.adjustment.value
    return value


@dataclass(frozen=True, slots=True)
class RunnerChainIdentity:
    """Explicit mapping from a runner chain index to public entity identity."""

    conf_chain_id: int
    entity_id: str
    entity_type: str
    chain_id: str
    ligand_molecule_id: str


def build_runner_chain_identities(system_obj: Any) -> tuple[RunnerChainIdentity, ...]:
    """Normalize system chains into one validated mapping used by runners and gather."""

    from cofolder.modules.input.system import iter_system_chains

    identities: list[RunnerChainIdentity] = []
    seen_chain_ids: set[str] = set()
    for conf_chain_id, chain in enumerate(iter_system_chains(system_obj)):
        if chain.chain_id in seen_chain_ids:
            raise ValueError(f"System declares duplicate chain ID {chain.chain_id!r}.")
        seen_chain_ids.add(chain.chain_id)
        molecule_id = (
            chain.entity_data.get("ccd")
            or chain.entity_data.get("smiles")
            or "UNKNOWN_LIGAND"
            if chain.entity_type == "ligand"
            else f"{chain.entity_type}_{chain.chain_id}"
        )
        identities.append(
            RunnerChainIdentity(
                conf_chain_id=conf_chain_id,
                entity_id=f"entity:{chain.sequence_index}",
                entity_type=chain.entity_type,
                chain_id=chain.chain_id,
                ligand_molecule_id=str(molecule_id),
            )
        )
    return tuple(identities)


def build_runner_public_records(
    system_df: Any,
    chain_df: Any,
    request: RunnerExecutionRequest,
) -> list[PublicRecord]:
    """Convert normalized runner frames to typed, in-memory public records."""

    if request.identity is None:
        return []
    from cofolder.modules.contracts.adapters import metric_records_from_frames

    normalized_chain_df = chain_df.copy()
    mapping = {item.conf_chain_id: item for item in request.chain_identities}
    if len(mapping) != len(request.chain_identities):
        raise ValueError("Runner request contains duplicate chain identity mappings.")
    if not normalized_chain_df.empty and not mapping:
        raise ValueError(
            "Typed runner records require an explicit chain identity mapping."
        )
    if not normalized_chain_df.empty:
        observed = {
            int(value) for value in normalized_chain_df["conf_chain_id"].dropna()
        }
        missing = sorted(observed - mapping.keys())
        if missing:
            raise ValueError(
                f"Runner records contain unmapped conf_chain_id values: {missing}."
            )
        normalized_chain_df["CHAIN_ID"] = normalized_chain_df["conf_chain_id"].map(
            lambda value: mapping[int(value)].chain_id
        )
        normalized_chain_df["ENTITY_ID"] = normalized_chain_df["conf_chain_id"].map(
            lambda value: mapping[int(value)].entity_id
        )
        normalized_chain_df["ENTITY_TYPE"] = normalized_chain_df["conf_chain_id"].map(
            lambda value: mapping[int(value)].entity_type
        )
        normalized_chain_df["ligand_molecule_id"] = normalized_chain_df[
            "conf_chain_id"
        ].map(lambda value: mapping[int(value)].ligand_molecule_id)
    return list(
        metric_records_from_frames(
            system_df,
            normalized_chain_df,
            base_identity=request.identity,
        )
    )


@dataclass(frozen=True, slots=True)
class RunnerInputCapabilities:
    """Molecular-system inputs a runner can faithfully consume."""

    entity_types: frozenset[str]
    constraint_types: frozenset[str]
    max_pocket_constraints: int | None = None
    required_pocket_distance: float | None = None
    supports_constraint_force: bool = True
    pocket_contacts_must_be_polymers: bool = False


@dataclass(slots=True)
class RunnerRuntime:
    """Shared runtime metadata that recipes may safely consume."""

    diffusion_samples: int = 1
    cache_path: str | None = None
    model_name: str | None = None


RunnerMetricState = Literal[
    "computed", "unsupported", "missing", "failed", "not_requested"
]


@dataclass(slots=True)
class RunnerMetricOutcome:
    """Explicit per-metric-group state produced at the runner boundary."""

    state: RunnerMetricState
    required_columns: tuple[str, ...] = ()
    required_artifacts: tuple[str, ...] = ()
    message: str | None = None


@dataclass(slots=True)
class RunnerCompanionArtifact:
    """Normalized non-tabular artifact exposed alongside CSV outputs."""

    label: str
    relative_path: str
    kind: str | None = None
    description: str | None = None


@dataclass(slots=True)
class RunnerNormalizedBundle:
    """Authoritative normalized output bundle consumed by shared validation."""

    runner_name: str
    raw_output_dir: Path
    normalized_dir: Path
    structures_dir: Path
    system_metrics_path: Path
    chain_metrics_path: Path
    manifest_path: Path
    diffusion_samples: int = 1
    capabilities: set[str] = field(default_factory=set)
    runtime: RunnerRuntime = field(default_factory=RunnerRuntime)
    sample_records: list[dict[str, Any]] = field(default_factory=list)
    metric_outcomes: dict[str, RunnerMetricOutcome] = field(default_factory=dict)
    companion_artifacts: list[RunnerCompanionArtifact] = field(default_factory=list)
    records: list[PublicRecord] = field(default_factory=list)
    chain_identities: tuple[RunnerChainIdentity, ...] = ()
    backend_identity: RunnerBackendIdentity | None = None
    seed_provenance: RepeatSeedProvenance | None = None


def merge_runner_runtime(*values: RunnerRuntime | None) -> RunnerRuntime:
    """Merge typed runtime metadata without falling back to ad hoc dict lookups."""

    merged = RunnerRuntime()
    for value in values:
        if value is None:
            continue
        if value.cache_path is not None:
            merged.cache_path = value.cache_path
        if value.model_name is not None:
            merged.model_name = value.model_name
        # Authoritative typed runtimes are merged in precedence order, so an
        # explicit single-sample value from a later runtime must still win.
        merged.diffusion_samples = value.diffusion_samples
    return merged


def copy_metric_outcomes(
    metric_outcomes: Mapping[str, RunnerMetricOutcome] | None,
) -> dict[str, RunnerMetricOutcome]:
    """Copy explicit metric outcomes into a mutable normalized-bundle mapping."""

    copied: dict[str, RunnerMetricOutcome] = {}
    if metric_outcomes is None:
        return copied

    for group_name, outcome in metric_outcomes.items():
        if not isinstance(outcome, RunnerMetricOutcome):
            raise TypeError(
                "metric_outcomes values must be RunnerMetricOutcome instances, "
                f"got {type(outcome)!r} for group {group_name!r}."
            )
        copied[str(group_name)] = RunnerMetricOutcome(
            state=outcome.state,
            required_columns=tuple(outcome.required_columns),
            required_artifacts=tuple(outcome.required_artifacts),
            message=outcome.message,
        )
    return copied


def copy_companion_artifacts(
    companion_artifacts: Sequence[RunnerCompanionArtifact] | None,
) -> list[RunnerCompanionArtifact]:
    """Copy normalized companion-artifact metadata into bundle-owned storage."""

    copied: list[RunnerCompanionArtifact] = []
    if companion_artifacts is None:
        return copied

    for artifact in companion_artifacts:
        if not isinstance(artifact, RunnerCompanionArtifact):
            raise TypeError(
                "companion_artifacts values must be RunnerCompanionArtifact instances, "
                f"got {type(artifact)!r}."
            )
        copied.append(
            RunnerCompanionArtifact(
                label=str(artifact.label),
                relative_path=str(artifact.relative_path),
                kind=None if artifact.kind is None else str(artifact.kind),
                description=None
                if artifact.description is None
                else str(artifact.description),
            )
        )
    return copied


@dataclass(slots=True, init=False)
class RunnerPreparationResult:
    system_obj: Any
    options_obj: Any
    warnings: list[str] = field(default_factory=list)
    runtime: RunnerRuntime = field(default_factory=RunnerRuntime)

    def __init__(
        self,
        system_obj: Any,
        options_obj: Any,
        warnings: list[str] | None = None,
        runtime: RunnerRuntime | None = None,
    ) -> None:
        self.system_obj = system_obj
        self.options_obj = options_obj
        self.warnings = list(warnings or [])
        self.runtime = runtime or RunnerRuntime()


@dataclass(slots=True)
class RunnerExecutionRequest:
    runner_name: str
    system_name: str
    system_path: Path
    system_obj: Any
    options_path: Path
    options_obj: Any
    repeat: int
    seed: int
    seed_provenance: RepeatSeedProvenance
    backend_identity: RunnerBackendIdentity
    repeat_dir: Path
    raw_dir: Path
    logger: logging.Logger | None
    timings: Any | None = None
    label_prefix: str | None = None
    runtime: RunnerRuntime = field(default_factory=RunnerRuntime)
    identity: OutputIdentity | None = None
    chain_identities: tuple[RunnerChainIdentity, ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.seed_provenance, RepeatSeedProvenance):
            raise SeedResolutionError(
                "Runner request seed_provenance must be RepeatSeedProvenance."
            )
        if not isinstance(self.backend_identity, RunnerBackendIdentity):
            raise RunnerProvenanceError(
                "Runner request backend_identity must be RunnerBackendIdentity."
            )
        if (
            isinstance(self.seed, bool)
            or not isinstance(self.seed, int)
            or not (0 <= self.seed <= 2**32 - 1)
        ):
            raise SeedResolutionError(
                "Runner request seed must be an integer between 0 and 4294967295."
            )
        if self.repeat != self.seed_provenance.repeat_id:
            raise SeedResolutionError(
                "Runner request repeat does not match seed provenance."
            )
        if self.seed != self.seed_provenance.effective_seed:
            raise SeedResolutionError(
                "Runner request seed does not match its effective seed."
            )
        if self.runner_name != self.backend_identity.runner_name:
            raise RunnerProvenanceError(
                "Runner request name does not match backend identity."
            )


@dataclass(slots=True, init=False)
class RunnerExecutionResult:
    runner_name: str
    raw_output_dir: Path
    normalized_dir: Path
    structures_dir: Path
    system_metrics_path: Path
    chain_metrics_path: Path
    manifest_path: Path
    diffusion_samples: int = 1
    capabilities: set[str] = field(default_factory=set)
    warnings: list[str] = field(default_factory=list)
    runtime: RunnerRuntime = field(default_factory=RunnerRuntime)
    sample_records: list[dict[str, Any]] = field(default_factory=list)
    records: list[PublicRecord] = field(default_factory=list)
    backend_identity: RunnerBackendIdentity | None = None
    seed_provenance: RepeatSeedProvenance | None = None
    normalized_bundle: RunnerNormalizedBundle = field(init=False)

    def __init__(
        self,
        runner_name: str,
        raw_output_dir: Path,
        normalized_dir: Path,
        structures_dir: Path,
        system_metrics_path: Path,
        chain_metrics_path: Path,
        manifest_path: Path,
        diffusion_samples: int = 1,
        capabilities: set[str] | None = None,
        warnings: list[str] | None = None,
        runtime: RunnerRuntime | None = None,
        sample_records: list[dict[str, Any]] | None = None,
        metric_outcomes: Mapping[str, RunnerMetricOutcome] | None = None,
        companion_artifacts: Sequence[RunnerCompanionArtifact] | None = None,
        chain_identities: Sequence[RunnerChainIdentity] | None = None,
        records: Sequence[PublicRecord] | None = None,
        backend_identity: RunnerBackendIdentity | None = None,
        seed_provenance: RepeatSeedProvenance | None = None,
    ) -> None:
        self.runner_name = runner_name
        self.raw_output_dir = raw_output_dir
        self.normalized_dir = normalized_dir
        self.structures_dir = structures_dir
        self.system_metrics_path = system_metrics_path
        self.chain_metrics_path = chain_metrics_path
        self.manifest_path = manifest_path
        self.diffusion_samples = int(diffusion_samples)
        self.capabilities = set(capabilities or set())
        self.warnings = list(warnings or [])
        self.runtime = runtime or RunnerRuntime(
            diffusion_samples=self.diffusion_samples
        )
        self.sample_records = list(sample_records or [])
        self.records = list(records or [])
        self.backend_identity = backend_identity
        self.seed_provenance = seed_provenance
        self.normalized_bundle = RunnerNormalizedBundle(
            runner_name=self.runner_name,
            raw_output_dir=self.raw_output_dir,
            normalized_dir=self.normalized_dir,
            structures_dir=self.structures_dir,
            system_metrics_path=self.system_metrics_path,
            chain_metrics_path=self.chain_metrics_path,
            manifest_path=self.manifest_path,
            diffusion_samples=self.diffusion_samples,
            capabilities=set(self.capabilities),
            runtime=self.runtime,
            sample_records=list(self.sample_records),
            metric_outcomes=copy_metric_outcomes(metric_outcomes),
            companion_artifacts=copy_companion_artifacts(companion_artifacts),
            chain_identities=tuple(chain_identities or ()),
            records=list(self.records),
            backend_identity=self.backend_identity,
            seed_provenance=self.seed_provenance,
        )

    @property
    def metric_outcomes(self) -> dict[str, RunnerMetricOutcome]:
        return self.normalized_bundle.metric_outcomes

    @property
    def companion_artifacts(self) -> list[RunnerCompanionArtifact]:
        return self.normalized_bundle.companion_artifacts
