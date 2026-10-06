"""Shared assembly and publication helpers for recipe-owned public results."""

from __future__ import annotations

from collections.abc import Collection, Iterable
from pathlib import Path

import pandas as pd

from cofolder.modules.contracts import (
    PUBLIC_SCHEMA_VERSION,
    ArtifactReference,
    EvidenceSource,
    FailureStage,
    OutputIdentity,
    PublicManifest,
    PublicOutputBundle,
    RunnerBackendIdentity,
    SeedPlan,
    WorkflowFailureRecord,
    bundle_from_frames,
    write_public_bundle,
)
from cofolder.recipes._diagnostics import workflow_failure


def manifest_provenance(
    backend: RunnerBackendIdentity | None,
    seed_plan: SeedPlan | None,
) -> dict[str, object]:
    """Return complete runner provenance or no partial provenance."""
    if backend is None or seed_plan is None:
        return {}
    return {"backend": backend, "seed_plan": seed_plan}


def write_failure_bundle(
    *,
    output_dir: Path,
    identity: OutputIdentity,
    exc: BaseException,
    stage: FailureStage,
    error_code: str | None,
    backend: RunnerBackendIdentity | None = None,
    seed_plan: SeedPlan | None = None,
) -> WorkflowFailureRecord:
    """Persist one terminal workflow failure and return its typed record."""
    failure = workflow_failure(
        exc,
        identity=identity,
        stage=stage,
        error_code=error_code,
    )
    write_public_bundle(
        PublicOutputBundle(
            manifest=PublicManifest(
                schema_version=PUBLIC_SCHEMA_VERSION,
                identity=identity,
                status="failed",
                **manifest_provenance(backend, seed_plan),
            ),
            records=(failure,),
        ),
        output_dir,
    )
    return failure


def standard_evidence(
    *,
    reference_path: Path | None = None,
    pocket_coverage_reference: object | None = None,
) -> list[EvidenceSource]:
    """Build the evidence shared by Validate-derived workflow outputs."""
    evidence: list[EvidenceSource] = []
    if reference_path is not None:
        evidence.append(
            EvidenceSource(
                kind="reference_structure",
                identifier=reference_path.name,
                path=str(reference_path),
            )
        )
    if pocket_coverage_reference:
        evidence.append(
            EvidenceSource(
                kind="custom_pocket",
                identifier="configured_custom_pocket",
            )
        )
    return evidence


def existing_directory_artifacts(
    root: Path,
    candidates: Iterable[tuple[str, str]],
) -> list[ArtifactReference]:
    """Describe only result directories that exist at publication time."""
    return [
        ArtifactReference(label=label, relative_path=relative_path, kind="directory")
        for label, relative_path in candidates
        if (root / relative_path).exists()
    ]


def write_frame_bundle(
    system_df: pd.DataFrame,
    chain_df: pd.DataFrame,
    *,
    output_dir: Path,
    identity: OutputIdentity,
    evidence: Iterable[EvidenceSource] = (),
    requested_metrics: Collection[str] = (),
    unavailable_groups: Collection[str] = (),
    failures: Iterable[WorkflowFailureRecord] = (),
    artifacts: Iterable[ArtifactReference] = (),
    robustness_df: pd.DataFrame | None = None,
    backend: RunnerBackendIdentity | None = None,
    seed_plan: SeedPlan | None = None,
) -> None:
    """Create and atomically publish a frame-backed public output bundle."""
    bundle = bundle_from_frames(
        system_df,
        chain_df,
        identity=identity,
        evidence=evidence,
        requested_metrics=requested_metrics,
        unavailable_groups=unavailable_groups,
        failures=failures,
        artifacts=artifacts,
        robustness_df=robustness_df,
        backend=backend,
        seed_plan=seed_plan,
    )
    write_public_bundle(bundle, output_dir)
