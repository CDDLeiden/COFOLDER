from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Collection

import pandas as pd

from cofolder.modules.contracts import RunnerProvenanceError
from cofolder.modules.runners.contracts import (
    RUNNER_PROVENANCE_COLUMNS,
    RunnerCompanionArtifact,
    RunnerExecutionResult,
    RunnerMetricOutcome,
    RunnerNormalizedBundle,
    backend_manifest_value,
    seed_manifest_value,
)

CANONICAL_SYSTEM_COLUMNS = (
    "cif_file",
    "model_name",
    "repeat",
    "diffusion_sample",
)
CANONICAL_CHAIN_COLUMNS = (
    "conf_chain_id",
    "cif_file",
    "model_name",
    "repeat",
    "diffusion_sample",
)
STRUCTURE_SUFFIXES = {".cif", ".mmcif", ".pdb"}
MetricColumnAliases = dict[str, tuple[str, ...]]
METRIC_PAYLOAD_SECTIONS = ("system_metrics", "chain_metrics")
CONFIDENCE_SYSTEM_MARKERS = (
    "ptm",
    "iptm",
    "confidence_score",
    "sample_ranking_score",
    "avg_plddt",
    "gpde",
    "disorder",
    "has_clash",
)
CONFIDENCE_CHAIN_MARKERS = (
    "chains_ptm",
    "chain_ptm",
    "prefix:chain_pair_iptm_",
    "prefix:bespoke_iptm_",
)
CONFIDENCE_ARTIFACT_LABELS = ("plddt", "pae", "pde")
DEFAULT_RUNNER_METRIC_REQUIREMENTS = {
    "confidence_metrics": {},
    "affinity_metrics": {
        "chain_metrics": {
            "affinity_pred_value": ("affinity_pred_value",),
            "affinity_probability_binary": ("affinity_probability_binary",),
        },
    },
    "affinity_metrics_ext": {
        "chain_metrics": {
            "pIC50": ("pIC50",),
            "IC50_M": ("IC50_M",),
            "pIC50_kcal_per_mol": ("pIC50_kcal_per_mol",),
        },
    },
}


class RunnerBundleValidationError(RuntimeError):
    """Raised when a runner hands invalid or unusable normalized outputs to shared code."""


def validate_runner_bundle(
    result: RunnerExecutionResult,
    *,
    requested_metric_groups: Collection[str] | None = None,
) -> RunnerNormalizedBundle:
    """Validate normalized runner outputs before shared gather or analytics proceed."""

    bundle = result.normalized_bundle
    requested_groups = sorted(set(requested_metric_groups or ()))
    _validate_bundle_layout(bundle)
    manifest = _validate_manifest(bundle)

    system_df = _read_csv(bundle, bundle.system_metrics_path, "system_metrics.csv")
    chain_df = _read_csv(bundle, bundle.chain_metrics_path, "chain_metrics.csv")
    _validate_required_columns(bundle, system_df, chain_df)
    _validate_execution_provenance(bundle, manifest, system_df, chain_df)
    _validate_chain_identity_mapping(bundle, chain_df)
    _validate_typed_records(bundle)
    structure_files = _validate_structures(bundle)
    _validate_referenced_structures(bundle, system_df, chain_df, structure_files)
    _validate_sample_records(bundle, structure_files)
    _validate_companion_artifacts(bundle, manifest)

    explicit_outcomes = dict(result.metric_outcomes)
    _validate_known_metric_groups(bundle, explicit_outcomes.keys(), requested_groups)

    validated_outcomes: dict[str, RunnerMetricOutcome] = {}
    for group_name, outcome in explicit_outcomes.items():
        validated_outcomes[group_name] = _validate_metric_outcome(
            bundle=bundle,
            group_name=group_name,
            outcome=outcome,
            system_df=system_df,
            chain_df=chain_df,
        )

    for group_name in requested_groups:
        outcome = explicit_outcomes.get(group_name)
        if outcome is None:
            # Story 2.4 retires the temporary "derive it from CSV shape" migration shim.
            raise RunnerBundleValidationError(
                f"{_bundle_prefix(bundle)} omits an explicit outcome for requested metric group "
                f"{group_name!r}."
            )
        validated_outcomes[group_name] = _validate_metric_outcome(
            bundle=bundle,
            group_name=group_name,
            outcome=outcome,
            system_df=system_df,
            chain_df=chain_df,
        )

    bundle.metric_outcomes = validated_outcomes
    _raise_for_non_continuable_outcomes(bundle, requested_groups=requested_groups)
    return bundle


def _validate_chain_identity_mapping(
    bundle: RunnerNormalizedBundle,
    chain_df: pd.DataFrame,
) -> None:
    if not bundle.chain_identities:
        return
    mapping = {item.conf_chain_id: item for item in bundle.chain_identities}
    if len(mapping) != len(bundle.chain_identities):
        raise RunnerBundleValidationError(
            f"Runner '{bundle.runner_name}' returned duplicate conf_chain_id values "
            "in its chain identity mapping."
        )
    if any(
        not item.entity_id.strip()
        or not item.entity_type.strip()
        or not item.chain_id.strip()
        for item in bundle.chain_identities
    ):
        raise RunnerBundleValidationError(
            f"Runner '{bundle.runner_name}' returned an incomplete chain identity mapping."
        )
    try:
        observed = {int(value) for value in chain_df["conf_chain_id"].dropna()}
    except (TypeError, ValueError) as exc:
        raise RunnerBundleValidationError(
            f"Runner '{bundle.runner_name}' emitted a non-integer conf_chain_id."
        ) from exc
    missing = sorted(observed - mapping.keys())
    if missing:
        raise RunnerBundleValidationError(
            f"Runner '{bundle.runner_name}' emitted chain indices without normalized "
            f"identity mappings: {missing}."
        )


def _validate_typed_records(bundle: RunnerNormalizedBundle) -> None:
    if not bundle.records:
        return
    from cofolder.modules.contracts import PublicContractError, validate_public_record

    try:
        for record in bundle.records:
            validate_public_record(record)
            if record.envelope.identity.runner_id != bundle.runner_name:
                raise RunnerBundleValidationError(
                    f"Runner '{bundle.runner_name}' emitted a typed record with runner_id "
                    f"{record.envelope.identity.runner_id!r}."
                )
            if (
                bundle.backend_identity is not None
                and bundle.seed_provenance is not None
            ):
                identity = record.envelope.identity
                if (
                    identity.runner_version != bundle.backend_identity.version
                    or identity.backend_name != bundle.backend_identity.backend_name
                    or identity.backend_version_status
                    != bundle.backend_identity.version_status
                    or identity.effective_seed
                    != bundle.seed_provenance.effective_seed
                ):
                    raise RunnerProvenanceError(
                        f"{_bundle_prefix(bundle)}: typed record provenance does not "
                        "match the bundle."
                    )
    except PublicContractError as exc:
        raise RunnerBundleValidationError(
            f"Runner '{bundle.runner_name}' emitted an invalid typed public record: {exc}"
        ) from exc


def _bundle_prefix(bundle: RunnerNormalizedBundle) -> str:
    return (
        f"Runner {bundle.runner_name!r} returned an invalid normalized bundle "
        f"at {str(bundle.normalized_dir)!r}"
    )


def _validate_bundle_layout(bundle: RunnerNormalizedBundle) -> None:
    if not bundle.normalized_dir.is_dir():
        raise RunnerBundleValidationError(
            f"{_bundle_prefix(bundle)}: normalized directory is missing."
        )

    path_expectations: tuple[tuple[str, Path, Path], ...] = (
        ("system_metrics.csv", bundle.system_metrics_path, bundle.normalized_dir),
        ("chain_metrics.csv", bundle.chain_metrics_path, bundle.normalized_dir),
        ("manifest.json", bundle.manifest_path, bundle.normalized_dir),
    )
    for label, path, expected_parent in path_expectations:
        if path.parent != expected_parent:
            raise RunnerBundleValidationError(
                f"{_bundle_prefix(bundle)}: {label} must live directly under "
                f"{str(bundle.normalized_dir)!r}."
            )
        if not path.is_file():
            raise RunnerBundleValidationError(
                f"{_bundle_prefix(bundle)}: required file {label!r} is missing."
            )

    if bundle.structures_dir.parent != bundle.normalized_dir:
        raise RunnerBundleValidationError(
            f"{_bundle_prefix(bundle)}: structures directory must live under "
            f"{str(bundle.normalized_dir)!r}."
        )
    if not bundle.structures_dir.is_dir():
        raise RunnerBundleValidationError(
            f"{_bundle_prefix(bundle)}: structures directory is missing."
        )


def _validate_manifest(bundle: RunnerNormalizedBundle) -> dict[str, Any]:
    try:
        manifest = json.loads(bundle.manifest_path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise RunnerBundleValidationError(
            f"{_bundle_prefix(bundle)}: manifest.json is not valid JSON: {exc}."
        ) from exc

    if not isinstance(manifest, dict):
        raise RunnerBundleValidationError(
            f"{_bundle_prefix(bundle)}: manifest.json must contain a JSON object."
        )

    manifest_runner = manifest.get("runner")
    if manifest_runner is not None and str(manifest_runner) != bundle.runner_name:
        raise RunnerBundleValidationError(
            f"{_bundle_prefix(bundle)}: manifest runner {manifest_runner!r} does not match "
            f"result runner {bundle.runner_name!r}."
        )

    manifest_capabilities = manifest.get("capabilities")
    if manifest_capabilities is not None:
        manifest_capability_set = {str(value) for value in manifest_capabilities}
        if manifest_capability_set != set(bundle.capabilities):
            raise RunnerBundleValidationError(
                f"{_bundle_prefix(bundle)}: manifest capabilities "
                f"{sorted(manifest_capability_set)!r} do not match result capabilities "
                f"{sorted(bundle.capabilities)!r}."
            )
    return manifest


def _validate_execution_provenance(
    bundle: RunnerNormalizedBundle,
    manifest: dict[str, Any],
    system_df: pd.DataFrame,
    chain_df: pd.DataFrame,
) -> None:
    backend = bundle.backend_identity
    seed = bundle.seed_provenance
    if (backend is None) != (seed is None):
        raise RunnerProvenanceError(
            f"{_bundle_prefix(bundle)}: backend and seed provenance must be "
            "declared together."
        )
    if backend is None or seed is None:
        return
    if manifest.get("backend") != backend_manifest_value(backend):
        raise RunnerProvenanceError(
            f"{_bundle_prefix(bundle)}: manifest backend provenance does not "
            "match the bundle."
        )
    if manifest.get("seed") != seed_manifest_value(seed):
        raise RunnerProvenanceError(
            f"{_bundle_prefix(bundle)}: manifest seed provenance does not match "
            "the bundle."
        )

    expected = {
        "runner_id": backend.runner_name,
        "backend_name": backend.backend_name,
        "runner_version": backend.version,
        "backend_version_status": backend.version_status.value,
        "effective_seed": seed.effective_seed,
    }
    for label, frame in (
        ("system_metrics.csv", system_df),
        ("chain_metrics.csv", chain_df),
    ):
        missing = sorted(set(RUNNER_PROVENANCE_COLUMNS) - set(frame.columns))
        if missing:
            raise RunnerProvenanceError(
                f"{_bundle_prefix(bundle)}: {label} is missing provenance columns "
                f"{missing!r}."
            )
        for column, value in expected.items():
            observed = frame[column]
            if value is None:
                matches = observed.isna()
            elif column == "effective_seed":
                numeric = pd.to_numeric(observed, errors="coerce")
                matches = numeric.notna() & (numeric == value) & (numeric % 1 == 0)
            else:
                matches = observed.astype(str) == str(value)
            if not bool(matches.all()):
                raise RunnerProvenanceError(
                    f"{_bundle_prefix(bundle)}: {label} {column!r} does not match "
                    "provenance."
                )

    for record in bundle.sample_records:
        for column, value in expected.items():
            if record.get(column) != value:
                raise RunnerProvenanceError(
                    f"{_bundle_prefix(bundle)}: sample record {column!r} does not "
                    "match provenance."
                )


def _validate_companion_artifacts(
    bundle: RunnerNormalizedBundle,
    manifest: dict[str, Any],
) -> None:
    manifest_entries = _normalize_manifest_artifacts(bundle, manifest)
    manifest_by_label = _artifact_entries_by_label(
        bundle=bundle,
        artifacts=manifest_entries,
        source_label="manifest companion_artifacts",
    )
    bundle_by_label = _artifact_entries_by_label(
        bundle=bundle,
        artifacts=bundle.companion_artifacts,
        source_label="bundle companion_artifacts",
    )
    if manifest_by_label != bundle_by_label:
        raise RunnerBundleValidationError(
            f"{_bundle_prefix(bundle)}: manifest companion_artifacts do not match the "
            "companion_artifacts declared on the normalized bundle."
        )

    for artifact in bundle.companion_artifacts:
        artifact_path = (bundle.normalized_dir / artifact.relative_path).resolve()
        normalized_root = bundle.normalized_dir.resolve()
        if not _path_is_within(artifact_path, normalized_root):
            raise RunnerBundleValidationError(
                f"{_bundle_prefix(bundle)}: companion artifact {artifact.label!r} must stay "
                f"inside {str(bundle.normalized_dir)!r}."
            )
        if not artifact_path.exists():
            raise RunnerBundleValidationError(
                f"{_bundle_prefix(bundle)}: companion artifact {artifact.label!r} points to "
                f"missing path {artifact.relative_path!r}."
            )
        if not artifact_path.is_file() and not artifact_path.is_dir():
            raise RunnerBundleValidationError(
                f"{_bundle_prefix(bundle)}: companion artifact {artifact.label!r} must point "
                "to a file or directory."
            )


def _normalize_manifest_artifacts(
    bundle: RunnerNormalizedBundle,
    manifest: dict[str, Any],
) -> list[RunnerCompanionArtifact]:
    manifest_artifacts = manifest.get("companion_artifacts", [])
    if manifest_artifacts is None:
        return []
    if not isinstance(manifest_artifacts, list):
        raise RunnerBundleValidationError(
            f"{_bundle_prefix(bundle)}: manifest companion_artifacts must be a JSON array."
        )

    normalized_entries: list[RunnerCompanionArtifact] = []
    for entry in manifest_artifacts:
        if not isinstance(entry, dict):
            raise RunnerBundleValidationError(
                f"{_bundle_prefix(bundle)}: each companion_artifacts entry must be a JSON object."
            )
        label = entry.get("label")
        relative_path = entry.get("relative_path")
        if not isinstance(label, str) or not label.strip():
            raise RunnerBundleValidationError(
                f"{_bundle_prefix(bundle)}: companion_artifacts entries require a non-empty "
                "'label' string."
            )
        if not isinstance(relative_path, str) or not relative_path.strip():
            raise RunnerBundleValidationError(
                f"{_bundle_prefix(bundle)}: companion_artifacts entries require a non-empty "
                "'relative_path' string."
            )
        normalized_entries.append(
            RunnerCompanionArtifact(
                label=label.strip(),
                relative_path=relative_path.strip(),
                kind=_optional_manifest_text(bundle, entry.get("kind")),
                description=_optional_manifest_text(bundle, entry.get("description")),
            )
        )
    return normalized_entries


def _artifact_entries_by_label(
    *,
    bundle: RunnerNormalizedBundle,
    artifacts: Collection[RunnerCompanionArtifact],
    source_label: str,
) -> dict[str, dict[str, str]]:
    artifact_entries: dict[str, dict[str, str]] = {}
    for artifact in artifacts:
        if artifact.label in artifact_entries:
            raise RunnerBundleValidationError(
                f"{_bundle_prefix(bundle)}: {source_label} labels must be unique; duplicate "
                f"label {artifact.label!r} found."
            )
        artifact_entries[artifact.label] = _artifact_to_manifest_entry(artifact)
    return artifact_entries


def _artifact_to_manifest_entry(artifact: RunnerCompanionArtifact) -> dict[str, str]:
    entry = {
        "label": artifact.label,
        "relative_path": artifact.relative_path,
    }
    if artifact.kind is not None:
        entry["kind"] = artifact.kind
    if artifact.description is not None:
        entry["description"] = artifact.description
    return entry


def _optional_manifest_text(
    bundle: RunnerNormalizedBundle,
    value: Any,
) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str):
        raise RunnerBundleValidationError(
            f"{_bundle_prefix(bundle)}: manifest companion_artifacts optional fields "
            "must be strings when provided."
        )
    stripped = value.strip()
    return stripped or None


def _path_is_within(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
    except ValueError:
        return False
    return True


def _read_csv(
    bundle: RunnerNormalizedBundle,
    path: Path,
    label: str,
) -> pd.DataFrame:
    try:
        return pd.read_csv(path)
    except Exception as exc:
        raise RunnerBundleValidationError(
            f"{_bundle_prefix(bundle)}: could not read {label}: {exc}."
        ) from exc


def _validate_required_columns(
    bundle: RunnerNormalizedBundle,
    system_df: pd.DataFrame,
    chain_df: pd.DataFrame,
) -> None:
    missing_system = [
        column for column in CANONICAL_SYSTEM_COLUMNS if column not in system_df.columns
    ]
    if missing_system:
        raise RunnerBundleValidationError(
            f"{_bundle_prefix(bundle)}: system_metrics.csv is missing required columns "
            f"{missing_system!r}."
        )

    missing_chain = [
        column for column in CANONICAL_CHAIN_COLUMNS if column not in chain_df.columns
    ]
    if missing_chain:
        raise RunnerBundleValidationError(
            f"{_bundle_prefix(bundle)}: chain_metrics.csv is missing required columns "
            f"{missing_chain!r}."
        )


def _validate_structures(bundle: RunnerNormalizedBundle) -> set[str]:
    structure_files = {
        path.name
        for path in bundle.structures_dir.iterdir()
        if path.is_file() and path.suffix.lower() in STRUCTURE_SUFFIXES
    }
    if not structure_files:
        raise RunnerBundleValidationError(
            f"{_bundle_prefix(bundle)}: structures directory does not contain any canonical "
            "structure files."
        )
    return structure_files


def _validate_referenced_structures(
    bundle: RunnerNormalizedBundle,
    system_df: pd.DataFrame,
    chain_df: pd.DataFrame,
    structure_files: set[str],
) -> None:
    referenced_files = {
        str(value)
        for frame in (system_df, chain_df)
        for value in frame.get("cif_file", pd.Series(dtype=object)).dropna().tolist()
    }
    missing_files = sorted(referenced_files - structure_files)
    if missing_files:
        raise RunnerBundleValidationError(
            f"{_bundle_prefix(bundle)}: normalized CSV rows reference missing structure files "
            f"{missing_files!r}."
        )


def _validate_sample_records(
    bundle: RunnerNormalizedBundle,
    structure_files: set[str],
) -> None:
    for record in bundle.sample_records:
        if (
            "repeat" not in record
            or "diffusion_sample" not in record
            or "cif_file" not in record
        ):
            raise RunnerBundleValidationError(
                f"{_bundle_prefix(bundle)}: sample_records entries must include "
                "'repeat', 'diffusion_sample', and 'cif_file'."
            )
        cif_file = str(record["cif_file"])
        if cif_file not in structure_files:
            raise RunnerBundleValidationError(
                f"{_bundle_prefix(bundle)}: sample_records references missing structure "
                f"{cif_file!r}."
            )


def _validate_known_metric_groups(
    bundle: RunnerNormalizedBundle,
    explicit_groups: Collection[str],
    requested_groups: Collection[str],
) -> None:
    unknown_groups = sorted(
        {
            str(group_name)
            for group_name in set(explicit_groups) | set(requested_groups)
            if str(group_name) not in DEFAULT_RUNNER_METRIC_REQUIREMENTS
        }
    )
    if unknown_groups:
        raise RunnerBundleValidationError(
            f"{_bundle_prefix(bundle)}: unknown runner metric groups {unknown_groups!r}."
        )


def _metric_requirement(group_name: str) -> dict[str, Any]:
    return DEFAULT_RUNNER_METRIC_REQUIREMENTS[group_name]


def _required_metric_columns(
    group_name: str,
) -> list[str]:
    requirement = _metric_requirement(group_name)
    required_columns: list[str] = []
    for file_label in METRIC_PAYLOAD_SECTIONS:
        columns = requirement.get(file_label, {})
        required_columns.extend(f"{file_label}.{column}" for column in columns.keys())
    return required_columns


def _missing_metric_columns(
    group_name: str,
    system_df: pd.DataFrame,
    chain_df: pd.DataFrame,
) -> list[str]:
    requirement = _metric_requirement(group_name)
    missing_columns: list[str] = []
    missing_columns.extend(
        _missing_metric_columns_for_frame(
            file_label="system_metrics",
            columns=requirement.get("system_metrics", {}),
            frame=system_df,
        )
    )
    missing_columns.extend(
        _missing_metric_columns_for_frame(
            file_label="chain_metrics",
            columns=requirement.get("chain_metrics", {}),
            frame=chain_df,
        )
    )
    return missing_columns


def _missing_metric_columns_for_frame(
    *,
    file_label: str,
    columns: MetricColumnAliases,
    frame: pd.DataFrame,
) -> list[str]:
    missing_columns: list[str] = []
    for canonical_name, accepted_columns in columns.items():
        if not any(
            _matches_metric_column(frame.columns, candidate)
            for candidate in accepted_columns
        ):
            missing_columns.append(f"{file_label}.{canonical_name}")
    return missing_columns


def _present_metric_columns(
    group_name: str,
    system_df: pd.DataFrame,
    chain_df: pd.DataFrame,
) -> tuple[str, ...]:
    requirement = _metric_requirement(group_name)
    present_columns = _present_metric_columns_for_frame(
        file_label="system_metrics",
        columns=requirement.get("system_metrics", {}),
        frame=system_df,
    )
    present_columns.extend(
        _present_metric_columns_for_frame(
            file_label="chain_metrics",
            columns=requirement.get("chain_metrics", {}),
            frame=chain_df,
        )
    )
    return tuple(present_columns)


def _present_metric_columns_for_frame(
    *,
    file_label: str,
    columns: MetricColumnAliases,
    frame: pd.DataFrame,
) -> list[str]:
    present_columns: list[str] = []
    for accepted_columns in columns.values():
        for column_name in accepted_columns:
            present_columns.extend(
                f"{file_label}.{matched_name}"
                for matched_name in _matched_metric_columns(frame.columns, column_name)
            )
    return present_columns


def _required_metric_artifacts(
    group_name: str,
) -> tuple[str, ...]:
    requirement = _metric_requirement(group_name)
    return tuple(
        f"companion_artifacts.{label}"
        for label in requirement.get("companion_artifacts", ())
    )


def _missing_metric_artifacts(
    bundle: RunnerNormalizedBundle,
    group_name: str,
) -> list[str]:
    required_labels = set(
        _metric_requirement(group_name).get("companion_artifacts", ())
    )
    present_labels = {artifact.label for artifact in bundle.companion_artifacts}
    return sorted(
        f"companion_artifacts.{label}" for label in required_labels - present_labels
    )


def _present_metric_artifacts(
    bundle: RunnerNormalizedBundle,
    group_name: str,
) -> tuple[str, ...]:
    required_labels = _metric_requirement(group_name).get("companion_artifacts", ())
    present_labels = {artifact.label for artifact in bundle.companion_artifacts}
    return tuple(
        f"companion_artifacts.{label}"
        for label in required_labels
        if label in present_labels
    )


def _uses_runner_authored_payload(group_name: str) -> bool:
    return group_name == "confidence_metrics"


def _matches_metric_column(
    available_columns: Collection[str],
    candidate: str,
) -> bool:
    if candidate.startswith("prefix:"):
        prefix = candidate.removeprefix("prefix:")
        return any(column_name.startswith(prefix) for column_name in available_columns)
    return candidate in available_columns


def _matched_metric_columns(
    available_columns: Collection[str],
    candidate: str,
) -> tuple[str, ...]:
    if candidate.startswith("prefix:"):
        prefix = candidate.removeprefix("prefix:")
        return tuple(
            column_name
            for column_name in available_columns
            if column_name.startswith(prefix)
        )
    if candidate in available_columns:
        return (candidate,)
    return ()


def _validate_metric_outcome(
    *,
    bundle: RunnerNormalizedBundle,
    group_name: str,
    outcome: RunnerMetricOutcome,
    system_df: pd.DataFrame,
    chain_df: pd.DataFrame,
) -> RunnerMetricOutcome:
    if _uses_runner_authored_payload(group_name):
        return _validate_runner_authored_metric_outcome(
            bundle=bundle,
            group_name=group_name,
            outcome=outcome,
            system_df=system_df,
            chain_df=chain_df,
        )

    missing_columns = _missing_metric_columns(group_name, system_df, chain_df)
    present_columns = _present_metric_columns(group_name, system_df, chain_df)
    missing_artifacts = _missing_metric_artifacts(bundle, group_name)
    present_artifacts = _present_metric_artifacts(bundle, group_name)

    if outcome.state == "unsupported":
        if group_name in bundle.capabilities:
            raise RunnerBundleValidationError(
                f"{_bundle_prefix(bundle)}: metric group {group_name!r} cannot be marked "
                "'unsupported' because the runner declares that capability."
            )
        if present_columns or present_artifacts:
            raise RunnerBundleValidationError(
                f"{_bundle_prefix(bundle)}: metric group {group_name!r} cannot be marked "
                f"'unsupported' while normalized payload is present: "
                f"{tuple(present_columns) + tuple(present_artifacts)!r}."
            )
        return RunnerMetricOutcome(state="unsupported")

    if group_name not in bundle.capabilities:
        raise RunnerBundleValidationError(
            f"{_bundle_prefix(bundle)}: metric group {group_name!r} cannot be marked "
            f"{outcome.state!r} because the runner does not declare that capability."
        )

    if outcome.state == "computed":
        if missing_columns or missing_artifacts:
            raise RunnerBundleValidationError(
                f"{_bundle_prefix(bundle)}: metric group {group_name!r} is marked 'computed' "
                "but required normalized payload entries are missing: "
                f"{missing_columns + missing_artifacts!r}."
            )
        return RunnerMetricOutcome(
            state="computed",
            required_columns=tuple(_required_metric_columns(group_name)),
            required_artifacts=_required_metric_artifacts(group_name),
        )

    if outcome.state == "missing":
        if not missing_columns and not missing_artifacts:
            raise RunnerBundleValidationError(
                f"{_bundle_prefix(bundle)}: metric group {group_name!r} is marked 'missing' "
                "but all required normalized payload entries are present."
            )
        return RunnerMetricOutcome(
            state="missing",
            required_columns=tuple(missing_columns),
            required_artifacts=tuple(missing_artifacts),
            message=outcome.message,
        )

    if outcome.state == "failed":
        if not outcome.message:
            raise RunnerBundleValidationError(
                f"{_bundle_prefix(bundle)}: metric group {group_name!r} is marked 'failed' "
                "but does not include an actionable message."
            )
        if not missing_columns and not missing_artifacts:
            raise RunnerBundleValidationError(
                f"{_bundle_prefix(bundle)}: metric group {group_name!r} is marked 'failed' "
                "even though the required normalized payload is complete."
            )
        return RunnerMetricOutcome(
            state="failed",
            required_columns=tuple(missing_columns),
            required_artifacts=tuple(missing_artifacts),
            message=outcome.message,
        )

    raise RunnerBundleValidationError(
        f"{_bundle_prefix(bundle)}: metric group {group_name!r} uses unknown state "
        f"{outcome.state!r}."
    )


def _validate_runner_authored_metric_outcome(
    *,
    bundle: RunnerNormalizedBundle,
    group_name: str,
    outcome: RunnerMetricOutcome,
    system_df: pd.DataFrame,
    chain_df: pd.DataFrame,
) -> RunnerMetricOutcome:
    present_columns, missing_columns = _runner_authored_payload_columns(
        outcome=outcome,
        system_df=system_df,
        chain_df=chain_df,
    )
    present_artifacts, missing_artifacts = _runner_authored_payload_artifacts(
        outcome=outcome,
        bundle=bundle,
    )

    if outcome.state == "unsupported":
        if group_name in bundle.capabilities:
            raise RunnerBundleValidationError(
                f"{_bundle_prefix(bundle)}: metric group {group_name!r} cannot be marked "
                "'unsupported' because the runner declares that capability."
            )
        if (
            outcome.required_columns
            or outcome.required_artifacts
            or _confidence_payload_present(bundle, system_df, chain_df)
        ):
            raise RunnerBundleValidationError(
                f"{_bundle_prefix(bundle)}: metric group {group_name!r} cannot be marked "
                "'unsupported' while confidence payload is present or declared."
            )
        return RunnerMetricOutcome(state="unsupported")

    if group_name not in bundle.capabilities:
        raise RunnerBundleValidationError(
            f"{_bundle_prefix(bundle)}: metric group {group_name!r} cannot be marked "
            f"{outcome.state!r} because the runner does not declare that capability."
        )

    if outcome.state == "computed":
        if missing_columns or missing_artifacts:
            raise RunnerBundleValidationError(
                f"{_bundle_prefix(bundle)}: metric group {group_name!r} is marked 'computed' "
                "but runner-declared normalized payload entries are missing: "
                f"{missing_columns + missing_artifacts!r}."
            )
        return RunnerMetricOutcome(
            state="computed",
            required_columns=tuple(outcome.required_columns),
            required_artifacts=tuple(outcome.required_artifacts),
            message=outcome.message,
        )

    if outcome.state == "missing":
        if present_columns or present_artifacts:
            raise RunnerBundleValidationError(
                f"{_bundle_prefix(bundle)}: metric group {group_name!r} is marked 'missing' "
                "but runner-declared missing payload entries are present: "
                f"{present_columns + present_artifacts!r}."
            )
        return RunnerMetricOutcome(
            state="missing",
            required_columns=tuple(outcome.required_columns),
            required_artifacts=tuple(outcome.required_artifacts),
            message=outcome.message,
        )

    if outcome.state == "failed":
        if not outcome.message:
            raise RunnerBundleValidationError(
                f"{_bundle_prefix(bundle)}: metric group {group_name!r} is marked 'failed' "
                "but does not include an actionable message."
            )
        if present_columns or present_artifacts:
            raise RunnerBundleValidationError(
                f"{_bundle_prefix(bundle)}: metric group {group_name!r} is marked 'failed' "
                "but runner-declared missing payload entries are present: "
                f"{present_columns + present_artifacts!r}."
            )
        return RunnerMetricOutcome(
            state="failed",
            required_columns=tuple(outcome.required_columns),
            required_artifacts=tuple(outcome.required_artifacts),
            message=outcome.message,
        )

    raise RunnerBundleValidationError(
        f"{_bundle_prefix(bundle)}: metric group {group_name!r} uses unknown state "
        f"{outcome.state!r}."
    )


def _runner_authored_payload_columns(
    *,
    outcome: RunnerMetricOutcome,
    system_df: pd.DataFrame,
    chain_df: pd.DataFrame,
) -> tuple[list[str], list[str]]:
    present_columns: list[str] = []
    missing_columns: list[str] = []
    frame_lookup = {
        "system_metrics": system_df,
        "chain_metrics": chain_df,
    }
    for column_ref in outcome.required_columns:
        section, column_name = _split_runner_authored_column_ref(column_ref)
        frame = frame_lookup[section]
        if _matches_metric_column(frame.columns, column_name):
            present_columns.append(column_ref)
        else:
            missing_columns.append(column_ref)
    return present_columns, missing_columns


def _split_runner_authored_column_ref(column_ref: str) -> tuple[str, str]:
    section, separator, column_name = str(column_ref).partition(".")
    if separator != "." or section not in METRIC_PAYLOAD_SECTIONS or not column_name:
        raise RunnerBundleValidationError(
            "Runner-authored metric outcome columns must use 'system_metrics.<column>' or "
            f"'chain_metrics.<column>' references, got {column_ref!r}."
        )
    return section, column_name


def _runner_authored_payload_artifacts(
    *,
    outcome: RunnerMetricOutcome,
    bundle: RunnerNormalizedBundle,
) -> tuple[list[str], list[str]]:
    present_labels = {artifact.label for artifact in bundle.companion_artifacts}
    present_artifacts: list[str] = []
    missing_artifacts: list[str] = []
    for artifact_ref in outcome.required_artifacts:
        label = _split_runner_authored_artifact_ref(artifact_ref)
        if label in present_labels:
            present_artifacts.append(artifact_ref)
        else:
            missing_artifacts.append(artifact_ref)
    return present_artifacts, missing_artifacts


def _split_runner_authored_artifact_ref(artifact_ref: str) -> str:
    section, separator, label = str(artifact_ref).partition(".")
    if separator != "." or section != "companion_artifacts" or not label:
        raise RunnerBundleValidationError(
            "Runner-authored metric outcome artifacts must use "
            f"'companion_artifacts.<label>' references, got {artifact_ref!r}."
        )
    return label


def _confidence_payload_present(
    bundle: RunnerNormalizedBundle,
    system_df: pd.DataFrame,
    chain_df: pd.DataFrame,
) -> bool:
    if any(
        _matches_metric_column(system_df.columns, column)
        for column in CONFIDENCE_SYSTEM_MARKERS
    ):
        return True
    if any(
        _matches_metric_column(chain_df.columns, column)
        for column in CONFIDENCE_CHAIN_MARKERS
    ):
        return True
    artifact_labels = {artifact.label for artifact in bundle.companion_artifacts}
    return any(label in artifact_labels for label in CONFIDENCE_ARTIFACT_LABELS)


def _raise_for_non_continuable_outcomes(
    bundle: RunnerNormalizedBundle,
    *,
    requested_groups: Collection[str],
) -> None:
    blocking_messages: list[str] = []
    for group_name in requested_groups:
        outcome = bundle.metric_outcomes.get(group_name)
        if outcome is None:
            continue
        if outcome.state == "missing":
            blocking_messages.append(
                f"{group_name}: missing required normalized payload "
                f"{list(outcome.required_columns) + list(outcome.required_artifacts)!r}"
            )
        elif outcome.state == "failed":
            blocking_messages.append(
                f"{group_name}: {outcome.message or 'runner reported failed outcome'}"
            )

    if blocking_messages:
        raise RunnerBundleValidationError(
            f"Runner {bundle.runner_name!r} stopped after execution because the normalized bundle "
            f"cannot satisfy requested metric groups: {'; '.join(blocking_messages)}. "
            "Partial runner outputs remain in the raw repeat directory for debugging."
        )
