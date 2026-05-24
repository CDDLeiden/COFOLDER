from __future__ import annotations

import json
from pathlib import Path
from typing import Collection

import pandas as pd

from cofolder.modules.runners.contracts import (
    RunnerExecutionResult,
    RunnerMetricOutcome,
    RunnerNormalizedBundle,
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
RUNNER_METRIC_REQUIREMENTS = {
    "confidence_metrics": {
        "system_metrics": ("ptm", "iptm", "confidence_score"),
        "chain_metrics": ("chains_ptm",),
    },
    "affinity_metrics": {
        "chain_metrics": ("affinity_pred_value", "affinity_probability_binary"),
    },
    "affinity_metrics_ext": {
        "chain_metrics": ("pIC50", "IC50_M", "pIC50_kcal_per_mol"),
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
    _validate_manifest(bundle)

    system_df = _read_csv(bundle, bundle.system_metrics_path, "system_metrics.csv")
    chain_df = _read_csv(bundle, bundle.chain_metrics_path, "chain_metrics.csv")
    _validate_required_columns(bundle, system_df, chain_df)
    structure_files = _validate_structures(bundle)
    _validate_referenced_structures(bundle, system_df, chain_df, structure_files)
    _validate_sample_records(bundle, structure_files)

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


def _validate_manifest(bundle: RunnerNormalizedBundle) -> None:
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
    missing_system = [column for column in CANONICAL_SYSTEM_COLUMNS if column not in system_df.columns]
    if missing_system:
        raise RunnerBundleValidationError(
            f"{_bundle_prefix(bundle)}: system_metrics.csv is missing required columns "
            f"{missing_system!r}."
        )

    missing_chain = [column for column in CANONICAL_CHAIN_COLUMNS if column not in chain_df.columns]
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
        if "repeat" not in record or "diffusion_sample" not in record or "cif_file" not in record:
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
            if str(group_name) not in RUNNER_METRIC_REQUIREMENTS
        }
    )
    if unknown_groups:
        raise RunnerBundleValidationError(
            f"{_bundle_prefix(bundle)}: unknown runner metric groups {unknown_groups!r}."
        )


def _required_metric_columns(group_name: str) -> list[str]:
    requirement = RUNNER_METRIC_REQUIREMENTS[group_name]
    required_columns: list[str] = []
    for file_label, columns in requirement.items():
        required_columns.extend(f"{file_label}.{column}" for column in columns)
    return required_columns


def _missing_metric_columns(
    group_name: str,
    system_df: pd.DataFrame,
    chain_df: pd.DataFrame,
) -> list[str]:
    requirement = RUNNER_METRIC_REQUIREMENTS[group_name]
    missing_columns: list[str] = []
    for column in requirement.get("system_metrics", ()):
        if column not in system_df.columns:
            missing_columns.append(f"system_metrics.{column}")
    for column in requirement.get("chain_metrics", ()):
        if column not in chain_df.columns:
            missing_columns.append(f"chain_metrics.{column}")
    return missing_columns


def _validate_metric_outcome(
    *,
    bundle: RunnerNormalizedBundle,
    group_name: str,
    outcome: RunnerMetricOutcome,
    system_df: pd.DataFrame,
    chain_df: pd.DataFrame,
) -> RunnerMetricOutcome:
    missing_columns = _missing_metric_columns(group_name, system_df, chain_df)
    present_columns = tuple(
        column
        for column in _required_metric_columns(group_name)
        if column not in missing_columns
    )

    if outcome.state == "unsupported":
        if group_name in bundle.capabilities:
            raise RunnerBundleValidationError(
                f"{_bundle_prefix(bundle)}: metric group {group_name!r} cannot be marked "
                "'unsupported' because the runner declares that capability."
            )
        if present_columns:
            raise RunnerBundleValidationError(
                f"{_bundle_prefix(bundle)}: metric group {group_name!r} cannot be marked "
                f"'unsupported' while normalized payload columns are present: {present_columns!r}."
            )
        return RunnerMetricOutcome(state="unsupported")

    if group_name not in bundle.capabilities:
        raise RunnerBundleValidationError(
            f"{_bundle_prefix(bundle)}: metric group {group_name!r} cannot be marked "
            f"{outcome.state!r} because the runner does not declare that capability."
        )

    if outcome.state == "computed":
        if missing_columns:
            raise RunnerBundleValidationError(
                f"{_bundle_prefix(bundle)}: metric group {group_name!r} is marked 'computed' "
                f"but required normalized columns are missing: {missing_columns!r}."
            )
        return RunnerMetricOutcome(
            state="computed",
            required_columns=tuple(_required_metric_columns(group_name)),
        )

    if outcome.state == "missing":
        if not missing_columns:
            raise RunnerBundleValidationError(
                f"{_bundle_prefix(bundle)}: metric group {group_name!r} is marked 'missing' "
                "but all required normalized columns are present."
            )
        return RunnerMetricOutcome(
            state="missing",
            required_columns=tuple(missing_columns),
            message=outcome.message,
        )

    if outcome.state == "failed":
        if not outcome.message:
            raise RunnerBundleValidationError(
                f"{_bundle_prefix(bundle)}: metric group {group_name!r} is marked 'failed' "
                "but does not include an actionable message."
            )
        if not missing_columns:
            raise RunnerBundleValidationError(
                f"{_bundle_prefix(bundle)}: metric group {group_name!r} is marked 'failed' "
                "even though the required normalized payload is complete."
            )
        return RunnerMetricOutcome(
            state="failed",
            required_columns=tuple(missing_columns),
            message=outcome.message,
        )

    raise RunnerBundleValidationError(
        f"{_bundle_prefix(bundle)}: metric group {group_name!r} uses unknown state "
        f"{outcome.state!r}."
    )


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
                f"{group_name}: missing required normalized columns {list(outcome.required_columns)!r}"
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
