"""Workflow filtering policy for normalized interaction fingerprints."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from typing import Literal

from .reference_ifp import (
    IFPComparison,
    IFPComparisonStatus,
    IFPSimilarityMetric,
    InteractionFingerprint,
    InteractionKey,
    ReferenceIFPInputError,
)


class IFPFilterMode(StrEnum):
    SIMILARITY = "similarity"
    REQUIRED = "required"


@dataclass(frozen=True, slots=True)
class ReferenceIFPFilterPolicy:
    mode: IFPFilterMode
    similarity_metric: IFPSimilarityMetric = IFPSimilarityMetric.JACCARD
    threshold: float | None = None
    required_interactions: frozenset[InteractionKey] | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "mode", IFPFilterMode(self.mode))
        object.__setattr__(self, "similarity_metric", IFPSimilarityMetric(self.similarity_metric))
        if self.mode is IFPFilterMode.SIMILARITY:
            if self.threshold is None or not 0 <= float(self.threshold) <= 1:
                raise ReferenceIFPInputError("Similarity filtering requires a threshold in [0, 1].")
        elif self.threshold is not None:
            raise ReferenceIFPInputError("Required-interaction filtering does not accept a threshold.")


@dataclass(frozen=True, slots=True)
class ReferenceIFPFilterOutcome:
    passed: bool | None
    status: Literal["accepted", "rejected", "not_evaluable", "not_applied"]
    reason: str
    comparison: IFPComparison | None


def evaluate_reference_ifp_filter(
    reference: InteractionFingerprint,
    comparison: IFPComparison,
    policy: ReferenceIFPFilterPolicy,
) -> ReferenceIFPFilterOutcome:
    if comparison.status is not IFPComparisonStatus.COMPARABLE:
        reason = comparison.mapping.failures[0] if comparison.mapping.failures else "empty_reference"
        return ReferenceIFPFilterOutcome(None, "not_evaluable", reason, comparison)

    if policy.mode is IFPFilterMode.SIMILARITY:
        score = comparison.similarities[policy.similarity_metric]
        passed = score >= float(policy.threshold)
        return ReferenceIFPFilterOutcome(
            passed, "accepted" if passed else "rejected",
            "threshold_met" if passed else "below_threshold", comparison,
        )

    required = policy.required_interactions or reference.interactions
    unknown = required - reference.interactions
    if unknown:
        raise ReferenceIFPInputError(
            "Required interactions are absent from the reference fingerprint: "
            + ", ".join(str(item) for item in sorted(unknown))
        )
    missing = required & set(comparison.missing_interactions)
    passed = not missing
    return ReferenceIFPFilterOutcome(
        passed, "accepted" if passed else "rejected",
        "required_interactions_met" if passed else "required_interactions_missing",
        comparison,
    )

