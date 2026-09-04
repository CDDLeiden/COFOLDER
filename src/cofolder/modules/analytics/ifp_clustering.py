"""Deterministic clustering and summaries for binary interaction fingerprints."""

from __future__ import annotations

import json
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import cdist, pdist


@dataclass(frozen=True)
class IFPClusteringResult:
    """Stable cluster assignments and their audit-ready summary."""

    cluster_ids: list[str]
    summary: pd.DataFrame


SUMMARY_COLUMNS = (
    "ifp_cluster_id",
    "size",
    "member_ids",
    "medoid_compound_id",
    "consensus_ifp",
    "mean_within_cluster_jaccard_similarity",
)


def cluster_binary_ifps(
    fingerprints: Sequence[Sequence[int]],
    member_ids: Sequence[str],
    *,
    similarity_threshold: float,
) -> IFPClusteringResult:
    """Cluster compatible binary IFPs using average-linkage Jaccard distance.

    Cluster IDs are normalized by the earliest input position, so callers never
    observe the implementation-specific integer labels returned by SciPy.
    """

    if len(fingerprints) != len(member_ids):
        raise ValueError("fingerprints and member_ids must have equal length")
    if not fingerprints:
        return IFPClusteringResult([], pd.DataFrame(columns=SUMMARY_COLUMNS))

    width = len(fingerprints[0])
    if width == 0 or any(len(value) != width for value in fingerprints):
        raise ValueError("fingerprints must be non-empty and have equal lengths")
    features = np.asarray(fingerprints, dtype=np.bool_)
    if features.ndim != 2 or features.shape != (len(fingerprints), width):
        raise ValueError("fingerprints must form a two-dimensional binary matrix")

    if len(features) == 1:
        raw_labels = np.array([1], dtype=int)
    else:
        distances = pdist(features, metric="jaccard")
        # Explicitly define the empty-set Jaccard distance as zero. This keeps
        # all-zero fingerprints as one reproducible "no contacts" pattern.
        distances = np.nan_to_num(distances, nan=0.0)
        tree = linkage(distances, method="average", optimal_ordering=False)
        raw_labels = fcluster(
            tree,
            t=1.0 - float(similarity_threshold),
            criterion="distance",
        ).astype(int)

    ordered_raw_labels = sorted(
        np.unique(raw_labels),
        key=lambda label: int(np.flatnonzero(raw_labels == label)[0]),
    )
    stable_names = {
        label: f"IFP{position:03d}"
        for position, label in enumerate(ordered_raw_labels, start=1)
    }
    cluster_ids = [stable_names[label] for label in raw_labels]

    rows: list[dict[str, object]] = []
    member_ids_as_text = [str(value) for value in member_ids]
    for raw_label in ordered_raw_labels:
        positions = np.flatnonzero(raw_labels == raw_label)
        cluster_features = features[positions]
        cluster_member_ids = [member_ids_as_text[position] for position in positions]
        medoid_offset, mean_similarity = _cluster_medoid_and_mean_similarity(
            cluster_features
        )
        consensus = (cluster_features.mean(axis=0) >= 0.5).astype(int).tolist()
        rows.append(
            {
                "ifp_cluster_id": stable_names[raw_label],
                "size": len(positions),
                "member_ids": json.dumps(cluster_member_ids),
                "medoid_compound_id": cluster_member_ids[medoid_offset],
                "consensus_ifp": json.dumps(consensus),
                "mean_within_cluster_jaccard_similarity": mean_similarity,
            }
        )

    return IFPClusteringResult(
        cluster_ids=cluster_ids,
        summary=pd.DataFrame(rows, columns=SUMMARY_COLUMNS),
    )


def _cluster_medoid_and_mean_similarity(features: np.ndarray) -> tuple[int, float]:
    """Return earliest medoid on ties and mean similarity over unique pairs."""

    size = len(features)
    if size == 1:
        return 0, 1.0

    distance_sums = np.zeros(size, dtype=float)
    block_size = 256
    for start in range(0, size, block_size):
        stop = min(start + block_size, size)
        block = cdist(features[start:stop], features, metric="jaccard")
        block = np.nan_to_num(block, nan=0.0)
        distance_sums[start:stop] = block.sum(axis=1)

    # np.argmin returns the earliest member when medoid scores tie.
    medoid_offset = int(np.argmin(distance_sums))
    mean_distance = float(distance_sums.sum() / (size * (size - 1)))
    return medoid_offset, 1.0 - mean_distance
