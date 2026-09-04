"""Tests for deterministic distance-IFP clustering."""

import json

import pytest

from cofolder.modules.analytics.ifp_clustering import cluster_binary_ifps


def test_fixed_matrix_has_stable_assignments_and_summary():
    result = cluster_binary_ifps(
        [
            [1, 1, 0, 0],
            [1, 0, 0, 0],
            [0, 0, 1, 1],
            [0, 0, 1, 1],
            [0, 0, 0, 0],
            [0, 0, 0, 0],
        ],
        ["A", "B", "C", "D", "E", "F"],
        similarity_threshold=0.5,
    )

    assert result.cluster_ids == [
        "IFP001",
        "IFP001",
        "IFP002",
        "IFP002",
        "IFP003",
        "IFP003",
    ]
    summary = result.summary.set_index("ifp_cluster_id")
    assert summary["size"].to_dict() == {"IFP001": 2, "IFP002": 2, "IFP003": 2}
    assert json.loads(summary.at["IFP001", "member_ids"]) == ["A", "B"]
    assert summary.at["IFP001", "medoid_compound_id"] == "A"
    assert json.loads(summary.at["IFP001", "consensus_ifp"]) == [1, 1, 0, 0]
    assert summary.at[
        "IFP001", "mean_within_cluster_jaccard_similarity"
    ] == pytest.approx(0.5)
    assert summary.at[
        "IFP003", "mean_within_cluster_jaccard_similarity"
    ] == pytest.approx(1.0)


def test_similarity_threshold_changes_assignments_at_inclusive_boundary():
    fingerprints = [[1, 1, 0], [1, 0, 0], [0, 0, 1]]
    member_ids = ["first", "second", "third"]

    inclusive = cluster_binary_ifps(
        fingerprints, member_ids, similarity_threshold=0.5
    )
    strict = cluster_binary_ifps(
        fingerprints, member_ids, similarity_threshold=0.500001
    )

    assert inclusive.cluster_ids == ["IFP001", "IFP001", "IFP002"]
    assert strict.cluster_ids == ["IFP001", "IFP002", "IFP003"]


def test_single_fingerprint_is_a_complete_cluster():
    result = cluster_binary_ifps([[0, 0]], ["empty"], similarity_threshold=1.0)

    assert result.cluster_ids == ["IFP001"]
    assert result.summary.loc[0, "medoid_compound_id"] == "empty"
    assert result.summary.loc[0, "mean_within_cluster_jaccard_similarity"] == 1.0
