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
    assert result.member_ids == ("A", "B", "C", "D", "E", "F")
    assert result.linkage_matrix == (
        (2.0, 3.0, 0.0, 2.0),
        (4.0, 5.0, 0.0, 2.0),
        (0.0, 1.0, 0.5, 2.0),
        (6.0, 8.0, 1.0, 4.0),
        (7.0, 9.0, 1.0, 6.0),
    )
    assert result.leaf_indices == (4, 5, 2, 3, 0, 1)
    assert result.leaf_member_ids == ("E", "F", "C", "D", "A", "B")
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

    inclusive = cluster_binary_ifps(fingerprints, member_ids, similarity_threshold=0.5)
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
    assert result.linkage_matrix == ()
    assert result.leaf_indices == (0,)
    assert result.leaf_member_ids == ("empty",)


def test_empty_input_and_ambiguous_member_ids_are_explicit():
    empty = cluster_binary_ifps([], [], similarity_threshold=0.5)
    assert empty.member_ids == ()
    assert empty.linkage_matrix == ()
    assert empty.leaf_indices == ()
    assert empty.leaf_member_ids == ()

    with pytest.raises(ValueError, match="unique"):
        cluster_binary_ifps([[1], [0]], ["same", "same"], similarity_threshold=0.5)
    with pytest.raises(ValueError, match="non-empty"):
        cluster_binary_ifps([[1]], ["  "], similarity_threshold=0.5)
