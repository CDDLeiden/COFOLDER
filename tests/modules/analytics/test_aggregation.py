"""Canonical analytics aggregation boundary and compatibility coverage."""

import logging
import math

import pandas as pd

from cofolder.modules.analytics import aggregation
from cofolder.modules.utils import gather


AGGREGATION_FUNCTIONS = (
    "gather_structures",
    "merge_runner_results",
    "add_chain_info",
    "assess_numeric_variance",
    "assess_bitstring_similarity",
    "gather_robustness_results",
    "is_metadata_column",
    "is_numeric_metric",
    "is_ifp_metric",
)


def test_legacy_gather_facade_aliases_all_canonical_functions():
    for name in AGGREGATION_FUNCTIONS:
        assert getattr(gather, name) is getattr(aggregation, name)

    assert gather.align is aggregation.align


def test_canonical_aggregation_classifies_and_summarizes_metrics():
    assert aggregation.assess_numeric_variance([1.0, 3.0], "score") == {
        "score_mean": 2.0,
        "score_std": math.sqrt(2.0),
    }
    assert aggregation.is_metadata_column("repeat") is True
    assert aggregation.is_metadata_column("confidence_score") is False
    assert bool(aggregation.is_numeric_metric(pd.Series(["1.0", "", None]))) is True
    assert bool(aggregation.is_numeric_metric(pd.Series(["value"]))) is False
    assert aggregation.is_ifp_metric(
        "ifp_distance", pd.Series([None, "101"])
    ) is True
    assert aggregation.is_ifp_metric("ifp_distance", pd.Series([None])) is False


def test_canonical_gather_structures_preserves_structure_inventory(temp_dir, caplog):
    caplog.set_level(logging.INFO)
    source = temp_dir / "raw" / "repeat_1" / "normalized" / "structures"
    source.mkdir(parents=True)
    (source / "model.cif").write_text("data_model\n", encoding="utf-8")
    (source / "ignored.txt").write_text("ignored\n", encoding="utf-8")

    aggregation.gather_structures(temp_dir, "system", repeats=1)

    target = temp_dir / "results" / "structures"
    assert [path.name for path in target.iterdir()] == ["model.cif"]
    assert {record.name for record in caplog.records} == {
        "cofolder.modules.utils.gather"
    }
