"""Compatibility facade for the canonical analytics aggregation module.

This import path remains supported until its separately approved R09 removal.
"""

from cofolder.modules.analytics import aggregation

align = aggregation.align
gather_structures = aggregation.gather_structures
merge_runner_results = aggregation.merge_runner_results
add_chain_info = aggregation.add_chain_info
assess_numeric_variance = aggregation.assess_numeric_variance
assess_bitstring_similarity = aggregation.assess_bitstring_similarity
gather_robustness_results = aggregation.gather_robustness_results
is_metadata_column = aggregation.is_metadata_column
is_numeric_metric = aggregation.is_numeric_metric
is_ifp_metric = aggregation.is_ifp_metric

__all__ = [
    "add_chain_info",
    "align",
    "assess_bitstring_similarity",
    "assess_numeric_variance",
    "gather_robustness_results",
    "gather_structures",
    "is_ifp_metric",
    "is_metadata_column",
    "is_numeric_metric",
    "merge_runner_results",
]
