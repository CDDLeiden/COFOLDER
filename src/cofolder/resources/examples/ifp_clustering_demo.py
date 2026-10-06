#!/usr/bin/env python3
"""Run a small, deterministic distance-IFP clustering example."""

from __future__ import annotations

import json

import pandas as pd

from cofolder.modules.analytics.ifp_clustering import cluster_binary_ifps


def main() -> None:
    compounds = pd.DataFrame(
        {
            "compound_id": [
                "hinge-a",
                "hinge-b",
                "back-pocket-a",
                "back-pocket-b",
                "no-contacts-a",
                "no-contacts-b",
            ],
            "ifp_distance": [
                [1, 1, 0, 0],
                [1, 0, 0, 0],
                [0, 0, 1, 1],
                [0, 0, 1, 1],
                [0, 0, 0, 0],
                [0, 0, 0, 0],
            ],
        }
    )

    result = cluster_binary_ifps(
        compounds["ifp_distance"].tolist(),
        compounds["compound_id"].tolist(),
        similarity_threshold=0.5,
    )
    compounds["ifp_cluster_id"] = result.cluster_ids

    print("Row assignments")
    print(compounds.to_string(index=False))
    print("\nCluster summary")
    print(result.summary.to_string(index=False))
    print("\nNative linkage matrix")
    print(result.linkage_matrix)
    print("\nDeterministic leaf labels")
    print(result.leaf_member_ids)

    print("\nDecoded JSON fields")
    for row in result.summary.itertuples(index=False):
        print(
            row.ifp_cluster_id,
            "members=",
            json.loads(row.member_ids),
            "consensus=",
            json.loads(row.consensus_ifp),
        )


if __name__ == "__main__":
    main()
