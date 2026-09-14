from __future__ import annotations

import json
from pathlib import Path

import yaml

from cofolder.tools import run_bias_matrix


def test_runs_the_four_required_cutoff_and_custom_scenarios(monkeypatch, temp_dir):
    calls = []

    class FakeBias:
        def __init__(self, **kwargs):
            self.kwargs = kwargs

        def run(self):
            calls.append(self.kwargs)

    monkeypatch.setattr("cofolder.recipes.bias.Bias", FakeBias)
    config = temp_dir / "matrix.yaml"
    config.write_text(yaml.safe_dump({
        "schema_version": 1,
        "workflow": "bias",
        "assessment_root": "assessment",
        "system_path": "system.yaml",
        "bias_training_data_protein_path": "protein",
        "bias_training_data_ligand_path": "ligand",
        "bias_query_cache_path": "cache",
        "custom_bias_reference_path": "custom",
        "release_cutoff": "2023-06-01",
        "bias_protein_similarity_threshold": 0.2,
        "bias_ligand_similarity_threshold": 0.4,
    }), encoding="utf-8")
    assert run_bias_matrix.main([str(config)]) == 0
    assert len(calls) == 4
    names = {Path(call["wrk_dir"]).name for call in calls}
    assert names == {
        "pre-2023-06-01",
        "whole-snapshot",
        "pre-2023-06-01-plus-custom",
        "custom-only",
    }
    by_name = {Path(call["wrk_dir"]).name: call for call in calls}
    assert by_name["pre-2023-06-01"]["bias_release_cutoff"] == "2023-06-01"
    assert by_name["whole-snapshot"]["bias_release_cutoff"] == "whole"
    assert by_name["pre-2023-06-01-plus-custom"]["custom_bias_reference_path"].endswith("custom")
    assert by_name["custom-only"]["bias_training_data_protein_path"] is None
    assert by_name["custom-only"]["bias_training_data_ligand_path"] is None
    assert by_name["custom-only"]["custom_bias_reference_path"].endswith("custom")
    assert {call["bias_protein_similarity_threshold"] for call in calls} == {0.2}
    assert {call["bias_ligand_similarity_threshold"] for call in calls} == {0.4}
    manifest = json.loads((temp_dir / "assessment/matrix_manifest.json").read_text())
    assert len(manifest["scenarios"]) == 4
    assert {row["name"] for row in manifest["scenarios"]} == names
