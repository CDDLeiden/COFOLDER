"""IFP clustering and end-to-end Screen output tests."""

from __future__ import annotations

import json
from unittest.mock import patch

import pandas as pd
import pytest
import yaml

from cofolder.modules.contracts import WorkflowExecutionError
from cofolder.recipes.screen import Screen


def _screen(system_path, options_path, csv_path, work_dir, **kwargs):
    ligand_chain = kwargs.pop("ligand_chain", "B")
    return Screen(
        wrk_dir=str(work_dir),
        system_path=str(system_path),
        options_path=str(options_path),
        ligand_chain=ligand_chain,
        library=str(csv_path),
        smiles_column="smiles",
        col_id="compound_id",
        **kwargs,
    )


def _write_ifp(validator, value):
    results_dir = validator.wrk_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [
            {"CHAIN_ID": "A", "ENTITY_TYPE": "protein", "ifp_distance": None},
            {"CHAIN_ID": "B", "ENTITY_TYPE": "ligand", "ifp_distance": value},
        ]
    ).to_csv(results_dir / "chain_metrics.csv", index=False)


@patch("cofolder.recipes.screen.Validate.run", autospec=True)
def test_screen_clusters_valid_ifps_and_excludes_invalid_or_incompatible_rows(
    mock_validate_run,
    sample_system_yaml,
    sample_options_yaml,
    temp_dir,
):
    csv_path = temp_dir / "cluster_library.csv"
    pd.DataFrame(
        {
            "compound_id": ["A", "B", "C", "D", "E", "F"],
            "smiles": ["CC", "CCC", "CCCC", "CCO", "CCN", "CO"],
        }
    ).to_csv(csv_path, index=False)
    values = {
        "compound_000001": "[1, 1, 0, 0]",
        "compound_000002": "[1, 0, 0, 0]",
        "compound_000003": "[0, 0, 1, 1]",
        "compound_000004": "[0, 0, 0, 0]",
        "compound_000005": "not-json",
        "compound_000006": "[1, 0, 0]",
    }

    def write_metrics(validator):
        _write_ifp(validator, values[validator.wrk_dir.name])

    mock_validate_run.side_effect = write_metrics
    results = _screen(
        sample_system_yaml,
        sample_options_yaml,
        csv_path,
        temp_dir / "screen",
        cluster_ifps=True,
        ifp_cluster_similarity_threshold=0.5,
    ).run()

    assert results["ifp_cluster_id"].iloc[:4].tolist() == [
        "IFP001",
        "IFP001",
        "IFP002",
        "IFP003",
    ]
    assert results["ifp_cluster_id"].iloc[4:].isna().all()
    assert results["ifp_cluster_status"].tolist() == [
        "clustered",
        "clustered",
        "clustered",
        "clustered",
        "not_evaluable",
        "not_evaluable",
    ]
    assert (temp_dir / "screen" / "results" / "records.jsonl").is_file()

    summary = pd.read_csv(temp_dir / "screen" / "results" / "ifp_cluster_summary.csv")
    assert summary["ifp_cluster_id"].tolist() == ["IFP001", "IFP002", "IFP003"]
    assert summary["size"].tolist() == [2, 1, 1]
    assert json.loads(summary.loc[0, "member_ids"]) == ["A", "B"]


@patch("cofolder.recipes.screen.Validate.run", autospec=True)
def test_acceptance_screen_combines_filtering_clustering_and_returned_dataframe(
    mock_validate_run,
    sample_system_yaml,
    sample_options_yaml,
    temp_dir,
):
    csv_path = temp_dir / "acceptance_library.csv"
    pd.DataFrame(
        {
            "compound_id": ["hit-1", "hit-2", "miss"],
            "smiles": ["CC", "CCC", "CCO"],
            "series": ["one", "one", "two"],
        }
    ).to_csv(csv_path, index=False)
    values = ["[1, 1, 0]", "[1, 0, 0]", "[0, 0, 1]"]

    def write_metrics(validator):
        position = int(validator.wrk_dir.name.removeprefix("compound_")) - 1
        _write_ifp(validator, values[position])

    mock_validate_run.side_effect = write_metrics
    work_dir = temp_dir / "acceptance_screen"
    results = _screen(
        sample_system_yaml,
        sample_options_yaml,
        csv_path,
        work_dir,
        merge_data="series",
        ifp_filter_threshold=0.5,
        pocket_coverage_reference="110",
        cluster_ifps=True,
    ).run()

    assert isinstance(results, pd.DataFrame)
    assert results["ifp_filter_status"].tolist() == [
        "accepted",
        "accepted",
        "rejected",
    ]
    assert results["ifp_cluster_id"].tolist() == ["IFP001", "IFP001", "IFP002"]
    assert len(results) == 3
    assert (work_dir / "results" / "records.jsonl").is_file()
    assert (work_dir / "results" / "ifp_cluster_summary.csv").is_file()


@patch("cofolder.recipes.screen.Validate.run")
def test_clustering_disabled_keeps_stable_columns_without_summary(
    mock_validate_run,
    sample_system_yaml,
    sample_options_yaml,
    sample_csv_file,
    temp_dir,
):
    results = _screen(
        sample_system_yaml,
        sample_options_yaml,
        sample_csv_file,
        temp_dir,
    ).run()

    assert results["ifp_cluster_status"].tolist() == ["not_applied", "not_applied"]
    assert results["ifp_cluster_id"].isna().all()
    assert not (temp_dir / "results" / "ifp_cluster_summary.csv").exists()


@pytest.mark.parametrize("threshold", [-0.01, 1.01, float("nan")])
def test_invalid_cluster_threshold_is_rejected_before_predictions(
    threshold,
    sample_system_yaml,
    sample_options_yaml,
    sample_csv_file,
    temp_dir,
):
    with pytest.raises(WorkflowExecutionError, match=r"\[0, 1\]"):
        _screen(
            sample_system_yaml,
            sample_options_yaml,
            sample_csv_file,
            temp_dir,
            cluster_ifps=True,
            ifp_cluster_similarity_threshold=threshold,
        ).run()


def test_clustering_requires_distance_ifp_and_valid_ligand_selector(
    sample_system_yaml,
    sample_options_yaml,
    sample_csv_file,
    temp_dir,
):
    with pytest.raises(WorkflowExecutionError, match="requires distance IFP scoring"):
        _screen(
            sample_system_yaml,
            sample_options_yaml,
            sample_csv_file,
            temp_dir,
            cluster_ifps=True,
            scoring_functions=["sasa"],
        ).run()

    system = yaml.safe_load(sample_system_yaml.read_text(encoding="utf-8"))
    system["sequences"].append({"ligand": {"id": "C", "smiles": "CC"}})
    sample_system_yaml.write_text(yaml.safe_dump(system), encoding="utf-8")
    with pytest.raises(WorkflowExecutionError, match="unknown chain"):
        _screen(
            sample_system_yaml,
            sample_options_yaml,
            sample_csv_file,
            temp_dir,
            cluster_ifps=True,
            ligand_chain="Z",
        ).run()
