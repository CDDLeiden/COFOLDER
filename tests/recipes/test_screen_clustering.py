"""IFP clustering and end-to-end Screen output tests."""

from __future__ import annotations

import json
import subprocess
from unittest.mock import patch

import gemmi
import pandas as pd
import pytest
import yaml

from cofolder.modules.contracts import (
    OutputIdentity,
    WorkflowExecutionError,
    WorkflowKind,
)
from cofolder.recipes._results import write_frame_bundle
from cofolder.recipes.screen import Screen
from tests.modules.analytics.test_reproduction import _write_predicted_pdb


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


@pytest.mark.parametrize("structure_suffix", [".pdb", ".cif"])
@pytest.mark.parametrize("metric_source", ["public_records", "csv_fallback"])
@pytest.mark.parametrize("taxonomy", ["distance", "prolif"])
def test_supported_screen_metrics_produce_identity_bearing_fingerprints(
    structure_suffix,
    metric_source,
    taxonomy,
    sample_system_yaml,
    sample_options_yaml,
    sample_csv_file,
    temp_dir,
    monkeypatch,
):
    run_dir = temp_dir / (
        f"{taxonomy}-{metric_source}-{structure_suffix.removeprefix('.')}"
    )
    structures_dir = run_dir / "results" / "structures"
    structures_dir.mkdir(parents=True)
    pdb_path = structures_dir / "prediction.pdb"
    _write_predicted_pdb(
        pdb_path,
        ((1.3, 1.2, 0.0), (2.7, 1.2, 0.0)),
    )
    structure_path = pdb_path
    if structure_suffix == ".cif":
        structure_path = structures_dir / "prediction.cif"
        gemmi.read_structure(str(pdb_path)).make_mmcif_document().write_file(
            str(structure_path)
        )
        pdb_path.unlink()

    system_df = pd.DataFrame(
        [
            {
                "cif_file": structure_path.name,
                "model_name": "model-a",
                "repeat": 2,
                "diffusion_sample": 1,
                "confidence_score": 0.9,
            }
        ]
    )
    chain_df = pd.DataFrame(
        [
            {
                "CHAIN_ID": "A",
                "ENTITY_TYPE": "protein",
                "cif_file": structure_path.name,
                "model_name": "model-a",
                "repeat": 2,
                "diffusion_sample": 1,
                "chains_ptm": 0.8,
            },
            {
                "CHAIN_ID": "Z",
                "ENTITY_TYPE": "ligand",
                "cif_file": structure_path.name,
                "model_name": "model-a",
                "repeat": 2,
                "diffusion_sample": 1,
                "chains_ptm": 0.7,
            },
        ]
    )
    if metric_source == "public_records":
        write_frame_bundle(
            system_df,
            chain_df,
            output_dir=run_dir / "results",
            identity=OutputIdentity(
                workflow=WorkflowKind.VALIDATE,
                run_id="identity-proof",
                system_id="system",
                runner_id="boltz2",
            ),
        )
    else:
        system_df.to_csv(run_dir / "results" / "system_metrics.csv", index=False)
        chain_df.to_csv(run_dir / "results" / "chain_metrics.csv", index=False)

    if taxonomy == "prolif":
        worker_payload = {
            "ligand": {
                "chain_id": "Z",
                "residue_number": 1,
                "insertion_code": "",
                "residue_name": "LIG",
            },
            "receptor_chains": ["A"],
            "interactions": ["A:1:HBAcceptor"],
            "events": [
                {
                    "interaction_key": "A:1:HBAcceptor",
                    "ligand_role": "acceptor",
                    "protein_role": "donor",
                    "ligand_atoms": [
                        {
                            "chain_id": "Z",
                            "residue_number": 1,
                            "insertion_code": "",
                            "residue_name": "LIG",
                            "atom_name": "O1",
                            "element": "O",
                            "atom_serial": 3,
                            "source_index": 2,
                        }
                    ],
                    "protein_atoms": [
                        {
                            "chain_id": "A",
                            "residue_number": 1,
                            "insertion_code": "",
                            "residue_name": "ALA",
                            "atom_name": "N",
                            "element": "N",
                            "atom_serial": 1,
                            "source_index": 0,
                        }
                    ],
                    "geometry": [
                        {"name": "distance", "value": 3.0, "unit": "angstrom"}
                    ],
                }
            ],
        }
        monkeypatch.setattr(
            subprocess,
            "run",
            lambda *args, **kwargs: subprocess.CompletedProcess(
                args[0], 0, json.dumps(worker_payload), ""
            ),
        )

    screen = _screen(
        sample_system_yaml,
        sample_options_yaml,
        sample_csv_file,
        temp_dir / "screen",
        ligand_chain="Z",
        cluster_ifps=True,
        ifp_taxonomy=taxonomy,
    )
    (
        fingerprint,
        reason,
    ) = screen._screen_postprocessor()._load_selected_interaction_fingerprint(
        run_dir,
        repeat_id=2,
        sample_id=1,
    )

    assert reason == ""
    assert fingerprint is not None
    assert fingerprint.ligand.chain_id == "Z"
    assert fingerprint.ligand.residue_number == 1
    assert fingerprint.receptor_chains == ("A",)
    assert fingerprint.interactions
    assert {item.receptor.chain_id for item in fingerprint.interactions} == {"A"}

    clustered_rows = pd.DataFrame(
        [
            {
                "compound_id": "identity-proof",
                "execution_key": "identity-proof|repeat=2|model=model-a|sample=1",
                "run_dir": str(run_dir),
                "repeat_id": 2,
                "sample_id": 1,
                "model_id": "model-a",
                "status": "success",
                **screen._default_cluster_result(),
            }
        ]
    )
    detailed_rows = clustered_rows.copy()
    screen._apply_ifp_clustering(clustered_rows, detailed_rows)
    screen._publish_prolif_events(clustered_rows)
    assert clustered_rows["ifp_cluster_status"].tolist() == ["clustered"]
    assert detailed_rows["ifp_cluster_id"].tolist() == ["IFP001"]
    assert (temp_dir / "screen" / "results" / "ifp_cluster_summary.csv").is_file()
    linkage = pd.read_csv(temp_dir / "screen" / "results" / "ifp_cluster_linkage.csv")
    leaves = pd.read_csv(temp_dir / "screen" / "results" / "ifp_cluster_leaf_order.csv")
    assert linkage.empty
    assert leaves["member_id"].tolist() == [
        "identity-proof|repeat=2|model=model-a|sample=1"
    ]
    if taxonomy == "prolif":
        event_path = temp_dir / "screen" / "results" / "ifp_interaction_events.jsonl"
        event = json.loads(event_path.read_text(encoding="utf-8"))
        assert event["compound_id"] == "identity-proof"
        assert event["protein_atoms"][0]["atom_name"] == "N"
    else:
        assert not (
            temp_dir / "screen" / "results" / "ifp_interaction_events.jsonl"
        ).exists()


@patch("cofolder.recipes.screen.Validate.run", autospec=True)
def test_screen_does_not_cluster_vector_only_ifps_without_residue_identities(
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

    assert results["ifp_cluster_id"].isna().all()
    assert results["ifp_cluster_status"].tolist() == ["not_evaluable"] * 6
    assert (temp_dir / "screen" / "results" / "records.jsonl").is_file()

    summary = pd.read_csv(temp_dir / "screen" / "results" / "ifp_cluster_summary.csv")
    assert summary.empty
    assert pd.read_csv(
        temp_dir / "screen" / "results" / "ifp_cluster_linkage.csv"
    ).empty
    assert pd.read_csv(
        temp_dir / "screen" / "results" / "ifp_cluster_leaf_order.csv"
    ).empty
    manifest = json.loads(
        (temp_dir / "screen" / "results" / "manifest.json").read_text(encoding="utf-8")
    )
    assert {artifact["label"] for artifact in manifest["artifacts"]} >= {
        "ifp_cluster_summary",
        "ifp_cluster_linkage",
        "ifp_cluster_leaf_order",
    }


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
    assert results["ifp_cluster_id"].isna().all()
    assert results["ifp_cluster_status"].tolist() == ["not_evaluable"] * 3
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
    assert not (temp_dir / "results" / "ifp_cluster_linkage.csv").exists()
    assert not (temp_dir / "results" / "ifp_cluster_leaf_order.csv").exists()


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
    # Clustering now requests its fingerprint extraction independently of the
    # Validate scoring-function selection.
    _screen(
        sample_system_yaml,
        sample_options_yaml,
        sample_csv_file,
        temp_dir,
        cluster_ifps=True,
        scoring_functions=["sasa"],
    )._validate_config()

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
