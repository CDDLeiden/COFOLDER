"""Reference-overlap filtering tests for the Screen recipe."""

from __future__ import annotations

import json
from unittest.mock import patch

import pandas as pd
import pytest
import yaml

from cofolder.recipes.screen import Screen


def _screen(
    sample_system_yaml,
    sample_options_yaml,
    sample_csv_file,
    temp_dir,
    **kwargs,
) -> Screen:
    return Screen(
        wrk_dir=str(temp_dir),
        system_path=str(sample_system_yaml),
        options_path=str(sample_options_yaml),
        variable=["sequences,1,ligand,smiles"],
        variable_csv=str(sample_csv_file),
        col_variable=["smiles"],
        col_id="compound_id",
        **kwargs,
    )


def _write_chain_metrics(validator, ifp, *, ligand_chain="B") -> None:
    results_dir = validator.wrk_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [
            {
                "CHAIN_ID": "A",
                "ENTITY_TYPE": "protein",
                "repeat": 1,
                "diffusion_sample": 0,
            },
            {
                "CHAIN_ID": ligand_chain,
                "ENTITY_TYPE": "ligand",
                "repeat": 1,
                "diffusion_sample": 0,
                "ifp_distance": ifp,
            },
        ]
    ).to_csv(results_dir / "chain_metrics.csv", index=False)


@patch("cofolder.recipes.screen.Validate.run", autospec=True)
def test_filter_accepts_rejects_returns_and_preserves_all_rows(
    mock_validate_run,
    sample_system_yaml,
    sample_options_yaml,
    sample_csv_file,
    temp_dir,
):
    def write_metrics(validator):
        ifp = [1, 1, 0] if validator.wrk_dir.name.startswith("1_") else [0, 0, 1]
        _write_chain_metrics(validator, json.dumps(ifp))

    mock_validate_run.side_effect = write_metrics
    reference_path = temp_dir / "reference.ifp"
    reference_path.write_text("110\n", encoding="utf-8")
    results = _screen(
        sample_system_yaml,
        sample_options_yaml,
        sample_csv_file,
        temp_dir,
        ifp_filter_threshold=0.5,
        pocket_coverage_reference=str(reference_path),
    ).run()

    assert isinstance(results, pd.DataFrame)
    assert results["ifp_filter_status"].tolist() == ["accepted", "rejected"]
    assert results["ifp_filter_pass"].tolist() == [True, False]
    assert results["ifp_filter_overlap"].tolist() == [1.0, 0.0]
    assert len(results) == 2

    expected_columns = {
        "ifp_filter_pass",
        "ifp_filter_status",
        "ifp_filter_reason",
        "ifp_filter_overlap",
        "ifp_filter_threshold",
        "ifp_filter_reference",
    }
    assert expected_columns.issubset(results.columns)
    assert (temp_dir / "results" / "records.jsonl").is_file()
    assert not (temp_dir / "screen_results.csv").exists()
    assert not (temp_dir / "screen_results_with_scores.csv").exists()


@patch("cofolder.recipes.screen.Validate.run", autospec=True)
def test_filter_threshold_boundary_is_inclusive(
    mock_validate_run,
    sample_system_yaml,
    sample_options_yaml,
    sample_csv_file,
    temp_dir,
):
    mock_validate_run.side_effect = lambda validator: _write_chain_metrics(
        validator, "[1, 0]"
    )
    results = _screen(
        sample_system_yaml,
        sample_options_yaml,
        sample_csv_file,
        temp_dir,
        ifp_filter_threshold=0.5,
        pocket_coverage_reference="11",
    ).run()

    assert results["ifp_filter_status"].tolist() == ["accepted", "accepted"]
    assert results["ifp_filter_overlap"].tolist() == [0.5, 0.5]


@patch("cofolder.recipes.screen.Validate.run", autospec=True)
def test_residue_reference_is_mapped_against_predicted_protein_numbering(
    mock_validate_run,
    sample_system_yaml,
    sample_options_yaml,
    sample_csv_file,
    temp_dir,
):
    def write_metrics(validator):
        results_dir = validator.wrk_dir / "results"
        structures_dir = results_dir / "structures"
        structures_dir.mkdir(parents=True, exist_ok=True)
        structure_name = "model.pdb"
        (structures_dir / structure_name).write_text(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\n"
            "ATOM      2  CA  GLY A   2       3.000   0.000   0.000  1.00 20.00           C\n"
            "TER\nEND\n",
            encoding="utf-8",
        )
        pd.DataFrame(
            [
                {
                    "CHAIN_ID": "A",
                    "ENTITY_TYPE": "protein",
                    "cif_file": structure_name,
                    "ifp_distance": None,
                },
                {
                    "CHAIN_ID": "B",
                    "ENTITY_TYPE": "ligand",
                    "cif_file": structure_name,
                    "ifp_distance": "[1, 0]",
                },
            ]
        ).to_csv(results_dir / "chain_metrics.csv", index=False)

    mock_validate_run.side_effect = write_metrics
    results = _screen(
        sample_system_yaml,
        sample_options_yaml,
        sample_csv_file,
        temp_dir,
        ifp_filter_threshold=1.0,
        pocket_coverage_reference="A1",
    ).run()

    assert results["ifp_filter_status"].tolist() == ["accepted", "accepted"]
    assert results["ifp_filter_overlap"].tolist() == [1.0, 1.0]


@pytest.mark.parametrize(
    ("ifp", "reference", "reason"),
    [
        (None, "10", "missing_ifp"),
        ("not-json", "10", "malformed_ifp"),
        ("[1, 0]", "00", "empty_reference"),
        ("[1, 0, 1]", "10", "incompatible_vector_lengths"),
    ],
)
@patch("cofolder.recipes.screen.Validate.run", autospec=True)
def test_filter_not_evaluable_states(
    mock_validate_run,
    ifp,
    reference,
    reason,
    sample_system_yaml,
    sample_options_yaml,
    sample_csv_file,
    temp_dir,
):
    mock_validate_run.side_effect = lambda validator: _write_chain_metrics(
        validator, ifp
    )
    results = _screen(
        sample_system_yaml,
        sample_options_yaml,
        sample_csv_file,
        temp_dir,
        ifp_filter_threshold=0.5,
        pocket_coverage_reference=reference,
    ).run()

    assert results["ifp_filter_status"].tolist() == ["not_evaluable"] * 2
    assert results["ifp_filter_reason"].tolist() == [reason] * 2
    assert results["ifp_filter_pass"].isna().all()


@patch("cofolder.recipes.screen.Validate.run")
def test_disabled_filter_still_has_audit_columns(
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

    assert results["ifp_filter_status"].tolist() == ["not_applied"] * 2
    assert results["ifp_filter_pass"].isna().all()


@pytest.mark.parametrize("threshold", [-0.01, 1.01, float("nan")])
def test_invalid_filter_threshold_is_rejected_before_predictions(
    threshold,
    sample_system_yaml,
    sample_options_yaml,
    sample_csv_file,
    temp_dir,
):
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        _screen(
            sample_system_yaml,
            sample_options_yaml,
            sample_csv_file,
            temp_dir,
            ifp_filter_threshold=threshold,
            pocket_coverage_reference="10",
        )


def test_filter_requires_distance_ifp_and_reference(
    sample_system_yaml,
    sample_options_yaml,
    sample_csv_file,
    temp_dir,
):
    with pytest.raises(ValueError, match="requires distance IFP scoring"):
        _screen(
            sample_system_yaml,
            sample_options_yaml,
            sample_csv_file,
            temp_dir,
            scoring_functions=["sasa"],
            ifp_filter_threshold=0.5,
            pocket_coverage_reference="10",
        )
    with pytest.raises(ValueError, match="requires --pocket_coverage_reference"):
        _screen(
            sample_system_yaml,
            sample_options_yaml,
            sample_csv_file,
            temp_dir,
            ifp_filter_threshold=0.5,
        )


def test_malformed_and_nonexistent_filter_references_are_rejected(
    sample_system_yaml,
    sample_options_yaml,
    sample_csv_file,
    temp_dir,
):
    with pytest.raises(ValueError, match="Invalid --pocket_coverage_reference"):
        _screen(
            sample_system_yaml,
            sample_options_yaml,
            sample_csv_file,
            temp_dir,
            ifp_filter_threshold=0.5,
            pocket_coverage_reference="not-a-reference",
        )
    with pytest.raises(ValueError, match="file does not exist"):
        _screen(
            sample_system_yaml,
            sample_options_yaml,
            sample_csv_file,
            temp_dir,
            ifp_filter_threshold=0.5,
            pocket_coverage_reference=str(temp_dir / "missing.ifp"),
        )
    empty_reference = temp_dir / "empty.ifp"
    empty_reference.write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match="file is empty"):
        _screen(
            sample_system_yaml,
            sample_options_yaml,
            sample_csv_file,
            temp_dir,
            ifp_filter_threshold=0.5,
            pocket_coverage_reference=str(empty_reference),
        )


def test_multiple_ligands_require_explicit_filter_chain(
    sample_system_yaml,
    sample_options_yaml,
    sample_csv_file,
    temp_dir,
):
    data = yaml.safe_load(sample_system_yaml.read_text(encoding="utf-8"))
    data["sequences"].append({"ligand": {"id": "C", "smiles": "CC"}})
    sample_system_yaml.write_text(yaml.safe_dump(data), encoding="utf-8")

    with pytest.raises(ValueError, match="--ifp_ligand_chain is required"):
        _screen(
            sample_system_yaml,
            sample_options_yaml,
            sample_csv_file,
            temp_dir,
            ifp_filter_threshold=0.5,
            pocket_coverage_reference="10",
        )


@patch("cofolder.recipes.screen.Validate.run", autospec=True)
def test_explicit_filter_chain_selects_one_of_multiple_ligands(
    mock_validate_run,
    sample_options_yaml,
    sample_csv_file,
    temp_dir,
):
    system_path = temp_dir / "multi_ligand.yaml"
    system_path.write_text(
        yaml.safe_dump(
            {
                "sequences": [
                    {"protein": {"id": "A", "sequence": "MKRAAT"}},
                    {"ligand": {"id": "B", "smiles": "CCO"}},
                    {"ligand": {"id": "C", "smiles": "CCN"}},
                ]
            }
        ),
        encoding="utf-8",
    )

    def write_metrics(validator):
        results_dir = validator.wrk_dir / "results"
        results_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(
            [
                {"CHAIN_ID": "A", "ENTITY_TYPE": "protein", "ifp_distance": None},
                {"CHAIN_ID": "B", "ENTITY_TYPE": "ligand", "ifp_distance": "[1, 0]"},
                {"CHAIN_ID": "C", "ENTITY_TYPE": "ligand", "ifp_distance": "[0, 1]"},
            ]
        ).to_csv(results_dir / "chain_metrics.csv", index=False)

    mock_validate_run.side_effect = write_metrics
    results = Screen(
        wrk_dir=str(temp_dir / "screen"),
        system_path=str(system_path),
        options_path=str(sample_options_yaml),
        variable=["sequences,1,ligand,smiles"],
        variable_csv=str(sample_csv_file),
        col_variable=["smiles"],
        col_id="compound_id",
        ifp_filter_threshold=0.5,
        ifp_ligand_chain="C",
        pocket_coverage_reference="10",
    ).run()

    assert results["ifp_filter_status"].tolist() == ["rejected", "rejected"]
