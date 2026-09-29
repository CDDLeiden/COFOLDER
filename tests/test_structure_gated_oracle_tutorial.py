"""Scientific and integration checks for the MAPK14 structure-gated tutorial."""

from __future__ import annotations

import json
from pathlib import Path
import sys
from types import SimpleNamespace
from unittest.mock import patch

import pandas as pd
import pytest
from rdkit import Chem


REPO_ROOT = Path(__file__).resolve().parents[1]
TUTORIALS_DIR = REPO_ROOT / "tutorials"
sys.path.insert(0, str(TUTORIALS_DIR))

from _structure_gated_oracle import (  # noqa: E402
    CONFIG_PATH,
    ILLUSTRATIVE_RESULTS_PATH,
    analyze_pose,
    bounded_score,
    load_config,
    load_ligands,
    make_oracle_scoring_function,
    match_interaction_features,
    match_key_interactions,
    rank_candidates,
    structure_gated_reward,
    write_audit_bundle,
)


def test_literature_panel_has_six_valid_unique_molecules_and_matched_kd_series():
    ligands = load_ligands()
    assert set(ligands["candidate_id"]) == {
        "sb203580",
        "birb43",
        "birb48",
        "birb796",
        "sb202474",
        "um101",
    }
    assert all(
        Chem.MolFromSmiles(value) is not None for value in ligands["canonical_smiles"]
    )
    assert not ligands["canonical_smiles"].duplicated().any()

    series = ligands.set_index("candidate_id").loc[["birb43", "birb48", "birb796"]]
    assert series["literature_endpoint"].tolist() == ["Kd", "Kd", "Kd"]
    assert series["literature_source"].nunique() == 1
    assert series["literature_value_nm"].tolist() == pytest.approx([14.0, 0.52, 0.046])


def test_reward_enforces_type_then_interactions_then_score():
    assert structure_gated_reward(
        is_type_ii=True,
        matched_interactions=0,
        raw_score=-1000,
        key_interaction_count=3,
    ) > structure_gated_reward(
        is_type_ii=False,
        matched_interactions=3,
        raw_score=1000,
        key_interaction_count=3,
    )
    assert structure_gated_reward(
        is_type_ii=False,
        matched_interactions=2,
        raw_score=-1000,
        key_interaction_count=3,
    ) > structure_gated_reward(
        is_type_ii=False,
        matched_interactions=1,
        raw_score=1000,
        key_interaction_count=3,
    )
    assert bounded_score(-2) < bounded_score(0) < bounded_score(2)


def test_interaction_matching_deduplicates_atom_level_occurrences():
    configured = ("A:71:hb_donor", "A:109:hb_acceptor")
    matched, missing = match_key_interactions(
        ["A:71:hb_donor", "A:71:hb_donor", "A:71:hb_donor"], configured
    )
    assert matched == ("A:71:hb_donor",)
    assert missing == ("A:109:hb_acceptor",)


def test_atom_filtered_features_deduplicate_and_reject_asp_sidechain():
    features = load_config()["interaction_policy"]["features"]
    event = {
        "interaction_key": "A:168:hb_acceptor",
        "ligand_role": "acceptor",
        "protein_role": "donor",
        "protein_atoms": [{"atom_name": "N"}],
    }
    sidechain = {**event, "protein_atoms": [{"atom_name": "OD1"}]}
    matched, _, evidence = match_interaction_features(
        [event, event, sidechain], features
    )
    assert matched == ("A:168:hb_acceptor",)
    assert len(evidence["asp168_backbone_hbond"]) == 2


def test_illustrative_ranking_matches_the_presented_rule_and_is_label_independent():
    frame = pd.read_csv(ILLUSTRATIVE_RESULTS_PATH)
    ranked = rank_candidates(frame).sort_values("structure_gated_rank")
    assert ranked["candidate_id"].tolist() == [
        "birb796",
        "birb48",
        "birb43",
        "sb203580",
        "um101",
        "sb202474",
    ]
    changed = frame.copy()
    changed["expected_role"] = list(reversed(frame["expected_role"].tolist()))
    assert (
        rank_candidates(changed)
        .sort_values("structure_gated_rank")["candidate_id"]
        .tolist()
        == ranked["candidate_id"].tolist()
    )


def test_equal_scores_have_a_stable_candidate_id_tie_break_and_missing_rows_are_unranked():
    frame = pd.DataFrame(
        [
            {
                "candidate_id": "b",
                "raw_score": 1.0,
                "matched_interactions": 1,
                "structural_classification": "type_II",
                "evaluation_status": "evaluable",
            },
            {
                "candidate_id": "a",
                "raw_score": 1.0,
                "matched_interactions": 1,
                "structural_classification": "type_II",
                "evaluation_status": "evaluable",
            },
            {
                "candidate_id": "missing",
                "raw_score": None,
                "matched_interactions": None,
                "structural_classification": None,
                "evaluation_status": "not_evaluable",
            },
        ]
    )
    ranked = rank_candidates(frame).set_index("candidate_id")
    assert ranked.at["a", "structure_gated_rank"] == 1
    assert ranked.at["b", "structure_gated_rank"] == 2
    assert pd.isna(ranked.at["missing", "structure_gated_rank"])


def test_frozen_dfg_rule_brackets_type_i_and_type_ii_references():
    definition = load_config()["structural_predicate"]["dfg_conformation"]
    d1 = definition["d1"]
    d2 = definition["d2"]
    assert d1["type_ii_reference_angstrom"] <= d1["threshold_angstrom"]
    assert d2["type_ii_reference_angstrom"] >= d2["threshold_angstrom"]
    assert d1["type_i_reference_angstrom"] > d1["threshold_angstrom"]
    assert d2["type_i_reference_angstrom"] < d2["threshold_angstrom"]


def test_configured_interactions_use_serialized_prolif_taxonomy():
    features = load_config()["interaction_policy"]["features"]
    assert [feature["key"] for feature in features] == [
        "A:71:hb_donor",
        "A:168:hb_acceptor",
        "A:109:hb_acceptor",
    ]
    assert [feature["receptor_atoms"] for feature in features] == [
        ["OE1", "OE2"],
        ["N"],
        ["N"],
    ]


def test_pose_analysis_pairs_geometry_interactions_and_score_from_one_structure(
    tmp_path: Path,
):
    structure = tmp_path / "prediction.pdb"
    structure.write_text(
        """\
ATOM      1  CA  ALA A  51       0.000   0.000   0.000  1.00 20.00           C
ATOM      2  CA  GLU A  71       0.000   0.000   0.000  1.00 20.00           C
ATOM      3  CA  ASN A 155       5.000   0.000   0.000  1.00 20.00           C
ATOM      4  CA  PHE A 169      10.000   0.000   0.000  1.00 20.00           C
HETATM    5  C1  LIG B   1       0.500   0.000   0.000  1.00 20.00           C
END
""",
        encoding="utf-8",
    )
    policy = load_config()
    policy["structural_predicate"]["atp_site_residues"] = [51]
    policy["structural_predicate"]["back_pocket_residues"] = [71]
    policy["structural_predicate"]["remote_um101_region"] = []
    prepared = {
        "events": [
            {
                "interaction_key": "A:71:hb_donor",
                "ligand_role": "donor",
                "protein_role": "acceptor",
                "protein_atoms": [{"atom_name": "OE1"}],
            }
        ],
        "preparation": {"ligand": {"formal_charge": 0}},
    }
    with patch("_structure_gated_oracle._prepared_interactions", return_value=prepared):
        row = analyze_pose(
            structure,
            raw_score=2.5,
            candidate_id="candidate",
            prediction_id="prediction-1",
            canonical_smiles="C",
            config=policy,
        )
    assert row["evaluation_status"] == "evaluable"
    assert row["structural_classification"] == "type_II"
    assert row["matched_interactions"] == 1
    assert row["raw_score"] == pytest.approx(2.5)


def test_preparation_failure_is_not_evaluable(tmp_path: Path):
    structure = tmp_path / "prediction.pdb"
    structure.write_text(
        """\
ATOM      1  CA  ALA A  51       0.000   0.000   0.000  1.00 20.00           C
ATOM      2  CA  GLU A  71       0.000   0.000   0.000  1.00 20.00           C
ATOM      3  CA  ASN A 155       5.000   0.000   0.000  1.00 20.00           C
ATOM      4  CA  PHE A 169      10.000   0.000   0.000  1.00 20.00           C
HETATM    5  C1  LIG B   1       0.500   0.000   0.000  1.00 20.00           C
END
""",
        encoding="utf-8",
    )
    policy = load_config()
    policy["structural_predicate"]["atp_site_residues"] = [51]
    policy["structural_predicate"]["back_pocket_residues"] = [71]
    policy["structural_predicate"]["remote_um101_region"] = []
    with patch(
        "_structure_gated_oracle._prepared_interactions",
        side_effect=ValueError("ligand_atom_mapping_failed"),
    ):
        row = analyze_pose(
            structure,
            raw_score=1.0,
            candidate_id="candidate",
            prediction_id="prediction-1",
            canonical_smiles="C",
            config=policy,
        )
    assert row["evaluation_status"] == "not_evaluable"
    assert "ligand_atom_mapping_failed" in row["evaluation_reason"]
    assert pd.isna(row["glu71_sidechain_hbond"])


def test_oracle_callback_rejects_multiple_prediction_paths(tmp_path: Path):
    context = SimpleNamespace(
        aggregated_metrics={"ligand_B__affinity_pred_value": 1.0},
        system_metrics=pd.DataFrame(
            {"cif_file": ["one.cif", "two.cif"], "model_name": ["one", "two"]}
        ),
        chain_metrics=pd.DataFrame(),
        run_dir=tmp_path,
    )
    sink = []
    callback = make_oracle_scoring_function(
        candidate_id="candidate", expected_role="unknown", audit_sink=sink
    )
    with pytest.raises(ValueError, match="exactly one predicted structure"):
        callback(context)
    assert sink == []


def test_oracle_callback_resolves_normalized_structure_and_orients_affinity(
    tmp_path: Path,
):
    structure = tmp_path / "results" / "structures" / "prediction.cif"
    structure.parent.mkdir(parents=True)
    structure.touch()
    context = SimpleNamespace(
        aggregated_metrics={"ligand_B__affinity_pred_value": -1.25},
        system_metrics=pd.DataFrame(
            {"cif_file": [structure.name], "model_name": ["prediction"]}
        ),
        chain_metrics=pd.DataFrame(),
        run_dir=tmp_path,
    )
    sink = []
    callback = make_oracle_scoring_function(
        candidate_id="candidate", expected_role="unknown", audit_sink=sink
    )
    evaluable = {
        "evaluation_status": "evaluable",
        "scalar_reward": 2.0,
    }
    with patch(
        "_structure_gated_oracle.analyze_pose", return_value=evaluable
    ) as analyze:
        assert callback(context) == pytest.approx(2.0)
    assert analyze.call_args.args[0] == structure
    assert analyze.call_args.kwargs["raw_score"] == pytest.approx(1.25)
    assert sink == [evaluable]


def test_audit_bundle_records_policy_and_input_hashes(tmp_path: Path):
    frame = pd.read_csv(ILLUSTRATIVE_RESULTS_PATH)
    csv_path, evidence_path, provenance_path = write_audit_bundle(
        frame,
        tmp_path,
        inputs=(CONFIG_PATH, ILLUSTRATIVE_RESULTS_PATH),
        metadata={"seed": 7},
    )
    assert csv_path.is_file()
    assert evidence_path.is_file()
    evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
    assert evidence["schema_version"] == 1
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    assert provenance["seed"] == 7
    assert provenance["policy_sha256"]
    assert len(provenance["input_sha256"]) == 2


def test_versioned_acceptance_records_reference_recovery_and_asset_hashes():
    asset_dir = CONFIG_PATH.parent
    acceptance = json.loads((asset_dir / "acceptance.json").read_text())
    reference = acceptance["reference_checks"]["1KV2"]
    assert set(reference) >= {
        "glu71_sidechain_hbond",
        "asp168_backbone_hbond",
        "met109_backbone_hbond",
    }
    for key in ("summary_csv", "atom_evidence", "provenance", "reference_evidence"):
        record = acceptance[key]
        assert record["sha256"] == _sha256(asset_dir / record["path"])


def test_acceptance_preparation_audit_preserves_curated_graphs_and_coordinates():
    asset_dir = CONFIG_PATH.parent
    evidence = json.loads((asset_dir / "acceptance-v2-interactions.json").read_text())
    ligands = load_ligands().set_index("candidate_id")
    assert len(evidence["candidates"]) == 6
    for candidate in evidence["candidates"]:
        candidate_id = candidate["candidate_id"]
        preparation = candidate["preparation"]
        ligand = preparation["ligand"]
        expected = str(ligands.at[candidate_id, "canonical_smiles"])
        assert ligand["canonical_isomeric_smiles"] == expected
        assert len(ligand["atom_mapping"]) == ligand["heavy_atom_count"]
        assert ligand["max_heavy_atom_displacement_angstrom"] == 0
        assert preparation["protein"]["max_heavy_atom_displacement_angstrom"] == 0
        assert isinstance(ligand["formal_charge"], int)


def _sha256(path: Path) -> str:
    import hashlib

    return hashlib.sha256(path.read_bytes()).hexdigest()
