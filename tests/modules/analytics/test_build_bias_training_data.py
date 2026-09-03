"""Regression tests for MMseqs provenance in the bias-training builder."""

from __future__ import annotations

from datetime import date
import sys

import pandas as pd

from cofolder.modules.analytics import build_bias_training_data as builder


def _run_builder(
    monkeypatch,
    temp_dir,
    *,
    skip_mmseqs: bool = False,
    reuse_protein: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    system_path = temp_dir / "system.yaml"
    system_path.write_text(
        "sequences:\n"
        "  - protein:\n"
        "      id: A\n"
        "      sequence: AAAA\n"
        "  - ligand:\n"
        "      id: B\n"
        "      smiles: CCO\n",
        encoding="utf-8",
    )
    components_path = temp_dir / "components.cif"
    components_path.touch()
    mmseqs_db_path = temp_dir / "mmseqs-db"
    mmseqs_db_path.touch()
    protein_path = temp_dir / "protein.csv"
    ligand_path = temp_dir / "ligand.csv"
    bias_path = temp_dir / "bias.csv"
    if reuse_protein:
        protein_path.write_text(
            "query_chain_id,pdb_id,release_date,sequence_similarity,sequence\n"
            "A,HIGH,2020-01-01,30.0,AAAA\n",
            encoding="utf-8",
        )

    monkeypatch.setattr(
        builder,
        "_load_ccd_smiles",
        lambda _path: pd.DataFrame([{"ligand_id": "LIG", "smiles": "CCO"}]),
    )
    monkeypatch.setattr(
        builder,
        "_collect_ccd_hits",
        lambda **_kwargs: pd.DataFrame(
            [
                {
                    "query_chain_id": "B",
                    "ligand_pdbid": "LIG",
                    "hit_smiles": "CCO",
                    "ecfp_similarity": 0.4,
                }
            ]
        ),
    )
    monkeypatch.setattr(
        builder,
        "_search_entries_for_ccd",
        lambda _ccd_id, _timeout: ["LOW1", "NONE", "EDGE", "HIGH"],
    )
    monkeypatch.setattr(
        builder,
        "_entry_release_date",
        lambda _pdb_id, _timeout: date(2020, 1, 1),
    )
    mmseqs_output = pd.DataFrame(
        [
            {"target": "low1_A", "pident": 24.0, "tseq": "AAAT"},
            {"target": "edge_A", "pident": 25.0, "tseq": "AATT"},
            {"target": "high_A", "pident": 30.0, "tseq": "AAAA"},
        ]
    )
    monkeypatch.setattr(builder, "_run_mmseqs", lambda **_kwargs: mmseqs_output.copy())

    argv = [
        "build_bias_training_data",
        "--system_path",
        str(system_path),
        "--components_cif",
        str(components_path),
        "--output_protein_csv",
        str(protein_path),
        "--output_ligand_csv",
        str(ligand_path),
        "--output_bias_csv",
        str(bias_path),
        "--mmseqs_db_path",
        str(mmseqs_db_path),
        "--protein_similarity_threshold",
        "25.0",
        "--workers",
        "1",
    ]
    if not reuse_protein:
        argv.append("--overwrite")
    if skip_mmseqs:
        argv.append("--skip_protein_mmseqs")
    monkeypatch.setattr(sys, "argv", argv)

    assert builder.main() == 0
    return pd.read_csv(protein_path), pd.read_csv(bias_path)


def test_below_threshold_mmseqs_hit_is_kept_for_combined_lookup(monkeypatch, temp_dir):
    protein, combined = _run_builder(monkeypatch, temp_dir)

    assert protein["pdb_id"].tolist() == ["HIGH"]
    low = combined.loc[combined["pdb_id"].eq("LOW1")].iloc[0]
    assert low["sequence_similarity"] == 24.0
    assert low["sequence_similarity_method"] == builder.MMSEQS_PIDENT_METHOD
    assert not (combined["sequence_similarity"].dropna() > 30.0).any()


def test_exact_threshold_mmseqs_hit_is_excluded_from_filtered_protein_view(
    monkeypatch,
    temp_dir,
):
    protein, combined = _run_builder(monkeypatch, temp_dir)

    assert "EDGE" not in set(protein["pdb_id"])
    edge = combined.loc[combined["pdb_id"].eq("EDGE")].iloc[0]
    assert edge["sequence_similarity"] == 25.0
    assert edge["sequence_similarity_method"] == builder.MMSEQS_PIDENT_METHOD


def test_find_ccd_hits_excludes_exact_ligand_threshold(monkeypatch):
    monkeypatch.setattr(builder, "_morgan_fp", lambda _smiles: object())
    monkeypatch.setattr(
        builder.DataStructs,
        "TanimotoSimilarity",
        lambda _query, _reference: 0.35,
    )

    hits = builder._find_ccd_hits(
        query_smiles="QUERY",
        ccd_df=pd.DataFrame([{"ligand_id": "EDGE", "smiles": "REFERENCE"}]),
        threshold=0.35,
        top_k=10,
    )

    assert hits == []


def test_missing_mmseqs_hit_is_unavailable(monkeypatch, temp_dir):
    _, combined = _run_builder(monkeypatch, temp_dir)

    missing = combined.loc[combined["pdb_id"].eq("NONE")].iloc[0]
    assert pd.isna(missing["sequence_similarity"])
    assert (
        missing["sequence_similarity_method"] == builder.UNAVAILABLE_SIMILARITY_METHOD
    )


def test_reused_filtered_protein_csv_still_builds_complete_mmseqs_lookup(
    monkeypatch, temp_dir
):
    protein, combined = _run_builder(monkeypatch, temp_dir, reuse_protein=True)

    assert protein["pdb_id"].tolist() == ["HIGH"]
    low = combined.loc[combined["pdb_id"].eq("LOW1")].iloc[0]
    assert low["sequence_similarity"] == 24.0
    assert low["sequence_similarity_method"] == builder.MMSEQS_PIDENT_METHOD


def test_skip_mmseqs_marks_all_protein_similarities_unavailable(monkeypatch, temp_dir):
    protein, combined = _run_builder(monkeypatch, temp_dir, skip_mmseqs=True)

    assert protein.empty
    assert combined["sequence_similarity"].isna().all()
    assert set(combined["sequence_similarity_method"]) == {
        builder.UNAVAILABLE_SIMILARITY_METHOD
    }


def test_legacy_combined_csv_without_provenance_is_not_reused(temp_dir):
    legacy_path = temp_dir / "bias.csv"
    legacy_path.write_text(
        "pdb_id,sequence_similarity,ecfp_similarity\nLOW1,35.0,0.4\n",
        encoding="utf-8",
    )

    assert builder._read_reusable_bias_csv(legacy_path) is None


def test_combined_csv_with_non_mmseqs_similarity_is_not_reused(temp_dir):
    stale_path = temp_dir / "bias.csv"
    stale_path.write_text(
        "pdb_id,sequence_similarity,sequence_similarity_method,ecfp_similarity\n"
        "LOW1,35.0,pairwise_aligner,0.4\n",
        encoding="utf-8",
    )

    assert builder._read_reusable_bias_csv(stale_path) is None


def test_normalize_mmseqs_hits_does_not_apply_threshold():
    normalized = builder._normalize_mmseqs_hits(
        pd.DataFrame(
            [
                {"target": "4kao_A", "pident": "24.0", "tseq": "AAAA"},
                {"target": "bad_A", "pident": "not-a-number", "tseq": "BBBB"},
                {"target": "none_A", "pident": "20.0", "tseq": None},
            ]
        )
    )

    assert normalized.to_dict("records") == [
        {"pdb_id": "4KAO", "sequence_similarity": 24.0, "sequence": "AAAA"}
    ]
