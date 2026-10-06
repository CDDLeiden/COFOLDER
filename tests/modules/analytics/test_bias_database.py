from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest
from unittest.mock import Mock
from concurrent.futures import ThreadPoolExecutor
import time

from cofolder.modules.analytics.bias_database import (
    BIAS_DATABASE_SCHEMA_VERSION,
    LIGAND_TABLE_NAME,
    PROTEIN_METADATA_NAME,
    PROTEIN_SEQUENCE_INDEX_NAME,
    materialize_bias_references,
    parse_bias_release_policy,
    validate_bias_database_bundle,
)
from cofolder.modules.analytics import bias_database
from cofolder.modules.input.system import System


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _bundle(root: Path, kind: str) -> Path:
    root.mkdir(parents=True)
    if kind == "protein":
        (root / "mmseqs").mkdir()
        (root / "mmseqs/db").write_text("database", encoding="utf-8")
        pd.DataFrame(
            [
                {"pdb_id": "1AAA", "release_date": "2023-05-31"},
                {"pdb_id": "2BBB", "release_date": "2023-06-01"},
            ]
        ).to_csv(root / PROTEIN_METADATA_NAME, index=False)
        pd.DataFrame(
            [
                {"pdb_id": "1AAA", "target_id": "1AAA_A", "sequence": "AAAA"},
                {"pdb_id": "2BBB", "target_id": "2BBB_A", "sequence": "BBBB"},
            ]
        ).to_csv(root / PROTEIN_SEQUENCE_INDEX_NAME, index=False)
    else:
        pd.DataFrame(
            [
                {"pdb_id": "1AAA", "release_date": "2023-05-31", "ligand_id": "ETH", "smiles": "CCO"},
                {"pdb_id": "2BBB", "release_date": "2023-06-01", "ligand_id": "ETH", "smiles": "CCO"},
            ]
        ).to_csv(root / LIGAND_TABLE_NAME, index=False)
    files = {
        str(path.relative_to(root)): _sha(path)
        for path in root.rglob("*")
        if path.is_file()
    }
    (root / "manifest.json").write_text(
        json.dumps({"schema_version": BIAS_DATABASE_SCHEMA_VERSION, "kind": kind, "files": files}),
        encoding="utf-8",
    )
    return root


def test_bundle_validation_is_relocatable_and_detects_checksum_changes(temp_dir):
    root = _bundle(temp_dir / "ligand", "ligand")
    assert validate_bias_database_bundle(root, "ligand").root == root.resolve()
    moved = temp_dir / "moved"
    root.rename(moved)
    validate_bias_database_bundle(moved, "ligand")
    with (moved / LIGAND_TABLE_NAME).open("a", encoding="utf-8") as handle:
        handle.write("corrupt")
    with pytest.raises(ValueError, match="checksum mismatch"):
        validate_bias_database_bundle(moved, "ligand")


def test_materializes_distinct_multi_entity_references_and_strict_cutoff(
    monkeypatch, temp_dir
):
    protein = validate_bias_database_bundle(
        _bundle(temp_dir / "protein", "protein"), "protein"
    )
    ligand = validate_bias_database_bundle(
        _bundle(temp_dir / "ligand", "ligand"), "ligand"
    )
    calls = []

    def fake_mmseqs(binary, fasta, target, workers, tmp, max_seqs):
        sequence = "".join(
            line for line in Path(fasta).read_text().splitlines() if not line.startswith(">")
        )
        calls.append(sequence)
        return pd.DataFrame(
            [{"pdb_id": "1AAA", "sequence_similarity": 90.0, "sequence": sequence},
             {"pdb_id": "2BBB", "sequence_similarity": 80.0, "sequence": sequence}]
        )

    monkeypatch.setattr(
        "cofolder.tools.build_bias_training_data._run_mmseqs", fake_mmseqs
    )
    monkeypatch.setattr(
        "cofolder.modules.analytics.bias_training._resolve_mmseqs_bin", lambda: "mmseqs"
    )
    system = System(
        system={
            "sequences": [
                {"protein": {"id": ["A", "C"], "sequence": "AAAA"}},
                {"protein": {"id": "D", "sequence": "DDDD"}},
                {"ligand": {"id": ["B", "E"], "smiles": "CCO"}},
            ]
        }
    )
    artifacts = materialize_bias_references(
        system_obj=system,
        protein_bundle=protein,
        ligand_bundle=ligand,
        output_dir=temp_dir / "results",
        release_cutoff="2023-06-01",
        ligand_similarity_threshold=0.35,
        selected_chains=None,
        query_cache_path=temp_dir / "cache",
    )
    protein_rows = pd.read_csv(artifacts.protein_path)
    ligand_rows = pd.read_csv(artifacts.ligand_path)
    assert calls == ["AAAA", "DDDD"]
    assert set(protein_rows["query_chain_id"]) == {"A", "C", "D"}
    assert set(ligand_rows["query_chain_id"]) == {"B", "E"}
    assert set(protein_rows["pdb_id"]) == {"1AAA"}
    assert set(ligand_rows["pdb_id"]) == {"1AAA"}
    again = materialize_bias_references(
        system_obj=system,
        protein_bundle=protein,
        ligand_bundle=ligand,
        output_dir=temp_dir / "results",
        release_cutoff="2023-06-01",
        ligand_similarity_threshold=0.35,
        selected_chains=None,
        query_cache_path=temp_dir / "cache",
    )
    assert again.reused is True
    assert calls == ["AAAA", "DDDD"]


def test_whole_policy_and_shared_component_cache_across_work_directories(
    monkeypatch, temp_dir
):
    protein = validate_bias_database_bundle(_bundle(temp_dir / "protein", "protein"), "protein")
    ligand = validate_bias_database_bundle(_bundle(temp_dir / "ligand", "ligand"), "ligand")
    calls = []

    def fake_mmseqs(binary, fasta, target, workers, tmp, max_seqs):
        sequence = "".join(line for line in Path(fasta).read_text().splitlines() if not line.startswith(">"))
        calls.append(sequence)
        return pd.DataFrame(
            [
                {"target": "1AAA_A", "pident": 90.0, "tseq": "AAAA"},
                {"target": "2BBB:B", "pident": 30.0, "tseq": "BBBB"},
            ]
        )

    monkeypatch.setattr("cofolder.tools.build_bias_training_data._run_mmseqs", fake_mmseqs)
    monkeypatch.setattr("cofolder.modules.analytics.bias_training._resolve_mmseqs_bin", lambda: "mmseqs")
    monkeypatch.setattr(bias_database, "_mmseqs_version", lambda binary: "test-version")
    monkeypatch.setattr(bias_database, "_run_mmseqs_backfill", lambda *args: pd.DataFrame())
    system = System(system={"sequences": [
        {"protein": {"id": "A", "sequence": "AAAA"}},
        {"protein": {"id": "C", "sequence": "CCCC"}},
        {"ligand": {"id": "B", "smiles": "CCO"}},
        {"ligand": {"id": "D", "smiles": "CCN"}},
    ]})
    cache = temp_dir / "shared-cache"
    dated = materialize_bias_references(
        system_obj=system, protein_bundle=protein, ligand_bundle=ligand,
        output_dir=temp_dir / "dated", release_cutoff="2023-06-01",
        ligand_similarity_threshold=0.35, selected_chains=None, query_cache_path=cache,
    )
    whole = materialize_bias_references(
        system_obj=system, protein_bundle=protein, ligand_bundle=ligand,
        output_dir=temp_dir / "whole", release_cutoff="whole",
        ligand_similarity_threshold=0.35, selected_chains=None, query_cache_path=cache,
    )
    assert calls == ["AAAA", "CCCC"]
    assert set(pd.read_csv(dated.protein_path)["pdb_id"]) == {"1AAA"}
    assert set(pd.read_csv(whole.protein_path)["pdb_id"]) == {"1AAA", "2BBB"}
    manifest = json.loads(whole.manifest_path.read_text())
    assert manifest["request"]["release_policy"]["mode"] == "whole"
    assert {event["status"] for event in manifest["cache"]["protein"]} == {"hit"}
    assert parse_bias_release_policy("whole").includes("2999-01-01")


def test_component_thresholds_filter_after_shared_raw_cache(monkeypatch, temp_dir):
    protein = validate_bias_database_bundle(
        _bundle(temp_dir / "protein", "protein"), "protein"
    )
    ligand = validate_bias_database_bundle(
        _bundle(temp_dir / "ligand", "ligand"), "ligand"
    )
    calls = []

    def fake_mmseqs(*args):
        calls.append(1)
        return pd.DataFrame(
            [
                {"target": "1AAA_A", "pident": 40.0, "tseq": "AAAA"},
                {"target": "2BBB_A", "pident": 60.0, "tseq": "BBBB"},
            ]
        )

    monkeypatch.setattr(
        "cofolder.tools.build_bias_training_data._run_mmseqs",
        fake_mmseqs,
    )
    monkeypatch.setattr(
        "cofolder.modules.analytics.bias_training._resolve_mmseqs_bin",
        lambda: "mmseqs",
    )
    monkeypatch.setattr(bias_database, "_mmseqs_version", lambda binary: "fixture")
    cache = temp_dir / "cache"
    protein_system = System(
        system={"sequences": [{"protein": {"id": "A", "sequence": "AAAA"}}]}
    )
    low = materialize_bias_references(
        system_obj=protein_system,
        protein_bundle=protein,
        ligand_bundle=None,
        output_dir=temp_dir / "protein-low",
        release_cutoff="whole",
        protein_similarity_threshold=0.25,
        selected_chains=None,
        query_cache_path=cache,
    )
    high = materialize_bias_references(
        system_obj=protein_system,
        protein_bundle=protein,
        ligand_bundle=None,
        output_dir=temp_dir / "protein-high",
        release_cutoff="whole",
        protein_similarity_threshold=0.5,
        selected_chains=None,
        query_cache_path=cache,
    )
    assert calls == [1]
    assert set(pd.read_csv(low.protein_path)["pdb_id"]) == {"1AAA", "2BBB"}
    assert set(pd.read_csv(high.protein_path)["pdb_id"]) == {"2BBB"}

    ligand_system = System(
        system={"sequences": [{"ligand": {"id": "B", "smiles": "CCO"}}]}
    )
    ligand_low = materialize_bias_references(
        system_obj=ligand_system,
        protein_bundle=None,
        ligand_bundle=ligand,
        output_dir=temp_dir / "ligand-low",
        release_cutoff="whole",
        ligand_similarity_threshold=0.35,
        selected_chains=None,
        query_cache_path=cache,
    )
    ligand_high = materialize_bias_references(
        system_obj=ligand_system,
        protein_bundle=None,
        ligand_bundle=ligand,
        output_dir=temp_dir / "ligand-high",
        release_cutoff="whole",
        ligand_similarity_threshold=1.0,
        selected_chains=None,
        query_cache_path=cache,
    )
    assert len(pd.read_csv(ligand_low.ligand_path)) == 2
    assert pd.read_csv(ligand_high.ligand_path).empty
    high_manifest = json.loads(ligand_high.manifest_path.read_text())
    assert high_manifest["cache"]["ligand"][0]["status"] == "hit"


def test_second_system_hits_shared_components_and_misses_only_changed_components(
    monkeypatch, temp_dir
):
    protein = validate_bias_database_bundle(
        _bundle(temp_dir / "protein", "protein"), "protein"
    )
    ligand = validate_bias_database_bundle(
        _bundle(temp_dir / "ligand", "ligand"), "ligand"
    )
    calls = []

    def fake_mmseqs(binary, fasta, *args):
        sequence = "".join(
            line
            for line in Path(fasta).read_text().splitlines()
            if not line.startswith(">")
        )
        calls.append(sequence)
        return pd.DataFrame(
            [{"target": "1AAA_A", "pident": 80.0, "tseq": sequence}]
        )

    monkeypatch.setattr(
        "cofolder.tools.build_bias_training_data._run_mmseqs",
        fake_mmseqs,
    )
    monkeypatch.setattr(
        "cofolder.modules.analytics.bias_training._resolve_mmseqs_bin",
        lambda: "mmseqs",
    )
    monkeypatch.setattr(bias_database, "_mmseqs_version", lambda binary: "fixture")
    monkeypatch.setattr(
        bias_database, "_run_mmseqs_backfill", lambda *args: pd.DataFrame()
    )
    shared = System(
        system={
            "sequences": [
                {"protein": {"id": "A", "sequence": "AAAA"}},
                {"ligand": {"id": "B", "smiles": "CCO"}},
            ]
        }
    )
    changed = System(
        system={
            "sequences": [
                {"protein": {"id": "A", "sequence": "AAAA"}},
                {"protein": {"id": "C", "sequence": "CCCC"}},
                {"ligand": {"id": "B", "smiles": "CCO"}},
                {"ligand": {"id": "D", "smiles": "CCN"}},
            ]
        }
    )
    cache = temp_dir / "shared-cache"
    materialize_bias_references(
        system_obj=shared,
        protein_bundle=protein,
        ligand_bundle=ligand,
        output_dir=temp_dir / "first-system",
        release_cutoff="whole",
        selected_chains=None,
        query_cache_path=cache,
    )
    second = materialize_bias_references(
        system_obj=changed,
        protein_bundle=protein,
        ligand_bundle=ligand,
        output_dir=temp_dir / "second-system",
        release_cutoff="whole",
        selected_chains=None,
        query_cache_path=cache,
    )
    assert calls == ["AAAA", "CCCC"]
    manifest = json.loads(second.manifest_path.read_text())
    assert [event["status"] for event in manifest["cache"]["protein"]] == [
        "hit",
        "miss",
    ]
    assert [event["status"] for event in manifest["cache"]["ligand"]] == [
        "hit",
        "miss",
    ]


def test_mmseqs_backfill_adapter_and_plot_scale(monkeypatch, temp_dir):
    commands = []

    def fake_run(command, **kwargs):
        commands.append(command)
        if command[1] == "convertalis":
            Path(command[5]).write_text("target_0\t12.5\n", encoding="utf-8")
        return Mock(stdout="", stderr="")

    monkeypatch.setattr(bias_database.subprocess, "run", fake_run)
    result = bias_database._run_mmseqs_backfill("mmseqs", "QUERYSEQ", ["TARGETSEQ"])
    assert [command[1] for command in commands] == ["createdb", "createdb", "search", "convertalis"]
    assert result["sequence_similarity"].tolist() == [12.5]
    assert result["sequence_similarity_pairwise"].isna().all()
    assert result["sequence_similarity_method"].tolist() == ["mmseqs_pident"]
    assert result["sequence_similarity"].iloc[0] / 100.0 == 0.125 < 0.25
    assert commands[-1][-1] == "target,pident"


def test_materialization_backfills_ligand_only_pdb_offline_and_reuses_cache(
    monkeypatch, temp_dir
):
    protein = validate_bias_database_bundle(_bundle(temp_dir / "protein", "protein"), "protein")
    ligand = validate_bias_database_bundle(_bundle(temp_dir / "ligand", "ligand"), "ligand")
    bulk_calls = []
    backfill_calls = []

    def bulk(*args):
        bulk_calls.append(1)
        return pd.DataFrame([{"target": "1AAA_A", "pident": 90.0, "tseq": "AAAA"}])

    def backfill(binary, query, targets):
        backfill_calls.append((query, tuple(targets)))
        return pd.DataFrame([{
            "sequence": targets[0], "sequence_similarity": 12.5,
            "sequence_similarity_pairwise": pd.NA,
            "sequence_similarity_method": "mmseqs_pident",
        }])

    monkeypatch.setattr("cofolder.tools.build_bias_training_data._run_mmseqs", bulk)
    monkeypatch.setattr("cofolder.modules.analytics.bias_training._resolve_mmseqs_bin", lambda: "mmseqs")
    monkeypatch.setattr(bias_database, "_mmseqs_version", lambda binary: "test-version")
    monkeypatch.setattr(bias_database, "_run_mmseqs_backfill", backfill)
    monkeypatch.setattr(
        "cofolder.modules.analytics.bias._enrichment._entry_fasta_sequences",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("network fallback must not run")),
    )
    system = System(system={"sequences": [
        {"protein": {"id": "A", "sequence": "AAAA"}},
        {"ligand": {"id": "B", "smiles": "CCO"}},
    ]})
    cache = temp_dir / "cache"
    first = materialize_bias_references(
        system_obj=system, protein_bundle=protein, ligand_bundle=ligand,
        output_dir=temp_dir / "first", release_cutoff="whole",
        ligand_similarity_threshold=0.35, selected_chains=None, query_cache_path=cache,
    )
    second = materialize_bias_references(
        system_obj=system, protein_bundle=protein, ligand_bundle=ligand,
        output_dir=temp_dir / "second", release_cutoff="whole",
        ligand_similarity_threshold=0.35, selected_chains=None, query_cache_path=cache,
    )
    rows = pd.read_csv(first.protein_path)
    recovered = rows[rows["pdb_id"] == "2BBB"].iloc[0]
    assert recovered["sequence_similarity"] == 12.5
    assert recovered["sequence_similarity_method"] == "mmseqs_pident"
    assert pd.isna(recovered["sequence_similarity_pairwise"])
    assert bulk_calls == [1]
    assert backfill_calls == [("AAAA", ("BBBB",))]
    manifest = json.loads(second.manifest_path.read_text())
    assert manifest["backfill"][0]["status"] == "hit"


@pytest.mark.parametrize(
    ("similarity", "warns"),
    [(24.9, False), (25.0, True), (25.1, True)],
)
def test_backfill_boundary_warning_uses_configured_protein_threshold(
    monkeypatch, temp_dir, caplog, similarity, warns
):
    protein = validate_bias_database_bundle(
        _bundle(temp_dir / "protein", "protein"), "protein"
    )
    ligand = validate_bias_database_bundle(
        _bundle(temp_dir / "ligand", "ligand"), "ligand"
    )
    monkeypatch.setattr(
        "cofolder.tools.build_bias_training_data._run_mmseqs",
        lambda *args: pd.DataFrame(
            [{"target": "1AAA_A", "pident": 90.0, "tseq": "AAAA"}]
        ),
    )
    monkeypatch.setattr(
        "cofolder.modules.analytics.bias_training._resolve_mmseqs_bin",
        lambda: "mmseqs",
    )
    monkeypatch.setattr(bias_database, "_mmseqs_version", lambda binary: "fixture")
    monkeypatch.setattr(
        bias_database,
        "_run_mmseqs_backfill",
        lambda *args: pd.DataFrame(
            [
                {
                    "sequence": "BBBB",
                    "sequence_similarity": similarity,
                    "sequence_similarity_pairwise": pd.NA,
                    "sequence_similarity_method": "mmseqs_pident",
                }
            ]
        ),
    )
    system = System(
        system={
            "sequences": [
                {"protein": {"id": "A", "sequence": "AAAA"}},
                {"ligand": {"id": "B", "smiles": "CCO"}},
            ]
        }
    )
    with caplog.at_level("WARNING"):
        artifacts = materialize_bias_references(
            system_obj=system,
            protein_bundle=protein,
            ligand_bundle=ligand,
            output_dir=temp_dir / "result",
            release_cutoff="whole",
            protein_similarity_threshold=0.25,
            selected_chains=None,
            query_cache_path=temp_dir / "cache",
        )
    recovered = pd.read_csv(artifacts.protein_path)
    recovered = recovered[recovered["pdb_id"] == "2BBB"].iloc[0]
    assert recovered["sequence_similarity"] == similarity
    assert recovered["sequence_similarity_method"] == "mmseqs_pident"
    assert ("reached or exceeded 25.0%" in caplog.text) is warns


def test_cache_rebuilds_corruption_and_serializes_concurrent_writers(temp_dir):
    root = temp_dir / "cache"
    identity = {"bundle": "fixture", "query": "same"}
    key = bias_database._cache_key("protein", identity)
    calls = []

    def builder():
        calls.append(1)
        time.sleep(0.02)
        return pd.DataFrame([{"value": 1}])

    with ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(lambda _: bias_database._cached_frame(
            root=root, kind="protein", key=key, identity=identity, builder=builder
        ), range(2)))
    assert calls == [1]
    assert sorted(status for _, status in results) == ["hit", "miss"]
    payload = root / "protein" / key[:2] / key / "payload.csv.gz"
    payload.write_bytes(b"corrupt")
    frame, status = bias_database._cached_frame(
        root=root, kind="protein", key=key, identity=identity, builder=builder
    )
    assert status == "rebuilt"
    assert frame["value"].tolist() == [1]
    assert calls == [1, 1]
