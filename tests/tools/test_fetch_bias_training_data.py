from __future__ import annotations

from gzip import compress
import json

import pandas as pd

from cofolder.tools import fetch_bias_training_data as fetcher


def test_fetches_ccd_and_uses_existing_mmseqs_defaults(monkeypatch, temp_dir):
    payload = compress(b"data_atp\n#\ndata_HEM\n")
    monkeypatch.setattr(fetcher, "_fetch_bytes", lambda url, timeout: payload)
    captured = {}

    def fake_mmseqs(db_name, output_db_path, tmp_dir, mmseqs_bin=None):
        captured.update(
            db_name=db_name,
            output_db_path=output_db_path,
            tmp_dir=tmp_dir,
            mmseqs_bin=mmseqs_bin,
        )
        output_db_path.parent.mkdir(parents=True, exist_ok=True)
        output_db_path.write_text("database", encoding="utf-8")

    monkeypatch.setattr(fetcher, "_run_mmseqs_databases", fake_mmseqs)
    monkeypatch.setattr(
        fetcher,
        "_extract_mmseqs_sequence_index",
        lambda *args: pd.DataFrame(
            [{"pdb_id": "1ABC", "target_id": "1ABC_A", "sequence": "AAAA"}]
        ),
    )
    monkeypatch.setattr(
        fetcher,
        "_collect_pdb_metadata",
        lambda timeout, workers, checkpoint_path=None: (
            pd.DataFrame([{"pdb_id": "1ABC", "release_date": "2020-01-01"}]),
            pd.DataFrame([{"pdb_id": "1ABC", "release_date": "2020-01-01", "ligand_id": "ATP", "smiles": "CCO"}]),
        ),
    )
    output_root = temp_dir / "training data"

    assert fetcher.main(["--output_root", str(output_root)]) == 0

    assert (output_root / "ccd/components.cif").read_bytes() == b"data_atp\n#\ndata_HEM\n"
    assert (output_root / "ccd/ccd_ids.txt").read_text(encoding="utf-8") == "ATP\nHEM\n"
    assert captured == {
        "db_name": "PDB",
        "output_db_path": output_root / "mmseqs/pdb/db",
        "tmp_dir": output_root / "mmseqs/tmp",
        "mmseqs_bin": None,
    }
    assert (output_root / "protein/manifest.json").is_file()
    assert (output_root / "ligand/manifest.json").is_file()
    assert fetcher.main(["--output_root", str(output_root), "--validate_only"]) == 0
    monkeypatch.setattr(
        fetcher,
        "_fetch_bytes",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("unexpected download")),
    )
    assert fetcher.main(["--output_root", str(output_root)]) == 0


def test_skip_mmseqs_does_not_resolve_or_run_it(monkeypatch, temp_dir):
    payload = compress(b"data_atp\n")
    monkeypatch.setattr(fetcher, "_fetch_bytes", lambda url, timeout: payload)
    monkeypatch.setattr(
        fetcher,
        "_run_mmseqs_databases",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("unexpected")),
    )
    monkeypatch.setattr(
        fetcher,
        "_collect_pdb_metadata",
        lambda timeout, workers, checkpoint_path=None: (
            pd.DataFrame(columns=["pdb_id", "release_date"]),
            pd.DataFrame(columns=["pdb_id", "release_date", "ligand_id", "smiles"]),
        ),
    )

    assert fetcher.main(
        ["--output_root", str(temp_dir / "data"), "--skip_mmseqs", "--skip_ligand"]
    ) == 0


def test_fetch_json_retries_transient_failures(monkeypatch):
    attempts = iter([OSError("temporary"), OSError("temporary"), b'{"ok": true}'])

    def fake_fetch(*args, **kwargs):
        result = next(attempts)
        if isinstance(result, Exception):
            raise result
        return result

    monkeypatch.setattr(fetcher, "_fetch_bytes", fake_fetch)
    monkeypatch.setattr(fetcher.time, "sleep", lambda delay: None)

    assert fetcher._fetch_json("https://example.test/data", 1) == {"ok": True}


def test_metadata_collection_resumes_from_checkpoint(monkeypatch, temp_dir):
    checkpoint = temp_dir / "metadata.jsonl"
    saved = {
        "protein": {"pdb_id": "1ABC", "release_date": "2020-01-01"},
        "ligands": [],
    }
    checkpoint.write_text(json.dumps(saved) + "\n", encoding="utf-8")
    calls = []
    monkeypatch.setattr(fetcher, "_fetch_json", lambda url, timeout: ["1ABC", "2DEF"])

    def fake_entry(pdb_id, timeout):
        calls.append(pdb_id)
        return {"pdb_id": pdb_id, "release_date": "2021-02-03"}, [
            {
                "pdb_id": pdb_id,
                "release_date": "2021-02-03",
                "ligand_id": "ATP",
                "smiles": "CCO",
            }
        ]

    monkeypatch.setattr(fetcher, "_entry_record", fake_entry)
    proteins, ligands = fetcher._collect_pdb_metadata(1, 2, checkpoint)

    assert calls == ["2DEF"]
    assert proteins["pdb_id"].tolist() == ["1ABC", "2DEF"]
    assert ligands["pdb_id"].tolist() == ["2DEF"]
    assert len(checkpoint.read_text(encoding="utf-8").splitlines()) == 2


def test_packages_existing_custom_database_sources(monkeypatch, temp_dir):
    output_root = temp_dir / "custom"
    imported_db = temp_dir / "imports" / "expanded" / "db"
    imported_db.parent.mkdir(parents=True)
    imported_db.write_text("expanded database", encoding="utf-8")
    protein_metadata = temp_dir / "imports" / "releases.csv"
    ligand_table = temp_dir / "imports" / "ligands.csv"
    pd.DataFrame(
        [{"pdb_id": "9XYZ", "release_date": "2025-01-02"}]
    ).to_csv(protein_metadata, index=False)
    pd.DataFrame(
        [
            {
                "pdb_id": "9XYZ",
                "release_date": "2025-01-02",
                "ligand_id": "NEW",
                "smiles": "CCN",
            }
        ]
    ).to_csv(ligand_table, index=False)
    monkeypatch.setattr(fetcher, "_fetch_bytes", lambda url, timeout: compress(b"data_NEW\n"))
    monkeypatch.setattr(
        fetcher,
        "_collect_pdb_metadata",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("unexpected RCSB fetch")),
    )
    monkeypatch.setattr(
        fetcher,
        "_run_mmseqs_databases",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("unexpected MMseqs fetch")),
    )
    monkeypatch.setattr(
        fetcher,
        "_extract_mmseqs_sequence_index",
        lambda *args: pd.DataFrame(
            [{"pdb_id": "9XYZ", "target_id": "9XYZ_A", "sequence": "AAAA"}]
        ),
    )

    assert fetcher.main(
        [
            "--output_root", str(output_root),
            "--mmseqs_output_db", str(imported_db),
            "--protein_release_metadata", str(protein_metadata),
            "--ligand_occurrence_table", str(ligand_table),
            "--overwrite",
        ]
    ) == 0

    assert (output_root / "protein/mmseqs/db").read_text(encoding="utf-8") == "expanded database"
    assert pd.read_csv(output_root / "ligand/pdb_ligands.csv.gz")["ligand_id"].tolist() == ["NEW"]
