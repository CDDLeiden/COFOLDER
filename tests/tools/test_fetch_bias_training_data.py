from __future__ import annotations

from gzip import compress

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

    monkeypatch.setattr(fetcher, "_run_mmseqs_databases", fake_mmseqs)
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


def test_skip_mmseqs_does_not_resolve_or_run_it(monkeypatch, temp_dir):
    payload = compress(b"data_atp\n")
    monkeypatch.setattr(fetcher, "_fetch_bytes", lambda url, timeout: payload)
    monkeypatch.setattr(
        fetcher,
        "_run_mmseqs_databases",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("unexpected")),
    )

    assert fetcher.main(
        ["--output_root", str(temp_dir / "data"), "--skip_mmseqs"]
    ) == 0
