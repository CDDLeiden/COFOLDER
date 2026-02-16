#!/usr/bin/env python3
"""Fetch CCD data and MMseqs2 sequence databases for bias precomputation.

This script intentionally downloads only the Chemical Component Dictionary (CCD)
entries. The downloaded CCD can be used to compute ligand similarity first.
A later step can map top-hit CCD IDs back to PDB entries/structures.

Sources:
- https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz
- MMseqs2 databases module (e.g. PDB, UniRef*):
  mmseqs databases <name> <o:sequenceDB> <tmpDir>
"""

from __future__ import annotations

import argparse
from gzip import decompress
from pathlib import Path
import shutil
import subprocess
import sys
from urllib.request import urlopen

CCD_GZ_URL = "https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz"


def _fetch_bytes(url: str, timeout: int) -> bytes:
    with urlopen(url, timeout=timeout) as resp:
        return resp.read()


def _extract_ccd_ids(components_cif: Path) -> list[str]:
    """Extract CCD IDs by scanning data blocks: data_<CCD_ID>."""
    ids: list[str] = []
    with components_cif.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith("data_"):
                continue
            ccd_id = line[5:].strip().upper()
            if ccd_id:
                ids.append(ccd_id)
    return ids


def _run_mmseqs_databases(
    db_name: str,
    output_db_path: Path,
    tmp_dir: Path,
    mmseqs_bin: str | None = None,
) -> None:
    mmseqs_cmd = mmseqs_bin or shutil.which("mmseqs")
    if mmseqs_cmd is None:
        raise RuntimeError(
            "mmseqs binary not found in PATH. Install MMseqs2, pass --mmseqs_bin, or disable with --skip_mmseqs."
        )

    tmp_dir.mkdir(parents=True, exist_ok=True)
    output_db_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        mmseqs_cmd,
        "databases",
        db_name,
        str(output_db_path),
        str(tmp_dir),
    ]
    print(f"[info] running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output_root",
        type=Path,
        required=True,
        help="Output directory for CCD files.",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=60,
        help="HTTP timeout in seconds.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Re-download/rewrite files if they already exist.",
    )
    parser.add_argument(
        "--skip_mmseqs",
        action="store_true",
        help="Skip MMseqs2 database download step.",
    )
    parser.add_argument(
        "--mmseqs_db_name",
        type=str,
        default="PDB",
        help="MMseqs2 database name to download (e.g. PDB, UniRef50, UniRef90).",
    )
    parser.add_argument(
        "--mmseqs_bin",
        type=str,
        default=None,
        help="Path to mmseqs executable (optional). If omitted, resolve from PATH.",
    )
    parser.add_argument(
        "--mmseqs_output_db",
        type=Path,
        default=None,
        help="Output DB path for `mmseqs databases`. Defaults to <output_root>/mmseqs/<db_name_lower>/db.",
    )
    parser.add_argument(
        "--mmseqs_tmp_dir",
        type=Path,
        default=None,
        help="Temporary directory for MMseqs2. Defaults to <output_root>/mmseqs/tmp.",
    )
    args = parser.parse_args()

    output_root = args.output_root
    ccd_dir = output_root / "ccd"
    ccd_dir.mkdir(parents=True, exist_ok=True)

    components_gz = ccd_dir / "components.cif.gz"
    components_cif = ccd_dir / "components.cif"
    ccd_ids_txt = ccd_dir / "ccd_ids.txt"

    if components_gz.exists() and components_cif.exists() and not args.overwrite:
        print(f"[info] Reusing existing CCD files in {ccd_dir}")
    else:
        payload = _fetch_bytes(CCD_GZ_URL, timeout=args.timeout)
        components_gz.write_bytes(payload)
        components_cif.write_bytes(decompress(payload))
        print(f"[done] downloaded {components_gz}")
        print(f"[done] extracted {components_cif}")

    ccd_ids = _extract_ccd_ids(components_cif)
    ccd_ids_txt.write_text("\n".join(ccd_ids) + "\n", encoding="utf-8")

    if not args.skip_mmseqs:
        mmseqs_root = output_root / "mmseqs"
        mmseqs_db_name = args.mmseqs_db_name.strip()
        mmseqs_db_key = mmseqs_db_name.replace("/", "_").lower()
        mmseqs_output_db = args.mmseqs_output_db or (mmseqs_root / mmseqs_db_key / "db")
        mmseqs_tmp_dir = args.mmseqs_tmp_dir or (mmseqs_root / "tmp")
        try:
            _run_mmseqs_databases(
                db_name=mmseqs_db_name,
                output_db_path=mmseqs_output_db,
                tmp_dir=mmseqs_tmp_dir,
                mmseqs_bin=args.mmseqs_bin,
            )
            print(f"[done] mmseqs_db={mmseqs_db_name} at {mmseqs_output_db}")
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"mmseqs databases failed for {mmseqs_db_name} (exit {e.returncode})"
            ) from e

    print(f"[done] wrote {ccd_ids_txt}")
    print(f"[done] ccd_entry_count={len(ccd_ids)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
