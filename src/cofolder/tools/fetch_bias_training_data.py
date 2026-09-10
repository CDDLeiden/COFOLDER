"""Fetch CCD data and optional MMseqs2 databases for bias preparation."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from gzip import decompress
import os
from pathlib import Path
import shutil
import subprocess
import sys
from urllib.request import urlopen


CCD_GZ_URL = "https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz"


def _resolve_mmseqs_bin(mmseqs_bin: str | None = None) -> str | None:
    candidates: list[str] = []
    if mmseqs_bin:
        candidates.append(str(mmseqs_bin))
    candidates.extend(
        [
            str(Path.home() / ".cofolder/vendor/mmseqs/bin/mmseqs"),
            str(Path(__file__).resolve().parents[3] / "vendor/mmseqs/bin/mmseqs"),
            str(Path(__file__).resolve().parents[3] / "vendor/mmseqs/mmseqs"),
        ]
    )
    from_env = shutil.which("mmseqs")
    if from_env:
        candidates.append(from_env)
    explicit_env = os.environ.get("COFOLDER_MMSEQS_BIN")
    if explicit_env:
        candidates.insert(0, explicit_env)

    for candidate in candidates:
        if not candidate:
            continue
        path = Path(candidate).expanduser()
        if path.exists() and path.is_file():
            return str(path)
        found = shutil.which(candidate)
        if found:
            return found
    return None


def _fetch_bytes(url: str, timeout: int) -> bytes:
    with urlopen(url, timeout=timeout) as response:
        return response.read()


def _extract_ccd_ids(components_cif: Path) -> list[str]:
    ids: list[str] = []
    with components_cif.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
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
    mmseqs_command = _resolve_mmseqs_bin(mmseqs_bin)
    if mmseqs_command is None:
        raise RuntimeError(
            "mmseqs binary not found in PATH. Install MMseqs2, pass --mmseqs_bin, "
            "or disable with --skip_mmseqs."
        )

    tmp_dir.mkdir(parents=True, exist_ok=True)
    output_db_path.parent.mkdir(parents=True, exist_ok=True)
    command = [
        mmseqs_command,
        "databases",
        db_name,
        str(output_db_path),
        str(tmp_dir),
    ]
    print(f"[info] running: {' '.join(command)}")
    subprocess.run(command, check=True)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cofolder-tools fetch-bias-training-data", description=__doc__
    )
    parser.add_argument("--output_root", type=Path, required=True, help="Output directory for CCD files.")
    parser.add_argument("--timeout", type=int, default=60, help="HTTP timeout in seconds.")
    parser.add_argument("--overwrite", action="store_true", help="Re-download/rewrite files if they already exist.")
    parser.add_argument("--skip_mmseqs", action="store_true", help="Skip MMseqs2 database download step.")
    parser.add_argument("--mmseqs_db_name", type=str, default="PDB", help="MMseqs2 database name to download (e.g. PDB, UniRef50, UniRef90).")
    parser.add_argument("--mmseqs_bin", type=str, default=None, help="Path to mmseqs executable (optional). If omitted, resolve from PATH.")
    parser.add_argument("--mmseqs_output_db", type=Path, default=None, help="Output DB path for `mmseqs databases`. Defaults to <output_root>/mmseqs/<db_name_lower>/db.")
    parser.add_argument("--mmseqs_tmp_dir", type=Path, default=None, help="Temporary directory for MMseqs2. Defaults to <output_root>/mmseqs/tmp.")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
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
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(
                f"mmseqs databases failed for {mmseqs_db_name} (exit {exc.returncode})"
            ) from exc

    print(f"[done] wrote {ccd_ids_txt}")
    print(f"[done] ccd_entry_count={len(ccd_ids)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
