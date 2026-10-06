"""Prepare validated protein and ligand database bundles for bias assessment."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from gzip import decompress
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
from urllib.request import urlopen

import pandas as pd

from cofolder import __version__
from cofolder.modules.analytics.bias_database import (
    BIAS_DATABASE_SCHEMA_VERSION,
    LIGAND_TABLE_NAME,
    PROTEIN_METADATA_NAME,
    PROTEIN_SEQUENCE_INDEX_NAME,
    normalize_pdb_id,
    validate_bias_database_bundle,
)
from cofolder.modules.utils.executables import resolve_mmseqs_executable


CCD_GZ_URL = "https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz"
CURRENT_ENTRY_IDS_URL = "https://data.rcsb.org/rest/v1/holdings/current/entry_ids"
CORE_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
CORE_NONPOLY_URL = "https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{pdb_id}/{entity_id}"


def _resolve_mmseqs_bin(mmseqs_bin: str | None = None) -> str | None:
    return resolve_mmseqs_executable(mmseqs_bin).path


def _fetch_bytes(url: str, timeout: int) -> bytes:
    with urlopen(url, timeout=timeout) as response:
        return response.read()


def _fetch_json(url: str, timeout: int, attempts: int = 4):
    last_error: Exception | None = None
    for attempt in range(attempts):
        try:
            return json.loads(_fetch_bytes(url, timeout).decode("utf-8"))
        except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
            last_error = exc
            if attempt + 1 < attempts:
                time.sleep(0.25 * (2**attempt))
    raise RuntimeError(f"Failed to retrieve valid JSON from {url}: {last_error}") from last_error


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _entry_record(pdb_id: str, timeout: int) -> tuple[dict, list[dict]]:
    entry = _fetch_json(CORE_ENTRY_URL.format(pdb_id=pdb_id), timeout)
    release = str(entry.get("rcsb_accession_info", {}).get("initial_release_date", ""))[:10]
    entity_ids = entry.get("rcsb_entry_container_identifiers", {}).get(
        "non_polymer_entity_ids", []
    )
    ligand_rows: list[dict] = []
    for entity_id in entity_ids:
        payload = _fetch_json(
            CORE_NONPOLY_URL.format(pdb_id=pdb_id, entity_id=entity_id), timeout
        )
        identifiers = payload.get("rcsb_nonpolymer_entity_container_identifiers", {})
        ligand_id = identifiers.get("nonpolymer_comp_id")
        descriptors = payload.get("rcsb_nonpolymer_entity", {}).get(
            "pdbx_description", ""
        )
        descriptor_rows = payload.get("pdbx_chem_comp_descriptor", [])
        if isinstance(descriptor_rows, dict):
            descriptor_rows = [descriptor_rows]
        smiles = ""
        for descriptor in descriptor_rows:
            if "SMILES" in str(descriptor.get("type", "")).upper():
                smiles = str(descriptor.get("descriptor", "")).strip()
                if smiles:
                    break
        if ligand_id and release:
            ligand_rows.append(
                {
                    "pdb_id": pdb_id,
                    "release_date": release,
                    "ligand_id": str(ligand_id).upper(),
                    "smiles": smiles,
                    "description": descriptors,
                }
            )
    return {"pdb_id": pdb_id, "release_date": release}, ligand_rows


def _load_checkpoint(path: Path | None) -> dict[str, tuple[dict, list[dict]]]:
    records: dict[str, tuple[dict, list[dict]]] = {}
    if path is None or not path.is_file():
        return records
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                payload = json.loads(line)
                protein = payload["protein"]
                ligands = payload["ligands"]
                pdb_id = str(protein["pdb_id"]).upper()
            except (KeyError, TypeError, json.JSONDecodeError) as exc:
                raise RuntimeError(
                    f"Invalid RCSB checkpoint {path} at line {line_number}. "
                    "Remove it and rerun preparation."
                ) from exc
            records[pdb_id] = (protein, ligands)
    return records


def _collect_pdb_metadata(
    timeout: int, workers: int, checkpoint_path: Path | None = None
) -> tuple[pd.DataFrame, pd.DataFrame]:
    identifiers = _fetch_json(CURRENT_ENTRY_IDS_URL, timeout)
    if not isinstance(identifiers, list):
        raise RuntimeError("RCSB holdings endpoint did not return an entry-ID list.")
    identifiers = [str(identifier).upper() for identifier in identifiers]
    records = _load_checkpoint(checkpoint_path)
    remaining = [identifier for identifier in identifiers if identifier not in records]
    checkpoint = None
    if checkpoint_path is not None:
        checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
        checkpoint = checkpoint_path.open("a", encoding="utf-8")
    with ThreadPoolExecutor(max_workers=max(1, workers)) as executor:
        try:
            batch_size = max(1, workers) * 4
            for offset in range(0, len(remaining), batch_size):
                batch = remaining[offset : offset + batch_size]
                futures = {
                    executor.submit(_entry_record, identifier, timeout): identifier
                    for identifier in batch
                }
                for future in as_completed(futures):
                    protein, ligands = future.result()
                    records[str(protein["pdb_id"]).upper()] = (protein, ligands)
                    if checkpoint is not None:
                        checkpoint.write(
                            json.dumps({"protein": protein, "ligands": ligands}) + "\n"
                        )
                        checkpoint.flush()
        finally:
            if checkpoint is not None:
                checkpoint.close()
    protein_rows = [records[pdb_id][0] for pdb_id in identifiers if pdb_id in records]
    ligand_rows = [
        ligand
        for pdb_id in identifiers
        if pdb_id in records
        for ligand in records[pdb_id][1]
    ]
    protein_rows = [row for row in protein_rows if row["release_date"]]
    return (
        pd.DataFrame(protein_rows, columns=["pdb_id", "release_date"]),
        pd.DataFrame(
            ligand_rows,
            columns=["pdb_id", "release_date", "ligand_id", "smiles", "description"],
        ),
    )


def _write_manifest(
    root: Path, *, kind: str, source: str, files: list[Path], record_count: int
) -> None:
    inventory = {
        str(path.relative_to(root)): _sha256(path)
        for path in sorted(files)
        if path.is_file()
    }
    manifest = {
        "schema_version": BIAS_DATABASE_SCHEMA_VERSION,
        "kind": kind,
        "source": source,
        "source_snapshot": datetime.now(timezone.utc).isoformat(),
        "created_at": datetime.now(timezone.utc).isoformat(),
        "record_count": int(record_count),
        "files": inventory,
        "tools": {"cofolder": __version__},
    }
    (root / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def _publish_bundle(staged: Path, destination: Path, overwrite: bool) -> None:
    if destination.exists():
        if not overwrite:
            validate_bias_database_bundle(destination, staged.name)
            print(f"[info] Reusing existing {staged.name} bias database: {destination}")
            return
        backup = destination.with_name(f".{destination.name}.previous")
        if backup.exists():
            shutil.rmtree(backup)
        destination.replace(backup)
        try:
            staged.replace(destination)
        except Exception:
            backup.replace(destination)
            raise
        shutil.rmtree(backup)
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    staged.replace(destination)


def _needs_preparation(destination: Path, kind: str, overwrite: bool) -> bool:
    if overwrite or not destination.exists():
        return True
    validate_bias_database_bundle(destination, kind)
    print(f"[info] Reusing existing {kind} bias database: {destination}")
    return False


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
        diagnostics = resolve_mmseqs_executable(mmseqs_bin).diagnostic_text()
        raise RuntimeError(
            "mmseqs binary not found in PATH. Install MMseqs2, pass --mmseqs_bin, "
            f"or disable with --skip_mmseqs. Attempted: {diagnostics}."
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


def _extract_mmseqs_sequence_index(
    mmseqs_db: Path, mmseqs_bin: str | None, temporary_root: Path
) -> pd.DataFrame:
    """Export a relocatable PDB sequence index from a prepared MMseqs database."""
    executable = _resolve_mmseqs_bin(mmseqs_bin)
    if executable is None:
        diagnostics = resolve_mmseqs_executable(mmseqs_bin).diagnostic_text()
        raise RuntimeError(
            "MMseqs2 is required to build the offline PDB sequence index. "
            f"Attempted: {diagnostics}."
        )
    fasta = temporary_root / "pdb_sequences.fasta"
    subprocess.run(
        [executable, "convert2fasta", str(mmseqs_db), str(fasta)], check=True
    )
    rows: list[dict[str, str]] = []
    header: str | None = None
    sequence: list[str] = []

    def append_record() -> None:
        if header is None:
            return
        normalized = "".join(sequence).strip().upper()
        pdb_id = normalize_pdb_id(header)
        if pdb_id and normalized:
            rows.append(
                {"pdb_id": pdb_id, "target_id": header.split()[0], "sequence": normalized}
            )

    with fasta.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                append_record()
                header = line[1:].strip()
                sequence = []
            else:
                sequence.append(line.strip())
    append_record()
    result = pd.DataFrame(rows, columns=["pdb_id", "target_id", "sequence"])
    return result.drop_duplicates().sort_values(["pdb_id", "target_id"])


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cofolder-tools fetch-bias-training-data", description=__doc__
    )
    parser.add_argument("--output_root", type=Path, required=True, help="Output root for database bundles and preparation files.")
    parser.add_argument("--timeout", type=int, default=60, help="HTTP timeout in seconds.")
    parser.add_argument("--overwrite", action="store_true", help="Re-download/rewrite files if they already exist.")
    parser.add_argument("--skip_mmseqs", action="store_true", help="Skip MMseqs2 database download step.")
    parser.add_argument("--mmseqs_db_name", type=str, default="PDB", help="MMseqs2 database name to download (e.g. PDB, UniRef50, UniRef90).")
    parser.add_argument(
        "--mmseqs_bin",
        type=str,
        default=None,
        help=(
            "MMseqs executable path or command name (highest precedence; otherwise "
            "use COFOLDER_MMSEQS_BIN, managed vendor paths, then PATH)."
        ),
    )
    parser.add_argument("--mmseqs_output_db", type=Path, default=None, help="Output DB path for `mmseqs databases`. Defaults to <output_root>/mmseqs/<db_name_lower>/db.")
    parser.add_argument("--mmseqs_tmp_dir", type=Path, default=None, help="Temporary directory for MMseqs2. Defaults to <output_root>/mmseqs/tmp.")
    parser.add_argument(
        "--protein_release_metadata",
        type=Path,
        default=None,
        help="Existing CSV/CSV.gz PDB release mapping to package instead of fetching it.",
    )
    parser.add_argument(
        "--protein_sequence_table",
        type=Path,
        default=None,
        help="Existing pdb_id/sequence CSV to package instead of exporting the MMseqs target.",
    )
    parser.add_argument(
        "--ligand_occurrence_table",
        type=Path,
        default=None,
        help="Existing conforming ligand CSV/CSV.gz to package instead of fetching it.",
    )
    parser.add_argument(
        "--bias_training_data_protein_path",
        type=Path,
        default=None,
        help="Protein source-bundle destination (default: <output_root>/protein).",
    )
    parser.add_argument(
        "--bias_training_data_ligand_path",
        type=Path,
        default=None,
        help="Ligand source-bundle destination (default: <output_root>/ligand).",
    )
    parser.add_argument("--skip_ligand", action="store_true", help="Skip the offline PDB ligand bundle.")
    parser.add_argument("--workers", type=int, default=16, help="Concurrent RCSB metadata requests.")
    parser.add_argument(
        "--validate_only",
        action="store_true",
        help="Validate existing protein and ligand bundles without downloading data.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    output_root = args.output_root
    protein_root = args.bias_training_data_protein_path or output_root / "protein"
    ligand_root = args.bias_training_data_ligand_path or output_root / "ligand"
    if args.validate_only:
        if not args.skip_mmseqs:
            validate_bias_database_bundle(protein_root, "protein")
        if not args.skip_ligand:
            validate_bias_database_bundle(ligand_root, "ligand")
        print("[done] bias database bundles are valid")
        return 0
    prepare_protein = not args.skip_mmseqs and _needs_preparation(
        protein_root, "protein", args.overwrite
    )
    prepare_ligand = not args.skip_ligand and _needs_preparation(
        ligand_root, "ligand", args.overwrite
    )
    if not prepare_protein and not prepare_ligand:
        print("[done] requested bias database bundles already exist and are valid")
        return 0
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

    protein_metadata = (
        pd.read_csv(args.protein_release_metadata)
        if prepare_protein and args.protein_release_metadata
        else pd.DataFrame(columns=["pdb_id", "release_date"])
    )
    protein_sequences = (
        pd.read_csv(args.protein_sequence_table)
        if prepare_protein and args.protein_sequence_table
        else None
    )
    ligand_database = (
        pd.read_csv(args.ligand_occurrence_table)
        if prepare_ligand and args.ligand_occurrence_table
        else pd.DataFrame(columns=["pdb_id", "release_date", "ligand_id", "smiles"])
    )
    checkpoint_path = output_root / ".rcsb-bias-metadata.checkpoint.jsonl"
    fetch_protein_metadata = prepare_protein and not args.protein_release_metadata
    fetch_ligand_metadata = prepare_ligand and not args.ligand_occurrence_table
    if fetch_protein_metadata or fetch_ligand_metadata:
        fetched_proteins, fetched_ligands = _collect_pdb_metadata(
            args.timeout, args.workers, checkpoint_path
        )
        if fetch_protein_metadata:
            protein_metadata = fetched_proteins
        if fetch_ligand_metadata:
            ligand_database = fetched_ligands
    if prepare_ligand:
        missing = (
            ligand_database["smiles"].fillna("").astype(str).str.strip().eq("")
            if not ligand_database.empty and "smiles" in ligand_database
            else pd.Series(dtype=bool)
        )
        if prepare_ligand and missing.any():
            from cofolder.tools.build_bias_training_data import _load_ccd_smiles

            ccd = _load_ccd_smiles(components_cif)
            ccd_smiles = dict(zip(ccd["ligand_id"].astype(str), ccd["smiles"].astype(str)))
            ligand_database.loc[missing, "smiles"] = ligand_database.loc[
                missing, "ligand_id"
            ].astype(str).map(ccd_smiles)
            ligand_database = ligand_database.dropna(subset=["smiles"])

    if prepare_protein:
        mmseqs_root = output_root / "mmseqs"
        mmseqs_db_name = args.mmseqs_db_name.strip()
        mmseqs_db_key = mmseqs_db_name.replace("/", "_").lower()
        mmseqs_output_db = args.mmseqs_output_db or (mmseqs_root / mmseqs_db_key / "db")
        mmseqs_tmp_dir = args.mmseqs_tmp_dir or (mmseqs_root / "tmp")
        existing_mmseqs = any(
            mmseqs_output_db.parent.glob(f"{mmseqs_output_db.name}*")
        )
        if existing_mmseqs and (args.mmseqs_output_db is not None or not args.overwrite):
            print(f"[info] Reusing existing MMseqs database at {mmseqs_output_db}")
        else:
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

        protein_root.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=".protein-bias-", dir=protein_root.parent) as raw:
            staged = Path(raw) / "protein"
            (staged / "mmseqs").mkdir(parents=True)
            for source in sorted(mmseqs_output_db.parent.glob(f"{mmseqs_output_db.name}*")):
                if source.is_file():
                    suffix = source.name[len(mmseqs_output_db.name) :]
                    shutil.copy2(source, staged / "mmseqs" / f"db{suffix}")
            protein_metadata.to_csv(staged / PROTEIN_METADATA_NAME, index=False)
            if protein_sequences is None:
                protein_sequences = _extract_mmseqs_sequence_index(
                    mmseqs_output_db, args.mmseqs_bin, Path(raw)
                )
            missing_sequence_columns = {"pdb_id", "sequence"} - set(protein_sequences.columns)
            if missing_sequence_columns:
                raise ValueError(
                    "Protein sequence table is missing columns: "
                    f"{sorted(missing_sequence_columns)}"
                )
            protein_sequences = protein_sequences.copy()
            protein_sequences["pdb_id"] = protein_sequences["pdb_id"].map(normalize_pdb_id)
            protein_sequences.to_csv(staged / PROTEIN_SEQUENCE_INDEX_NAME, index=False)
            files = [path for path in staged.rglob("*") if path.is_file()]
            _write_manifest(
                staged,
                kind="protein",
                source=(
                    f"{mmseqs_db_name} + {args.protein_release_metadata}"
                    if args.protein_release_metadata
                    else mmseqs_db_name
                ),
                files=files,
                record_count=len(protein_metadata),
            )
            validate_bias_database_bundle(staged, "protein")
            _publish_bundle(staged, protein_root, args.overwrite)

    if prepare_ligand:
        ligand_root.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=".ligand-bias-", dir=ligand_root.parent) as raw:
            staged = Path(raw) / "ligand"
            staged.mkdir(parents=True)
            ligand_database.to_csv(staged / LIGAND_TABLE_NAME, index=False)
            _write_manifest(
                staged,
                kind="ligand",
                source=(
                    str(args.ligand_occurrence_table)
                    if args.ligand_occurrence_table
                    else "wwPDB CCD + RCSB PDB Data API"
                ),
                files=[staged / LIGAND_TABLE_NAME],
                record_count=len(ligand_database),
            )
            validate_bias_database_bundle(staged, "ligand")
            _publish_bundle(staged, ligand_root, args.overwrite)

    checkpoint_path.unlink(missing_ok=True)

    print(f"[done] wrote {ccd_ids_txt}")
    print(f"[done] ccd_entry_count={len(ccd_ids)}")
    if prepare_protein:
        print(f"[done] protein_bias_database={protein_root}")
    if prepare_ligand:
        print(f"[done] ligand_bias_database={ligand_root}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
