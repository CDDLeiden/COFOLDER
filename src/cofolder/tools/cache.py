"""Explicit Boltz2 cache setup and CCD population tools."""

from __future__ import annotations

import argparse
from collections import Counter
from collections.abc import Sequence
import fcntl
import logging
from pathlib import Path
import tempfile

from cofolder.modules.utils import read, write

logger = logging.getLogger(__name__)

_BUILTIN_CCD_MARKERS = ("ALA.pkl", "GLY.pkl")


def missing_boltz2_cache_components(cache_path: str | Path) -> tuple[str, ...]:
    """Return required Boltz2 cache components that are absent or empty."""
    root = Path(cache_path).expanduser()
    missing = []
    for name in ("boltz2_conf.ckpt", "mols.tar"):
        candidate = root / name
        if not candidate.is_file() or candidate.stat().st_size == 0:
            missing.append(name)
    mols = root / "mols"
    if not mols.is_dir():
        missing.append("mols/")
    else:
        for name in _BUILTIN_CCD_MARKERS:
            candidate = mols / name
            if not candidate.is_file() or candidate.stat().st_size == 0:
                missing.append(f"mols/{name}")
    return tuple(missing)


def _download_boltz2_cache(cache_path: Path) -> None:
    with tempfile.TemporaryDirectory(prefix="cofolder-boltz-cache-") as raw:
        workspace = Path(raw)
        fasta = workspace / "setup.fasta"
        write.write_fasta(fasta, sequence="A", header="A|protein|")
        argv = [
            "boltz", "predict", str(fasta), "--cache", str(cache_path),
            "--out_dir", str(workspace), "--use_msa_server",
            "--recycling_steps", "1", "--sampling_steps", "10",
        ]
        from cofolder.modules.runners.boltz_runner import run_boltz
        run_boltz(argv)


def setup_boltz2_cache(cache_path: str | Path) -> Path:
    """Explicitly initialize an incomplete Boltz2 cache and return its path."""
    root = Path(cache_path).expanduser()
    root.mkdir(parents=True, exist_ok=True)
    with (root / ".cofolder_setup.lock").open("a+b") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        missing = missing_boltz2_cache_components(root)
        if not missing:
            return root
        mols = root / "mols"
        if mols.is_dir() and not any(mols.iterdir()):
            mols.rmdir()
        _download_boltz2_cache(root)
        missing = missing_boltz2_cache_components(root)
        if missing:
            raise RuntimeError(
                f"Boltz2 cache setup did not produce required components at {root}: "
                + ", ".join(missing)
            )
    return root


def populate_ccd_cache_from_sdf(
    file_path: str | Path,
    property_id: str,
    *,
    cache_path: str | Path,
    on_conflict: str = "overwrite",
) -> None:
    """Validate an SDF and populate an already initialized CCD cache."""
    source = Path(file_path).expanduser()
    root = Path(cache_path).expanduser()
    if on_conflict not in {"use_cache", "overwrite"}:
        raise ValueError(f"Invalid on_conflict mode: {on_conflict}")
    if not source.is_file():
        raise FileNotFoundError(f"Input file not found: {source}")
    molecules = read.read_sdf(str(source))
    if not molecules:
        raise ValueError(f"SDF contains no valid molecule records: {source}")
    missing_ids = [
        index for index, molecule in enumerate(molecules)
        if not molecule.HasProp(property_id) or not molecule.GetProp(property_id).strip()
    ]
    if missing_ids:
        raise ValueError(
            f"SDF records {missing_ids} are missing required property '{property_id}'"
        )
    identifiers = [molecule.GetProp(property_id) for molecule in molecules]
    duplicates = [key for key, count in Counter(identifiers).items() if count > 1]
    if duplicates:
        raise ValueError(f"Duplicate molecule IDs found in SDF: {duplicates}")
    missing = missing_boltz2_cache_components(root)
    if missing:
        raise RuntimeError(
            f"Boltz2 cache at {root} is incomplete ({', '.join(missing)}). Run "
            f"`cofolder-tools setup-boltz2-cache --cache-path {root}` first."
        )
    from cofolder.modules.entities.ligand import mol_to_ccd
    mols = root / "mols"
    for molecule in molecules:
        identifier = molecule.GetProp(property_id)
        destination = mols / f"{identifier}.pkl"
        if destination.exists() and on_conflict == "use_cache":
            logger.info("Using cached CCD for %s", identifier)
            continue
        try:
            mol_to_ccd(identifier, molecule, boltz_path=root)
        except Exception as exc:
            logger.error("Failed to process ID %s: %s", identifier, exc)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="cofolder-tools")
    subparsers = parser.add_subparsers(dest="operation", required=True)
    setup = subparsers.add_parser("setup-boltz2-cache")
    setup.add_argument("--cache-path", required=True, type=Path)
    populate = subparsers.add_parser("populate-ccd-cache")
    populate.add_argument("--sdf", required=True, type=Path)
    populate.add_argument("--property-id", required=True)
    populate.add_argument("--cache-path", required=True, type=Path)
    populate.add_argument("--on-conflict", choices=("use_cache", "overwrite"), default="overwrite")
    args = parser.parse_args(argv)
    if args.operation == "setup-boltz2-cache":
        setup_boltz2_cache(args.cache_path)
    else:
        populate_ccd_cache_from_sdf(
            args.sdf, args.property_id, cache_path=args.cache_path,
            on_conflict=args.on_conflict,
        )
    return 0
