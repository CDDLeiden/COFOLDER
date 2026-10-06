"""Validated source bundles and query-specific bias reference materialization."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import date
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import tempfile
from typing import Any

import pandas as pd
from rdkit import Chem, DataStructs
from rdkit import rdBase
from rdkit.Chem import rdFingerprintGenerator

from cofolder.modules.input.system import System, iter_system_chains
from cofolder.modules.utils.executables import resolve_mmseqs_executable


BIAS_DATABASE_SCHEMA_VERSION = 2
BIAS_REFERENCE_SCHEMA_VERSION = 2
BIAS_QUERY_CACHE_SCHEMA_VERSION = 1
DEFAULT_BIAS_ROOT = Path("~/.cofolder/data/bias").expanduser()
DEFAULT_BIAS_QUERY_CACHE = Path("~/.cofolder/cache/bias_queries").expanduser()
DEFAULT_PROTEIN_DATABASE = DEFAULT_BIAS_ROOT / "protein"
DEFAULT_LIGAND_DATABASE = DEFAULT_BIAS_ROOT / "ligand"
PROTEIN_METADATA_NAME = "pdb_release_metadata.csv.gz"
PROTEIN_SEQUENCE_INDEX_NAME = "pdb_sequences.csv.gz"
LIGAND_TABLE_NAME = "pdb_ligands.csv.gz"
_MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
PROTEIN_SIMILARITY_CUTOFF = 0.25
PROTEIN_SIMILARITY_PERCENT_CUTOFF = PROTEIN_SIMILARITY_CUTOFF * 100.0
LIGAND_SIMILARITY_CUTOFF = 0.35


@dataclass(frozen=True, slots=True)
class BiasReleasePolicy:
    mode: str
    cutoff: date | None = None

    @property
    def label(self) -> str:
        return "whole" if self.mode == "whole" else self.cutoff.isoformat()

    def includes(self, value: object) -> bool:
        try:
            released = date.fromisoformat(str(value)[:10])
        except (TypeError, ValueError):
            return False
        return self.mode == "whole" or released < self.cutoff


def parse_bias_release_policy(value: str | BiasReleasePolicy) -> BiasReleasePolicy:
    if isinstance(value, BiasReleasePolicy):
        return value
    text = str(value).strip().lower()
    if text == "whole":
        return BiasReleasePolicy("whole")
    try:
        return BiasReleasePolicy("before", date.fromisoformat(text))
    except ValueError as exc:
        raise ValueError(
            "--bias_release_cutoff must be an ISO date (YYYY-MM-DD) or 'whole'."
        ) from exc


@dataclass(frozen=True, slots=True)
class BiasDatabaseBundle:
    kind: str
    root: Path
    manifest: dict[str, Any]
    fingerprint: str
    data_path: Path
    metadata_path: Path | None = None
    sequence_index_path: Path | None = None


@dataclass(frozen=True, slots=True)
class BiasReferenceArtifacts:
    protein_path: Path | None
    ligand_path: Path | None
    manifest_path: Path
    reused: bool


@dataclass(frozen=True, slots=True)
class CustomBiasReferenceBundle:
    root: Path
    manifest: dict[str, Any]
    fingerprint: str
    protein_path: Path
    ligand_path: Path


def validate_custom_bias_reference_bundle(path: str | Path) -> CustomBiasReferenceBundle:
    root = Path(path).expanduser().resolve()
    manifest_path = root / "manifest.json"
    if not root.is_dir() or not manifest_path.is_file():
        raise ValueError(f"Custom bias reference bundle is incomplete: {root}")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"Invalid custom bias reference manifest: {manifest_path}") from exc
    if manifest.get("schema_version") != 1 or manifest.get("kind") != "custom_complexes":
        raise ValueError(f"Unsupported custom bias reference bundle: {manifest_path}")
    files = manifest.get("files", {})
    for relative, checksum in files.items():
        candidate = (root / relative).resolve()
        try:
            candidate.relative_to(root)
        except ValueError as exc:
            raise ValueError(f"Unsafe custom reference path: {relative}") from exc
        if not candidate.is_file() or _sha256(candidate) != checksum:
            raise ValueError(f"Custom reference checksum mismatch: {candidate}")
    protein = root / "custom_protein_references.csv"
    ligand = root / "custom_ligand_references.csv"
    if not protein.is_file() or not ligand.is_file():
        raise ValueError(f"Custom bundle lacks paired reference tables: {root}")
    return CustomBiasReferenceBundle(
        root, manifest, _manifest_fingerprint(manifest), protein, ligand
    )


def normalize_pdb_id(value: object) -> str:
    """Normalize common PDB/MMseqs identifiers to their four-character PDB ID."""
    text = str(value).strip().upper().lstrip(">")
    for prefix in ("PDB|", "PDB:"):
        if text.startswith(prefix):
            text = text[len(prefix):]
    candidate = text.split()[0].split("|")[0].split("_")[0].split(":")[0]
    return candidate[:4] if len(candidate) >= 4 else candidate


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _manifest_fingerprint(manifest: dict[str, Any]) -> str:
    raw = json.dumps(manifest, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def validate_bias_database_bundle(path: str | Path, expected_kind: str) -> BiasDatabaseBundle:
    root = Path(path).expanduser()
    path_option = f"--bias_training_data_{expected_kind}_path"
    skip_option = "--skip_ligand" if expected_kind == "protein" else "--skip_mmseqs"
    preparation = "Prepare it with `cofolder-tools fetch-bias-training-data " + " ".join(
        (
            f"--output_root {shlex.quote(str(root.parent))}",
            f"{path_option} {shlex.quote(str(root))}",
            skip_option,
        )
    ) + "`."

    def invalid(message: str) -> ValueError:
        return ValueError(f"{message} {preparation}")

    if not root.is_dir():
        raise invalid(f"Bias {expected_kind} database directory not found: {root}.")
    manifest_path = root / "manifest.json"
    if not manifest_path.is_file():
        raise invalid(f"Bias {expected_kind} database manifest not found: {manifest_path}.")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise invalid(f"Invalid bias database manifest {manifest_path}: {exc}.") from exc
    if manifest.get("schema_version") != BIAS_DATABASE_SCHEMA_VERSION:
        raise invalid(
            f"Unsupported bias database schema in {manifest_path}: "
            f"{manifest.get('schema_version')!r}; expected {BIAS_DATABASE_SCHEMA_VERSION}."
        )
    if manifest.get("kind") != expected_kind:
        raise invalid(
            f"Bias database {root} has kind {manifest.get('kind')!r}; expected {expected_kind!r}."
        )
    files = manifest.get("files")
    if not isinstance(files, dict) or not files:
        raise invalid(f"Bias database manifest {manifest_path} has no file inventory.")
    resolved: dict[str, Path] = {}
    for relative, expected_hash in files.items():
        candidate = (root / str(relative)).resolve()
        try:
            candidate.relative_to(root.resolve())
        except ValueError as exc:
            raise invalid(f"Bias database manifest contains an unsafe path: {relative}.") from exc
        if not candidate.is_file():
            raise invalid(f"Bias database file is missing: {candidate}.")
        if _sha256(candidate) != expected_hash:
            raise invalid(f"Bias database checksum mismatch: {candidate}.")
        resolved[str(relative)] = candidate
    if expected_kind == "protein":
        data_path = root / "mmseqs" / "db"
        metadata_path = root / PROTEIN_METADATA_NAME
        sequence_index_path = root / PROTEIN_SEQUENCE_INDEX_NAME
        if (
            not metadata_path.is_file()
            or not sequence_index_path.is_file()
            or not any(data_path.parent.glob("db*"))
        ):
            raise invalid(
                f"Protein bias bundle {root} lacks its MMseqs database, release metadata, "
                "or offline PDB sequence index."
            )
        columns = set(pd.read_csv(metadata_path, nrows=0).columns)
        required = {"pdb_id", "release_date"}
        if missing := required - columns:
            raise invalid(
                f"Protein bias bundle metadata is missing columns: {sorted(missing)}."
            )
        sequence_columns = set(pd.read_csv(sequence_index_path, nrows=0).columns)
        if missing := {"pdb_id", "sequence"} - sequence_columns:
            raise invalid(
                f"Protein bias sequence index is missing columns: {sorted(missing)}."
            )
        return BiasDatabaseBundle(
            expected_kind, root.resolve(), manifest, _manifest_fingerprint(manifest),
            data_path.resolve(), metadata_path.resolve(), sequence_index_path.resolve(),
        )
    data_path = root / LIGAND_TABLE_NAME
    if not data_path.is_file():
        raise invalid(f"Ligand bias bundle {root} lacks {LIGAND_TABLE_NAME}.")
    columns = set(pd.read_csv(data_path, nrows=0).columns)
    required = {"pdb_id", "release_date", "ligand_id", "smiles"}
    if missing := required - columns:
        raise invalid(f"Ligand bias bundle table is missing columns: {sorted(missing)}.")
    return BiasDatabaseBundle(
        expected_kind, root.resolve(), manifest, _manifest_fingerprint(manifest), data_path.resolve()
    )


def _selected_queries(
    system_obj: System,
    selected_chains: set[str] | None,
    ligand_smiles_by_id: dict[str, str] | None = None,
) -> tuple[dict[str, str], dict[str, str]]:
    proteins: dict[str, str] = {}
    ligands: dict[str, str] = {}
    sequences = system_obj.find_value(key="sequences") or []
    for chain in iter_system_chains(system_obj):
        chain_id = chain.chain_id.strip().upper()
        if selected_chains and chain_id not in selected_chains:
            continue
        payload = sequences[chain.sequence_index].get(chain.entity_type, {})
        if chain.entity_type == "protein":
            sequence = str(payload.get("sequence") or payload.get("fasta") or "").strip()
            if sequence:
                proteins[chain_id] = sequence
        elif chain.entity_type == "ligand":
            smiles = str(payload.get("smiles") or "").strip()
            if not smiles:
                raw_ccd = payload.get("ccd") or payload.get("ccd_codes") or []
                ccd_ids = raw_ccd if isinstance(raw_ccd, list) else [raw_ccd]
                smiles = next(
                    (
                        str((ligand_smiles_by_id or {}).get(str(ccd).upper()) or "").strip()
                        for ccd in ccd_ids
                        if (ligand_smiles_by_id or {}).get(str(ccd).upper())
                    ),
                    "",
                )
            if smiles:
                ligands[chain_id] = smiles
    return proteins, ligands


def _normalized_sequence(sequence: str) -> str:
    return "".join(str(sequence).split()).upper()


def _canonical_smiles(smiles: str) -> str:
    molecule = Chem.MolFromSmiles(str(smiles))
    if molecule is None:
        raise ValueError(f"Could not parse ligand SMILES for bias query: {smiles}")
    return Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)


def _mmseqs_version(binary: str) -> str:
    try:
        completed = subprocess.run(
            [binary, "version"], check=True, capture_output=True, text=True, timeout=30
        )
        return (completed.stdout or completed.stderr).strip().splitlines()[0]
    except (OSError, subprocess.SubprocessError, IndexError):
        return "unknown"


def _cache_key(kind: str, material: dict[str, Any]) -> str:
    encoded = json.dumps(
        {"schema_version": BIAS_QUERY_CACHE_SCHEMA_VERSION, "kind": kind, **material},
        sort_keys=True,
        separators=(",", ":"),
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


def _cache_marker(root: Path) -> Path:
    return root / "cache.json"


def bias_query_cache_readiness(path: str | Path | None) -> tuple[Path, bool, str]:
    root = Path(path or DEFAULT_BIAS_QUERY_CACHE).expanduser().resolve()
    marker = _cache_marker(root)
    if root.exists():
        if not root.is_dir():
            return root, False, "cache path is not a directory"
        if marker.exists():
            try:
                payload = json.loads(marker.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                return root, False, "cache marker is invalid"
            if payload.get("schema_version") != BIAS_QUERY_CACHE_SCHEMA_VERSION:
                return root, False, "cache schema is unsupported"
        return root, os.access(root, os.R_OK | os.W_OK), "existing"
    ancestor = root.parent
    while not ancestor.exists() and ancestor != ancestor.parent:
        ancestor = ancestor.parent
    ready = ancestor.is_dir() and os.access(ancestor, os.W_OK)
    return root, ready, "will_create" if ready else "parent is not writable"


def _ensure_cache_root(root: Path) -> None:
    root.mkdir(parents=True, exist_ok=True)
    marker = _cache_marker(root)
    if marker.exists():
        payload = json.loads(marker.read_text(encoding="utf-8"))
        if payload.get("schema_version") != BIAS_QUERY_CACHE_SCHEMA_VERSION:
            raise ValueError(f"Unsupported bias query cache schema at {marker}.")
    else:
        with tempfile.NamedTemporaryFile("w", dir=root, delete=False, encoding="utf-8") as handle:
            json.dump({"schema_version": BIAS_QUERY_CACHE_SCHEMA_VERSION}, handle)
            temporary = Path(handle.name)
        temporary.replace(marker)
    (root / ".locks").mkdir(exist_ok=True)


def _cached_frame(
    *,
    root: Path,
    kind: str,
    key: str,
    identity: dict[str, Any],
    builder,
) -> tuple[pd.DataFrame, str]:
    _ensure_cache_root(root)
    global_lock_path = root / ".invalidate.lock"
    with global_lock_path.open("a+", encoding="utf-8") as global_lock:
        fcntl.flock(global_lock.fileno(), fcntl.LOCK_SH)
        return _cached_frame_entry(
            root=root, kind=kind, key=key, identity=identity, builder=builder
        )


def _cached_frame_entry(
    *,
    root: Path,
    kind: str,
    key: str,
    identity: dict[str, Any],
    builder,
) -> tuple[pd.DataFrame, str]:
    destination = root / kind / key[:2] / key
    lock_path = root / ".locks" / f"{kind}-{key}.lock"
    lock_path.parent.mkdir(parents=True, exist_ok=True)

    def read_entry() -> pd.DataFrame | None:
        manifest_path = destination / "manifest.json"
        payload_path = destination / "payload.csv.gz"
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            if (
                manifest.get("schema_version") != BIAS_QUERY_CACHE_SCHEMA_VERSION
                or manifest.get("key") != key
                or manifest.get("identity") != identity
                or manifest.get("payload_sha256") != _sha256(payload_path)
            ):
                return None
            return pd.read_csv(payload_path)
        except (OSError, ValueError, json.JSONDecodeError):
            return None

    existing = read_entry() if destination.exists() else None
    if existing is not None:
        return existing, "hit"
    with lock_path.open("a+", encoding="utf-8") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        existing = read_entry() if destination.exists() else None
        if existing is not None:
            return existing, "hit"
        status = "rebuilt" if destination.exists() else "miss"
        frame = builder()
        destination.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=f".{key}-", dir=destination.parent) as raw:
            staged = Path(raw) / key
            staged.mkdir()
            payload_path = staged / "payload.csv.gz"
            frame.to_csv(payload_path, index=False)
            (staged / "manifest.json").write_text(
                json.dumps(
                    {
                        "schema_version": BIAS_QUERY_CACHE_SCHEMA_VERSION,
                        "kind": kind,
                        "key": key,
                        "identity": identity,
                        "payload_sha256": _sha256(payload_path),
                    },
                    indent=2,
                ),
                encoding="utf-8",
            )
            if destination.exists():
                shutil.rmtree(destination)
            staged.replace(destination)
        return frame, status


def _protein_references(
    bundle: BiasDatabaseBundle,
    queries: dict[str, str],
    policy: BiasReleasePolicy,
    output_dir: Path,
    cache_root: Path,
    threshold: float,
) -> tuple[pd.DataFrame, list[dict[str, str]]]:
    from cofolder.tools.build_bias_training_data import _run_mmseqs
    from cofolder.modules.analytics.bias_training import _resolve_mmseqs_bin

    metadata = pd.read_csv(bundle.metadata_path)
    required = {"pdb_id", "release_date"}
    if required - set(metadata.columns):
        raise ValueError(f"Protein database metadata missing columns: {sorted(required - set(metadata.columns))}")
    metadata["pdb_id"] = metadata["pdb_id"].map(normalize_pdb_id)
    metadata = metadata[metadata["release_date"].map(policy.includes)]
    mmseqs = _resolve_mmseqs_bin()
    if mmseqs is None:
        diagnostics = resolve_mmseqs_executable().diagnostic_text()
        raise RuntimeError(
            "MMseqs2 is required to search the protein bias database. "
            "Run `cofolder-tools install-mmseqs` or set COFOLDER_MMSEQS_BIN. "
            f"Attempted: {diagnostics}."
        )
    rows: list[pd.DataFrame] = []
    local: dict[str, tuple[pd.DataFrame, str, str]] = {}
    provenance: list[dict[str, str]] = []
    mmseqs_version = _mmseqs_version(mmseqs)
    for chain_id, sequence in queries.items():
        normalized = _normalized_sequence(sequence)
        sequence_hash = hashlib.sha256(normalized.encode()).hexdigest()
        identity = {
            "bundle": bundle.fingerprint,
            "sequence_sha256": sequence_hash,
            "mmseqs_version": mmseqs_version,
            "parameters": {"workers": 1, "max_seqs": 10000},
        }
        key = _cache_key("protein", identity)
        if sequence_hash not in local:
            def build(sequence=normalized, sequence_hash=sequence_hash):
                with tempfile.TemporaryDirectory(prefix="cofolder-protein-query-") as raw:
                    query_root = Path(raw)
                    fasta = query_root / "query.fasta"
                    fasta.write_text(f">query\n{sequence}\n", encoding="utf-8")
                    hits = _run_mmseqs(
                        mmseqs, fasta, bundle.data_path, 1, query_root / "tmp", 10000
                    )
                if "target" in hits:
                    hits = hits.rename(columns={"target": "pdb_id", "pident": "sequence_similarity", "tseq": "sequence"})
                hits["pdb_id"] = hits["pdb_id"].map(normalize_pdb_id)
                return hits
            raw_hits, status = _cached_frame(
                root=cache_root, kind="protein", key=key, identity=identity, builder=build
            )
            local[sequence_hash] = (raw_hits, status, key)
        raw_hits, status, key = local[sequence_hash]
        frame = raw_hits.merge(metadata, on="pdb_id", how="inner")
        frame = frame[
            pd.to_numeric(frame["sequence_similarity"], errors="coerce")
            > float(threshold) * 100.0
        ].copy()
        frame["query_chain_id"] = chain_id
        rows.append(frame)
        provenance.append({"chain_id": chain_id, "key": key, "status": status})
    columns = ["query_chain_id", "pdb_id", "release_date", "sequence_similarity", "sequence"]
    if not rows:
        return pd.DataFrame(columns=columns), provenance
    result = pd.concat(rows, ignore_index=True, sort=False)
    return result.loc[:, [name for name in columns if name in result.columns]], provenance


def _ligand_references(
    bundle: BiasDatabaseBundle,
    queries: dict[str, str],
    policy: BiasReleasePolicy,
    threshold: float,
    cache_root: Path,
) -> tuple[pd.DataFrame, list[dict[str, str]]]:
    database = pd.read_csv(bundle.data_path)
    required = {"pdb_id", "release_date", "ligand_id", "smiles"}
    if required - set(database.columns):
        raise ValueError(f"Ligand database missing columns: {sorted(required - set(database.columns))}")
    reference_fps = [
        _MORGAN.GetFingerprint(mol) if (mol := Chem.MolFromSmiles(str(smiles))) is not None else None
        for smiles in database["smiles"]
    ]
    rows: list[pd.DataFrame] = []
    local: dict[str, tuple[pd.DataFrame, str, str]] = {}
    provenance: list[dict[str, str]] = []
    for chain_id, smiles in queries.items():
        canonical = _canonical_smiles(smiles)
        identity = {
            "bundle": bundle.fingerprint,
            "canonical_smiles_sha256": hashlib.sha256(canonical.encode()).hexdigest(),
            "fingerprint": {"algorithm": "Morgan", "radius": 2, "bits": 2048},
            "rdkit_version": rdBase.rdkitVersion,
        }
        key = _cache_key("ligand", identity)
        if canonical not in local:
            def build(canonical=canonical):
                query_fp = _MORGAN.GetFingerprint(Chem.MolFromSmiles(canonical))
                frame = database.copy()
                frame["ecfp_similarity"] = [
                    float(DataStructs.TanimotoSimilarity(query_fp, fp)) if fp is not None else float("nan")
                    for fp in reference_fps
                ]
                return frame
            raw, status = _cached_frame(
                root=cache_root, kind="ligand", key=key, identity=identity, builder=build
            )
            local[canonical] = (raw, status, key)
        raw, status, key = local[canonical]
        frame = raw[raw["release_date"].map(policy.includes)].copy()
        frame = frame[frame["ecfp_similarity"] > float(threshold)].copy()
        frame["query_chain_id"] = chain_id
        rows.append(frame)
        provenance.append({"chain_id": chain_id, "key": key, "status": status})
    columns = ["query_chain_id", "pdb_id", "release_date", "ligand_id", "ecfp_similarity", "smiles"]
    if not rows:
        return pd.DataFrame(columns=columns), provenance
    return pd.concat(rows, ignore_index=True, sort=False).loc[:, columns], provenance


def _run_mmseqs_backfill(
    mmseqs: str, query_sequence: str, target_sequences: list[str]
) -> pd.DataFrame:
    with tempfile.TemporaryDirectory(prefix="cofolder-bias-backfill-") as raw:
        root = Path(raw)
        query_fasta = root / "query.fasta"
        target_fasta = root / "target.fasta"
        query_fasta.write_text(f">query\n{query_sequence}\n", encoding="utf-8")
        target_fasta.write_text(
            "".join(f">target_{index}\n{sequence}\n" for index, sequence in enumerate(target_sequences)),
            encoding="utf-8",
        )
        qdb, tdb, rdb = root / "qdb", root / "tdb", root / "rdb"
        result = root / "result.tsv"
        temporary = root / "tmp"
        commands = [
            [mmseqs, "createdb", str(query_fasta), str(qdb)],
            [mmseqs, "createdb", str(target_fasta), str(tdb)],
            [mmseqs, "search", str(qdb), str(tdb), str(rdb), str(temporary), "--threads", "1", "--max-seqs", str(max(1, len(target_sequences)))],
            [mmseqs, "convertalis", str(qdb), str(tdb), str(rdb), str(result), "--format-output", "target,pident"],
        ]
        for command in commands:
            subprocess.run(command, check=True, capture_output=True, text=True, timeout=120)
        hits = (
            pd.read_csv(result, sep="\t", header=None, names=["target", "sequence_similarity"])
            if result.is_file() and result.stat().st_size
            else pd.DataFrame(columns=["target", "sequence_similarity"])
        )
    best: dict[int, float] = {}
    for _, row in hits.iterrows():
        match = re.search(r"(?:^|[|:_])target_(\d+)(?:$|[|:_])", str(row["target"]), re.I)
        if match is None:
            match = re.search(r"target_(\d+)", str(row["target"]), re.I)
        if match is None:
            continue
        index = int(match.group(1))
        similarity = float(row["sequence_similarity"])
        best[index] = max(similarity, best.get(index, 0.0))
    return pd.DataFrame(
        [
            {
                "sequence": sequence,
                "sequence_similarity": best.get(index, 0.0),
                "sequence_similarity_pairwise": pd.NA,
                "sequence_similarity_method": "mmseqs_pident",
            }
            for index, sequence in enumerate(target_sequences)
        ]
    )


def _backfill_protein_references(
    *,
    bundle: BiasDatabaseBundle,
    queries: dict[str, str],
    protein_df: pd.DataFrame,
    ligand_df: pd.DataFrame,
    policy: BiasReleasePolicy,
    cache_root: Path,
    protein_similarity_threshold: float,
) -> tuple[pd.DataFrame, list[dict[str, str]]]:
    if ligand_df.empty or bundle.sequence_index_path is None:
        return protein_df, []
    from cofolder.modules.analytics.bias_training import _resolve_mmseqs_bin

    mmseqs = _resolve_mmseqs_bin()
    if mmseqs is None:
        return protein_df, []
    version = _mmseqs_version(mmseqs)
    sequences = pd.read_csv(bundle.sequence_index_path)
    sequences["pdb_id"] = sequences["pdb_id"].map(normalize_pdb_id)
    metadata = pd.read_csv(bundle.metadata_path)
    metadata["pdb_id"] = metadata["pdb_id"].map(normalize_pdb_id)
    metadata = metadata[metadata["release_date"].map(policy.includes)]
    target_pdbs = sorted(set(ligand_df["pdb_id"].map(normalize_pdb_id)))
    rows = [protein_df]
    events: list[dict[str, str]] = []
    for chain_id, raw_query in queries.items():
        query = _normalized_sequence(raw_query)
        present = set(
            protein_df.loc[protein_df["query_chain_id"] == chain_id, "pdb_id"].map(normalize_pdb_id)
        ) if not protein_df.empty else set()
        for pdb_id in target_pdbs:
            if pdb_id in present:
                continue
            targets = sequences.loc[sequences["pdb_id"] == pdb_id, "sequence"].dropna().astype(str).tolist()
            if not targets:
                continue
            identity = {
                "mode": "pdb_backfill",
                "bundle": bundle.fingerprint,
                "pdb_id": pdb_id,
                "sequence_sha256": hashlib.sha256(query.encode()).hexdigest(),
                "mmseqs_version": version,
                "parameters": {"threads": 1, "max_seqs": len(targets)},
            }
            key = _cache_key("protein", identity)
            frame, status = _cached_frame(
                root=cache_root,
                kind="protein",
                key=key,
                identity=identity,
                builder=lambda query=query, targets=targets: _run_mmseqs_backfill(mmseqs, query, targets),
            )
            if frame.empty:
                events.append(
                    {
                        "chain_id": chain_id,
                        "pdb_id": pdb_id,
                        "key": key,
                        "status": status,
                    }
                )
                continue
            release = metadata.loc[metadata["pdb_id"] == pdb_id, "release_date"]
            frame = frame.copy()
            frame["query_chain_id"] = chain_id
            frame["pdb_id"] = pdb_id
            frame["release_date"] = release.iloc[0] if len(release) else pd.NA
            threshold_percent = float(protein_similarity_threshold) * 100.0
            if (pd.to_numeric(frame["sequence_similarity"], errors="coerce") >= threshold_percent).any():
                import logging
                logging.getLogger(__name__).warning(
                    "Protein MMseqs fallback reached or exceeded %.1f%% for query chain %s and PDB %s; check reference threshold/top-N consistency.",
                    threshold_percent, chain_id, pdb_id,
                )
            rows.append(frame)
            events.append({"chain_id": chain_id, "pdb_id": pdb_id, "key": key, "status": status})
    combined = pd.concat(rows, ignore_index=True, sort=False)
    columns = ["query_chain_id", "pdb_id", "release_date", "sequence_similarity", "sequence", "sequence_similarity_pairwise", "sequence_similarity_method"]
    return combined.loc[:, [column for column in columns if column in combined]], events


def materialize_bias_references(
    *,
    system_obj: System,
    protein_bundle: BiasDatabaseBundle | None,
    ligand_bundle: BiasDatabaseBundle | None,
    output_dir: Path,
    release_cutoff: str,
    protein_similarity_threshold: float = PROTEIN_SIMILARITY_CUTOFF,
    ligand_similarity_threshold: float = LIGAND_SIMILARITY_CUTOFF,
    selected_chains: set[str] | None,
    query_cache_path: str | Path | None = None,
    custom_fingerprint: str | None = None,
) -> BiasReferenceArtifacts:
    policy = parse_bias_release_policy(release_cutoff)
    cache_root = Path(query_cache_path or DEFAULT_BIAS_QUERY_CACHE).expanduser().resolve()
    ligand_smiles_by_id: dict[str, str] = {}
    if ligand_bundle is not None:
        ligand_source = pd.read_csv(ligand_bundle.data_path)
        if {"ligand_id", "smiles"} <= set(ligand_source.columns):
            ligand_smiles_by_id = {
                str(row["ligand_id"]).upper(): str(row["smiles"])
                for _, row in ligand_source.dropna(subset=["ligand_id", "smiles"]).iterrows()
            }
    proteins, ligands = _selected_queries(
        system_obj, selected_chains, ligand_smiles_by_id
    )
    expected_ligands = {
        chain.chain_id.strip().upper()
        for chain in iter_system_chains(system_obj)
        if chain.entity_type == "ligand"
        and (not selected_chains or chain.chain_id.strip().upper() in selected_chains)
    }
    unresolved_ligands = sorted(expected_ligands - set(ligands))
    if unresolved_ligands:
        raise ValueError(
            "The ligand bias database cannot resolve SMILES for query chain(s): "
            + ", ".join(unresolved_ligands)
        )
    query_hashes = {
        "protein": {key: hashlib.sha256(_normalized_sequence(value).encode()).hexdigest() for key, value in proteins.items()},
        "ligand": {key: hashlib.sha256(_canonical_smiles(value).encode()).hexdigest() for key, value in ligands.items()},
    }
    request = {
        "schema_version": BIAS_REFERENCE_SCHEMA_VERSION,
        "sources": {
            "protein": _bundle_provenance(protein_bundle),
            "ligand": _bundle_provenance(ligand_bundle),
            "custom_fingerprint": custom_fingerprint,
        },
        "queries": query_hashes,
        "release_policy": _policy_provenance(policy, protein_bundle, ligand_bundle),
        "protein_similarity_threshold": float(protein_similarity_threshold),
        "ligand_similarity_threshold": float(ligand_similarity_threshold),
        "selected_chains": sorted(selected_chains or ()),
        "query_cache_path": str(cache_root),
    }
    manifest_path = output_dir / "reference_manifest.json"
    protein_path = output_dir / "protein_references.csv" if proteins else None
    ligand_path = output_dir / "ligand_references.csv" if ligands else None
    if manifest_path.is_file():
        try:
            previous = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            previous = None
        previous_request = previous.get("request") if isinstance(previous, dict) else None
        outputs = previous.get("outputs", {}) if isinstance(previous, dict) else {}
        output_valid = all(
            path is None or (
                path.is_file() and outputs.get(path.name) == _sha256(path)
            )
            for path in (protein_path, ligand_path)
        )
        if previous_request == request and output_valid:
            return BiasReferenceArtifacts(protein_path, ligand_path, manifest_path, True)
    output_dir.mkdir(parents=True, exist_ok=True)
    if proteins and protein_bundle is None:
        raise ValueError("Protein queries require a protein bias database bundle.")
    if ligands and ligand_bundle is None:
        raise ValueError("Ligand queries require a ligand bias database bundle.")
    if proteins:
        protein_df, protein_cache = _protein_references(
            protein_bundle,
            proteins,
            policy,
            output_dir,
            cache_root,
            protein_similarity_threshold,
        )
    else:
        protein_df, protein_cache = None, []
    if ligands:
        ligand_df, ligand_cache = _ligand_references(
            ligand_bundle, ligands, policy, ligand_similarity_threshold, cache_root
        )
    else:
        ligand_df, ligand_cache = None, []
    if protein_df is not None and ligand_df is not None:
        protein_df, backfill_cache = _backfill_protein_references(
            bundle=protein_bundle,
            queries=proteins,
            protein_df=protein_df,
            ligand_df=ligand_df,
            policy=policy,
            cache_root=cache_root,
            protein_similarity_threshold=protein_similarity_threshold,
        )
    else:
        backfill_cache = []
    with tempfile.TemporaryDirectory(prefix=".bias-references-", dir=output_dir) as temporary:
        temp = Path(temporary)
        if protein_df is not None:
            protein_df.to_csv(temp / protein_path.name, index=False)
        if ligand_df is not None:
            ligand_df.to_csv(temp / ligand_path.name, index=False)
            for chain_id, frame in ligand_df.groupby("query_chain_id"):
                frame.to_csv(temp / f"ligand_references_{chain_id}.csv", index=False)
        output_hashes = {
            candidate.name: _sha256(candidate)
            for candidate in temp.iterdir()
            if candidate.is_file()
        }
        manifest = {
            "schema_version": BIAS_REFERENCE_SCHEMA_VERSION,
            "request": request,
            "cache": {
                "root": str(cache_root),
                "protein": protein_cache,
                "ligand": ligand_cache,
            },
            "backfill": backfill_cache,
            "outputs": output_hashes,
        }
        (temp / manifest_path.name).write_text(json.dumps(manifest, indent=2), encoding="utf-8")
        for candidate in temp.iterdir():
            candidate.replace(output_dir / candidate.name)
    return BiasReferenceArtifacts(protein_path, ligand_path, manifest_path, False)


def _bundle_provenance(bundle: BiasDatabaseBundle | None) -> dict[str, Any] | None:
    if bundle is None:
        return None
    return {
        "path": str(bundle.root),
        "fingerprint": bundle.fingerprint,
        "source": bundle.manifest.get("source"),
        "source_snapshot": bundle.manifest.get("source_snapshot"),
    }


def _maximum_release(bundle: BiasDatabaseBundle | None) -> str | None:
    if bundle is None:
        return None
    source = bundle.metadata_path or bundle.data_path
    frame = pd.read_csv(source, usecols=["release_date"])
    values = []
    for value in frame["release_date"]:
        candidate = str(value)[:10]
        try:
            date.fromisoformat(candidate)
        except ValueError:
            continue
        values.append(candidate)
    return max(values) if values else None


def _policy_provenance(
    policy: BiasReleasePolicy,
    protein_bundle: BiasDatabaseBundle | None,
    ligand_bundle: BiasDatabaseBundle | None,
) -> dict[str, Any]:
    return {
        "mode": policy.mode,
        "cutoff": policy.cutoff.isoformat() if policy.cutoff else None,
        "strict_before": policy.mode == "before",
        "source_max_release_dates": {
            "protein": _maximum_release(protein_bundle),
            "ligand": _maximum_release(ligand_bundle),
        },
    }


__all__ = [
    "BIAS_DATABASE_SCHEMA_VERSION",
    "BIAS_QUERY_CACHE_SCHEMA_VERSION",
    "BiasReleasePolicy",
    "BiasDatabaseBundle",
    "CustomBiasReferenceBundle",
    "BiasReferenceArtifacts",
    "DEFAULT_LIGAND_DATABASE",
    "DEFAULT_PROTEIN_DATABASE",
    "DEFAULT_BIAS_QUERY_CACHE",
    "LIGAND_TABLE_NAME",
    "PROTEIN_METADATA_NAME",
    "PROTEIN_SEQUENCE_INDEX_NAME",
    "PROTEIN_SIMILARITY_CUTOFF",
    "PROTEIN_SIMILARITY_PERCENT_CUTOFF",
    "LIGAND_SIMILARITY_CUTOFF",
    "bias_query_cache_readiness",
    "materialize_bias_references",
    "normalize_pdb_id",
    "parse_bias_release_policy",
    "validate_bias_database_bundle",
    "validate_custom_bias_reference_bundle",
]
