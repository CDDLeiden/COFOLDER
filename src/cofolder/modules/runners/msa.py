"""Helpers for persisting and reusing runner-generated protein MSAs."""

from __future__ import annotations

import csv
import hashlib
import json
import re
import shutil
from pathlib import Path
from typing import Any

from cofolder.modules.input.system import iter_system_chains


MANIFEST_NAME = "manifest.json"


def protein_sequence(payload: dict[str, Any]) -> str:
    """Return a normalized protein sequence from a system entity payload."""

    raw = payload.get("sequence") or payload.get("fasta") or ""
    return re.sub(r"[^A-Za-z]", "", str(raw)).upper()


def sequence_key(sequence: str) -> str:
    """Return the stable cache key for a normalized protein sequence."""

    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()


def protein_payloads(system_obj: Any) -> list[dict[str, Any]]:
    """Return unique protein entity payloads in system order."""

    payloads: list[dict[str, Any]] = []
    seen_indices: set[int] = set()
    for chain in iter_system_chains(system_obj):
        if chain.entity_type != "protein" or chain.sequence_index in seen_indices:
            continue
        seen_indices.add(chain.sequence_index)
        payloads.append(chain.entity_data)
    return payloads


def resolve_declared_msa_paths(system_obj: Any, *, base_dir: Path) -> int:
    """Resolve relative declared MSA paths before writing relocated system YAMLs."""

    resolved = 0
    for payload in protein_payloads(system_obj):
        raw_path = payload.get("msa")
        if raw_path is None or not str(raw_path).strip():
            continue
        if str(raw_path).strip().lower() == "empty":
            continue
        path = Path(str(raw_path)).expanduser()
        if not path.is_absolute():
            path = (base_dir / path).resolve()
        payload["msa"] = str(path)
        resolved += 1
    return resolved


def unresolved_protein_sequences(system_obj: Any) -> list[str]:
    """Return normalized sequences for protein entities without a declared MSA."""

    unresolved: list[str] = []
    for payload in protein_payloads(system_obj):
        if payload.get("msa") is not None and str(payload.get("msa")).strip():
            continue
        sequence = protein_sequence(payload)
        if sequence and sequence not in unresolved:
            unresolved.append(sequence)
    return unresolved


def inject_cached_msas(system_obj: Any, cache_dir: Path) -> int:
    """Inject cached MSA paths into matching unresolved protein entities."""

    manifest = _read_manifest(cache_dir)
    entries = manifest.get("proteins", {})
    injected = 0
    for payload in protein_payloads(system_obj):
        if payload.get("msa") is not None and str(payload.get("msa")).strip():
            continue
        sequence = protein_sequence(payload)
        key = sequence_key(sequence) if sequence else None
        if key is None or key not in entries:
            continue
        entry = entries[key]
        if not isinstance(entry, dict):
            raise ValueError(
                f"Reusable MSA manifest entry for sequence {key} has an invalid schema."
            )
        relative_path = entry.get("path")
        if not relative_path:
            raise ValueError(
                f"Reusable MSA manifest entry for sequence {key} "
                "does not define an artifact path. Refusing to recalculate silently."
            )
        path = (cache_dir / str(relative_path)).resolve()
        if not path.is_file():
            raise ValueError(
                f"Reusable MSA artifact is missing: {path}. "
                "Refusing to recalculate silently."
            )
        if _read_msa_query(path) != sequence:
            raise ValueError(
                f"Reusable MSA artifact does not match its protein sequence: {path}. "
                "Refusing to recalculate silently."
            )
        payload["msa"] = str(path)
        injected += 1
    return injected


def capture_generated_msas(
    system_obj: Any,
    *,
    generated_dir: Path,
    cache_dir: Path,
    runner_name: str,
) -> int:
    """Copy generated MSAs into a stable screen-level sequence-keyed cache."""

    unresolved = set(unresolved_protein_sequences(system_obj))
    if not unresolved or not generated_dir.is_dir():
        return 0

    candidates = sorted(
        (
            path
            for path in generated_dir.rglob("*")
            if path.is_file() and path.suffix.lower() in {".csv", ".a3m"}
        ),
        key=lambda path: (
            len(path.relative_to(generated_dir).parts),
            path.suffix != ".csv",
            str(path),
        ),
    )
    matched: dict[str, Path] = {}
    for candidate in candidates:
        query = _read_msa_query(candidate)
        if query in unresolved and query not in matched:
            matched[query] = candidate

    if not matched:
        return 0

    cache_dir.mkdir(parents=True, exist_ok=True)
    manifest = _read_manifest(cache_dir)
    manifest["version"] = 1
    manifest["runner"] = str(runner_name)
    proteins = manifest.setdefault("proteins", {})
    for sequence, source in matched.items():
        key = sequence_key(sequence)
        destination = cache_dir / f"protein_{key[:16]}{source.suffix.lower()}"
        if source.resolve() != destination.resolve():
            shutil.copy2(source, destination)
        proteins[key] = {
            "path": destination.name,
            "sequence_sha256": key,
            "source": str(source),
        }
    _write_manifest(cache_dir, manifest)
    return len(matched)


def _read_msa_query(path: Path) -> str | None:
    try:
        if path.suffix.lower() == ".csv":
            with path.open(encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle)
                first = next(reader, None)
            if not first or "sequence" not in first:
                return None
            return re.sub(r"[^A-Za-z]", "", str(first["sequence"])).upper()

        with path.open(encoding="utf-8") as handle:
            sequence_parts: list[str] = []
            found_header = False
            for raw_line in handle:
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if found_header and sequence_parts:
                        break
                    found_header = True
                    continue
                if found_header:
                    sequence_parts.append(line)
            if not sequence_parts:
                return None
            return re.sub(r"[^A-Za-z]", "", "".join(sequence_parts)).upper()
    except (OSError, UnicodeError, csv.Error):
        return None


def _read_manifest(cache_dir: Path) -> dict[str, Any]:
    path = cache_dir / MANIFEST_NAME
    if not path.is_file():
        return {"version": 1, "proteins": {}}
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"Reusable MSA manifest is unreadable: {path}") from exc
    if not isinstance(value, dict) or not isinstance(value.get("proteins", {}), dict):
        raise ValueError(f"Reusable MSA manifest has an invalid schema: {path}")
    return value


def _write_manifest(cache_dir: Path, manifest: dict[str, Any]) -> None:
    path = cache_dir / MANIFEST_NAME
    temporary = cache_dir / f".{MANIFEST_NAME}.tmp"
    temporary.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    temporary.replace(path)
