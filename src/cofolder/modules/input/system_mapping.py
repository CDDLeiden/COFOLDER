"""CSV-to-system mappings for row-oriented parameter screens."""

from __future__ import annotations

import copy
import csv
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from cofolder.modules.input.compound_library import (
    CompoundLibraryFormat,
    CompoundLibrarySchemaError,
    CompoundSourceRecord,
    DuplicateIdPolicy,
    _execution_ids,
)
from cofolder.modules.input.config import InputValidationError
from cofolder.modules.input.system import System


@dataclass(frozen=True, slots=True)
class SystemMapping:
    column: str
    path: tuple[str | int, ...]
    path_text: str


@dataclass(frozen=True, slots=True)
class MappedSystemMember:
    source: CompoundSourceRecord
    execution_id: str
    execution_directory: str
    system: System
    coordinate_mode: None = None
    conformer_molblock: None = None


@dataclass(frozen=True, slots=True)
class MappedSystemMemberFailure:
    source: CompoundSourceRecord
    execution_id: str
    execution_directory: str
    exception: InputValidationError
    coordinate_mode: None = None
    conformer_molblock: None = None


@dataclass(frozen=True, slots=True)
class MappedSystemLibrary:
    source_path: Path
    source_format: CompoundLibraryFormat
    outcomes: tuple[MappedSystemMember | MappedSystemMemberFailure, ...]
    mappings: tuple[SystemMapping, ...]


class SystemMappingError(InputValidationError):
    error_code = "system_mapping_failed"


class InheritedMsaConflictError(SystemMappingError):
    """A changed protein sequence would inherit a template alignment."""

    error_code = "inherited_msa_conflict"


def _split_path(text: str) -> tuple[str | int, ...]:
    tokens: list[str] = []
    token: list[str] = []
    escaped = False
    for character in text:
        if escaped:
            if character not in {".", "\\"}:
                raise SystemMappingError(
                    f"Invalid escape in mapping path {text!r}; only dots and backslashes may be escaped."
                )
            token.append(character)
            escaped = False
        elif character == "\\":
            escaped = True
        elif character == ".":
            if not token:
                raise SystemMappingError(f"Mapping path {text!r} contains an empty segment.")
            tokens.append("".join(token))
            token = []
        else:
            token.append(character)
    if escaped or not token:
        raise SystemMappingError(f"Invalid mapping path {text!r}.")
    tokens.append("".join(token))
    return tuple(int(value) if value.isdigit() else value for value in tokens)


def parse_system_mappings(values: list[str] | tuple[str, ...]) -> tuple[SystemMapping, ...]:
    mappings: list[SystemMapping] = []
    for raw in values:
        column, separator, path_text = str(raw).partition("=")
        column = column.strip()
        path_text = path_text.strip()
        if not separator or not column or not path_text:
            raise SystemMappingError(
                f"Invalid mapping {raw!r}; expected COLUMN=YAML_PATH."
            )
        mappings.append(SystemMapping(column, _split_path(path_text), path_text))
    paths = [mapping.path for mapping in mappings]
    for index, path in enumerate(paths):
        for other in paths[index + 1 :]:
            common = min(len(path), len(other))
            if path[:common] == other[:common]:
                relation = "duplicate" if path == other else "overlapping"
                raise SystemMappingError(
                    f"Mappings have {relation} destinations: "
                    f"{_display_path(path)!r} and {_display_path(other)!r}."
                )
    columns = [mapping.column for mapping in mappings]
    duplicates = sorted({value for value in columns if columns.count(value) > 1})
    if duplicates:
        raise SystemMappingError(
            "CSV columns may be mapped only once: " + ", ".join(duplicates)
        )
    chemistry_destinations: dict[tuple[str | int, ...], list[str]] = {}
    for mapping in mappings:
        chemistry = _ligand_chemistry_parent(mapping.path)
        if chemistry:
            parent, representation = chemistry
            chemistry_destinations.setdefault(parent, []).append(representation)
    conflicts = [
        (_display_path(parent), values)
        for parent, values in chemistry_destinations.items()
        if len(values) > 1
    ]
    if conflicts:
        parent, values = conflicts[0]
        raise SystemMappingError(
            f"Mappings select competing ligand representations at {parent!r}: "
            + ", ".join(values)
            + "."
        )
    return tuple(mappings)


def _display_path(path: tuple[str | int, ...]) -> str:
    return ".".join(str(value) for value in path)


def _target(root: Any, path: tuple[str | int, ...]) -> Any:
    current = root
    for part in path:
        if isinstance(current, dict) and isinstance(part, str) and part in current:
            current = current[part]
        elif isinstance(current, list) and isinstance(part, int) and 0 <= part < len(current):
            current = current[part]
        else:
            raise SystemMappingError(
                f"Mapping destination {_display_path(path)!r} does not exist in the system template."
            )
    return current


def _convert_cell(raw: str, template: Any, *, column: str, row: int) -> Any:
    if raw.strip() == "":
        raise SystemMappingError(f"Mapped column {column!r} is blank in CSV row {row}.")
    if isinstance(template, str):
        return raw
    try:
        value = json.loads(raw)
    except json.JSONDecodeError as exc:
        raise SystemMappingError(
            f"Mapped column {column!r} in CSV row {row} must contain valid JSON: {exc.msg}."
        ) from exc
    if template is None:
        return value
    expected = type(template)
    valid = (
        (expected is bool and type(value) is bool)
        or (expected is int and type(value) is int)
        or (expected is float and type(value) in {int, float})
        or (expected is list and isinstance(value, list))
        or (expected is dict and isinstance(value, dict))
    )
    if not valid:
        raise SystemMappingError(
            f"Mapped column {column!r} in CSV row {row} has type "
            f"{type(value).__name__}; destination requires {expected.__name__}."
        )
    return float(value) if expected is float else value


def _set_target(root: Any, path: tuple[str | int, ...], value: Any) -> None:
    parent = root
    for part in path[:-1]:
        parent = parent[part]
    parent[path[-1]] = value


def _ligand_chemistry_parent(path: tuple[str | int, ...]) -> tuple[tuple[str | int, ...], str] | None:
    if len(path) >= 2 and path[-1] in {"smiles", "ccd", "ccd_codes"}:
        return path[:-1], str(path[-1])
    return None


def _normalized_sequence(payload: dict[str, Any]) -> str:
    raw = payload.get("sequence") or payload.get("fasta") or ""
    return re.sub(r"[^A-Za-z]", "", str(raw)).upper()


def _mapping_covers(
    mapping_path: tuple[str | int, ...], target_path: tuple[str | int, ...]
) -> bool:
    return len(mapping_path) <= len(target_path) and (
        target_path[: len(mapping_path)] == mapping_path
    )


def _reject_inherited_msa_conflicts(
    root: dict[str, Any],
    base: dict[str, Any],
    mappings: tuple[SystemMapping, ...],
    *,
    row_number: int,
    row_id: str | None,
) -> None:
    sequences = root.get("sequences", [])
    base_sequences = base.get("sequences", [])
    for index, entry in enumerate(sequences):
        if (
            not isinstance(entry, dict)
            or "protein" not in entry
            or index >= len(base_sequences)
            or not isinstance(base_sequences[index], dict)
            or "protein" not in base_sequences[index]
        ):
            continue
        payload = entry["protein"]
        old = base_sequences[index].get("protein", {})
        if not isinstance(payload, dict) or not isinstance(old, dict):
            continue
        sequence_changed = _normalized_sequence(payload) != _normalized_sequence(old)
        inherited_msa = old.get("msa")
        if (
            not sequence_changed
            or inherited_msa is None
            or not str(inherited_msa).strip()
            or str(inherited_msa).strip().lower() == "empty"
        ):
            continue
        msa_path = ("sequences", index, "protein", "msa")
        msa_is_mapped = any(
            _mapping_covers(mapping.path, msa_path) for mapping in mappings
        )
        if msa_is_mapped or payload.get("msa") != inherited_msa:
            continue
        sequence_paths = tuple(
            mapping.path_text
            for mapping in mappings
            if any(
                _mapping_covers(
                    mapping.path,
                    ("sequences", index, "protein", sequence_key),
                )
                for sequence_key in ("sequence", "fasta")
            )
        )
        raw_ids = payload.get("id")
        chain_ids = raw_ids if isinstance(raw_ids, list) else [raw_ids]
        chains = ", ".join(str(value) for value in chain_ids if value is not None)
        identity = row_id or f"CSV row {row_number}"
        sequence_mapping = ", ".join(sequence_paths) or "a structured sequence mapping"
        raise InheritedMsaConflictError(
            f"Screen setup rejected: row {identity!r} changes protein chain(s) "
            f"{chains or f'entity:{index}'} through {sequence_mapping!r} but inherits "
            f"the fixed MSA {str(inherited_msa)!r} from the system template. Map a "
            f"matching alignment to 'sequences.{index}.protein.msa', or remove the "
            "template MSA so COFOLDER can generate and cache one per unique sequence."
        )


def build_mapped_system(
    base: System,
    mappings: tuple[SystemMapping, ...],
    row: dict[str, str],
    *,
    row_number: int,
    row_id: str | None = None,
) -> System:
    result = copy.deepcopy(base.system)
    for mapping in mappings:
        template_value = _target(base.system, mapping.path)
        value = _convert_cell(row[mapping.column], template_value, column=mapping.column, row=row_number)
        chemistry = _ligand_chemistry_parent(mapping.path)
        if chemistry:
            parent_path, selected = chemistry
            parent = _target(result, parent_path)
            if not isinstance(parent, dict):
                raise SystemMappingError(f"Ligand chemistry parent {_display_path(parent_path)!r} is not a mapping.")
            for key in {"smiles", "ccd", "ccd_codes"} - {selected}:
                parent.pop(key, None)
        _set_target(result, mapping.path, value)
    _reject_inherited_msa_conflicts(
        result,
        base.system,
        mappings,
        row_number=row_number,
        row_id=row_id,
    )
    return System(system=result)


def load_mapped_system_library(
    path: Path,
    *,
    base_system: System,
    mapping_values: list[str] | tuple[str, ...],
    id_column: str,
    duplicate_policy: DuplicateIdPolicy,
) -> MappedSystemLibrary:
    mappings = parse_system_mappings(mapping_values)
    for mapping in mappings:
        _target(base_system.system, mapping.path)
    try:
        with Path(path).open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.reader(handle)
            rows = list(reader)
    except (OSError, csv.Error) as exc:
        raise CompoundLibrarySchemaError(f"Unable to read mapped CSV: {exc}", source_path=path) from exc
    if not rows:
        raise CompoundLibrarySchemaError("Mapped CSV must contain a header and at least one row.", source_path=path)
    header = rows[0]
    duplicates = sorted({name for name in header if header.count(name) > 1})
    if duplicates:
        raise CompoundLibrarySchemaError("CSV contains duplicate headers: " + ", ".join(duplicates), source_path=path)
    required = [id_column, *(mapping.column for mapping in mappings)]
    missing = [name for name in required if name not in header]
    if missing:
        raise CompoundLibrarySchemaError("CSV missing required columns: " + ", ".join(missing), source_path=path)
    if len(rows) == 1:
        raise CompoundLibrarySchemaError("Mapped CSV must contain at least one data row.", source_path=path)
    records: list[tuple[CompoundSourceRecord, dict[str, str]]] = []
    for index, values in enumerate(rows[1:], 1):
        if len(values) != len(header):
            raise CompoundLibrarySchemaError(
                f"CSV row {index} has {len(values)} fields; expected {len(header)}.", source_path=path
            )
        row = dict(zip(header, values, strict=True))
        original_id = row[id_column].strip() or None
        metadata = dict(row)
        records.append((CompoundSourceRecord(CompoundLibraryFormat.CSV, Path(path), index, f"record_{index:06d}", original_id, metadata), row))
    execution_ids = _execution_ids([record for record, _ in records], DuplicateIdPolicy(duplicate_policy))
    outcomes: list[MappedSystemMember | MappedSystemMemberFailure] = []
    msa_conflicts: list[str] = []
    for (record, row), execution_id in zip(records, execution_ids, strict=True):
        directory = f"compound_{record.source_record_index:06d}"
        try:
            row_system = build_mapped_system(
                base_system,
                mappings,
                row,
                row_number=record.source_record_index,
                row_id=record.original_id or record.source_record_id,
            )
            outcomes.append(MappedSystemMember(record, execution_id, directory, row_system))
        except InheritedMsaConflictError as exc:
            msa_conflicts.append(str(exc))
        except InputValidationError as exc:
            outcomes.append(MappedSystemMemberFailure(record, execution_id, directory, exc))
    if msa_conflicts:
        raise InheritedMsaConflictError("\n".join(msa_conflicts), source_path=Path(path))
    return MappedSystemLibrary(Path(path), CompoundLibraryFormat.CSV, tuple(outcomes), mappings)


__all__ = [
    "MappedSystemLibrary", "MappedSystemMember", "MappedSystemMemberFailure",
    "InheritedMsaConflictError", "SystemMapping", "SystemMappingError", "build_mapped_system",
    "load_mapped_system_library", "parse_system_mappings",
]
