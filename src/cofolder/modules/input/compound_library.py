"""Canonical, record-oriented compound-library ingestion."""

from __future__ import annotations

from collections import Counter
from collections.abc import Collection, Mapping
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import Literal, TypeAlias

import pandas as pd

from cofolder.modules.contracts.models import JSONValue
from cofolder.modules.entities.ligand import (
    RawMolblockRecord,
    iter_sdf_records,
    parse_molblock,
    read_molblock_record,
)
from cofolder.modules.input.config import InputValidationError, LigandValidationError
from cofolder.modules.input.ligand import (
    LigandSourceIdentity,
    LigandTarget,
    NormalizedLigand,
    validate_smiles,
)


class CompoundLibraryFormat(StrEnum):
    CSV = "csv"
    SDF = "sdf"
    MOL = "mol"


class DuplicateIdPolicy(StrEnum):
    REJECT = "reject"
    SUFFIX = "suffix"
    SOURCE_INDEX = "source_index"


class CompoundLibraryError(InputValidationError):
    error_code = "compound_library_failed"


class UnsupportedCompoundLibraryFormatError(CompoundLibraryError):
    error_code = "compound_library_format_unsupported"


class CompoundLibrarySchemaError(CompoundLibraryError):
    error_code = "compound_library_schema_failed"


class DuplicateCompoundIdError(CompoundLibraryError):
    error_code = "duplicate_compound_id"


class MolblockRecordParseError(LigandValidationError):
    error_code = "molblock_record_parse_failed"


@dataclass(frozen=True, slots=True)
class CompoundSourceRecord:
    source_format: CompoundLibraryFormat
    source_path: Path
    source_record_index: int
    source_record_id: str
    original_id: str | None
    metadata: Mapping[str, JSONValue]
    molblock: str | None = None
    first_line: int | None = None


@dataclass(frozen=True, slots=True)
class CompoundMember:
    source: CompoundSourceRecord
    execution_id: str
    execution_directory: str
    ligand: NormalizedLigand
    conformer_molblock: str | None
    coordinate_mode: Literal["source", "generated", "native_smiles"]


@dataclass(frozen=True, slots=True)
class CompoundMemberFailure:
    source: CompoundSourceRecord
    execution_id: str
    execution_directory: str
    exception: InputValidationError


CompoundLibraryOutcome: TypeAlias = CompoundMember | CompoundMemberFailure


@dataclass(frozen=True, slots=True)
class CompoundLibrary:
    source_path: Path
    source_format: CompoundLibraryFormat
    outcomes: tuple[CompoundLibraryOutcome, ...]


def _format_for(path: Path, requested: CompoundLibraryFormat | str | None) -> CompoundLibraryFormat:
    if requested is not None:
        try:
            return CompoundLibraryFormat(requested)
        except ValueError as exc:
            raise UnsupportedCompoundLibraryFormatError(
                f"Unsupported compound library format: {requested!r}.", source_path=path
            ) from exc
    suffix = path.suffix.lower()
    inferred = {
        ".csv": CompoundLibraryFormat.CSV,
        ".sdf": CompoundLibraryFormat.SDF,
        ".sd": CompoundLibraryFormat.SDF,
        ".mol": CompoundLibraryFormat.MOL,
    }.get(suffix)
    if inferred is None:
        raise UnsupportedCompoundLibraryFormatError(
            f"Cannot infer compound library format from {path.name!r}; use --library_format.",
            source_path=path,
        )
    return inferred


def _json_scalar(value: object) -> JSONValue:
    if value is None or pd.isna(value):
        return None
    if hasattr(value, "item"):
        value = value.item()  # numpy scalar
    if isinstance(value, (bool, int, float, str)):
        return value
    return str(value)


def _clean_id(value: object) -> str | None:
    if value is None or pd.isna(value):
        return None
    text = str(value).strip()
    return text or None


def _csv_records(
    path: Path,
    *,
    smiles_column: str | None,
    id_column: str | None,
    metadata_fields: Collection[str],
) -> tuple[tuple[CompoundSourceRecord, object], ...]:
    if not smiles_column or not id_column:
        raise CompoundLibrarySchemaError(
            "CSV libraries require --smiles_column and --col_id.", source_path=path
        )
    try:
        frame = pd.read_csv(path)
    except Exception as exc:
        raise CompoundLibrarySchemaError(
            f"Unable to read CSV compound library: {exc}", source_path=path
        ) from exc
    required = [id_column, smiles_column, *metadata_fields]
    missing = [name for name in required if name not in frame.columns]
    if missing:
        raise CompoundLibrarySchemaError(
            f"CSV missing required columns: {', '.join(missing)}", source_path=path
        )
    records = []
    for index, (_, row) in enumerate(frame.iterrows(), 1):
        record_id = f"record_{index:06d}"
        metadata = {name: _json_scalar(row[name]) for name in frame.columns}
        records.append(
            (
                CompoundSourceRecord(
                    CompoundLibraryFormat.CSV,
                    path,
                    index,
                    record_id,
                    _clean_id(row[id_column]),
                    metadata,
                ),
                row[smiles_column],
            )
        )
    return tuple(records)


def _structure_records(
    path: Path,
    *,
    source_format: CompoundLibraryFormat,
    id_property: str,
    metadata_fields: Collection[str],
) -> tuple[tuple[CompoundSourceRecord, RawMolblockRecord], ...]:
    raw_records = (
        tuple(iter_sdf_records(path))
        if source_format is CompoundLibraryFormat.SDF
        else (read_molblock_record(path),)
    )
    records = []
    for raw in raw_records:
        original_id = raw.title if id_property == "_Name" else raw.properties.get(id_property)
        metadata = {name: raw.properties.get(name) for name in metadata_fields}
        records.append(
            (
                CompoundSourceRecord(
                    source_format,
                    path,
                    raw.index,
                    raw.source_record_id,
                    _clean_id(original_id),
                    metadata,
                    raw.molblock,
                    raw.first_line,
                ),
                raw,
            )
        )
    return tuple(records)


def _execution_ids(
    records: Collection[CompoundSourceRecord], policy: DuplicateIdPolicy
) -> tuple[str, ...]:
    base_ids = [record.original_id or record.source_record_id for record in records]
    if policy is DuplicateIdPolicy.SOURCE_INDEX:
        return tuple(record.source_record_id for record in records)
    if policy is DuplicateIdPolicy.REJECT:
        duplicates = sorted(name for name, count in Counter(base_ids).items() if count > 1)
        if duplicates:
            locations = {
                name: [
                    record.source_record_index
                    for record, base in zip(records, base_ids, strict=True)
                    if base == name
                ]
                for name in duplicates
            }
            raise DuplicateCompoundIdError(
                "Duplicate compound IDs: "
                + ", ".join(
                    f"{name} (records {', '.join(map(str, locations[name]))})"
                    for name in duplicates
                ),
                source_path=next(iter(records)).source_path if records else None,
            )
        return tuple(base_ids)
    allocated: set[str] = set()
    resolved: list[str] = []
    for base in base_ids:
        candidate = base
        occurrence = 2
        while candidate in allocated:
            candidate = f"{base}__{occurrence}"
            occurrence += 1
        allocated.add(candidate)
        resolved.append(candidate)
    return tuple(resolved)


def load_compound_library(
    path: Path,
    *,
    target: LigandTarget,
    source_format: CompoundLibraryFormat | None = None,
    smiles_column: str | None = None,
    id_column: str | None = None,
    id_property: str = "_Name",
    metadata_fields: Collection[str] = (),
    duplicate_policy: DuplicateIdPolicy = DuplicateIdPolicy.REJECT,
) -> CompoundLibrary:
    path = Path(path)
    if not path.is_file():
        raise CompoundLibrarySchemaError(
            f"Compound library is not a readable file: {path}", source_path=path
        )
    resolved_format = _format_for(path, source_format)
    try:
        policy = DuplicateIdPolicy(duplicate_policy)
    except ValueError as exc:
        raise CompoundLibrarySchemaError(
            f"Unsupported duplicate ID policy: {duplicate_policy!r}.", source_path=path
        ) from exc
    if resolved_format is CompoundLibraryFormat.CSV:
        pending = _csv_records(
            path,
            smiles_column=smiles_column,
            id_column=id_column,
            metadata_fields=metadata_fields,
        )
    else:
        pending = _structure_records(
            path,
            source_format=resolved_format,
            id_property=id_property,
            metadata_fields=metadata_fields,
        )
    if not pending:
        raise CompoundLibrarySchemaError(
            "Compound library must contain at least one source record.", source_path=path
        )
    execution_ids = _execution_ids([item[0] for item in pending], policy)
    outcomes: list[CompoundLibraryOutcome] = []
    for (record, payload), execution_id in zip(pending, execution_ids, strict=True):
        source = LigandSourceIdentity(
            target.entity_id, target.chain_ids, record.source_record_id
        )
        directory = f"compound_{record.source_record_index:06d}"
        try:
            if resolved_format is CompoundLibraryFormat.CSV:
                normalized = validate_smiles(
                    payload,
                    source=source,
                    source_path=path,
                    field_path=(record.source_record_index, str(smiles_column)),
                )
                conformer = None
                coordinate_mode = "native_smiles"
            else:
                normalized = parse_molblock(
                    payload, source=source, source_path=path
                )
                conformer = record.molblock
                coordinate_mode = "source"
            outcomes.append(
                CompoundMember(
                    record,
                    execution_id,
                    directory,
                    normalized,
                    conformer,
                    coordinate_mode,
                )
            )
        except InputValidationError as exc:
            outcomes.append(
                CompoundMemberFailure(record, execution_id, directory, exc)
            )
    return CompoundLibrary(path, resolved_format, tuple(outcomes))


__all__ = [
    "CompoundLibrary",
    "CompoundLibraryError",
    "CompoundLibraryFormat",
    "CompoundLibraryOutcome",
    "CompoundLibrarySchemaError",
    "CompoundMember",
    "CompoundMemberFailure",
    "CompoundSourceRecord",
    "DuplicateCompoundIdError",
    "DuplicateIdPolicy",
    "MolblockRecordParseError",
    "UnsupportedCompoundLibraryFormatError",
    "load_compound_library",
]
