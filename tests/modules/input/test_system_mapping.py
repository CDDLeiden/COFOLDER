"""Tests for CSV-driven system mappings."""

import pytest

from cofolder.modules.input.system import System
from cofolder.modules.input.system_mapping import (
    InheritedMsaConflictError,
    MappedSystemMember,
    MappedSystemMemberFailure,
    SystemMappingError,
    build_mapped_system,
    load_mapped_system_library,
    parse_system_mappings,
)
from cofolder.modules.input.compound_library import (
    CompoundLibrarySchemaError,
    DuplicateIdPolicy,
)


def _base_system(*, with_msa=True):
    return System(
        system={
            "sequences": [
                {
                    "protein": {
                        "id": "A",
                        "sequence": "AAAA",
                        **({"msa": "/template/a.a3m"} if with_msa else {}),
                    }
                },
                {"ligand": {"id": "B", "ccd": "EDO", "smiles": None}},
            ],
            "metadata": {"literal.key": "old", "enabled": False},
            "constraints": [],
        }
    )


def test_build_mapped_system_changes_multiple_fields():
    mappings = parse_system_mappings(
        [
            "sequence=sequences.0.protein.sequence",
            "smiles=sequences.1.ligand.smiles",
            "constraints=constraints",
        ]
    )

    result = build_mapped_system(
        _base_system(with_msa=False),
        mappings,
        {
            "sequence": "CCCC",
            "smiles": '"CCN"',
            "constraints": '[{"pocket": {"binder": "B"}}]',
        },
        row_number=1,
    ).system

    assert result["sequences"][0]["protein"]["sequence"] == "CCCC"
    assert "msa" not in result["sequences"][0]["protein"]
    assert result["sequences"][1]["ligand"]["smiles"] == "CCN"
    assert "ccd" not in result["sequences"][1]["ligand"]
    assert result["constraints"] == [{"pocket": {"binder": "B"}}]
    assert _base_system().system["sequences"][0]["protein"]["sequence"] == "AAAA"


def test_changed_sequence_rejects_inherited_fixed_msa():
    mappings = parse_system_mappings(["sequence=sequences.0.protein.sequence"])
    with pytest.raises(InheritedMsaConflictError, match="inherits the fixed MSA"):
        build_mapped_system(
            _base_system(),
            mappings,
            {"sequence": "CCCC"},
            row_number=2,
            row_id="variant_2",
        )


def test_normalized_unchanged_sequence_may_keep_fixed_msa():
    mappings = parse_system_mappings(["sequence=sequences.0.protein.sequence"])
    result = build_mapped_system(
        _base_system(),
        mappings,
        {"sequence": " a a a a "},
        row_number=1,
    )
    assert result.system["sequences"][0]["protein"]["msa"] == "/template/a.a3m"


def test_direct_msa_mapping_allows_sequence_change():
    mappings = parse_system_mappings(
        [
            "sequence=sequences.0.protein.sequence",
            "msa=sequences.0.protein.msa",
        ]
    )
    result = build_mapped_system(
        _base_system(),
        mappings,
        {"sequence": "CCCC", "msa": "new.a3m"},
        row_number=1,
    )
    assert result.system["sequences"][0]["protein"]["msa"] == "new.a3m"


def test_structured_protein_mapping_can_supply_sequence_and_msa_together():
    mappings = parse_system_mappings(["protein=sequences.0.protein"])
    result = build_mapped_system(
        _base_system(),
        mappings,
        {"protein": '{"id":"A","sequence":"CCCC","msa":"new.a3m"}'},
        row_number=1,
    )
    assert result.system["sequences"][0]["protein"] == {
        "id": "A",
        "sequence": "CCCC",
        "msa": "new.a3m",
    }


def test_empty_msa_does_not_conflict_with_sequence_iteration():
    base = _base_system()
    base.system["sequences"][0]["protein"]["msa"] = "empty"
    result = build_mapped_system(
        base,
        parse_system_mappings(["sequence=sequences.0.protein.sequence"]),
        {"sequence": "CCCC"},
        row_number=1,
    )
    assert result.system["sequences"][0]["protein"]["msa"] == "empty"


def test_mapping_supports_escaped_key_and_preserves_string_exactly():
    mappings = parse_system_mappings([r"value=metadata.literal\.key"])
    result = build_mapped_system(
        _base_system(), mappings, {"value": "  001  "}, row_number=1
    )
    assert result.system["metadata"]["literal.key"] == "  001  "


@pytest.mark.parametrize(
    "values, message",
    [
        (["a=sequences", "b=sequences.0"], "overlapping"),
        (
            ["a=sequences.1.ligand.smiles", "b=sequences.1.ligand.ccd"],
            "competing ligand representations",
        ),
    ],
)
def test_mapping_rejects_conflicting_destinations(values, message):
    with pytest.raises(SystemMappingError, match=message):
        parse_system_mappings(values)


def test_mapping_rejects_blank_and_wrong_json_type():
    mapping = parse_system_mappings(["enabled=metadata.enabled"])
    with pytest.raises(SystemMappingError, match="blank"):
        build_mapped_system(_base_system(), mapping, {"enabled": ""}, row_number=2)
    with pytest.raises(SystemMappingError, match="requires bool"):
        build_mapped_system(
            _base_system(), mapping, {"enabled": '"false"'}, row_number=2
        )


def test_load_mapped_library_preserves_rows_and_isolates_cell_failures(tmp_path):
    path = tmp_path / "variants.csv"
    path.write_text(
        "id,sequence,enabled\n"
        "first,CCCC,true\n"
        "second,DDDD,\n",
        encoding="utf-8",
    )

    library = load_mapped_system_library(
        path,
        base_system=_base_system(with_msa=False),
        mapping_values=[
            "sequence=sequences.0.protein.sequence",
            "enabled=metadata.enabled",
        ],
        id_column="id",
        duplicate_policy=DuplicateIdPolicy.REJECT,
    )

    assert isinstance(library.outcomes[0], MappedSystemMember)
    assert library.outcomes[0].system.system["metadata"]["enabled"] is True
    assert isinstance(library.outcomes[1], MappedSystemMemberFailure)
    assert "blank" in str(library.outcomes[1].exception)


def test_load_mapped_library_rejects_duplicate_headers(tmp_path):
    path = tmp_path / "variants.csv"
    path.write_text("id,sequence,sequence\none,AAAA,BBBB\n", encoding="utf-8")
    with pytest.raises(CompoundLibrarySchemaError, match="duplicate headers"):
        load_mapped_system_library(
            path,
            base_system=_base_system(),
            mapping_values=["sequence=sequences.0.protein.sequence"],
            id_column="id",
            duplicate_policy=DuplicateIdPolicy.REJECT,
        )
