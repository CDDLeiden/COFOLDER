from pathlib import Path

import pytest
from rdkit import Chem

from cofolder.modules.entities.ligand import iter_sdf_records
from cofolder.modules.input.compound_library import (
    CompoundMember,
    CompoundMemberFailure,
    DuplicateCompoundIdError,
    DuplicateIdPolicy,
    load_compound_library,
)
from cofolder.modules.input.ligand import LigandTarget

TARGET = LigandTarget("entity:1", 1, ("B",))


def _sdf_record(smiles: str, name: str = "") -> str:
    mol = Chem.MolFromSmiles(smiles)
    if name:
        mol.SetProp("_Name", name)
    return Chem.MolToMolBlock(mol)


def test_record_iterator_preserves_malformed_and_empty_records(tmp_path: Path):
    path = tmp_path / "mixed.sdf"
    path.write_text(
        _sdf_record("CCO", "valid")
        + "\n$$$$\n"
        + "\n$$$$\n"
        + "not a molblock\n$$$$\n",
        encoding="utf-8",
    )

    records = list(iter_sdf_records(path))

    assert [record.source_record_id for record in records] == [
        "record_000001",
        "record_000002",
        "record_000003",
    ]
    assert records[1].molblock == ""
    assert records[2].first_line > records[1].first_line


def test_mixed_sdf_returns_member_or_failure_for_every_record(tmp_path: Path):
    path = tmp_path / "mixed.sdf"
    path.write_text(
        _sdf_record("CCO", "valid")
        + "\n$$$$\n"
        + "broken\n$$$$\n",
        encoding="utf-8",
    )

    library = load_compound_library(path, target=TARGET)

    assert len(library.outcomes) == 2
    assert isinstance(library.outcomes[0], CompoundMember)
    assert isinstance(library.outcomes[1], CompoundMemberFailure)
    assert library.outcomes[1].exception.error_code == "molblock_record_parse_failed"
    assert library.outcomes[1].source.source_record_id == "record_000002"


def test_named_property_and_missing_id_fallback(tmp_path: Path):
    path = tmp_path / "ids.sdf"
    first = _sdf_record("CCO") + "\n>  <compound>\nA\n"
    second = _sdf_record("CCN")
    path.write_text(first + "\n$$$$\n" + second + "\n$$$$\n", encoding="utf-8")

    library = load_compound_library(path, target=TARGET, id_property="compound")

    assert [item.execution_id for item in library.outcomes] == ["A", "record_000002"]


def test_duplicate_policies_are_deterministic(tmp_path: Path):
    path = tmp_path / "compounds.csv"
    path.write_text("id,smiles\nA,CCO\nA__2,CCN\nA,CCC\n", encoding="utf-8")

    with pytest.raises(DuplicateCompoundIdError):
        load_compound_library(
            path, target=TARGET, smiles_column="smiles", id_column="id"
        )

    suffixed = load_compound_library(
        path,
        target=TARGET,
        smiles_column="smiles",
        id_column="id",
        duplicate_policy=DuplicateIdPolicy.SUFFIX,
    )
    indexed = load_compound_library(
        path,
        target=TARGET,
        smiles_column="smiles",
        id_column="id",
        duplicate_policy=DuplicateIdPolicy.SOURCE_INDEX,
    )

    assert [item.execution_id for item in suffixed.outcomes] == ["A", "A__2", "A__3"]
    assert [item.execution_id for item in indexed.outcomes] == [
        "record_000001",
        "record_000002",
        "record_000003",
    ]


def test_csv_sdf_and_mol_normalize_to_equivalent_chemistry(tmp_path: Path):
    csv_path = tmp_path / "library.csv"
    csv_path.write_text("id,smiles\nethanol,C(C)O\n", encoding="utf-8")
    molecule = Chem.MolFromSmiles("CCO")
    molecule.SetProp("_Name", "ethanol")
    molblock = Chem.MolToMolBlock(molecule)
    sdf_path = tmp_path / "library.sdf"
    sdf_path.write_text(molblock + "\n$$$$\n", encoding="utf-8")
    mol_path = tmp_path / "library.mol"
    mol_path.write_text(molblock, encoding="utf-8")

    csv = load_compound_library(
        csv_path, target=TARGET, smiles_column="smiles", id_column="id"
    )
    sdf = load_compound_library(sdf_path, target=TARGET)
    mol = load_compound_library(mol_path, target=TARGET)

    assert {
        library.outcomes[0].ligand.canonical_smiles
        for library in (csv, sdf, mol)
        if isinstance(library.outcomes[0], CompoundMember)
    } == {"CCO"}


def test_unsanitizable_record_is_retained_as_ligand_failure(tmp_path: Path):
    molecule = Chem.MolFromSmiles("[CH5]", sanitize=False)
    path = tmp_path / "invalid.sdf"
    path.write_text(Chem.MolToMolBlock(molecule) + "\n$$$$\n", encoding="utf-8")

    library = load_compound_library(path, target=TARGET)

    outcome = library.outcomes[0]
    assert isinstance(outcome, CompoundMemberFailure)
    assert outcome.exception.error_code == "ligand_validation_failed"
    assert outcome.exception.line == 1
    assert outcome.source.source_record_id == "record_000001"
