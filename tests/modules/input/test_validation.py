from __future__ import annotations

import pytest

from cofolder.modules.input import MsaValidationError, SequenceValidationError
from cofolder.modules.input.system import System
from cofolder.modules.input.validation import (
    SystemInputValidationError,
    WorkflowInputRequirements,
)
from cofolder.modules.runners.boltz1_runner import Boltz1Runner
from cofolder.modules.runners.boltz2_runner import Boltz2Runner
from cofolder.modules.runners.boltz_community_runner import BoltzCommunityRunner
from cofolder.modules.runners.openfold3_runner import OpenFold3Runner


def _system(constraints=None) -> System:
    return System(
        system={
            "sequences": [
                {"protein": {"id": "A", "sequence": "AC"}},
                {"dna": {"id": "D", "sequence": "AT"}},
                {"rna": {"id": "R", "sequence": "GU"}},
                {"ligand": {"id": "L", "smiles": "CCO"}},
            ],
            "constraints": constraints or [],
        }
    )


@pytest.mark.parametrize(
    "runner",
    [Boltz1Runner(), Boltz2Runner(), BoltzCommunityRunner(), OpenFold3Runner()],
)
def test_all_runners_declare_nucleic_acid_support(runner):
    assert runner.input_capabilities.entity_types == {"protein", "ligand", "dna", "rna"}
    runner.validate_system(_system(), {}, check_atom_names=True)


def test_boltz2_accepts_bond_pocket_and_contact_constraints():
    constraints = [
        {"bond": {"atom1": ["A", 2, "SG"], "atom2": ["L", 1, "O1"]}},
        {"pocket": {"binder": "L", "contacts": [["A", 1]], "max_distance": 7}},
        {"contact": {"token1": ["D", 1], "token2": ["R", 2], "max_distance": 5}},
    ]

    Boltz2Runner().validate_system(_system(constraints), {}, check_atom_names=False)
    BoltzCommunityRunner().validate_system(
        _system(constraints), {}, check_atom_names=False
    )


def test_boltz1_rejects_contact_and_non_default_pocket_distance():
    with pytest.raises(SystemInputValidationError, match="contact.*unsupported"):
        Boltz1Runner().validate_system(
            _system([{"contact": {"token1": ["A", 1], "token2": ["D", 1]}}]),
            {},
        )
    with pytest.raises(SystemInputValidationError, match="must be 6"):
        Boltz1Runner().validate_system(
            _system(
                [{"pocket": {"binder": "L", "contacts": [["A", 1]], "max_distance": 7}}]
            ),
            {},
        )


@pytest.mark.parametrize("kind", ["bond", "contact"])
def test_openfold3_rejects_constraints_it_cannot_consume(kind):
    payload = (
        {"atom1": ["A", 1, "N"], "atom2": ["L", 1, "C1"]}
        if kind == "bond"
        else {"token1": ["A", 1], "token2": ["D", 1]}
    )
    with pytest.raises(SystemInputValidationError, match=rf"{kind}.*unsupported"):
        OpenFold3Runner().validate_system(
            _system([{kind: payload}]), {}, check_atom_names=False
        )


def test_openfold3_rejects_multiple_pockets_and_force_key():
    pocket = {"binder": "L", "contacts": [["A", 1]]}
    with pytest.raises(SystemInputValidationError, match="at most 1"):
        OpenFold3Runner().validate_system(
            _system([{"pocket": pocket}, {"pocket": pocket}]), {}
        )
    with pytest.raises(SystemInputValidationError, match="uses 'force'"):
        OpenFold3Runner().validate_system(
            _system([{"pocket": {**pocket, "force": False}}]), {}
        )


@pytest.mark.parametrize(
    ("constraints", "message"),
    [
        (
            [{"bond": {"atom1": ["Z", 1, "N"], "atom2": ["L", 1, "C1"]}}],
            "unknown chain",
        ),
        (
            [{"bond": {"atom1": ["A", 9, "N"], "atom2": ["L", 1, "C1"]}}],
            "outside chain",
        ),
        (
            [{"bond": {"atom1": ["A", 1, "SG"], "atom2": ["L", 1, "C1"]}}],
            "does not exist",
        ),
    ],
)
def test_constraint_references_are_actionable(constraints, message):
    with pytest.raises(SystemInputValidationError, match=message):
        Boltz2Runner().validate_system(_system(constraints), {}, check_atom_names=True)


def test_duplicate_chain_ids_are_rejected():
    value = _system()
    value.system["sequences"][1]["dna"]["id"] = "A"
    with pytest.raises(SystemInputValidationError, match="declared more than once"):
        Boltz2Runner().validate_system(value, {}, check_atom_names=False)


@pytest.mark.parametrize(
    "representation",
    [
        {"smiles": "CCO"},
        {"ccd": "ATP"},
        {"ccd_codes": ["ATP", "MG"]},
    ],
)
def test_ligand_accepts_exactly_one_nonempty_representation(representation):
    value = System(
        system={"sequences": [{"ligand": {"id": "L", **representation}}]}
    )

    Boltz2Runner().validate_system(value, {}, check_atom_names=False)


@pytest.mark.parametrize(
    ("representation", "message"),
    [
        ({}, "requires exactly one"),
        ({"smiles": ""}, "requires exactly one"),
        ({"smiles": "   "}, "requires exactly one"),
        ({"ccd": ""}, "requires exactly one"),
        ({"ccd": "   "}, "requires exactly one"),
        ({"ccd_codes": []}, "requires exactly one"),
        ({"ccd_codes": [""]}, "requires exactly one"),
        ({"smiles": "CCO", "ccd": "ETH"}, "'smiles', 'ccd'"),
        ({"smiles": "CCO", "ccd_codes": ["ETH"]}, "'smiles', 'ccd_codes'"),
        ({"ccd": "ETH", "ccd_codes": ["ETH"]}, "'ccd', 'ccd_codes'"),
    ],
)
def test_ligand_rejects_missing_or_conflicting_representations(
    representation, message
):
    value = System(
        system={"sequences": [{"ligand": {"id": "L", **representation}}]}
    )

    with pytest.raises(SystemInputValidationError, match=message):
        Boltz2Runner().validate_system(value, {}, check_atom_names=False)


@pytest.mark.parametrize(
    ("entity_type", "sequence", "normalized"),
    [
        ("protein", "ac dX\n", "ACDX"),
        ("dna", "acg tn", "ACGTN"),
        ("rna", "acg un", "ACGUN"),
    ],
)
def test_sequence_content_is_normalized_by_entity(entity_type, sequence, normalized):
    value = System(
        system={"sequences": [{entity_type: {"id": "A", "sequence": sequence}}]}
    )

    Boltz2Runner().validate_system(value, {}, check_atom_names=False)

    assert value.system["sequences"][0][entity_type]["sequence"] == normalized


@pytest.mark.parametrize(
    ("entity_type", "sequence", "symbol", "position"),
    [("protein", "ACB", "B", 3), ("dna", "ACU", "U", 3), ("rna", "ACT", "T", 3)],
)
def test_invalid_sequence_reports_chain_position_and_symbol(
    entity_type, sequence, symbol, position, temp_dir
):
    value = System(
        system={"sequences": [{entity_type: {"id": "A", "sequence": sequence}}]}
    )

    with pytest.raises(SequenceValidationError) as caught:
        Boltz2Runner().validate_system(
            value,
            {},
            check_atom_names=False,
            source_path=temp_dir / "system.yaml",
        )

    assert caught.value.chain_id == "A"
    assert caught.value.entity_id == "entity:0"
    assert caught.value.details["residue_position"] == position
    assert caught.value.details["character"] == symbol
    assert symbol in str(caught.value)
    assert f"residue {position}" in str(caught.value)


def test_empty_normalized_sequence_is_rejected(temp_dir):
    value = System(
        system={"sequences": [{"protein": {"id": "A", "sequence": " \n\t"}}]}
    )

    with pytest.raises(SequenceValidationError, match="sequence is empty"):
        Boltz2Runner().validate_system(
            value,
            {},
            check_atom_names=False,
            source_path=temp_dir / "system.yaml",
        )


def test_workflow_composition_requirements_are_distinct():
    protein_only = System(
        system={"sequences": [{"protein": {"id": "A", "sequence": "AC"}}]}
    )
    ligand_only = System(
        system={"sequences": [{"ligand": {"id": "L", "smiles": "CCO"}}]}
    )
    requirements = WorkflowInputRequirements(require_protein=True, require_ligand=True)

    with pytest.raises(SystemInputValidationError, match="ligand entity"):
        Boltz2Runner().validate_system(
            protein_only, {}, requirements=requirements, check_atom_names=False
        )
    with pytest.raises(SystemInputValidationError, match="protein entity"):
        Boltz2Runner().validate_system(
            ligand_only, {}, requirements=requirements, check_atom_names=False
        )


@pytest.mark.parametrize("suffix", [".txt", ".fasta"])
def test_declared_msa_rejects_unsupported_formats(temp_dir, suffix):
    msa = temp_dir / f"query{suffix}"
    msa.write_text(">query\nAC\n", encoding="utf-8")
    value = System(
        system={
            "sequences": [
                {"protein": {"id": "A", "sequence": "AC", "msa": str(msa)}}
            ]
        }
    )

    with pytest.raises(MsaValidationError, match="Unsupported MSA format"):
        Boltz2Runner().validate_system(value, {}, check_atom_names=False)


@pytest.mark.parametrize(
    ("suffix", "content"),
    [(".a3m", ">query\nA-cC\n"), (".csv", "key,sequence\n-1,AC\n")],
)
def test_declared_msa_query_must_match_owning_protein(temp_dir, suffix, content):
    msa = temp_dir / f"query{suffix}"
    msa.write_text(content, encoding="utf-8")
    value = System(
        system={
            "sequences": [
                {"protein": {"id": "A", "sequence": "AC", "msa": str(msa)}}
            ]
        }
    )

    validated = Boltz2Runner().validate_system(value, {}, check_atom_names=False)
    assert validated.system.system["sequences"][0]["protein"]["msa"] == str(msa)

    value.system["sequences"][0]["protein"]["sequence"] = "AA"
    with pytest.raises(MsaValidationError, match="does not match"):
        Boltz2Runner().validate_system(value, {}, check_atom_names=False)


@pytest.mark.parametrize(
    ("suffix", "content"),
    [
        (".csv", "sequence\nAC\n"),
        (".csv", "key,sequence\n,AC\n"),
        (".a3m", "AC\n"),
        (".a3m", ">query\nAC*\n"),
    ],
)
def test_declared_msa_rejects_malformed_structure(temp_dir, suffix, content):
    msa = temp_dir / f"query{suffix}"
    msa.write_text(content, encoding="utf-8")
    value = System(
        system={
            "sequences": [
                {"protein": {"id": "A", "sequence": "AC", "msa": str(msa)}}
            ]
        }
    )

    with pytest.raises(MsaValidationError, match="unreadable or malformed"):
        Boltz2Runner().validate_system(value, {}, check_atom_names=False)
