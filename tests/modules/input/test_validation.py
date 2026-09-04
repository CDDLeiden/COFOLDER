from __future__ import annotations

import pytest

from cofolder.modules.input.system import System
from cofolder.modules.input.validation import SystemInputValidationError
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
