"""Tests for structure analytics compatibility helpers."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

from cofolder.modules.analytics.reference_ifp import (
    AtomIdentity,
    IFPTaxonomy,
    InteractionEvent,
    InteractionFingerprint,
    InteractionGeometry,
    InteractionKey,
    LigandIdentity,
    ResidueIdentity,
)
from cofolder.modules.analytics.structure import Structure, _run_pdb2pqr


def test_run_pdb2pqr_uses_legacy_entrypoint(monkeypatch):
    calls: list[list[str]] = []

    def legacy_runner(args):
        calls.append(args)

    monkeypatch.setattr(
        "importlib.import_module",
        lambda name: SimpleNamespace(run_pdb2pqr=legacy_runner),
    )

    _run_pdb2pqr(["input.pdb", "output.pdb"])

    assert calls == [["input.pdb", "output.pdb"]]


def test_run_pdb2pqr_uses_parser_and_main_driver_when_legacy_entrypoint_is_absent(
    monkeypatch,
):
    parser_calls: list[list[str]] = []
    driver_calls: list[object] = []

    class _Parser:
        def parse_args(self, args):
            parser_calls.append(args)
            return SimpleNamespace(parsed_args=args)

    def parser_factory():
        return _Parser()

    def main_driver(namespace):
        driver_calls.append(namespace)

    monkeypatch.setattr(
        "importlib.import_module",
        lambda name: SimpleNamespace(
            build_main_parser=parser_factory,
            main_driver=main_driver,
        ),
    )

    _run_pdb2pqr(["input.pdb", "output.pdb", "--ff", "PARSE"])

    assert parser_calls == [["input.pdb", "output.pdb", "--ff", "PARSE"]]
    assert len(driver_calls) == 1
    assert driver_calls[0].parsed_args == ["input.pdb", "output.pdb", "--ff", "PARSE"]


def test_add_ifp_prolif_publishes_compact_and_event_sidecars(temp_dir, monkeypatch):
    structures = temp_dir / "results" / "structures"
    structures.mkdir(parents=True)
    cif_path = structures / "prediction.cif"
    cif_path.write_text("placeholder", encoding="utf-8")
    key = InteractionKey(ResidueIdentity("A", 168, residue_name="ASP"), "hb_acceptor")
    fingerprint = InteractionFingerprint(
        taxonomy=IFPTaxonomy.PROLIF,
        ligand=LigandIdentity("B", 1, residue_name="LIG"),
        receptor_chains=("A",),
        interactions=frozenset({key}),
        source_path=cif_path,
        events=(
            InteractionEvent(
                interaction=key,
                ligand_atoms=(AtomIdentity("B", 1, "", "LIG", "O1", "O"),),
                protein_atoms=(AtomIdentity("A", 168, "", "ASP", "N", "N"),),
                ligand_role="acceptor",
                protein_role="donor",
                geometry=(InteractionGeometry("distance", 2.9, "angstrom"),),
            ),
        ),
    )
    monkeypatch.setattr(
        "cofolder.modules.analytics.structure.extract_interaction_fingerprint",
        lambda *args, **kwargs: fingerprint,
    )
    chain_df = pd.DataFrame(
        [
            {"CHAIN_ID": "A", "ENTITY_TYPE": "protein", "cif_file": "prediction.cif"},
            {"CHAIN_ID": "B", "ENTITY_TYPE": "ligand", "cif_file": "prediction.cif"},
        ]
    )

    result = Structure(temp_dir, chain_df, structures).add_ifp_prolif()

    compact_name = result.loc[1, "ifp_prolif"]
    event_name = result.loc[1, "ifp_prolif_events"]
    output_dir = temp_dir / "results" / "ifp" / "prolif"
    assert json.loads((output_dir / compact_name).read_text()) == ["A:168:hb_acceptor"]
    events = json.loads((temp_dir / "results" / event_name).read_text())
    assert events[0]["protein_atoms"][0]["atom_name"] == "N"
    assert Path(event_name).parts[:2] == ("ifp", "prolif")
