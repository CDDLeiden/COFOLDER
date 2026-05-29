"""Tests for structure analytics compatibility helpers."""

from __future__ import annotations

from types import SimpleNamespace

from cofolder.modules.analytics.structure import _run_pdb2pqr


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


def test_run_pdb2pqr_uses_parser_and_main_driver_when_legacy_entrypoint_is_absent(monkeypatch):
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
