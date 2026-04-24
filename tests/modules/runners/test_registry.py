"""Tests for runner discovery and lookup."""

from types import SimpleNamespace

import pytest

from cofolder.modules.runners import discover_runners, get_runner, list_runner_names


class _FakeRunner:
    def __init__(self, name):
        self.name = name


def test_discover_runners_finds_runner_modules(monkeypatch):
    monkeypatch.setattr(
        "cofolder.modules.runners.pkgutil.iter_modules",
        lambda paths: [
            SimpleNamespace(name="alpha_runner"),
            SimpleNamespace(name="base"),
            SimpleNamespace(name="_private"),
            SimpleNamespace(name="not_a_runner"),
        ],
    )

    def _fake_import(name):
        if name.endswith("alpha_runner"):
            return SimpleNamespace(RUNNER=_FakeRunner("alpha"))
        return SimpleNamespace()

    monkeypatch.setattr("cofolder.modules.runners.importlib.import_module", _fake_import)

    runners = discover_runners()

    assert list(runners) == ["alpha"]


def test_list_runner_names_returns_sorted_names(monkeypatch):
    monkeypatch.setattr(
        "cofolder.modules.runners.discover_runners",
        lambda: {"boltz2": _FakeRunner("boltz2"), "alpha": _FakeRunner("alpha")},
    )

    assert list_runner_names() == ["alpha", "boltz2"]


def test_get_runner_raises_with_available_names(monkeypatch):
    monkeypatch.setattr(
        "cofolder.modules.runners.discover_runners",
        lambda: {"boltz2": _FakeRunner("boltz2")},
    )

    with pytest.raises(ValueError, match="Available runners: boltz2"):
        get_runner("missing")


def test_discover_real_runners_include_boltz1_boltz2_and_boltz_community():
    runners = discover_runners()

    assert "boltz1" in runners
    assert "boltz2" in runners
    assert "boltz-community" in runners
