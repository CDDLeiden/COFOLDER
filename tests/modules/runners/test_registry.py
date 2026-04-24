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
        lambda: {"boltz": _FakeRunner("boltz"), "alpha": _FakeRunner("alpha")},
    )

    assert list_runner_names() == ["alpha", "boltz"]


def test_get_runner_raises_with_available_names(monkeypatch):
    monkeypatch.setattr(
        "cofolder.modules.runners.discover_runners",
        lambda: {"boltz": _FakeRunner("boltz")},
    )

    with pytest.raises(ValueError, match="Available runners: boltz"):
        get_runner("missing")
