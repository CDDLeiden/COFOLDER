"""Tests for runner discovery and lookup."""

from types import SimpleNamespace

import pytest

import cofolder.modules.runners as runner_exports
from cofolder.modules.runners import discover_runners, get_runner, list_runner_names
from cofolder.modules.runners.contracts import (
    RunnerExecutionRequest,
    RunnerExecutionResult,
    RunnerPreparationResult,
)


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


def test_runner_module_keeps_legacy_aliases_as_compatibility_exports():
    assert runner_exports.RunnerExecutionRequest is RunnerExecutionRequest
    assert runner_exports.RunnerExecutionResult is RunnerExecutionResult
    assert runner_exports.RunnerPreparationResult is RunnerPreparationResult
    assert runner_exports.RunnerRequest is RunnerExecutionRequest
    assert runner_exports.RunnerResult is RunnerExecutionResult
    assert runner_exports.RunnerPreparation is RunnerPreparationResult
