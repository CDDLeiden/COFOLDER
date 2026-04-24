"""Tests for the boltz1 runner."""

from importlib.metadata import PackageNotFoundError
from unittest.mock import patch

import yaml

from cofolder.modules.runners.boltz1_runner import Boltz1Runner


def test_boltz1_runner_requires_exact_100_package_line():
    runner = Boltz1Runner()

    def _fake_version(name):
        if name == "boltz":
            return "2.2.1"
        raise PackageNotFoundError

    with patch(
        "cofolder.modules.runners.base.metadata.version",
        side_effect=_fake_version,
    ):
        available, message = runner.check_availability()

    assert available is False
    assert "requires boltz==1.0.0" in message
    assert "cofolder[boltz1]" in message


def test_boltz1_runner_accepts_exact_100_package_line():
    runner = Boltz1Runner()

    def _fake_version(name):
        if name == "boltz":
            return "1.0.0"
        raise PackageNotFoundError

    with patch(
        "cofolder.modules.runners.base.metadata.version",
        side_effect=_fake_version,
    ):
        available, message = runner.check_availability()

    assert available is True
    assert message is None
    assert runner.capabilities == {"confidence_metrics"}


def test_boltz1_runner_load_options_removes_model_flag(temp_dir):
    options_path = temp_dir / "options.yaml"
    options_path.write_text(
        yaml.safe_dump({"options": [{"cache": "~/.boltz"}, {"model": "boltz2"}]}),
        encoding="utf-8",
    )

    command = Boltz1Runner().load_options(options_path)

    assert command.find_value(key="model") is None
