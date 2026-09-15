"""Tests for cofolder.modules.input.command module."""
from pathlib import Path

import pytest

from cofolder.modules.input import command
from cofolder.modules.input.command import Command
from cofolder.modules.input.system import System


def test_download_cache_uses_a_unique_cleaned_workspace(monkeypatch):
    workspaces = []

    def fake_run_boltz(argv):
        fasta_path = Path(argv[2])
        output_path = Path(argv[argv.index("--out_dir") + 1])
        assert fasta_path.read_text(encoding="utf-8") == ">A|protein|\nA\n"
        assert output_path == fasta_path.parent
        workspaces.append(fasta_path.parent)

    monkeypatch.setattr(
        "cofolder.modules.runners.boltz_runner.run_boltz", fake_run_boltz
    )

    command.download_cache("/cache with spaces")
    command.download_cache("/cache with spaces")

    assert len(set(workspaces)) == 2
    assert all(not workspace.exists() for workspace in workspaces)


def test_download_cache_cleans_workspace_after_failure(monkeypatch):
    workspaces = []

    def fail(argv):
        workspaces.append(Path(argv[2]).parent)
        raise RuntimeError("setup failed")

    monkeypatch.setattr("cofolder.modules.runners.boltz_runner.run_boltz", fail)

    with pytest.raises(RuntimeError, match="setup failed"):
        command.download_cache("/cache")

    assert len(workspaces) == 1
    assert not workspaces[0].exists()


class TestCommandInit:
    """Tests for Command initialization."""

    def test_init_with_dict(self):
        """Test initialization with a dictionary."""
        options_dict = {"options": [{"cache": "~/.boltz"}]}
        cmd = Command(options=options_dict)

        assert cmd.options == options_dict

    def test_init_with_yaml_path(self, sample_options_yaml):
        """Test initialization with YAML file path."""
        cmd = Command(options_path=str(sample_options_yaml))

        assert cmd.options is not None
        assert cmd.options["version"] == 1
        assert cmd.options["runner"]["recycling_steps"] == 3

    def test_init_with_both_raises_error(self, sample_options_yaml):
        """Test that providing both options and options_path raises error."""
        with pytest.raises(ValueError, match="Provide either"):
            Command(options={"test": "data"}, options_path=str(sample_options_yaml))

    def test_init_with_neither_raises_error(self):
        """Test that providing neither raises error."""
        with pytest.raises(ValueError, match="Either 'options' or 'options_path'"):
            Command()


class TestCommandUpdateOptions:
    """Tests for Command.update_options method."""

    def test_update_with_path(self):
        """Test updating options with path."""
        options_dict = {"options": [{"recycling_steps": 3}]}
        cmd = Command(options=options_dict)

        cmd.update_options(value=5, path=["options", 0, "recycling_steps"])

        assert cmd.options["options"][0]["recycling_steps"] == 5


class TestCommandFindValue:
    """Tests for Command.find_value method."""

    def test_find_by_key(self):
        """Test finding value by key."""
        options_dict = {"options": [{"cache": "~/.boltz"}]}
        cmd = Command(options=options_dict)

        result = cmd.find_value(key="cache")

        assert result == "~/.boltz"

    def test_find_by_path(self):
        """Test finding value by path."""
        options_dict = {"options": [{"recycling_steps": 3}]}
        cmd = Command(options=options_dict)

        result = cmd.find_value(path=["options", 0, "recycling_steps"])

        assert result == 3


class TestCommandSetCommand:
    """Tests for Command.set_command method."""

    def test_set_command_basic(self, sample_system_yaml, sample_options_yaml):
        """Test setting command with basic options."""
        cmd = Command(options_path=str(sample_options_yaml))
        sys = System(system_path=str(sample_system_yaml))

        cmd.out_dir = "/tmp/output"
        cmd.system_path = str(sample_system_yaml)

        command_list = cmd.set_command(system=sys)

        assert isinstance(command_list, list)
        assert "boltz" in command_list
        assert "predict" in command_list
        assert str(sample_system_yaml) in command_list

    def test_set_command_with_devices(self):
        """Test set_command with devices specified."""
        options_dict = {
            "options": [
                {"cache": "~/.boltz"},
                {"devices": [0, 1]}
            ]
        }
        cmd = Command(options=options_dict)
        sys = System(system={"sequences": []})

        cmd.out_dir = "/tmp/output"
        cmd.system_path = "/tmp/system.yaml"

        command_list = cmd.set_command(system=sys)

        assert "--devices" in command_list
        # Find index of --devices and check the next value
        idx = command_list.index("--devices")
        # The actual implementation converts list to string representation
        assert "[0, 1]" in command_list[idx + 1] or "0,1" in command_list[idx + 1]

    def test_set_command_does_not_treat_numeric_zero_or_one_as_booleans(self):
        cmd = Command(
            options={
                "options": [
                    {"recycling_steps": 0},
                    {"devices": 1},
                    {"write_full_pae": True},
                    {"override": False},
                ]
            }
        )
        sys = System(system={"sequences": []})
        cmd.out_dir = "/tmp/output"
        cmd.system_path = "/tmp/system.yaml"

        assert cmd.set_command(system=sys) == [
            "boltz",
            "predict",
            "/tmp/system.yaml",
            "--out_dir",
            "/tmp/output",
            "--seed",
            "0",
            "--recycling_steps",
            "0",
            "--devices",
            "1",
            "--write_full_pae",
        ]

    def test_set_command_requests_only_missing_protein_msas(self):
        cmd = Command(options={"options": []})
        sys = System(
            system={
                "sequences": [
                    {"protein": {"id": "A", "sequence": "AAAA", "msa": "/tmp/a.a3m"}},
                    {"protein": {"id": "B", "sequence": "BBBB"}},
                ]
            }
        )
        cmd.out_dir = "/tmp/output"
        cmd.system_path = "/tmp/system.yaml"

        assert "--use_msa_server" in cmd.set_command(system=sys)

        sys.system["sequences"][1]["protein"]["msa"] = "/tmp/b.a3m"
        assert "--use_msa_server" not in cmd.set_command(system=sys)
