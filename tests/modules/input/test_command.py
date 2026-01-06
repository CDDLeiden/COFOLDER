"""Tests for boltz_lab.modules.input.command module."""
import pytest

from boltz_lab.modules.input.command import Command
from boltz_lab.modules.input.system import System


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
        assert "options" in cmd.options

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
