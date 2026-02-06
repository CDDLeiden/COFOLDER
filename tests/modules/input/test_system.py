"""Tests for cofolder.modules.input.system module."""
import pytest
import yaml

from cofolder.modules.input.system import System


class TestSystemInit:
    """Tests for System initialization."""

    def test_init_with_dict(self):
        """Test initialization with a dictionary."""
        system_dict = {"sequences": [{"protein": {"id": "A"}}]}
        sys = System(system=system_dict)

        assert sys.system == system_dict

    def test_init_with_yaml_path(self, sample_system_yaml):
        """Test initialization with YAML file path."""
        sys = System(system_path=str(sample_system_yaml))

        assert sys.system is not None
        assert "sequences" in sys.system

    def test_init_with_both_raises_error(self, sample_system_yaml):
        """Test that providing both system and system_path raises error."""
        with pytest.raises(ValueError, match="Provide either"):
            System(system={"test": "data"}, system_path=str(sample_system_yaml))

    def test_init_with_neither_raises_error(self):
        """Test that providing neither raises error."""
        with pytest.raises(ValueError, match="Either 'system' or 'system_path'"):
            System()


class TestSystemUpdateSystem:
    """Tests for System.update_system method."""

    def test_update_with_path(self):
        """Test updating a nested value using path."""
        system_dict = {
            "sequences": [
                {"protein": {"id": "A", "fasta": "MKRAAT"}}
            ]
        }
        sys = System(system=system_dict)

        sys.update_system(value="NEWSEQ", path=["sequences", 0, "protein", "fasta"])

        assert sys.system["sequences"][0]["protein"]["fasta"] == "NEWSEQ"

    def test_update_creates_nested_structure(self):
        """Test that update creates nested structure if needed."""
        # Start with a minimal valid system
        sys = System(system={"existing": "data"})

        sys.update_system(value="test", path=["level1", "level2", "key"])

        assert sys.system["level1"]["level2"]["key"] == "test"
        assert sys.system["existing"] == "data"  # Original data preserved

    def test_update_list_index(self):
        """Test updating value at list index."""
        sys = System(system={"items": [1, 2, 3]})

        sys.update_system(value=99, path=["items", 1])

        assert sys.system["items"][1] == 99

    def test_update_with_parent_key_and_sub_key(self):
        """Test updating with parent_key and sub_key."""
        system_dict = {
            "sequences": [
                {"protein": {"id": "A", "msa": None}}
            ]
        }
        sys = System(system=system_dict)

        sys.update_system(value="/path/to/msa", parent_key="protein", sub_key="msa")

        # This should find the protein and update its msa
        assert sys.system["sequences"][0]["protein"]["msa"] == "/path/to/msa"


class TestSystemFindValue:
    """Tests for System.find_value method."""

    def test_find_by_key(self):
        """Test finding value by key name."""
        system_dict = {
            "sequences": [
                {"protein": {"id": "A", "msa": "/path/to/msa"}}
            ]
        }
        sys = System(system=system_dict)

        result = sys.find_value(key="msa")

        assert result == "/path/to/msa"

    def test_find_by_path(self):
        """Test finding value by path."""
        system_dict = {
            "sequences": [
                {"protein": {"id": "A", "fasta": "MKRAAT"}}
            ]
        }
        sys = System(system=system_dict)

        result = sys.find_value(path=["sequences", 0, "protein", "fasta"])

        assert result == "MKRAAT"

    def test_find_nonexistent_key_returns_none(self):
        """Test that finding nonexistent key returns None."""
        sys = System(system={"key": "value"})

        # The actual implementation returns None for nonexistent keys
        result = sys.find_value(key="nonexistent")
        assert result is None

    def test_find_without_key_or_path_returns_none(self):
        """Test that calling find_value without key or path returns None."""
        sys = System(system={"key": "value"})

        # The actual implementation returns None when neither key nor path provided
        result = sys.find_value()
        assert result is None


class TestSystemSaveSystemToYaml:
    """Tests for System.save_system_to_yaml method."""

    def test_save_to_yaml(self, temp_dir):
        """Test saving system to YAML file."""
        system_dict = {"sequences": [{"protein": {"id": "A"}}]}
        sys = System(system=system_dict)

        output_path = temp_dir / "output.yaml"
        sys.save_system_to_yaml(path=str(output_path))

        assert output_path.exists()

        # Verify content
        with open(output_path) as f:
            loaded_data = yaml.safe_load(f)

        assert loaded_data == system_dict
