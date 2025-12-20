"""Tests for boltz_lab.recipes.screen module."""
from unittest.mock import patch


from boltz_lab.recipes.screen import Screen


class TestScreenInit:
    """Tests for Screen initialization."""

    def test_init_basic(self, sample_system_yaml, sample_options_yaml, temp_dir):
        """Test basic initialization."""
        screener = Screen(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            variable="sequences,1,ligand,smiles"
        )

        assert screener.wrk_dir == str(temp_dir)
        assert screener.variable == ["sequences", 1, "ligand", "smiles"]

    def test_parse_list_integers(self, sample_system_yaml, sample_options_yaml, temp_dir):
        """Test parsing list with integers."""
        screener = Screen(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            variable="sequences,0,ligand,smiles"
        )

        assert screener.variable == ["sequences", 0, "ligand", "smiles"]

    def test_parse_merge_data(self, sample_system_yaml, sample_options_yaml, temp_dir):
        """Test parsing merge_data parameter."""
        screener = Screen(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            variable="sequences,0,ligand,smiles",
            merge_data="mw,logp,tpsa"
        )

        assert screener.merge_data == ["mw", "logp", "tpsa"]


class TestScreenParseList:
    """Tests for Screen._parse_list static method."""

    def test_parse_mixed_types(self):
        """Test parsing mixed strings and integers."""
        result = Screen._parse_list("sequences,0,ligand,smiles")
        assert result == ["sequences", 0, "ligand", "smiles"]

    def test_parse_none(self):
        """Test parsing None."""
        result = Screen._parse_list(None)
        assert result == []

    def test_parse_empty_string(self):
        """Test parsing empty string."""
        result = Screen._parse_list("")
        assert result == []


class TestScreenRun:
    """Tests for Screen.run method."""

    @patch.object(Screen, 'load_screen')
    @patch.object(Screen, 'iterate')
    def test_run_calls_methods(self, mock_iterate, mock_load_screen,
                               sample_system_yaml, sample_options_yaml, temp_dir):
        """Test that run calls load_screen and iterate."""
        screener = Screen(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            variable="sequences,0,ligand,smiles"
        )

        screener.run()

        assert mock_load_screen.called
        assert mock_iterate.called
