"""Tests for cofolder.cli module."""
from unittest.mock import patch

import pytest

from cofolder import cli


class TestCLIMain:
    """Tests for CLI main function."""

    @patch('cofolder.recipes.validate.Validate.run')
    @patch('cofolder.modules.utils.helpers.create_dir')
    def test_validate_command(self, mock_create_dir, mock_run, sample_system_yaml, sample_options_yaml):
        """Test validate command execution."""
        args = [
            "validate",
            "-s", str(sample_system_yaml),
            "-b", str(sample_options_yaml),
            "-w", "/tmp/test"
        ]

        cli.main(args)

        assert mock_run.called

    def test_main_no_args(self):
        """Test main with no arguments."""
        with pytest.raises(SystemExit):
            cli.main([])

    def test_main_help(self):
        """Test main with help flag."""
        with pytest.raises(SystemExit) as exc_info:
            cli.main(["-h"])
        assert exc_info.value.code == 0

    def test_version_flag(self):
        """Test version flag."""
        with pytest.raises(SystemExit) as exc_info:
            cli.main(["-v"])
        assert exc_info.value.code == 0


class TestValidateRecipe:
    """Tests for ValidateRecipe class."""

    def test_add_arguments(self):
        """Test that add_arguments adds correct arguments."""
        import argparse

        parser = argparse.ArgumentParser()
        subparser = parser.add_subparsers()
        validate_parser = subparser.add_parser("validate")

        cli.ValidateRecipe.add_arguments(validate_parser)

        args = validate_parser.parse_args([
            "-s", "system.yaml",
            "-b", "options.yaml",
            "--reference_path", "reference.pdb",
            "--pocket_coverage_reference", "A2 S8 T10",
            "--reproduction_metrics", "sucos",
        ])

        assert args.system_path == "system.yaml"
        assert args.options_path == "options.yaml"
        assert args.reference_path == "reference.pdb"
        assert args.pocket_coverage_reference == "A2 S8 T10"
        assert args.reproduction_metrics == ["sucos"]


class TestScreenRecipe:
    """Tests for ScreenRecipe class."""

    def test_add_arguments(self):
        """Test that add_arguments adds correct arguments."""
        import argparse

        parser = argparse.ArgumentParser()
        subparser = parser.add_subparsers()
        screen_parser = subparser.add_parser("screen")

        cli.ScreenRecipe.add_arguments(screen_parser)

        args = screen_parser.parse_args([
            "-s", "system.yaml",
            "-b", "options.yaml",
            "-c", "compounds.csv",
            "--col_id", "id",
            "-v", "sequences,0,ligand,smiles",
            "--col_variable", "smiles",
            "-v", "sequences,0,ligand,ccd",
            "--col_variable", "ccd",
        ])

        assert args.system_path == "system.yaml"
        assert args.variable == ["sequences,0,ligand,smiles", "sequences,0,ligand,ccd"]
        assert args.col_variable == ["smiles", "ccd"]

    def test_main_raises_on_mapping_count_mismatch(self, sample_system_yaml, sample_options_yaml, temp_dir):
        """Test that screen main rejects mismatched --variable/--col_variable counts."""
        args = [
            "screen",
            "-s", str(sample_system_yaml),
            "-b", str(sample_options_yaml),
            "-c", "compounds.csv",
            "--col_id", "id",
            "-v", "sequences,0,ligand,smiles",
            "-v", "sequences,0,ligand,ccd",
            "--col_variable", "smiles",
            "-w", str(temp_dir),
        ]

        with pytest.raises(ValueError, match="Number of --variable entries must match"):
            cli.main(args)


class TestOracleRecipe:
    """Tests for OracleRecipe class."""

    def test_add_arguments(self):
        """Test that add_arguments adds correct arguments."""
        import argparse

        parser = argparse.ArgumentParser()
        subparser = parser.add_subparsers()
        oracle_parser = subparser.add_parser("oracle")

        cli.OracleRecipe.add_arguments(oracle_parser)

        args = oracle_parser.parse_args([
            "-s", "system.yaml",
            "-b", "options.yaml",
        ])

        assert args.system_path == "system.yaml"
        assert args.options_path == "options.yaml"
