"""Tests for cofolder.cli module."""
import logging
from unittest.mock import MagicMock
from unittest.mock import Mock
from unittest.mock import patch

import pytest

from cofolder import cli


class TestCLIMain:
    """Tests for CLI main function."""

    @patch('cofolder.recipes.validate.Validate.run')
    @patch('cofolder.modules.utils.helpers.create_dir')
    def test_validate_command(self, mock_create_dir, mock_run, sample_system_yaml, sample_options_yaml, temp_dir):
        """Test validate command execution."""
        args = [
            "validate",
            "-s", str(sample_system_yaml),
            "-o", str(sample_options_yaml),
            "-w", str(temp_dir)
        ]

        cli.main(args)

        assert mock_run.called

    @patch('cofolder.recipes.validate.Validate.run')
    def test_validate_logs_cli_command_first(self, mock_run, sample_system_yaml, sample_options_yaml, temp_dir, caplog):
        """Test that the CLI command is logged before the logger-init message."""
        args = [
            "validate",
            "-s", str(sample_system_yaml),
            "-o", str(sample_options_yaml),
            "-w", str(temp_dir),
        ]

        with caplog.at_level(logging.INFO):
            cli.main(args)

        text = caplog.text
        assert "CLI command: cofolder validate" in text
        assert "Logger initialized." in text
        assert text.index("CLI command: cofolder validate") < text.index("Logger initialized.")

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

    def test_unavailable_runner_shows_install_hint(self, sample_system_yaml, sample_options_yaml, temp_dir):
        args = [
            "validate",
            "-s", str(sample_system_yaml),
            "-o", str(sample_options_yaml),
            "-w", str(temp_dir),
        ]
        fake_runner = MagicMock()
        fake_runner.ensure_available.side_effect = RuntimeError("install boltz")
        with patch("cofolder.cli.get_runner", return_value=fake_runner):
            with pytest.raises(RuntimeError, match="install boltz"):
                cli.main(args)

    def test_unavailable_runner_is_checked_before_input_path_validation(self, temp_dir):
        args = [
            "validate",
            "-s", str(temp_dir / "missing_system.yaml"),
            "-o", str(temp_dir / "missing_options.yaml"),
            "-w", str(temp_dir),
        ]
        fake_runner = MagicMock()
        fake_runner.ensure_available.side_effect = RuntimeError("install boltz first")

        with patch("cofolder.cli.get_runner", return_value=fake_runner):
            with pytest.raises(RuntimeError, match="install boltz first"):
                cli.main(args)


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
            "-o", "options.yaml",
            "--runner", "boltz2",
            "--reference_path", "reference.pdb",
            "--pocket_coverage_reference", "A2 S8 T10",
            "--reproduction_metrics", "sucos",
        ])

        assert args.system_path == "system.yaml"
        assert args.options_path == "options.yaml"
        assert args.runner == "boltz2"
        assert args.reference_path == "reference.pdb"
        assert args.pocket_coverage_reference == "A2 S8 T10"
        assert args.reproduction_metrics == ["sucos"]

    def test_runner_help_lists_available_runners(self):
        import argparse

        parser = argparse.ArgumentParser()
        subparser = parser.add_subparsers()
        validate_parser = subparser.add_parser("validate")

        cli.ValidateRecipe.add_arguments(validate_parser)

        help_text = validate_parser.format_help()

        assert "--runner" in help_text
        assert "Available runners:" in help_text
        assert "boltz1" in help_text
        assert "boltz2" in help_text
        assert "boltz-community" in help_text

    def test_runner_defaults_to_boltz2(self):
        import argparse

        parser = argparse.ArgumentParser()
        subparser = parser.add_subparsers()
        validate_parser = subparser.add_parser("validate")

        cli.ValidateRecipe.add_arguments(validate_parser)

        args = validate_parser.parse_args([
            "-s", "system.yaml",
            "-o", "options.yaml",
        ])

        assert args.runner == "boltz2"

    def test_main_uses_lazy_recipe_loader(self, sample_system_yaml, sample_options_yaml, temp_dir):
        validator = Mock()
        validate_cls = Mock(return_value=validator)

        args = [
            "validate",
            "-s", str(sample_system_yaml),
            "-o", str(sample_options_yaml),
            "-w", str(temp_dir),
        ]

        with patch("cofolder.cli._load_recipe_class", return_value=validate_cls):
            cli.main(args)

        validate_cls.assert_called_once()
        validator.run.assert_called_once()


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
            "-o", "options.yaml",
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
            "-o", str(sample_options_yaml),
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
            "-o", "options.yaml",
            "--input_smiles", "CCO",
            "--output_metric", "affinity_pred_value",
            "--aggregate", "first",
        ])

        assert args.system_path == "system.yaml"
        assert args.options_path == "options.yaml"
        assert args.input_smiles == "CCO"
        assert args.output_metric == "affinity_pred_value"
        assert args.aggregate == "first"
