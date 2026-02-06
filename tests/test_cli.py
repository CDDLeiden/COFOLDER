"""Tests for cofolder.cli module."""
from unittest.mock import patch

import pytest

from cofolder import cli


class TestCLIMain:
    """Tests for CLI main function."""

    @patch('cofolder.recipes.predict.Predict.run')
    @patch('cofolder.modules.utils.helpers.create_dir')
    def test_predict_command(self, mock_create_dir, mock_run, sample_system_yaml, sample_options_yaml):
        """Test predict command execution."""
        args = [
            "predict",
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
        # Help should exit with code 0
        assert exc_info.value.code == 0

    def test_version_flag(self):
        """Test version flag."""
        with pytest.raises(SystemExit) as exc_info:
            cli.main(["-v"])
        assert exc_info.value.code == 0


class TestPredictRecipe:
    """Tests for PredictRecipe class."""

    def test_add_arguments(self):
        """Test that add_arguments adds correct arguments."""
        import argparse
        parser = argparse.ArgumentParser()
        subparser = parser.add_subparsers()
        predict_parser = subparser.add_parser("predict")

        cli.PredictRecipe.add_arguments(predict_parser)

        # Parse test arguments
        args = predict_parser.parse_args([
            "-s", "system.yaml",
            "-b", "options.yaml"
        ])

        assert args.system_path == "system.yaml"
        assert args.options_path == "options.yaml"


class TestScreenRecipe:
    """Tests for ScreenRecipe class."""

    def test_add_arguments(self):
        """Test that add_arguments adds correct arguments."""
        import argparse
        parser = argparse.ArgumentParser()
        subparser = parser.add_subparsers()
        screen_parser = subparser.add_parser("screen")

        cli.ScreenRecipe.add_arguments(screen_parser)

        # Parse test arguments
        args = screen_parser.parse_args([
            "-s", "system.yaml",
            "-b", "options.yaml",
            "-v", "sequences,0,ligand,smiles",
            "-c", "compounds.csv",
            "--col_variable", "smiles",
            "--col_id", "id"
        ])

        assert args.system_path == "system.yaml"
        assert args.variable == "sequences,0,ligand,smiles"


class TestEvaluateRecipe:
    """Tests for EvaluateRecipe class."""

    def test_add_arguments(self):
        """Test that add_arguments adds correct arguments."""
        import argparse
        parser = argparse.ArgumentParser()
        subparser = parser.add_subparsers()
        eval_parser = subparser.add_parser("evaluate")

        cli.EvaluateRecipe.add_arguments(eval_parser)

        # Parse test arguments
        args = eval_parser.parse_args([
            "-s", "system.yaml",
            "-b", "options.yaml",
            "--repeats", "5"
        ])

        assert args.system_path == "system.yaml"
        assert args.repeats == 5


class TestOracleRecipe:
    """Tests for OracleRecipe class."""

    def test_add_arguments(self):
        """Test that add_arguments adds correct arguments."""
        import argparse
        parser = argparse.ArgumentParser()
        subparser = parser.add_subparsers()
        oracle_parser = subparser.add_parser("oracle")

        cli.OracleRecipe.add_arguments(oracle_parser)

        # Parse test arguments
        args = oracle_parser.parse_args([
            "-s", "system.yaml",
            "-b", "options.yaml"
        ])

        assert args.system_path == "system.yaml"
        assert args.options_path == "options.yaml"
