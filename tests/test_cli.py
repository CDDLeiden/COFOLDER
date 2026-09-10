"""Tests for cofolder.cli module."""
import logging
import runpy
from unittest.mock import MagicMock, Mock, patch

import pytest

from cofolder import cli


def test_module_entrypoint_propagates_main_exit_code():
    with patch("cofolder.cli.main", return_value=2):
        with pytest.raises(SystemExit) as caught:
            runpy.run_module("cofolder.__main__", run_name="__main__")

    assert caught.value.code == 2


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

    def test_version_flag(self, capsys):
        """Test version flag."""
        with pytest.raises(SystemExit) as exc_info:
            cli.main(["-v"])
        assert exc_info.value.code == 0
        assert capsys.readouterr().out == "cofolder 1.0.0\n"

    def test_unavailable_runner_shows_install_hint(self, sample_system_yaml, sample_options_yaml, temp_dir):
        from cofolder.modules.runners.boltz2_runner import Boltz2Runner

        args = [
            "validate",
            "-s", str(sample_system_yaml),
            "-o", str(sample_options_yaml),
            "-w", str(temp_dir),
        ]
        fake_runner = Boltz2Runner()
        fake_runner.ensure_available = MagicMock()
        fake_runner.ensure_available.side_effect = RuntimeError("install boltz")
        with patch("cofolder.recipes.validate.get_runner", return_value=fake_runner):
            assert cli.main(args) == 1

    def test_input_paths_are_checked_before_runner_availability(self, temp_dir):
        args = [
            "validate",
            "-s", str(temp_dir / "missing_system.yaml"),
            "-o", str(temp_dir / "missing_options.yaml"),
            "-w", str(temp_dir),
        ]
        assert cli.main(args) == 2

    def test_bias_command_does_not_require_options_or_runner(
        self,
        sample_system_yaml,
        temp_dir,
    ):
        protein_ref = temp_dir / "protein_training.csv"
        ligand_ref = temp_dir / "ligand_training.csv"
        protein_ref.write_text(
            "pdb_id,release_date,sequence\n"
            "1ABC,2022-01-01,MKRAAT\n",
            encoding="utf-8",
        )
        ligand_ref.write_text(
            "pdb_id,release_date,ligand_id,smiles\n"
            "1ABC,2022-01-01,ETH,CCO\n",
            encoding="utf-8",
        )

        bias_runner = Mock()
        bias_cls = Mock(return_value=bias_runner)

        args = [
            "bias",
            "-s", str(sample_system_yaml),
            "--protein_training_data_path", str(protein_ref),
            "--ligand_training_data_path", str(ligand_ref),
            "-w", str(temp_dir),
        ]

        with patch("cofolder.cli._load_recipe_class", return_value=bias_cls):
            cli.main(args)

        bias_cls.assert_called_once()
        bias_runner.run.assert_called_once()

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
        assert args.scoring_functions is None

    def test_scoring_functions_preserve_explicit_selection(self):
        import argparse

        parser = argparse.ArgumentParser()
        subparser = parser.add_subparsers()
        validate_parser = subparser.add_parser("validate")
        cli.ValidateRecipe.add_arguments(validate_parser)

        args = validate_parser.parse_args([
            "-s", "system.yaml",
            "-o", "options.yaml",
            "--scoring_functions", "confidence_metrics", "affinity_metrics",
        ])

        assert args.scoring_functions == ["confidence_metrics", "affinity_metrics"]

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
            "--ligand_chain", "B",
            "--smiles_column", "smiles",
            "--ifp_filter_threshold", "0.75",
            "--ifp_filter_source", "reference_complex",
            "--ifp_taxonomy", "prolif",
            "--ifp_similarity_metric", "jaccard",
            "--ifp_reference_ligand", "L:401",
            "--ifp_reference_receptor_chain", "A",
            "--cluster_ifps",
            "--ifp_cluster_similarity_threshold", "0.8",
        ])

        assert args.system_path == "system.yaml"
        assert args.ligand_chain == "B"
        assert args.smiles_column == "smiles"
        assert args.ifp_filter_threshold == 0.75
        assert args.ifp_filter_source == "reference_complex"
        assert args.ifp_taxonomy == "prolif"
        assert args.ifp_similarity_metric == "jaccard"
        assert args.ifp_reference_ligand == "L:401"
        assert args.ifp_reference_receptor_chains == ["A"]
        assert args.cluster_ifps is True
        assert args.ifp_cluster_similarity_threshold == 0.8
        assert args.scoring_functions is None

    def test_main_rejects_legacy_mutation_path_flags(self, sample_system_yaml, sample_options_yaml, temp_dir):
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

        with pytest.raises(SystemExit) as caught:
            cli.main(args)
        assert caught.value.code == 2


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


class TestBiasRecipe:
    """Tests for BiasRecipe class."""

    def test_add_arguments(self):
        import argparse

        parser = argparse.ArgumentParser()
        subparser = parser.add_subparsers()
        bias_parser = subparser.add_parser("bias")

        cli.BiasRecipe.add_arguments(bias_parser)

        args = bias_parser.parse_args([
            "-s", "system.yaml",
            "--protein_training_data_path", "protein_training.csv",
            "--ligand_training_data_path", "ligand_training.csv",
            "--custom_protein_reference_path", "custom_protein.csv",
            "--custom_ligand_reference_path", "custom_ligand.csv",
            "--bias_chains", "A", "B",
        ])

        assert args.system_path == "system.yaml"
        assert args.protein_training_data_path == "protein_training.csv"
        assert args.ligand_training_data_path == "ligand_training.csv"
        assert args.custom_protein_reference_path == "custom_protein.csv"
        assert args.custom_ligand_reference_path == "custom_ligand.csv"
        assert args.bias_chains == ["A", "B"]

    def test_main_uses_lazy_recipe_loader(self, sample_system_yaml, temp_dir):
        protein_ref = temp_dir / "protein_training.csv"
        ligand_ref = temp_dir / "ligand_training.csv"
        protein_ref.write_text(
            "pdb_id,release_date,sequence\n"
            "1ABC,2022-01-01,MKRAAT\n",
            encoding="utf-8",
        )
        ligand_ref.write_text(
            "pdb_id,release_date,ligand_id,smiles\n"
            "1ABC,2022-01-01,ETH,CCO\n",
            encoding="utf-8",
        )

        bias_runner = Mock()
        bias_cls = Mock(return_value=bias_runner)

        args = [
            "bias",
            "-s", str(sample_system_yaml),
            "--protein_training_data_path", str(protein_ref),
            "--ligand_training_data_path", str(ligand_ref),
            "-w", str(temp_dir),
        ]

        with patch("cofolder.cli._load_recipe_class", return_value=bias_cls):
            cli.main(args)

        bias_cls.assert_called_once()
        bias_runner.run.assert_called_once()

    def test_main_requires_ligand_training_data_without_build_flag(
        self,
        temp_dir,
    ):
        system_path = temp_dir / "system.yaml"
        system_path.write_text(
            "sequences:\n"
            "  - protein:\n"
            "      id: A\n"
            "      sequence: MKRAAT\n"
            "  - ligand:\n"
            "      id: B\n"
            "      smiles: CCO\n",
            encoding="utf-8",
        )
        protein_ref = temp_dir / "protein_training.csv"
        protein_ref.write_text(
            "pdb_id,release_date,sequence\n"
            "1ABC,2022-01-01,MKRAAT\n",
            encoding="utf-8",
        )

        args = [
            "bias",
            "-s", str(system_path),
            "--protein_training_data_path", str(protein_ref),
            "-w", str(temp_dir),
        ]

        assert cli.main(
            ["bias", "-s", str(system_path), "-w", str(temp_dir)]
        ) == 2

        assert cli.main(args) == 1

    def test_main_accepts_build_mode_with_fresh_output_paths(
        self,
        sample_system_yaml,
        temp_dir,
    ):
        components_cif = temp_dir / "components.cif"
        components_cif.write_text("data_components\n", encoding="utf-8")
        protein_output = temp_dir / "generated" / "protein_training.csv"

        bias_runner = Mock()
        bias_cls = Mock(return_value=bias_runner)

        args = [
            "bias",
            "-s", str(sample_system_yaml),
            "--build_bias_training_data",
            "--protein_training_data_path", str(protein_output),
            "--bias_training_components_cif", str(components_cif),
            "-w", str(temp_dir),
        ]

        with patch("cofolder.cli._load_recipe_class", return_value=bias_cls):
            cli.main(args)

        bias_cls.assert_called_once()
        assert bias_cls.call_args.kwargs["protein_training_data_path"] == str(protein_output)
        bias_runner.run.assert_called_once()

    def test_main_accepts_fresh_ligand_output_path_in_build_mode(
        self,
        sample_system_yaml,
        temp_dir,
    ):
        components_cif = temp_dir / "components.cif"
        components_cif.write_text("data_components\n", encoding="utf-8")
        protein_output = temp_dir / "generated" / "protein_training.csv"
        ligand_output = temp_dir / "generated" / "ligand_training.csv"
        bias_runner = Mock()
        bias_cls = Mock(return_value=bias_runner)

        args = [
            "bias",
            "-s", str(sample_system_yaml),
            "--build_bias_training_data",
            "--protein_training_data_path", str(protein_output),
            "--ligand_training_data_path", str(ligand_output),
            "--bias_training_components_cif", str(components_cif),
            "-w", str(temp_dir),
        ]

        with patch("cofolder.cli._load_recipe_class", return_value=bias_cls):
            cli.main(args)

        assert bias_cls.call_args.kwargs["ligand_training_data_path"] == str(ligand_output)
        bias_runner.run.assert_called_once()

    def test_main_rejects_non_numeric_bias_threshold_before_recipe_load(
        self,
        sample_system_yaml,
        temp_dir,
    ):
        protein_ref = temp_dir / "protein_training.csv"
        ligand_ref = temp_dir / "ligand_training.csv"
        protein_ref.write_text(
            "pdb_id,release_date,sequence\n"
            "1ABC,2022-01-01,MKRAAT\n",
            encoding="utf-8",
        )
        ligand_ref.write_text(
            "pdb_id,release_date,ligand_id,smiles\n"
            "1ABC,2022-01-01,ETH,CCO\n",
            encoding="utf-8",
        )

        args = [
            "bias",
            "-s", str(sample_system_yaml),
            "--protein_training_data_path", str(protein_ref),
            "--ligand_training_data_path", str(ligand_ref),
            "--bias_ligand_similarity_threshold", "not-a-number",
            "-w", str(temp_dir),
        ]

        with patch("cofolder.cli._load_recipe_class") as mock_loader:
            with pytest.raises(SystemExit) as exc_info:
                cli.main(args)

        assert exc_info.value.code == 2
        mock_loader.assert_not_called()

    @pytest.mark.parametrize("threshold", ["-0.01", "1.01"])
    def test_main_rejects_out_of_range_bias_threshold_before_recipe_load(
        self,
        sample_system_yaml,
        temp_dir,
        threshold,
    ):
        protein_ref = temp_dir / "protein_training.csv"
        ligand_ref = temp_dir / "ligand_training.csv"
        protein_ref.write_text(
            "pdb_id,release_date,sequence\n"
            "1ABC,2022-01-01,MKRAAT\n",
            encoding="utf-8",
        )
        ligand_ref.write_text(
            "pdb_id,release_date,ligand_id,smiles\n"
            "1ABC,2022-01-01,ETH,CCO\n",
            encoding="utf-8",
        )

        args = [
            "bias",
            "-s", str(sample_system_yaml),
            "--protein_training_data_path", str(protein_ref),
            "--ligand_training_data_path", str(ligand_ref),
            "--bias_ligand_similarity_threshold", threshold,
            "-w", str(temp_dir),
        ]

        with patch("cofolder.cli._load_recipe_class") as mock_loader:
            assert cli.main(args) == 2

        mock_loader.assert_not_called()

    def test_main_requires_protein_output_path_in_build_mode_before_recipe_load(
        self,
        sample_system_yaml,
        temp_dir,
    ):
        args = [
            "bias",
            "-s", str(sample_system_yaml),
            "--build_bias_training_data",
            "-w", str(temp_dir),
        ]

        with patch("cofolder.cli._load_recipe_class") as mock_loader:
            assert cli.main(args) == 2

        mock_loader.assert_not_called()

    def test_main_rejects_missing_components_cif_in_build_mode_before_recipe_load(
        self,
        sample_system_yaml,
        temp_dir,
    ):
        protein_output = temp_dir / "generated" / "protein_training.csv"
        missing_components = temp_dir / "missing_components.cif"

        args = [
            "bias",
            "-s", str(sample_system_yaml),
            "--build_bias_training_data",
            "--protein_training_data_path", str(protein_output),
            "--bias_training_components_cif", str(missing_components),
            "-w", str(temp_dir),
        ]

        with patch("cofolder.cli._load_recipe_class") as mock_loader:
            assert cli.main(args) == 2

        mock_loader.assert_not_called()

    def test_main_rejects_missing_default_components_cif_in_build_mode_before_recipe_load(
        self,
        sample_system_yaml,
        temp_dir,
    ):
        protein_output = temp_dir / "generated" / "protein_training.csv"

        args = [
            "bias",
            "-s", str(sample_system_yaml),
            "--build_bias_training_data",
            "--protein_training_data_path", str(protein_output),
            "-w", str(temp_dir),
        ]

        with patch("cofolder.cli._load_recipe_class") as mock_loader:
            assert cli.main(args) == 2

        mock_loader.assert_not_called()

    def test_main_rejects_missing_public_protein_reference_before_recipe_load(
        self,
        sample_system_yaml,
        temp_dir,
    ):
        missing_protein = temp_dir / "missing_protein.csv"
        ligand_ref = temp_dir / "ligand_training.csv"
        ligand_ref.write_text(
            "pdb_id,release_date,ligand_id,smiles\n"
            "1ABC,2022-01-01,ETH,CCO\n",
            encoding="utf-8",
        )

        args = [
            "bias",
            "-s", str(sample_system_yaml),
            "--protein_training_data_path", str(missing_protein),
            "--ligand_training_data_path", str(ligand_ref),
            "-w", str(temp_dir),
        ]

        with patch("cofolder.cli._load_recipe_class") as mock_loader:
            assert cli.main(args) == 2

        mock_loader.assert_not_called()

    def test_main_rejects_invalid_custom_protein_file_type_before_recipe_load(
        self,
        sample_system_yaml,
        temp_dir,
    ):
        custom_protein = temp_dir / "custom_protein.txt"
        custom_ligand = temp_dir / "custom_ligand.csv"
        custom_protein.write_text("sequence\nMKRAAT\n", encoding="utf-8")
        custom_ligand.write_text("smiles\nCCO\n", encoding="utf-8")

        with patch("cofolder.cli._load_recipe_class") as mock_loader:
            assert cli.main(
                [
                    "bias",
                    "-s", str(sample_system_yaml),
                    "--custom_protein_reference_path", str(custom_protein),
                    "--custom_ligand_reference_path", str(custom_ligand),
                    "-w", str(temp_dir),
                ]
            ) == 2

        mock_loader.assert_not_called()

    def test_main_accepts_custom_only_reference_inputs(
        self,
        sample_system_yaml,
        temp_dir,
    ):
        custom_protein = temp_dir / "custom_protein.csv"
        custom_ligand = temp_dir / "custom_ligand.csv"
        custom_protein.write_text("sequence,dataset_name\nMKRAAT,private_set\n", encoding="utf-8")
        custom_ligand.write_text("smiles,dataset_name\nCCO,private_set\n", encoding="utf-8")

        bias_runner = Mock()
        bias_cls = Mock(return_value=bias_runner)

        args = [
            "bias",
            "-s", str(sample_system_yaml),
            "--custom_protein_reference_path", str(custom_protein),
            "--custom_ligand_reference_path", str(custom_ligand),
            "-w", str(temp_dir),
        ]

        with patch("cofolder.cli._load_recipe_class", return_value=bias_cls):
            cli.main(args)

        assert bias_cls.call_args.kwargs["custom_protein_reference_path"] == str(custom_protein)
        assert bias_cls.call_args.kwargs["custom_ligand_reference_path"] == str(custom_ligand)
        bias_runner.run.assert_called_once()

    def test_main_rejects_invalid_custom_ligand_file_type(
        self,
        sample_system_yaml,
        temp_dir,
    ):
        custom_protein = temp_dir / "custom_protein.csv"
        custom_ligand = temp_dir / "custom_ligand.txt"
        custom_protein.write_text("sequence\nMKRAAT\n", encoding="utf-8")
        custom_ligand.write_text("smiles\nCCO\n", encoding="utf-8")

        assert cli.main(
            [
                "bias",
                "-s", str(sample_system_yaml),
                "--custom_protein_reference_path", str(custom_protein),
                "--custom_ligand_reference_path", str(custom_ligand),
                "-w", str(temp_dir),
            ]
        ) == 2
