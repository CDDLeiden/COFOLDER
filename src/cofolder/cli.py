import argparse
import os
import logging
import shlex
import sys
from pathlib import Path

from cofolder import __version__

from cofolder.recipes.validate import Validate
from cofolder.recipes.screen import Screen
from cofolder.recipes.oracle import Oracle

from cofolder.modules.utils import helpers
from cofolder.modules.utils.log import setup_root_logger

SCORING_FUNCTIONS = [
    "boltz_confidence_metrics",
    "boltz_affinity_metrics",
    "boltz_affinity_metrics_ext",
    "ifp_distance",
    "ifp_prolif",
    "sasa",
    "sasa_normalized",
]

DEFAULT_SCORING_FUNCTIONS = [
    f for f in SCORING_FUNCTIONS if f != "ifp_prolif"
]

REPRODUCTION_METRICS = [
    "protein_rmsd",
    "ligand_rmsd",
    "sucos",
    "pocket_coverage",
]

class BaseRecipe:
    LOGGER_NAME = "cofolder"

    @staticmethod
    def _format_cli_command(args) -> str:
        raw_argv = list(getattr(args, "_cli_argv", []) or [])
        return shlex.join(["cofolder", *raw_argv])

    @staticmethod
    def add_common_arguments(parser):
        parser.add_argument(
            '-w', '--wrk_dir',
            type=str,
            default=os.getcwd(),
            help='Set working directory (default: CWD).'
        )

        parser.add_argument(
            '-s', '--system_path',
            type=str,
            required=True,
            help='Path to system YAML file.'
        )

        parser.add_argument(
            '-b', '--boltz_options_path',
            dest='options_path',
            type=str,
            required=True,
            help='Path to Boltz options YAML file.'
        )

        parser.add_argument(
            '--repeats',
            type=int,
            default=1,
            help='Number of repeats.'
        )

        parser.add_argument(
            '--seed',
            type=int,
            default=None,
            help='Global seed used for predictions.'
        )

        parser.add_argument(
            '--scoring_functions',
            nargs='+',
            choices=SCORING_FUNCTIONS,
            metavar="SCORING_FUNCTION",
            default=DEFAULT_SCORING_FUNCTIONS,
            help=(
                "Scoring functions to compute. "
                f"Choices: {', '.join(SCORING_FUNCTIONS)}"
            ),
        )

        parser.add_argument(
            '--assess_robustness',
            dest='assess_robustness',
            action=argparse.BooleanOptionalAction,
            default=True,
            help='Assess robustness across repeats and diffusion samples (default: True).'
        )
        parser.add_argument(
            '--assess_bias',
            '--asess_bias',
            dest='assess_bias',
            action=argparse.BooleanOptionalAction,
            default=False,
            help='Assess bias against training data references (default: False).'
        )
        parser.add_argument(
            '--protein_training_data_path',
            type=str,
            default=None,
            help='Path to protein training reference CSV (requires release_date, pdb_id, sequence).'
        )
        parser.add_argument(
            '--ligand_training_data_path',
            type=str,
            default=None,
            help=(
                "Path to ligand training reference CSV/SDF "
                "(requires release_date, pdb_id, smiles) when --no-build_bias_training_data is used. "
                "When --build_bias_training_data is enabled, this is treated as output path for generated "
                "ligand training data. If omitted in build mode, defaults to "
                "<wrk_dir>/results/bias_train/ligand_training_data.csv."
            )
        )
        parser.add_argument(
            '--bias_release_cutoff',
            type=str,
            default='2023-06-01',
            help='Release-date cutoff (YYYY-MM-DD) for training data filtering.'
        )
        parser.add_argument(
            '--bias_ligand_similarity_threshold',
            type=float,
            default=0.35,
            help='Ligand ECFP/Tanimoto threshold used while building bias training data.'
        )
        parser.add_argument(
            '--bias_chains',
            nargs='+',
            default=None,
            help=(
                "Optional chain IDs to restrict bias assessment to. "
                "Example: --bias_chains A F"
            )
        )
        parser.add_argument(
            '--build_bias_training_data',
            action=argparse.BooleanOptionalAction,
            default=False,
            help='Build bias training CSVs inside validate() before assessing bias (default: False).'
        )
        parser.add_argument(
            '--bias_training_components_cif',
            type=str,
            default=None,
            help=(
                "Path to components.cif used by in-validate bias training build. "
                "Default: <protein_training_data_path parent>/ccd/components.cif"
            ),
        )

        parser.add_argument(
            '--conformers',
            choices=['2D', '3D', 'sdf'],
            default=None,
            help=(
                'Generate 2D or 3D conformers for CCD input, or use conformers '
                'from an existing SDF file. '
                'Only valid for SMILES-based inputs.'
            )
        )

        parser.add_argument(
            '--sdf_file',
            type=str,
            default=None,
            help=(
                'Path to an SDF file containing conformers to use. '
                'Required if --conformers is set to "sdf".'
            )
        )

        parser.add_argument(
            '--reference_path',
            type=str,
            default=None,
            help='Path to reference structure (PDB/CIF) used for model reproduction metrics.'
        )
        parser.add_argument(
            "--pocket_coverage_reference",
            type=str,
            default=None,
            help=(
                "Custom pocket-coverage reference (residue list or bitstring). "
                "Examples: '2 8 10', 'A2 S8 T10', '0100000101'. "
                "You may also provide a text file path containing one of these formats."
            ),
        )
        parser.add_argument(
            "--reproduction_metrics",
            nargs="+",
            choices=REPRODUCTION_METRICS,
            default=REPRODUCTION_METRICS,
            help=(
                "Reference-based reproduction metrics to compute when --reference_path is set. "
                f"Choices: {', '.join(REPRODUCTION_METRICS)}"
            ),
        )

    @staticmethod
    def add_final_arguments(parser):
        parser.add_argument(
            '--log_name',
            type=str,
            default='log',
            help='Base name for log file.'
        )

        parser.add_argument(
            '-d', '--debug',
            action='store_true',
            help='Enable debug logging.'
        )

    @staticmethod
    def common_kwargs(args):
        return dict(
            wrk_dir=args.wrk_dir,
            system_path=args.system_path,
            options_path=args.options_path,
            repeats=args.repeats,
            seed=args.seed,
            scoring_functions=args.scoring_functions,
            assess_robustness=args.assess_robustness,
            assess_bias=args.assess_bias,
            protein_training_data_path=args.protein_training_data_path,
            ligand_training_data_path=args.ligand_training_data_path,
            bias_release_cutoff=args.bias_release_cutoff,
            bias_ligand_similarity_threshold=args.bias_ligand_similarity_threshold,
            bias_chains=args.bias_chains,
            build_bias_training_data=args.build_bias_training_data,
            bias_training_components_cif=args.bias_training_components_cif,
            conformers=args.conformers,
            sdf_file=args.sdf_file,
            reference_path=args.reference_path,
            pocket_coverage_reference=args.pocket_coverage_reference,
            reproduction_metrics=args.reproduction_metrics,
        )

    @staticmethod
    def _validate_common_args(args):
        # ---- scoring functions validation ----
        if args.scoring_functions is not None:
            invalid = set(args.scoring_functions) - set(SCORING_FUNCTIONS)
            if invalid:
                raise ValueError(f"Invalid --scoring_functions: {sorted(invalid)}")
        
        # ---- path validation ----
        system_path = Path(args.system_path)
        if not system_path.exists():
            raise ValueError(f'--system_path does not exist: {system_path}')
        if not system_path.is_file():
            raise ValueError(f'--system_path is not a file: {system_path}')

        options_path = Path(args.options_path)
        if not options_path.exists():
            raise ValueError(f'--boltz_options_path does not exist: {options_path}')
        if not options_path.is_file():
            raise ValueError(f'--boltz_options_path is not a file: {options_path}')
        
        if args.sdf_file is not None:
            sdf_path = Path(args.sdf_file)
            if not sdf_path.exists():
                raise ValueError(f'--sdf_file does not exist: {sdf_path}')
            if not sdf_path.is_file():
                raise ValueError(f'--sdf_file is not a file: {sdf_path}')

        # ---- conformer/sdf validation ----
        if args.conformers == 'sdf' and args.sdf_file is None:
            raise ValueError(
                '--sdf_file must be provided when --conformers is "sdf"'
            )
        if args.assess_bias:
            # Runtime availability checks are handled in Validate.run so the
            # pipeline can continue while warning and skipping bias metrics.
            pass

    @classmethod
    def setup(cls, args):
        cls._validate_common_args(args)
        
        helpers.create_dir(args.wrk_dir)

        log_file = Path(args.wrk_dir) / f"{args.log_name}.log"

        setup_root_logger(
            level=logging.DEBUG if args.debug else logging.INFO,
            log_file=log_file
        )

        logger = logging.getLogger(cls.LOGGER_NAME)
        logger.info("CLI command: %s", cls._format_cli_command(args))
        logger.info(
            "Logger initialized. "
            f"Log file: {log_file}, debug={args.debug}"
        )

        return logger

class ValidateRecipe(BaseRecipe):
    LOGGER_NAME = "cofolder.validate"

    @staticmethod
    def add_arguments(parser):
        BaseRecipe.add_common_arguments(parser)
        BaseRecipe.add_final_arguments(parser)

    @staticmethod
    def main(args):
        logger = ValidateRecipe.setup(args)
        logger.info("Starting COFOLDER validation pipeline.")

        validator = Validate(**BaseRecipe.common_kwargs(args))

        validator.run()
        logger.info("Validation pipeline completed.")


class ScreenRecipe(BaseRecipe):
    LOGGER_NAME = "cofolder.screen"

    @staticmethod
    def add_arguments(parser):
        BaseRecipe.add_common_arguments(parser)

        parser.add_argument(
            '-v', '--variable',
            type=str,
            action='append',
            default=None,
            help='Repeatable comma-separated YAML path(s) to update.'
        )

        parser.add_argument('-c', '--variable_csv', type=str, required=True)
        parser.add_argument('--col_variable', type=str, action='append', default=None)
        parser.add_argument('--col_id', type=str, required=True)

        parser.add_argument(
            '--merge_data',
            type=str,
            help='Comma-separated list of metadata fields to merge.'
        )

        BaseRecipe.add_final_arguments(parser)

    @staticmethod
    def main(args):
        logger = ScreenRecipe.setup(args)
        logger.info("Starting COFOLDER screening pipeline.")

        if not args.variable or not args.col_variable:
            raise ValueError("At least one --variable/--col_variable pair is required.")
        if len(args.variable) != len(args.col_variable):
            raise ValueError(
                "Number of --variable entries must match number of --col_variable entries. "
                f"Got variable={len(args.variable)} col_variable={len(args.col_variable)}."
            )

        screener = Screen(
            variable=args.variable,
            variable_csv=args.variable_csv,
            col_variable=args.col_variable,
            col_id=args.col_id,
            merge_data=args.merge_data,
            **BaseRecipe.common_kwargs(args),
        )

        screener.run()
        logger.info("Screening pipeline completed.")


class OracleRecipe(BaseRecipe):
    LOGGER_NAME = "cofolder.oracle"

    @staticmethod
    def add_arguments(parser):
        BaseRecipe.add_common_arguments(parser)
        input_group = parser.add_mutually_exclusive_group(required=True)
        input_group.add_argument(
            "--input_smiles",
            type=str,
            default=None,
            help="Single query SMILES for oracle scoring."
        )
        input_group.add_argument(
            "--input_mol_file",
            type=str,
            default=None,
            help="Single query MOL/SDF file for oracle scoring."
        )
        parser.add_argument(
            "--output_metric",
            type=str,
            required=True,
            help="Metric column name to return as oracle value (e.g. affinity_pred_value)."
        )
        parser.add_argument(
            "--aggregate",
            type=str,
            choices=["first", "mean", "max", "min", "median"],
            default="first",
            help="Aggregation applied when multiple metric rows are present."
        )
        BaseRecipe.add_final_arguments(parser)

    @staticmethod
    def main(args):
        logger = OracleRecipe.setup(args)
        logger.info("Starting COFOLDER oracle pipeline.")

        oracle = Oracle(
            input_smiles=args.input_smiles,
            input_mol_file=args.input_mol_file,
            output_metric=args.output_metric,
            aggregate=args.aggregate,
            **BaseRecipe.common_kwargs(args),
        )

        value = oracle.run()
        logger.info("Oracle value (%s, aggregate=%s): %s", args.output_metric, args.aggregate, value)
        logger.info("Oracle pipeline completed.")


RECIPES = [
    ("validate", ValidateRecipe, "Basic protocol for co-folding and validating a single system using Boltz."),
    ("screen", ScreenRecipe, "Co-fold a library using Boltz for virtual screening."),
    ("oracle", OracleRecipe, "Use Boltz as an oracle function for single SMILES predictions."),
]


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog='cofolder',
        description='Collection of helpful tools for performing various co-folding tasks using Boltz.'
    )
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s {__version__}')
    
    subparsers = parser.add_subparsers(dest='command', required=True)

    for name, recipe, help_text in RECIPES:
        sp = subparsers.add_parser(name, help=help_text)
        recipe.add_arguments(sp)
        sp.set_defaults(func=recipe.main)

    raw_argv = list(argv) if argv is not None else sys.argv[1:]
    args = parser.parse_args(raw_argv)
    args._cli_argv = raw_argv
    args.func(args)
