import argparse
import os
import logging
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

class BaseRecipe:
    LOGGER_NAME = "cofolder"

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
            conformers=args.conformers,
            sdf_file=args.sdf_file,
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

        validator = Validate(
            **BaseRecipe.common_kwargs(args)
        )

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
            required=True,
            help='Comma-separated YAML path to variable.'
        )

        parser.add_argument('--variable_csv', type=str)
        parser.add_argument('--col_variable', type=str)
        parser.add_argument('--col_id', type=str)

        parser.add_argument('--variable_sdf', type=str)
        parser.add_argument('--property_id', type=str)

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

        screener = Screen(
            wrk_dir=args.wrk_dir,
            system_path=args.system_path,
            options_path=args.options_path,
            variable=args.variable,
            variable_csv=args.variable_csv,
            col_variable=args.col_variable,
            col_id=args.col_id,
            variable_sdf=args.variable_sdf,
            property_id=args.property_id,
            conformers=args.conformers,
            merge_data=args.merge_data
        )

        screener.run()
        logger.info("Screening pipeline completed.")


class OracleRecipe(BaseRecipe):
    LOGGER_NAME = "cofolder.oracle"

    @staticmethod
    def add_arguments(parser):
        BaseRecipe.add_common_arguments(parser)
        BaseRecipe.add_final_arguments(parser)

    @staticmethod
    def main(args):
        logger = OracleRecipe.setup(args)
        logger.info("Starting COFOLDER oracle pipeline.")

        oracle = Oracle(
            wrk_dir=args.wrk_dir,
            system_path=args.system_path,
            options_path=args.options_path,
            conformers=args.conformers
        )

        oracle.run()
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

    args = parser.parse_args(argv)
    args.func(args)

