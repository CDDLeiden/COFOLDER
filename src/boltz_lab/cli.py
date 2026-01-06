import argparse
import os
import logging

from boltz_lab import __version__

from boltz_lab.recipes.predict import Predict
from boltz_lab.recipes.screen import Screen
from boltz_lab.recipes.oracle import Oracle
from boltz_lab.recipes.evaluate import Evaluate

from boltz_lab.modules.utils import helpers
from boltz_lab.modules.utils.log import setup_root_logger

class PredictRecipe:
    @staticmethod
    def add_arguments(parser):
        """Add predict-specific CLI arguments."""
        parser.add_argument('-w', '--wrk_dir',
                            type=str,
                            dest='wrk_dir',
                            help='Set Working directory if different from CWD.',
                            default=os.getcwd())

        parser.add_argument('-s', '--system_path',
                            type=str,
                            dest='system_path',
                            help='Path to system YAML file.',
                            required=True)

        parser.add_argument('-b', '--boltz_options_path',
                            type=str,
                            dest='options_path',
                            help='Path to Boltz options YAML file.',
                            required=True)
        
        parser.add_argument('--generate_conformers',
                            choices=['2D', '3D'],
                            default=None,
                            dest='generate_conformers',
                            help='Generate 2D or 3D conformers for CCD input. If not specified, \
                                original SMILES (csv) or MolBlock (sdf) are used as system input. \
                                Note: This only works for SMILES, not any other variable type.')
        
        parser.add_argument('-d', '--debug',
                            action='store_true',
                            help='Enable debug logging')

    @staticmethod
    def main(args):
        """Run the predict tool."""
        helpers.create_dir(args.wrk_dir)  

        logger = logging.getLogger("boltz-lab.predict")

        logger.info("Starting Boltz-lab prediction pipeline.")

        predictor = Predict(
            wrk_dir=args.wrk_dir,
            system_path=args.system_path,
            options_path=args.options_path,
            generate_conformers=args.generate_conformers, # FIXME : Unexpected argument
        )

        predictor.run()
        logger.info("Prediction pipeline completed.")


class EvaluateRecipe:
    @staticmethod
    def add_arguments(parser):
        """Add evaluate-specific CLI arguments."""
        parser.add_argument('-w', '--wrk_dir',
                            type=str,
                            dest='wrk_dir',
                            help='Set Working directory if different from CWD.',
                            default=os.getcwd())

        parser.add_argument('-s', '--system_path',
                            type=str,
                            dest='system_path',
                            help='Path to system YAML file.',
                            required=True)

        parser.add_argument('-b', '--boltz_options_path',
                            type=str,
                            dest='options_path',
                            help='Path to Boltz options YAML file.',
                            required=True)

        parser.add_argument('--repeats',
                            type=int,
                            default=1,
                            dest='repeats',
                            help=('Number of repeats to run the calculation with different seeds. '
                                'Each repeat will use a different random seed if unspecified. Must be >=1.'))

        parser.add_argument('--seeds',
                            type=str,
                            default=None,
                            dest='seeds',
                            help=('Optional, A comma-seperated list of integer seeds to use for repeats. '
                                'Must be the same length as --repeats, e.g. --repeats 3 --seeds 42,123,999 '))
        
        parser.add_argument('-i', '--input_pdb',
                            type=str,
                            dest='input_pdb',
                            help='Path to reference CIF/PDB file containing a (ligand-)protein system. \
                                By adding a system, protein and ligand RMSD calculations will be \
                                performed between the co-folded and reference structures.')

        parser.add_argument('--ifp',
                            type=str,
                            dest='ifp',
                            help=('Interaction fingerprint (IFP) specification. Accepts: \n'
                                '  - "true": extract from crystal structure (requires input_pdb)\n'
                                '  - "false": skip IFP overlap calculation\n'
                                '  - A dictionary of interacting residues for manual specification, e.g.\n'
                                '      \'{"A:123":"ARG", "B:45":"TYR"}\''))

        parser.add_argument('--generate_conformers',
                            choices=['2D', '3D'],
                            default=None,
                            dest='generate_conformers',
                            help='Generate 2D or 3D conformers for CCD input. If not specified, \
                                original SMILES (csv) or MolBlock (sdf) are used as system input. \
                                Note: This only works for SMILES, not any other variable type.')
        
        parser.add_argument('--log_name',
                            type=str,
                            default='log',
                            help='Base name for log file. Defaults to "log".')

        parser.add_argument('-d', '--debug',
                            action='store_true',
                            help='Enable debug logging')
    
    @staticmethod
    def main(args):
        """Run the evaluate system recipe."""
        setup_root_logger(logging.DEBUG) if args.debug else setup_root_logger(logging.INFO)    
        logger = logging.getLogger(__name__)
        logger.info(f"Logger initialized. Log file: {os.path.join(args.wrk_dir, args.log_name)+'.log'}, debug={args.debug}")

        helpers.create_dir(args.wrk_dir) 
        
        logger.info("Starting Boltz-lab evaluation pipeline.")
        evaluator = Evaluate(
            wrk_dir=args.wrk_dir,
            system_path=args.system_path,
            options_path=args.options_path,
            repeats=args.repeats,
            seeds=args.seeds,
            input_pdb=args.input_pdb,
            ifp=args.ifp,
            generate_conformers=args.generate_conformers
        )

        evaluator.run()
        logger.info("Evaluation pipeline completed.")


class ScreenRecipe:
    @staticmethod
    def add_arguments(parser):
        """Add screen-specific CLI arguments."""
        parser.add_argument('-w', '--wrk_dir',
                            type=str,
                            dest='wrk_dir',
                            help='Set Working directory if different from CWD.',
                            default=os.getcwd())

        parser.add_argument('-s', '--system_path',
                            type=str,
                            dest='system_path',
                            help='Path to system YAML file.',
                            required=True)

        parser.add_argument('-b', '--boltz_options_path',
                            type=str,
                            dest='options_path',
                            help='Path to Boltz options YAML file.',
                            required=True)
        
        parser.add_argument('-v', '--variable',
                            type=str,
                            dest='variable',
                            help='Comma-seperated list of keys specifying the nested path in \
                                the system YAML to update. (e.g. "sequences,1,ligand,smiles")',
                            required=True)
        
        parser.add_argument('-c', '--variable_csv',
                            type=str,
                            default=None,
                            dest='variable_csv',
                            help='Path to CSV file containing variables. If provided, you must \
                                also specify --col_variable and --col_id.')
        
        parser.add_argument('--col_variable',
                            type=str,
                            default=None,
                            dest='col_variable',
                            help='Column containing variable (e.g. SMILES/CCD/FASTA) (required \
                                if --csv is used).')
        
        parser.add_argument('--col_id',
                            type=str,
                            default=None,
                            dest='col_id',
                            help='Column containing variable ID (required if --csv is used).')
        
        parser.add_argument('-s,', '--variable_sdf',
                            type=str,
                            default=None,
                            dest='variable_sdf',
                            help='Path to SDF file containing variables. If provided, you must \
                                also specify --propterty_id.')
        
        parser.add_argument('--property_id',
                            type=str,
                            default=None,
                            dest='property_id',
                            help='Property name for compound ID in SDF file. Required if \
                                --variable_sdf is used.')
        
        parser.add_argument('--generate_conformers',
                            choices=['2D', '3D'],
                            default=None,
                            dest='generate_conformers',
                            help='Generate 2D or 3D conformers for CCD input. If not specified, \
                                original SMILES (csv) or MolBlock (sdf) are used as system input. \
                                Note: This only works for SMILES, not any other variable type.')
        
        parser.add_argument('--merge_data',
                            type=str,
                            default=None,
                            help='Comma-separated list of CSV columns or SDF prperties from input \
                                to merge into output. (optional)')
        
        parser.add_argument('-d', '--debug',
                            action='store_true',
                            help='Enable debug logging')

    @staticmethod
    def main(args):
        """Run the virtual screening tool."""
        helpers.create_dir(args.wrk_dir)  

        logger = logging.getLogger("boltz-lab.evaluate")
        #initiate_logger(logger, debug=args.debug, wrk_dir=args.wrk_dir)

        logger.info("Starting Boltz-lab screening pipeline.")

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
            generate_conformers=args.generate_conformers,
            merge_data=args.merge_data,
        )

        screener.run()
        logger.info("Screening pipeline completed.")


class OracleRecipe:
    @staticmethod
    def add_arguments(parser):
        """Add oracle-specific CLI arguments."""
        parser.add_argument('-w', '--wrk_dir',
                            type=str,
                            dest='wrk_dir',
                            help='Set Working directory if different from CWD.',
                            default=os.getcwd())

        parser.add_argument('-s', '--system_path',
                            type=str,
                            dest='system_path',
                            help='Path to system YAML file.',
                            required=True)

        parser.add_argument('-b', '--boltz_options_path',
                            type=str,
                            dest='options_path',
                            help='Path to Boltz options YAML file.',
                            required=True)
        
        parser.add_argument('-d', '--debug',
                        action='store_true',
                        help='Enable debug logging')

    @staticmethod
    def main(args):
        """Run the oracle recipe."""
        helpers.create_dir(args.wrk_dir) 

        logger = logging.getLogger("boltz-lab.oracle")
        #initiate_logger(logger, debug=args.debug, wrk_dir=args.wrk_dir)

        logger.info("Starting Boltz-lab oracle pipeline.")

        oracle = Oracle(
            wrk_dir=args.wrk_dir,
            system_path=args.system_path,
            options_path=args.options_path
        )

        oracle.run()
        logger.info("Oracle pipeline completed.")


RECIPES = [
    ("predict", PredictRecipe, "Basic protocol for co-folding a single system using Boltz."),
    ("screen", ScreenRecipe, "Co-fold a library using Boltz for virtual screening."),
    ("oracle", OracleRecipe, "Use Boltz as an oracle function for single SMILES predictions."),
    ("evaluate", EvaluateRecipe, "Evaluate Boltz system configuration."),
]


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog='boltz-lab',
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


