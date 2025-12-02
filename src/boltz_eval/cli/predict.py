import os
import logging

from ..recipes.predict import Predict
from ..modules import utils
from ..modules.logger import initiate_logger

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

def main(args):
    """Run the predict tool."""
    utils.set_dir(args.wrk_dir)  

    logger = logging.getLogger("boltz-eval.predict")
    initiate_logger(logger, debug=args.debug, wrk_dir=args.wrk_dir)

    logger.info("Starting Boltz-eval prediction pipeline.")

    predictor = Predict(
        wrk_dir=args.wrk_dir,
        system_path=args.system_path,
        options_path=args.options_path,
        generate_conformers=args.generate_conformers,
    )

    predictor.run()
    logger.info("Prediction pipeline completed.")
