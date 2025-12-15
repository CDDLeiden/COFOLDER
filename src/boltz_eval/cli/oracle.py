import os
import logging

from ..recipes.oracle import Oracle
from ..modules import utils
from ..modules.logger import initiate_logger

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

def main(args):
    """Run the oracle recipe."""
    utils.create_dir(args.wrk_dir) 

    logger = logging.getLogger("boltz-eval.oracle")
    initiate_logger(logger, debug=args.debug, wrk_dir=args.wrk_dir)

    logger.info("Starting Boltz-eval oracle pipeline.")

    oracle = Oracle(
        wrk_dir=args.wrk_dir,
        system_path=args.system_path,
        options_path=args.options_path
    )

    oracle.run()
    logger.info("Oracle pipeline completed.")
