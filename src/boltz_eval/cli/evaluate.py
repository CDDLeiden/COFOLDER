import os
import logging

from ..recipes.evaluate import Evaluate
from ..modules import utils
from ..modules.logger import initiate_logger

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
    
    parser.add_argument('-d', '--debug',
                        action='store_true',
                        help='Enable debug logging')
    
def main(args):
    """Run the evaluate system recipe."""
    utils.create_dir(args.wrk_dir) 

    logger = logging.getLogger("evaluate")
    initiate_logger(logger, debug=args.debug, wrk_dir=args.wrk_dir)

    logger.info("Starting Boltz-eval evaluation pipeline.")

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
