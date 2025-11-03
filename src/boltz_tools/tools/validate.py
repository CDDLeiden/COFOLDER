# Script for validating Boltz system based on ability to recreate poses / find 
# correct binding pocket in reference PDB structure

import os
import logging

from ..helpers import command, conformers, system, utils

def add_arguments(parser):
    """Add validate-specific CLI arguments."""
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
                        type=list or int,
                        default=None,
                        dest='seeds',
                        help=('Optional, List of integer seeds to use for repeats. '
                              'Must be the same length as --repeats, e.g. "[42, 123, 999]".'))
    
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
    """Run the validate system tool."""
    utils.set_dir(args.wrk_dir) 

    logger = logging.getLogger('boltz-tools')
    initiate_logger(logger, debug=args.debug, wrk_dir=args.wrk_dir)
    logger.info("Boltz-tools validation started.")

    # Validate seeds length if provided
    if args.seeds is not None and len(args.seeds) != args.repeats:
        raise ValueError(f"Number of seeds ({len(args.seeds)}) must match --repeats ({args.repeats}).")

    validate = Validate(
         wrk_dir=args.wrk_dir,
         system_path=args.system_path,
         options_path=args.options_path,
         input_pdb=args.input_pdb,
         ifp=args.ifp,
         generate_conformers=args.generate_conformers,
         repeats=args.repeats,
         seeds=args.seeds
         )

def initiate_logger(logger, debug, wrk_dir):
    log_file = os.path.join(wrk_dir, 'boltz-tools.log')
    open(log_file, 'w+').close()

    fh = logging.FileHandler(log_file)
    ch = logging.StreamHandler()

    fh.setLevel(logging.DEBUG if debug else logging.INFO)
    ch.setLevel(logging.WARNING)

    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d,%H:%M:%S'
    )
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    logger.setLevel(logging.DEBUG if debug else logging.INFO)
    logger.addHandler(fh)
    logger.addHandler(ch)

class Validate(object):
    """
    Class to validate Boltz system based on ability to recreate poses / find 
    correct binding pocket in reference PDB structure.
    """
    def __init__(self, **kwargs):
        self.wrk_dir = kwargs.get("wrk_dir")
        self.system_path = kwargs.get("system_path")
        self.options_path = kwargs.get("options_path")
        self.repeats = kwargs.get("repeats")
        self.seeds = kwargs.get("seeds")
        self.input_pdb = kwargs.get("input_pdb")
        self.ifp = kwargs.get("ifp")
        self.generate_conformers = kwargs.get("generate_conformers")

        self.logger = logging.getLogger('boltz-tools.validate.Validate')
        self.logger.debug(f"Validate args: {kwargs}")

        # Set command and system objects
        self._options = utils.read_yaml(path=self.options_path)
        self.opt = command.Command(options=self._options)
        self._system = utils.read_yaml(path=self.system_path)
        self.sys = system.System(system=self._system)

        self.run()
