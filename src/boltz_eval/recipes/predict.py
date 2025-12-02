import os
import subprocess
import time
import logging

from ..helpers import command, system, utils

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

    logger = logging.getLogger('boltz-tools')
    initiate_logger(logger, debug=args.debug, wrk_dir=args.wrk_dir)
    logger.info("Boltz-tools predict started.")

    predict = Predict(
        wrk_dir=args.wrk_dir,
        system_path=args.system_path,
        options_path=args.options_path,
        generate_conformers=args.generate_conformers,
    )
    predict.run()

def initiate_logger(logger, debug, wrk_dir):
    log_file = os.path.join(wrk_dir, 'boltz-tools.log')
    open(log_file, 'w+').close()
    
    fh = logging.FileHandler(log_file)
    ch = logging.StreamHandler()
    
    fh.setLevel(logging.DEBUG if debug else logging.INFO)
    ch.setLevel(logging.WARNING)
    
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d,%H:%M:%S')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    
    logger.setLevel(logging.DEBUG if debug else logging.INFO)
    logger.addHandler(fh)
    logger.addHandler(ch)

class Predict(object):
    def __init__(self, **kwargs):
        self.wrk_dir = kwargs.get("wrk_dir")
        self.system_path = kwargs.get("system_path")
        self.options_path = kwargs.get("options_path")
        # TODO: Activate when conformer generation is implemented [in screen.py]
        #self.generate_conformers = kwargs.get("generate_conformers")
        self.logger = logging.getLogger('boltz-tools.screening.Screen')
        self.logger.debug(f"Screen args: {kwargs}")

        # Set command and system objects
        self._options = utils.read_yaml(path=self.options_path)
        self.opt = command.Command(options=self._options)
        self._system = utils.read_yaml(path=self.system_path)
        self.sys = system.System(system=self._system)

    def run(self):
        start_time = time.time()    

        # Set output directory and update system
        self.opt.out_dir = self.wrk_dir
        self.opt.system_path = self.system_path

        # Set and run command
        cmd = self.opt.set_command(system=self.sys)
        self.logger.info(f'Running: {" ".join(cmd)}')
        subprocess.run(cmd)

        self.logger.info(" pred time--- %.2f seconds ---" % (time.time() - start_time))