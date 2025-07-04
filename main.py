# Parser code for command: boltz-tools.

# Should call subcommands: validate, screen and oracle and initialize the
# respective command.

# TODO: make operational after all commands are working on their own.
#       (when calling __name__=="__main__")

import argparse
import os

from predict import Predict
from screening import Screen
#from oracle import Oracle
#from validate import Validate

import helpers

import logging

def main(): 
    parser = argparse.ArgumentParser(
        prog='boltz-tools',
        description=' == Collection of helpful tools for performing various \
            co-folding tasks using Boltz. =='
        )
    
    parser.add_argument('-v', '--version',
                        action='version',
                        version='%(prog)s 0.0'
        )
    
    subparsers = parser.add_subparsers(dest='tools')

    predict = subparsers.add_parser('predict', 
                                    help = 'Run Boltz Predict for a single \
                                        system.'
        )
    
    screen = subparsers.add_parser('screen',
                                   help = 'Run Boltz for a virtual screen.'

        )  
    
    parser.add_argument('-w', '--wrk_dir',
                        dest = 'wrk_dir',
                        help = 'Working dir if different from cwd.',
                        default = os.getcwd()
        )
    
    parser.add_argument('-y', '--yaml_system',                       
                        dest = 'yaml_system',
                        help = 'Path to system YAML file.',
        )
    
    parser.add_argument('-b', '--yaml_boltz',                       
                        dest = 'yaml_boltz',
                        help = 'Path to boltz options YAML file.',
        )
    
    screen.add_argument('-c', '--csv',
                        dest = 'csv',
                        help = 'Path to ligands CSV file'
        )
    
    screen.add_argument('-v', '--variable',
                        dest = 'variable',
                        help = 'Location to inject variable into YAML system',
                        default = 'INJECT')
    
    screen.add_argument('--col_smiles',
                        dest = 'col_smiles',
                        help = 'Column containing SMILES molecule.'
        )
    
    screen.add_argument('--col_id',
                        dest = 'col_id',
                        help = 'Column constaining ID molecule',
                        default = None
        )
    
    parser.add_argument('-d', '--debug',
                        action='store_true'
        )
    
    #TODO: add subparsers 'oracle' and 'validate'

    args = parser.parse_args()

    logger = logging.getLogger('boltz-tools')
    initiate_logger(logger, debug=args.debug)

    logger.info("Boltz-tools started.")
    helpers.set_dir(args.wrk_dir)

    run_tool(args)

def initiate_logger(logger, debug):
    # create file handler which logs to info messages normally and debug messages when in debug mode.
    with open('boltz-tools.log', 'w+'):
        pass
    fh = logging.FileHandler('boltz-tools.log')

    if debug:
        logger.setLevel(logging.DEBUG)
        fh.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
        fh.setLevel(logging.INFO)

    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.WARNING)
    
    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d,%H:%M:%S')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    
    # add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)

def run_tool(args):
    if args.tools == 'predict':
        Predict(wrk_dir = args.wrk_dir,
                yaml_system = args.yaml_system,
                yaml_options = args.yaml_options)
        
    if args.tools == 'screen':
        Screen(wrk_dir = args.wrk_dir,
               yaml_system = args.yaml_system,
               yaml_options = args.yaml_options,
               csv = args.csv,
               variable = args.variable,
               col_smiles = args.col_smiles,
               col_id = args.col_id
               )

    #TODO: add run_tools for 'oracle' and 'validate'

if __name__ == "__main__":
    main()