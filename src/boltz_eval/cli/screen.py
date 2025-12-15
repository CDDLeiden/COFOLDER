import os
import logging

from ..recipes.screen import Screen
from ..modules import utils
from ..modules.logger import initiate_logger

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

def main(args):
    """Run the virtual screening tool."""
    utils.create_dir(args.wrk_dir)  

    logger = logging.getLogger("boltz-eval.evaluate")
    initiate_logger(logger, debug=args.debug, wrk_dir=args.wrk_dir)

    logger.info("Starting Boltz-eval screening pipeline.")

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
