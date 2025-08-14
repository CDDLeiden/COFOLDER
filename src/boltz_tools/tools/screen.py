import os
import pandas as pd
import subprocess
import time
import json
import shutil
import logging

from .. import helpers
from ..helpers import command, system, utils

def add_arguments(parser):
    """Add screen-specific CLI arguments."""
    parser.add_argument('-w', '--wrk_dir',
                        dest='wrk_dir',
                        help='Working dir if different from cwd.',
                        default=os.getcwd())

    parser.add_argument('-y', '--yaml_system',
                        dest='yaml_system',
                        help='Path to system YAML file.')

    parser.add_argument('-b', '--yaml_boltz',
                        dest='yaml_boltz',
                        help='Path to boltz options YAML file.')

    parser.add_argument('-c', '--csv',
                        dest='csv',
                        help='Path to ligands CSV file')
    
    parser.add_argument('-v', '--variable',
                        dest='variable',
                        help='Location to inject variable into YAML system',
                        default='INJECT')
    
    parser.add_argument('--col_smiles',
                        dest='col_smiles',
                        help='Column containing SMILES molecule.')
    
    parser.add_argument('--col_id',
                        dest='col_id',
                        help='Column containing ID molecule',
                        default=None)
    
    parser.add_argument('--merge_columns',
                        type=str,
                        default=None,
                        help='Comma-separated list of CSV columns to merge into output.')
    
    parser.add_argument('-d', '--debug',
                    action='store_true',
                    help='Enable debug logging')

def main(args):
    """Run the virtual screening tool."""
    # Setup logging before running the tool
    logger = logging.getLogger('boltz-tools')
    initiate_logger(logger, debug=args.debug)

    logger.info("Boltz-tools started.")
    utils.set_dir(args.wrk_dir)

    
    merge_columns = [col.strip() for col in args.merge_columns.split(",")] if args.merge_columns else []

    screen = Screen(
        wrk_dir=args.wrk_dir,
        yaml_system=args.yaml_system,
        yaml_options=args.yaml_boltz,
        csv=args.csv,
        variable=args.variable,
        col_smiles=args.col_smiles,
        col_id=args.col_id,
        merge_columns=merge_columns
    )
    screen.run()

def initiate_logger(logger, debug):
    log_file = 'boltz-tools.log'
    with open(log_file, 'w+'):
        pass
    fh = logging.FileHandler(log_file)

    if debug:
        logger.setLevel(logging.DEBUG)
        fh.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
        fh.setLevel(logging.INFO)

    ch = logging.StreamHandler()
    ch.setLevel(logging.WARNING)

    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d,%H:%M:%S')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    logger.addHandler(fh)
    logger.addHandler(ch)

class Screen(object):
    def __init__(self, *args, **kwargs):
        self.wrk_dir = kwargs.get("wrk_dir")
        self.yaml_system = kwargs.get("yaml_system")
        self.yaml_options = kwargs.get("yaml_options")
        self.csv = kwargs.get("csv")
        self.variable = kwargs.get("variable")
        self.col_smiles = kwargs.get("col_smiles")
        self.col_id = kwargs.get("col_id")
        self.merge_columns = kwargs.get("merge_columns", [])

        self.logger = logging.getLogger('boltz-tools.screening.Screen')

        self.logger.debug('Screen arguments initialized')
        self.logger.debug(f'self.wrk_dir = {str(self.wrk_dir)}')
        self.logger.debug(f'self.yaml_system = {str(self.yaml_system)}')
        self.logger.debug(f'self.yaml_options = {str(self.yaml_options)}')
        self.logger.debug(f'self.csv = {str(self.csv)}')
        self.logger.debug(f'self.variable = {str(self.variable)}')
        self.logger.debug(f'self.col_smiles = {str(self.col_smiles)}')
        self.logger.debug(f'self.col_id = {str(self.col_id)}')

        # Read input files
        self.options = helpers.read_yaml(path=self.yaml_options)
        self.system = system.System(path=self.yaml_system)
        self.ligands = helpers.read_csv(path=self.csv, columns=[self.col_smiles, self.col_id] + self.merge_columns)

        # Validate required columns in YAML and CSV
        if not self.col_smiles or not self.col_id:
            raise ValueError("Both col_smiles and col_id must be specified in the YAML file.")
        missing_cols = [col for col in [self.col_smiles, self.col_id] + self.merge_columns if col not in self.ligands.columns]
        if self.ligands is None or missing_cols:
            raise ValueError(f"CSV file must contain columns: {', '.join(missing_cols)}")

        self.logger.debug(f'self.options = {str(self.options)}')
        self.logger.debug(f'self.system = {str(self.system)}')
        self.logger.debug(f'self.ligands = {str(self.ligands)}')

        self.run()

    def run(self):
        msa_path = None
        for i, (row_idx, row) in enumerate(self.ligands.iterrows(), start=1):
            start_time = time.time()

            name = str(row[self.col_id]) if self.col_id else str(i)
            basename = f'{i}_{name}'
            smiles = row[self.col_smiles]
            self.logger.debug(f'self.name = {name}')
            self.logger.debug(f'self.smiles = {smiles}')
            self.logger.info(f"({i}/{len(self.ligands)}) {name}: {smiles}")

            out_dir = os.path.join(self.wrk_dir, basename)
            helpers.set_dir(path=out_dir)

            sys_data = self.system.data.copy()
            # For the first ligand, calculate MSA; for others, reuse
            # current_msa_path = None if i == 1 else msa_path
            if i != 1 and (not msa_path or not os.path.exists(msa_path)):
                self.logger.error(f"MSA path {msa_path} does not exist for ligand {name}. Cannot reuse MSA.")
                raise FileNotFoundError(f"MSA path {msa_path} does not exist for ligand {name}. Cannot reuse MSA.")
            yaml_path, msa_path = system.System.set_yaml(
                out_dir=out_dir,
                sys_data=sys_data,
                i=i,
                row=row,
                id_col=self.col_id,
                msa_path=msa_path,
                smiles_col=self.col_smiles
            )

            cmd = command.set_command(yaml_path, self.options, i, row, out_dir, self.col_id)
            self.logger.info(f'running command: {" ".join(cmd)}')
            subprocess.run(cmd)

            # After first ligand, store and clean up MSA
            if i == 1:
                if msa_path and os.path.exists(msa_path):
                    self.logger.info(f'Cleaning up MSA file: {msa_path}')
                    helpers.delete_last_line(msa_path)

            self.gather_metrics(self.wrk_dir, out_dir, i, row, self.col_id)
            self.gather_structures(self.wrk_dir, out_dir, i, row, self.col_id)

            self.logger.info(" pred time--- %s seconds ---" % (time.time() - start_time))



    def gather_metrics(self, run_dir, out_dir, i, row, id_col):
        # function for gathering confidence and affinity metrics into csv
        basename = f'{i}_{row[id_col]}'
        csv_file = os.path.join(run_dir, "output.csv")

        conf_path = os.path.join(out_dir, f'boltz_results_{basename}/predictions/{basename}/confidence_{basename}_model_0.json')
        aff_path = os.path.join(out_dir, f'boltz_results_{basename}/predictions/{basename}/affinity_{basename}.json')

        if not os.path.exists(conf_path):
            self.logger.warning(f"Missing confidence file: {conf_path}")
            return
        try:
            with open(conf_path) as json_conf:
                data_conf = json.load(json_conf)
                df_conf = pd.json_normalize(data_conf)
        except Exception as e:
            self.logger.error(f"Error reading confidence file {conf_path}: {e}")
            return

        if os.path.exists(aff_path):
            try:
                with open(aff_path) as json_aff:
                    data_aff = json.load(json_aff)
                    df_aff = pd.json_normalize(data_aff)
                df_new = pd.concat([df_aff, df_conf], axis=1)
            except Exception as e:
                self.logger.error(f"Error reading affinity file {aff_path}: {e}")
                df_new = df_conf
        else:
            df_new = df_conf

        # Build a single output row as a dictionary
        output_row = {
            'index': i,
            'id': row[id_col] if id_col in row else None,
            'basename': basename,
            'smiles': row[self.col_smiles] if self.col_smiles in row else None
        }
        # Add merge_columns, skipping any already added
        for col in self.merge_columns:
            if col not in output_row:
                output_row[col] = row[col] if col in row else None
        # Add all metrics columns from df_new (flattened)
        for col in df_new.columns:
            output_row[col] = df_new.iloc[0][col]
        # Write to CSV
        output_df = pd.DataFrame([output_row])
        if not os.path.exists(csv_file):
            output_df.to_csv(csv_file, index=False)
        else:
            output_df.to_csv(csv_file, mode='a', header=False, index=False)

    def gather_structures(self, run_dir, out_dir, i , row, id_col):
        # function to gather cif or pdb files into single folder
        basename = f'{i}_{row[id_col]}'
        target_dir = os.path.join(run_dir, "structures")
        os.makedirs(target_dir, exist_ok=True)

        structure_dir = os.path.join(out_dir, f'boltz_results_{basename}/predictions/{basename}/')
        if not os.path.exists(structure_dir):
            self.logger.warning(f"Missing structure directory: {structure_dir}")
            return
        for file_name in os.listdir(structure_dir):
            if file_name.endswith(('.cif', '.pdb')):
                structure_path = os.path.join(structure_dir, file_name)
                target_path = os.path.join(target_dir, file_name)
                try:
                    shutil.copy2(structure_path, target_path)
                except Exception as e:
                    self.logger.error(f"Error copying {structure_path} to {target_path}: {e}")
