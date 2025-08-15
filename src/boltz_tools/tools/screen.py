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
    utils.set_dir(args.wrk_dir)  

    logger = logging.getLogger('boltz-tools')
    initiate_logger(logger, debug=args.debug, wrk_dir=args.wrk_dir)
    logger.info("Boltz-tools screen started.")

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

class Screen(object):
    def __init__(self, **kwargs):
        self.wrk_dir = kwargs.get("wrk_dir")
        self.yaml_system = kwargs.get("yaml_system")
        self.yaml_options = kwargs.get("yaml_options")
        self.csv = kwargs.get("csv")
        self.variable = kwargs.get("variable")
        self.col_smiles = kwargs.get("col_smiles")
        self.col_id = kwargs.get("col_id")
        self.merge_columns = kwargs.get("merge_columns", [])

        self.logger = logging.getLogger('boltz-tools.screening.Screen')
        self.logger.debug(f"Screen args: {kwargs}")

        # Read input files
        self.options = utils.read_yaml(path=self.yaml_options)
        self.system = system.System(path=self.yaml_system)
        self.ligands = utils.read_csv(path=self.csv, columns=[self.col_smiles, self.col_id] + self.merge_columns)

        if not self.col_smiles or not self.col_id:
            raise ValueError("Both col_smiles and col_id must be specified.")
        missing = [c for c in [self.col_smiles, self.col_id] + self.merge_columns if c not in self.ligands.columns]
        if missing:
            raise ValueError(f"CSV missing columns: {', '.join(missing)}")

    def run(self):
        msa_path = None
        for i, (_, row) in enumerate(self.ligands.iterrows(), 1):
            start_time = time.time()
            name = str(row[self.col_id])
            basename = f'{i}_{name}'
            smiles = row[self.col_smiles]
            self.logger.info(f"({i}/{len(self.ligands)}) {name}: {smiles}")

            out_dir = os.path.join(self.wrk_dir, basename)
            utils.set_dir(out_dir)

            sys_data = self.system.data.copy()
            if i != 1 and (not msa_path or not os.path.exists(msa_path)):
                raise FileNotFoundError(f"MSA path {msa_path} missing for initial system")

            yaml_path, msa_path = system.System.set_yaml(
                out_dir=out_dir,
                sys_data=sys_data,
                i=i,
                row=row,
                id_col=self.col_id,
                msa_path=msa_path,
                smiles_col=self.col_smiles
            )

            cmd = command.set_command(yaml_path, self.options, i, out_dir)
            self.logger.info(f'Running: {" ".join(cmd)}')
            subprocess.run(cmd)

            if i == 1 and msa_path and os.path.exists(msa_path):
                self.logger.info(f'Cleaning up MSA file: {msa_path}')
                utils.delete_last_line(msa_path)

            self.gather_metrics(out_dir, i, row)
            self.gather_structures(out_dir, i, row)
            self.logger.info(" pred time--- %.2f seconds ---" % (time.time() - start_time))

    def gather_metrics(self, out_dir, i, row):
        # function for gathering confidence and affinity metrics into csv
        basename = f'{i}_{row[self.col_id]}'
        csv_file = os.path.join(self.wrk_dir, "output.csv")
        conf_file = os.path.join(out_dir, f'boltz_results_{basename}/predictions/{basename}/confidence_{basename}_model_0.json')
        aff_file = os.path.join(out_dir, f'boltz_results_{basename}/predictions/{basename}/affinity_{basename}.json')

        if not os.path.exists(conf_file):
            self.logger.warning(f"Missing confidence file: {conf_file}")
            return

        try:
            df_conf = pd.json_normalize(json.load(open(conf_file)))
        except Exception as e:
            self.logger.error(f"Error reading {conf_file}: {e}")
            return

        df_new = df_conf
        if os.path.exists(aff_file):
            try:
                df_aff = pd.json_normalize(json.load(open(aff_file)))
                df_new = pd.concat([df_aff, df_conf], axis=1)
            except Exception as e:
                self.logger.error(f"Error reading {aff_file}: {e}")

        output = {**{'index': i, 'id': row[self.col_id], 'basename': basename, 'smiles': row[self.col_smiles]},
                  **{col: row[col] for col in self.merge_columns},
                  **df_new.iloc[0].to_dict()}

        pd.DataFrame([output]).to_csv(csv_file, mode='a', header=not os.path.exists(csv_file), index=False)

    def gather_structures(self, out_dir, i, row):
        # function to gather cif or pdb files into single folder
        basename = f'{i}_{row[self.col_id]}'
        target_dir = os.path.join(self.wrk_dir, "structures")
        os.makedirs(target_dir, exist_ok=True)
        structure_dir = os.path.join(out_dir, f'boltz_results_{basename}/predictions/{basename}/')
        if not os.path.exists(structure_dir):
            self.logger.warning(f"Missing structure directory: {structure_dir}")
            return
        for file_name in os.listdir(structure_dir):
            if file_name.endswith(('.cif', '.pdb')):
                try:
                    shutil.copy2(os.path.join(structure_dir, file_name), os.path.join(target_dir, file_name))
                except Exception as e:
                    self.logger.error(f"Error copying {file_name}: {e}")
