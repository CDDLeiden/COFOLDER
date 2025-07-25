# Script for performing virtual screening using Boltz.

# Input: 
#   - [OPTIONAL] SMILES or FASTA sequence to inject into YAML system
#   - [OPTIONAL] CSV file containing SMILES or FASTA sequences to inject into YAML system
#   - [OPTIONAL] SMILES column when iterating through CSV file
#   - [OPTIONAL] ID column when iterating through CSV file (if not specified, index will be used)
#   - YAML system file containing: 
#       - FASTA sequence of protein(s)
#       - [OPTIONAL] Co-factors
#       - [OPTIONAL] Boltz Constraints and/or Templates
# Output:
#   - output confidence metric (confidence_score, ptm, ...)    
#   - output affinity (affinity_pred_value, affinity_probability_binary, ...)
#   - RMSD of diffusion samples

# This script should (co-)fold a protein, protein-ligand complex or a virtual screening (when provided)
# with an CSV file. Results should be gathered into a singlular output file, and structures should be 
# copied to a single folder.

# EXTRA IDEAS:
#   - Toggle saving for saving either all data, only structures and output file, OR only output file
#   - Allow for Grid search of Boltz parameters, co-factors and Boltz Constraints/Templates
#   - Add IFP profiling for ligands

import os
import pandas as pd
import subprocess
import time

import json
import shutil

import command
import helpers
import system

import logging

import argparse

class Screen(object):
    def __init__(self, *args, **kwargs):
        self.wrk_dir = kwargs.get("wrk_dir")
        self.yaml_system = kwargs.get("yaml_system")
        self.yaml_options = kwargs.get("yaml_options")
        self.csv = kwargs.get("csv")
        self.variable = kwargs.get("variable")
        self.col_smiles = kwargs.get("col_smiles")
        self.col_id = kwargs.get("col_id")

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
        self.ligands = helpers.read_csv(path=self.csv, columns=[self.col_smiles, self.col_id])

        # Validate required columns in YAML and CSV
        if not self.col_smiles or not self.col_id:
            raise ValueError("Both col_smiles and col_id must be specified in the YAML file.")
        if self.ligands is None or self.col_smiles not in self.ligands.columns or self.col_id not in self.ligands.columns:
            raise ValueError(f"CSV file must contain columns: {self.col_smiles}, {self.col_id}")

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

        df_new.insert(0, 'id', f'{basename}')

        if not os.path.exists(csv_file):
            df_new.to_csv(csv_file, index=False)
        else:
            df_new.to_csv(csv_file, mode='a', header=False, index=False)

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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run virtual screening with Boltz.")
    parser.add_argument("options_path", type=str, help="Path to the options YAML file.")
    args = parser.parse_args()
    options_path = args.options_path

    logger = logging.getLogger(__name__)
    logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', level=logging.INFO)

    opt_data = helpers.read_yaml(path=options_path)
    sys_path = opt_data.get("wrapper")[1].get("system")
    sys_data = helpers.read_yaml(path=sys_path)

    lig_data = opt_data.get("wrapper")[2]["ligands"]
    smiles_col = lig_data[1].get("smiles_col")
    id_col = lig_data[2].get("id_col")
    lig_csv = lig_data[0].get("lig_csv")

    # Create and run the Screen instance
    screen = Screen(
        wrk_dir=opt_data.get("wrapper")[0].get("run_dir"),
        yaml_system=sys_path,
        yaml_options=options_path,
        csv=lig_csv,
        variable=None,
        col_smiles=smiles_col,
        col_id=id_col
    )
