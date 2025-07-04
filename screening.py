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
import pandas as pd
import json
import shutil

import command
import helpers
import system

import logging

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
        self.data_csv = helpers.read_csv(file_path=self.csv,
                                         columns=[self.col_smiles,self.col_id])

        ## current place ##

        lig_df = Screen.read_csv(file_path=lig_data[0].get("lig_csv"),smiles_col=smiles_col,id_col=smiles_col)
    
        # run boltz
        run_dir = helpers.set_dir(path=opt_data.get("wrapper")[0].get("run_dir"))
        
        Screen.run()

    def run(self):
  






    def run_boltz(self, run_dir, sys_data, opt_data, smiles_col, id_col, lig_df):
    msa_path = ''
    

    for i, (row_idx, row) in enumerate(lig_df.iterrows(), start=1):  
            self.logger.info(f"({i}/{len(lig_df)}) {row[id_col]}: {row[smiles_col]}")
            start_time = time.time()

            out_dir = helpers.set_dir(path=os.path.join(run_dir, f'{i}_{row[id_col]}'))
            yaml_path, msa_path = system.set_yaml(out_dir, sys_data, i, row, id_col, msa_path, smiles_col)    

            cmd = command.set_command(yaml_path, opt_data, i, row, run_dir, id_col)       
            self.logger.info(f'running command: {" ".join(cmd)}')

            subprocess.run(cmd)

            gather_metrics(run_dir, out_dir, i, row, id_col)
            gather_structures(run_dir, out_dir, i, row, id_col)

            self.logger.info(" pred time--- %s seconds ---" % (time.time() - start_time))
            
    def gather_metrics(self, run_dir, out_dir, i, row, id_col):
        # funtion for gathering confidence  and affinity metrics into csv
        csv_file = os.path.join(run_dir, "output.csv")
        
        conf_path = os.path.join(out_dir, f'boltz_results_{i}_{row[id_col]}/predictions/{i}_{row[id_col]}/confidence_{i}_{row[id_col]}_model_0.json')
        aff_path = os.path.join(out_dir, f'boltz_results_{i}_{row[id_col]}/predictions/{i}_{row[id_col]}/affinity_{i}_{row[id_col]}.json')
        
        with open(conf_path) as json_conf:
            data_conf = json.load(json_conf)
            df_conf = pd.json_normalize(data_conf)

        if os.path.exists(aff_path):
            with open(aff_path) as json_aff:
                data_aff = json.load(json_aff)
                df_aff = pd.json_normalize(data_aff)

            df_new = pd.concat([df_aff, df_conf], axis=1)
        else:
            df_new = df_conf
        
        df_new.insert(0, 'id', f'{i}_{row[id_col]}')
        
        if not os.path.exists(csv_file):
            df_new.to_csv(csv_file, index=False)
        else:
            df_new.to_csv(csv_file, mode='a', header=False, index=False)

    def gather_structures(self, run_dir, out_dir, i , row, id_col):
        # function to gather cif or pdb files into single folder
        target_dir = os.path.join(run_dir, "structures")
        os.makedirs(target_dir, exist_ok=True)
        
        stucture_dir = os.path.join(out_dir, f'boltz_results_{i}_{row[id_col]}/predictions/{i}_{row[id_col]}/')
        for file_name in os.listdir(stucture_dir):
            if file_name.endswith(('.cif', '.pdb')):
                stucture_path = os.path.join(stucture_dir, file_name)
                target_path = os.path.join(target_dir, file_name)
                shutil.copy2(stucture_path, target_path)

if __name__ == "__main__":

    # Load input files
    opt_data = system.read_yaml(file_path="options.yaml")
    sys_path = opt_data.get("wrapper")[1].get("system")
    sys_data = system.read_yaml(file_path=sys_path)

    lig_data = opt_data.get("wrapper")[2]["ligands"]
    smiles_col = lig_data[1].get("smiles_col")
    id_col = id_col=lig_data[2].get("id_col")
    
    lig_df = Screen.read_csv(file_path=lig_data[0].get("lig_csv"),smiles_col=smiles_col,id_col=smiles_col)
    
    # run boltz
    run_dir = helpers.set_dir(path=opt_data.get("wrapper")[0].get("run_dir"))
    Screen.run_boltz(
        run_dir,
        sys_data,
        opt_data,
        smiles_col,
        id_col,
        lig_df,
    )
