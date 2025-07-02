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

import logging
import os
import pandas as pd
import subprocess
import time

import command
import helpers
import system

def read_csv(file_path, smiles_col, id_col):
    lig_df = pd.read_csv(file_path)
    logging.info(f"Read lig_csv containing {len(lig_df)} ligands")

    if smiles_col in lig_df.columns:
        logging.info(f"SMILES col ({smiles_col}) in lig_csv")
    else:
        logging.warning(f"SMILES col not found ({smiles_col}) in lig_csv")

    if id_col in lig_df.columns:
        logging.info(f"ID col ({id_col}) in lig_csv")
    else:
        logging.warning(f"ID col not found ({id_col}) in lig_csv")

    return lig_df

def run_boltz(run_dir, sys_data, opt_data, smiles_col, id_col, lig_df):
   msa_path = ''
   for i, (row_idx, row) in enumerate(lig_df.iterrows(), start=1):  
        logging.info(f"({i}/{len(lig_df)}) {row[id_col]}: {row[smiles_col]}")
        start_time = time.time()

        out_dir = helpers.set_dir(path=os.path.join(run_dir, f'{i}_{row[id_col]}'))
        yaml_path, msa_path = system.set_yaml(out_dir, sys_data, i, row, id_col, msa_path, smiles_col)    

        cmd = command.set_command(yaml_path, opt_data, i, row, run_dir, id_col)       
        logging.info(f'running command: {" ".join(cmd)}')

        subprocess.run(cmd)
        logging.info(" pred time--- %s seconds ---" % (time.time() - start_time))
        

if __name__ == "__main__":
    logger = logging.getLogger(__name__)
    logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', level=logging.INFO)

    # Load input files
    opt_data = system.read_yaml(file_path="options.yaml")
    sys_path = opt_data.get("wrapper")[1].get("system")
    sys_data = system.read_yaml(file_path=sys_path)

    lig_data = opt_data.get("wrapper")[2]["ligands"]
    smiles_col = lig_data[1].get("smiles_col")
    id_col = id_col=lig_data[2].get("id_col")
    
    lig_df = read_csv(file_path=lig_data[0].get("lig_csv"),smiles_col=smiles_col,id_col=smiles_col)
    
    # run boltz
    run_dir = helpers.set_dir(path=opt_data.get("wrapper")[0].get("run_dir"))
    run_boltz(
        run_dir,
        sys_data,
        opt_data,
        smiles_col,
        id_col,
        lig_df,
    )
 