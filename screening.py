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
import multiprocessing
import yaml

def read_yaml(file_path):
    try:
        with open(file_path, 'r') as file:
            data = yaml.safe_load(file)
            logging.info(f"{file_path} loaded successfully. Contents:")
            for key, value in data.items():
                  logging.info(f"\t{key}: {value}")
            
            return data

    except FileNotFoundError:
        logging.error(f"File not found: {file_path}")
    except yaml.YAMLError as e:
        logging.error(f"Error parsing YAML: {e}")

    return None

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

def run_boltz(run_dir, sys_data, smiles_col, id_col, lig_df):
   msa_path = ''
   for i, (row_idx, row) in enumerate(lig_df.iterrows(), start=1):  
        logging.info(f"({i}/{len(lig_df)}) {row[id_col]}: {row[smiles_col]}")
        
        out_dir = set_dir(path=os.path.join(run_dir, f'{i}_{row[id_col]}'))
        yaml_path, msa_path = set_yaml(out_dir, sys_data, i, row, msa_path)    

        command = set_command(yaml_path, i, row)       
        logging.info(f'running command: {" ".join(command)}')

        subprocess.run(command)
        


def set_yaml(out_dir, sys_data, i, row, msa_path):
    out = os.path.join(out_dir, f'{i}_{row[id_col]}.yaml')
    sys_data['sequences'][1]['ligand']['smiles'] = row[smiles_col]
    
    if i == 1:
        base_file = os.path.basename(out)
        base_name = os.path.splitext(base_file)[0]
        msa_path = os.path.join(out_dir, f'boltz_results_{base_name}/msa/{base_name}_unpaired_tmp_env/uniref.a3m')
        logging.info(f'set msa: {msa_path}')
    if i != 1:
        sys_data['sequences'][0]['protein']['msa'] = msa_path
    
    logging.info(f'this is the new sys_data: {sys_data}')

    with open(out, "w") as file:
        yaml.dump(sys_data, file, sort_keys=False)

    return out, msa_path

def set_command(yaml_path, i, row):
    command = [
        "boltz",
        "predict",
        yaml_path
    ]

    if i == 1:
        command.extend([f"--use_msa_server"])

    options = opt_data.get("options")
    for n, item in enumerate(options):
        for key, value in item.items():
            value = str(value)

            if key == "out_dir":
                command.extend([f"--out_dir"])
                command.extend([os.path.join(run_dir, f'{i}_{row[id_col]}')])   
                continue     
            
            if key == "use_msa_server":
                continue

            if value == "None" or value == "False":
                continue

            if value == "True":
                command.extend([f"--{key}"])
                continue

            if value == "multiprocessing.cpu_count()":
                command.extend([f"--{key}"])
                command.extend([str(multiprocessing.cpu_count())])
                continue
            
            command.extend([f"--{key}"])
            command.extend([str(value)]) 

    return command


def set_dir(path):
    if path:
        os.makedirs(path, exist_ok=True)
        logging.info(f"Created or verified existence of directory: {path}")
    else:
        logging.warning("No path ound in wrapper entry")

    return path

if __name__ == "__main__":
    logger = logging.getLogger(__name__)
    logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', level=logging.INFO)

    # Load input files
    opt_data = read_yaml(file_path="options.yaml")
    sys_path = opt_data.get("wrapper")[1].get("system")
    sys_data = read_yaml(file_path=sys_path)

    lig_data = opt_data.get("wrapper")[2]["ligands"]
    smiles_col = lig_data[1].get("smiles_col")
    id_col = id_col=lig_data[2].get("id_col")
    
    lig_df = read_csv(file_path=lig_data[0].get("lig_csv"),smiles_col=smiles_col,id_col=smiles_col)
    
    # run boltz
    run_dir = set_dir(path=opt_data.get("wrapper")[0].get("run_dir"))
    run_boltz(
        run_dir,
        sys_data,
        smiles_col,
        id_col,
        lig_df,
    )
 