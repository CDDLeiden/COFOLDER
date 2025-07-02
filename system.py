# Script containing YAML functions, such injecting ligands into a system.

import logging
import os
import yaml

# TODO: make read YAML a class

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

def set_yaml(out_dir, sys_data, i, row, id_col, msa_path, smiles_col):
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