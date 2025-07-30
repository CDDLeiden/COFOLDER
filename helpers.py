# Script containing general functions
import os
import pandas as pd
import yaml
import numpy as np
from typing import Optional, List

import logging
helpers_logger = logging.getLogger('boltz-tools.helpers')

def set_dir(path):
    if not (os.path.isdir(path)):
        os.makedirs(path)
        helpers_logger.info("Created working dir {0}".format(path))
    # os.chdir(path)

def read_csv(path, columns):
    try:
        df = pd.read_csv(path)  
        helpers_logger.info(f"Read {path} containing {len(df)} entries")

        for col in columns:
            if col and col in df.columns:
                helpers_logger.info(f"Column {col} in {path}")
            else:
                helpers_logger.warning(f"Column {col} not in {path}")

        return df
    
    except FileNotFoundError:
        logging.error(f"File not found: {path}")

def read_yaml(path):
    try:
        with open(path, 'r') as file:
            data = yaml.safe_load(file)
            logging.info(f"{path} loaded successfully. Contents:")
            for key, value in data.items():
                logging.info(f"\t{key}: {value}")
            
            return data

    except FileNotFoundError:
        logging.error(f"File not found: {path}")
    except yaml.YAMLError as e:
        logging.error(f"Error parsing YAML: {e}")

def delete_last_line(file_path):
    """Delete the last line from a file (in-place)."""
    with open(file_path, 'r') as f:
        lines = f.readlines()
    if lines:
        with open(file_path, 'w') as f:
            f.writelines(lines[:-1])

def convert_affinity_to_ic50(
    csv_path: str,
    affinity_col: str = 'affinity_pred_value',
    output_path: Optional[str] = None
) -> pd.DataFrame:
    """
    Convert affinity predictions (log(IC50) in μM) to IC50 (μM) and pIC50 (kcal/mol).
    Adds two new columns: 'IC50_uM' and 'pIC50_kcal_per_mol'.
    Optionally saves the result to a new CSV file.

    Args:
        csv_path (str): Path to the input CSV file with affinity predictions.
        affinity_col (str): Column name for affinity predictions (default: 'affinity_pred_value').
        output_path (Optional[str]): If provided, save the new DataFrame to this path.

    Returns:
        pd.DataFrame: DataFrame with added columns.
    """
    df = pd.read_csv(csv_path)
    if affinity_col not in df.columns:
        raise ValueError(f"Column '{affinity_col}' not found in {csv_path}")
    df['IC50_uM'] = 10 ** df[affinity_col]
    df['pIC50_kcal_per_mol'] = (6 - df[affinity_col]) * 1.364
    if output_path:
        df.to_csv(output_path, index=False)
    return df
