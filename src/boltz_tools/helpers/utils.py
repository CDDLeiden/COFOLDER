# Script containing general functions
import os
import pandas as pd
import yaml
import numpy as np
from typing import Optional, List
from rdkit import Chem
import pickle

from scipy.stats import pearsonr, spearmanr, kendalltau
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error, root_mean_squared_error
import matplotlib.pyplot as plt
import seaborn as sns

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

def parse_censored_affinity(
    affinity_series: pd.Series,
    keep_sign: bool = True
) -> pd.DataFrame:
    """
    Parse affinity values with possible censoring signs (e.g., '>', '<', '>=', '<=') and separate them from the numeric part.

    Parameters
    ----------
    affinity_series : pd.Series
        Series of affinity values, possibly as strings with censoring signs.
    keep_sign : bool, default=True
        Whether to keep the censoring sign in the output DataFrame.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns 'affinity_value' (float) and 'affinity_sign' (str or None).
    """
    import re
    signs = ['>=', '<=', '>', '<']
    def split_sign(val):
        if pd.isnull(val):
            return (None, None)
        val = str(val).strip()
        for s in signs:
            if val.startswith(s):
                try:
                    return (float(val[len(s):].strip()), s)
                except ValueError:
                    return (None, s)
        try:
            return (float(val), None)
        except ValueError:
            return (None, None)
    parsed = affinity_series.apply(split_sign)
    df = pd.DataFrame(parsed.tolist(), columns=['affinity_value', 'affinity_sign'])
    if not keep_sign:
        df['affinity_sign'] = None
    return df

def remove_censored_affinity(
    df: pd.DataFrame,
    cols: list
) -> pd.DataFrame:
    """
    Remove rows where any of the specified columns contain censoring signs (>, <, >=, <=).

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to filter.
    cols : list of str
        List of column names to check for censoring signs.

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame with only rows where all specified columns are numeric.
    """
    import re
    censor_pattern = re.compile(r'^(>=|<=|>|<)')
    mask = pd.Series([True] * len(df))
    for col in cols:
        mask &= ~df[col].astype(str).str.strip().str.match(censor_pattern)
    return df[mask].copy()

def strip_censoring_signs(
    df: pd.DataFrame,
    cols: list
) -> pd.DataFrame:
    """
    Remove censoring signs (>, <, >=, <=) from the start of values in specified columns, converting them to floats.
    Non-numeric values after stripping will be set to NaN.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to process.
    cols : list of str
        List of column names to strip censoring signs from.

    Returns
    -------
    pd.DataFrame
        DataFrame with censoring signs removed and values converted to float in specified columns.
    """
    import re
    df = df.copy()
    censor_pattern = re.compile(r'^(>=|<=|>|<)')
    for col in cols:
        df[col] = df[col].astype(str).str.strip().str.replace(censor_pattern, '', regex=True)
        df[col] = pd.to_numeric(df[col], errors='coerce')
    return df

def prepare_affinity_dataframe(
    df: pd.DataFrame,
    cols: list,
    censoring: str = 'remove'  # options: 'remove', 'strip'
) -> pd.DataFrame:
    """
    Prepare a DataFrame for affinity correlation/plotting by handling censoring signs.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to process.
    cols : list of str
        List of column names to check/clean for censoring signs.
    censoring : {'remove', 'strip'}, default='remove'
        If 'remove', remove rows with censoring signs in any of the columns.
        If 'strip', remove censoring signs and use the numeric part.

    Returns
    -------
    pd.DataFrame
        Cleaned DataFrame ready for numeric analysis.
    """
    if censoring == 'remove':
        return remove_censored_affinity(df, cols)
    elif censoring == 'strip':
        return strip_censoring_signs(df, cols)
    else:
        raise ValueError("censoring must be 'remove' or 'strip'")

def convert_boltz_affinity_to_ic50(
    df: pd.DataFrame,
    affinity_col: str = 'affinity_pred_value',
    output_path: Optional[str] = None
) -> pd.DataFrame:
    """
    Convert affinity predictions (log(IC50) in μM) to IC50 (μM) and pIC50 (kcal/mol) from a DataFrame.
    Only works on numeric values (no censoring signs).
    Optionally saves the updated DataFrame to a CSV file.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing affinity predictions.
    affinity_col : str, default='affinity_pred_value'
        Column name for affinity predictions (must be numeric).
    output_path : str, optional
        If provided, save the updated DataFrame to this path as a CSV file.

    Returns
    -------
    pd.DataFrame
        DataFrame with added columns: 'IC50_uM', 'pIC50_kcal_per_mol'.
    """
    df = df.copy()
    df['affinity_value'] = df[affinity_col].astype(float)
    df['IC50_uM'] = 10 ** df['affinity_value']
    df['pIC50_kcal_per_mol'] = (6 - df['affinity_value']) * 1.364
    if output_path is not None:
        df.to_csv(output_path, index=False)
    return df

def drop_and_log_nans(df: pd.DataFrame, cols: list, context: str = "") -> pd.DataFrame:
    """
    Drop rows with NaN in any of the specified columns and log the number of dropped rows.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to process.
    cols : list of str
        Columns to check for NaN values.
    context : str, optional
        Context string to include in the log message.

    Returns
    -------
    pd.DataFrame
        DataFrame with rows containing NaN in specified columns removed.
    """
    nan_count = df[cols].isna().any(axis=1).sum()
    if nan_count > 0:
        logging.info(f"Removed {nan_count} rows with NaN in {cols} {f'for {context}' if context else ''}.")
    return df.dropna(subset=cols)

def calculate_affinity_correlations(
    df: pd.DataFrame,
    pred_col: str,
    exp_col: str,
    sample_size: Optional[int] = None,
    censoring: str = 'strip'  # options: 'remove', 'strip'
) -> dict:
    """
    Calculate correlation metrics between predicted and experimental affinities from a DataFrame.
    Handles censoring signs according to the 'censoring' argument.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing both predicted and experimental values.
    pred_col : str
        Column name for predicted affinity values.
    exp_col : str
        Column name for experimental (expected) affinity values.
    sample_size : int, optional
        If set, randomly sample this many rows for metrics/plots.
    censoring : {'remove', 'strip'}, default='strip'
        If 'remove', remove rows with censoring signs in either column.
        If 'strip', remove censoring signs and use the numeric part.

    Returns
    -------
    dict
        Dictionary of correlation metrics (R², Pearson, Spearman, Kendall, RMSE, MAE).
    """
    df = prepare_affinity_dataframe(df, [pred_col, exp_col], censoring=censoring)
    df = drop_and_log_nans(df, [pred_col, exp_col], context="correlation calculation")
    if sample_size is not None and sample_size < len(df):
        df = df.sample(n=sample_size, random_state=42)
    x = df[pred_col]
    y = df[exp_col]
    metrics = {
        'r2': float(r2_score(y, x)),
        'pearson': float(pearsonr(x, y)[0]),
        'spearman': float(spearmanr(x, y)[0]),
        'kendall': float(kendalltau(x, y)[0]),
        'rmse': float(root_mean_squared_error(y, x)),
        'mae': float(mean_absolute_error(y, x))
    }
    return metrics

def plot_affinity_correlation(
    df: pd.DataFrame,
    pred_col: str,
    exp_col: str,
    sample_size: Optional[int] = None,
    censoring: str = 'strip',
    output_path: Optional[str] = None
):
    """
    Plot predicted vs experimental affinities from a DataFrame with jointplot and correlation metrics.
    Handles censoring signs according to the 'censoring' argument.
    Saves the plot to the specified output_path or displays it if output_path is None.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing both predicted and experimental values.
    pred_col : str
        Column name for predicted affinity values.
    exp_col : str
        Column name for experimental (expected) affinity values.
    sample_size : int, optional
        If set, randomly sample this many rows for plotting.
    censoring : {'remove', 'strip'}, default='strip'
        If 'remove', remove rows with censoring signs in either column.
        If 'strip', remove censoring signs and use the numeric part.
    output_path : str, optional
        If provided, save the plot to this path. If None, display the plot interactively.
    """
    df = prepare_affinity_dataframe(df, [pred_col, exp_col], censoring=censoring)
    df = drop_and_log_nans(df, [pred_col, exp_col], context="correlation plotting")
    if sample_size is not None and sample_size < len(df):
        df = df.sample(n=sample_size, random_state=42)
    x = df[exp_col]
    y = df[pred_col]
    # Important to note here that we already sampled the DataFrame above,
    # so we don't need to sample again for metrics calculation.
    # Also censoring is handled in the prepare_affinity_dataframe function. so 'remove' or 'strip' is already applied.
    # And there is no need to pass
    metrics = calculate_affinity_correlations(df, pred_col, exp_col, sample_size=None)
    plt.figure(figsize=(7,7))
    g = sns.jointplot(x=x, y=y, kind='scatter', marginal_kws=dict(bins=30, fill=True))
    g.ax_joint.plot([x.min(), x.max()], [x.min(), x.max()], 'r--', alpha=0.5)
    legend = '\n'.join([
        f"R² = {metrics['r2']:.3f}",
        f"Pearson = {metrics['pearson']:.3f}",
        f"Spearman = {metrics['spearman']:.3f}",
        f"Kendall = {metrics['kendall']:.3f}",
        f"RMSE = {metrics['rmse']:.3f}",
        f"MAE = {metrics['mae']:.3f}"
    ])
    g.ax_joint.legend([legend], loc='upper left', fontsize=9, frameon=True)
    g.set_axis_labels('Experimental Affinity', 'Predicted Affinity' )
    plt.tight_layout()
    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        plt.savefig(output_path)
        plt.close()
    else:
        plt.show()

def cache_mols_from_file(file_path: str, id_col_or_prop: str, smiles_col: Optional[str] = None, boltz_cache_path: str = "~/.boltz/"):
    """
    Cache molecule objects from a CSV or SDF file into .pkl files named by compound ID.

    Parameters
    ----------
    file_path : str
        Path to the input file (CSV or SDF).
    id_col_or_prop : str
        Column name for compound ID in CSV or property name for compound ID in SDF.
    smiles_col : str, optional
        Column name for SMILES in CSV file. Required if file is CSV.
    boltz_cache_path : str, default='~/.boltz/'
        Path to the boltz cache directory. Defaults to '~/.boltz/'.

    Raises
    ------
    ValueError
        If required columns/properties are missing or file type is unsupported.
    FileNotFoundError
        If the input file does not exist.
    Exception
        For other errors during molecule creation or file writing.
    """
    file_path = os.path.expanduser(file_path)
    boltz_cache_path = os.path.expanduser(boltz_cache_path)
    mols_dir = os.path.join(boltz_cache_path, "mols")
    os.makedirs(mols_dir, exist_ok=True)

    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"Input file not found: {file_path}")

    ext = os.path.splitext(file_path)[1].lower()
    if ext == ".csv":
        if smiles_col is None:
            raise ValueError("For CSV, smiles_col must be provided.")
        df = pd.read_csv(file_path)
        if smiles_col not in df.columns or id_col_or_prop not in df.columns:
            raise ValueError(f"CSV missing required columns: {smiles_col}, {id_col_or_prop}")
        for idx, row in df.iterrows():
            smiles = row[smiles_col]
            mol_id = str(row[id_col_or_prop])
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    raise ValueError(f"Invalid SMILES: {smiles} (ID: {mol_id})")
                out_path = os.path.join(mols_dir, f"{mol_id}.pkl")
                with open(out_path, "wb") as f:
                    pickle.dump(mol, f)
            except Exception as e:
                helpers_logger.error(f"Failed to process ID {mol_id}: {e}")
    elif ext == ".sdf":
        if smiles_col is not None:
            helpers_logger.warning("smiles_col argument will not be used for SDF files.")
        suppl = Chem.SDMolSupplier(file_path)
        for mol in suppl:
            if mol is None:
                continue
            mol_id = mol.GetProp(id_col_or_prop) if mol.HasProp(id_col_or_prop) else None
            if not mol_id:
                helpers_logger.error(f"SDF molecule missing ID property '{id_col_or_prop}'")
                continue
            try:
                out_path = os.path.join(mols_dir, f"{mol_id}.pkl")
                with open(out_path, "wb") as f:
                    pickle.dump(mol, f)
            except Exception as e:
                helpers_logger.error(f"Failed to process ID {mol_id}: {e}")
    else:
        raise ValueError(f"Unsupported file type: {ext}. Only .csv and .sdf are supported.")
