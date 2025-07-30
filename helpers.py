# Script containing general functions
import os
import pandas as pd
import yaml
import numpy as np
from typing import Optional, List

from scipy.stats import pearsonr, spearmanr, kendalltau
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
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
    if sample_size is not None and sample_size < len(df):
        df = df.sample(n=sample_size, random_state=42)
    x = df[pred_col]
    y = df[exp_col]
    metrics = {
        'r2': r2_score(y, x),
        'pearson': pearsonr(x, y)[0],
        'spearman': spearmanr(x, y)[0],
        'kendall': kendalltau(x, y)[0],
        'rmse': mean_squared_error(y, x, squared=False),
        'mae': mean_absolute_error(y, x)
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
    if sample_size is not None and sample_size < len(df):
        df = df.sample(n=sample_size, random_state=42)
    x = df[pred_col]
    y = df[exp_col]
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
    g.set_axis_labels('Predicted Affinity', 'Experimental Affinity')
    plt.tight_layout()
    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        plt.savefig(output_path)
        plt.close()
    else:
        plt.show()
