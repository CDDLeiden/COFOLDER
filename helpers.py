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

def convert_affinity_to_ic50(
    df: pd.DataFrame,
    affinity_col: str = 'affinity_pred_value',
    keep_sign: bool = True
) -> pd.DataFrame:
    """
    Convert affinity predictions (log(IC50) in μM) to IC50 (μM) and pIC50 (kcal/mol) from a DataFrame.
    Handles censored values (e.g., '>', '<', '>=', '<=') and can optionally keep the sign in the output.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing affinity predictions.
    affinity_col : str, default='affinity_pred_value'
        Column name for affinity predictions (may contain censoring signs).
    keep_sign : bool, default=True
        Whether to keep the censoring sign in the output DataFrame.

    Returns
    -------
    pd.DataFrame
        DataFrame with added columns: 'IC50_uM', 'pIC50_kcal_per_mol', and optionally 'affinity_sign'.
    """
    parsed = parse_censored_affinity(df[affinity_col], keep_sign=keep_sign)
    df = df.copy()
    df['affinity_value'] = parsed['affinity_value']
    if keep_sign:
        df['affinity_sign'] = parsed['affinity_sign']
    df['IC50_uM'] = 10 ** df['affinity_value']
    df['pIC50_kcal_per_mol'] = (6 - df['affinity_value']) * 1.364
    if keep_sign:
        # Merge sign as prefix if present
        def prefix_sign(val, sign):
            if pd.isnull(sign) or sign is None:
                return str(val) if not pd.isnull(val) else None
            return f"{sign}{val}" if not pd.isnull(val) else None
        df['IC50_uM'] = [prefix_sign(v, s) for v, s in zip(df['IC50_uM'], df['affinity_sign'])]
        df['pIC50_kcal_per_mol'] = [prefix_sign(v, s) for v, s in zip(df['pIC50_kcal_per_mol'], df['affinity_sign'])]
    return df

def calculate_affinity_correlations(
    df: pd.DataFrame,
    pred_col: str,
    exp_col: str,
    sample_size: Optional[int] = None
) -> dict:
    """
    Calculate correlation metrics between predicted and experimental affinities from a DataFrame.

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

    Returns
    -------
    dict
        Dictionary of correlation metrics (R², Pearson, Spearman, Kendall, RMSE, MAE).
    """
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
    csv_path: str,
    pred_col: str,
    exp_col: str,
    sample_size: Optional[int] = None,
    outdir: str = 'figures',
    outname: str = 'affinity_correlation.png'
):
    """
    Plot predicted vs experimental affinities from a single CSV file with jointplot and correlation metrics.
    Saves the plot to the specified directory.

    Parameters
    ----------
    csv_path : str
        Path to the CSV file containing both predicted and experimental values.
    pred_col : str
        Column name for predicted affinity values.
    exp_col : str
        Column name for experimental (expected) affinity values.
    sample_size : int, optional
        If set, randomly sample this many rows for plotting.
    outdir : str, default='figures'
        Directory to save the plot.
    outname : str, default='affinity_correlation.png'
        Filename for the saved plot.
    """
    df = pd.read_csv(csv_path)
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
    if not os.path.exists(outdir):
        os.makedirs(outdir)
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
    plt.savefig(os.path.join(outdir, outname))
    plt.close()
