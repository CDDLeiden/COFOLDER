# Script containing general functions
import os
import pandas as pd
import yaml
import random
from rdkit import Chem
from typing import Optional, List


from scipy.stats import pearsonr, spearmanr, kendalltau
from sklearn.metrics import r2_score, mean_absolute_error, root_mean_squared_error
import matplotlib.pyplot as plt
import seaborn as sns

import logging

logger = logging.getLogger(__name__)

def create_dir(path: str):
    """Ensure that a directory exists, creating it if necessary.

    Parameters
    ----------
    path : str
        The directory path to check and/or create.

    Examples
    --------
    >>> create_dir("./output")
    >>> create_dir("/tmp/my_project/results")
    """
    logger.debug(f"Checking directory: '{path}'")

    if not os.path.isdir(path):
        logger.debug(f"Directory does not exist. Creating: '{path}'")
        os.makedirs(path, exist_ok=True)
        logger.info(f"Created working directory '{path}'")
    else:
        logger.debug(f"Directory already exists: '{path}'")

def parse_list_as_str(
    list_as_str: str,
    separator: str = ",",
    item_type=str,
    expected_length: int = None
):
    """Convert a separated string into a typed Python list.

    Parses a command-line string value containing delimited items and
    converts each item to the specified type. Useful for parsing CLI
    arguments like device lists, seeds, or other numeric sequences.

    Parameters
    ----------
    list_as_str : str
        The raw string containing a list (e.g., "1,2,3").
    separator : str, default=","
        Delimiter used to separate items in the string.
    item_type : type, default=str
        The expected type of each element (str, int, float, etc.).
    expected_length : int, optional
        If provided, the parsed list must match this exact length.

    Returns
    -------
    list
        The processed and type-casted list.

    Raises
    ------
    ValueError
        If type conversion fails or if the parsed list length doesn't
        match expected_length.

    Examples
    --------
    >>> parse_list_as_str("1,2,3", item_type=int)
    [1, 2, 3]

    >>> parse_list_as_str("a;b;c", separator=";")
    ['a', 'b', 'c']

    >>> parse_list_as_str("1.5,2.5", item_type=float, expected_length=2)
    [1.5, 2.5]
    """
    logger.debug(f"Parsing list: raw='{list_as_str}', "
                 f"separator='{separator}', item_type={item_type.__name__}, "
                 f"expected_length={expected_length}")

    # Split into items
    items = [item.strip() for item in list_as_str.split(separator) if item.strip()]

    # Convert each item to specified type
    try:
        typed_items = [item_type(item) for item in items]
    except Exception:
        raise ValueError(
            f"Failed to convert items to type {item_type.__name__}: {items}"
        )

    # Check list length (if included)
    if expected_length is not None and len(typed_items) != expected_length:
        raise ValueError(
            f"List length must be {expected_length}, but got {len(typed_items)}"
        )
    
    logger.debug(f"Finished parsing. Result: {typed_items}")

    return typed_items



def set_dir(path):
    """Create a directory if it doesn't exist.

    Parameters
    ----------
    path : str
        The directory path to create.

    Notes
    -----
    This is a legacy function. Consider using `create_dir()` instead,
    which has more detailed logging and type hints.

    Examples
    --------
    >>> set_dir("./experiments")
    """
    if not (os.path.isdir(path)):
        os.makedirs(path)
        logger.info("Created working dir {0}".format(path))

def read_csv(path, columns):
    """Read a CSV file and validate expected columns.

    Loads a CSV file into a pandas DataFrame and checks for the presence
    of specified columns, logging warnings for any missing columns.

    Parameters
    ----------
    path : str
        Path to the CSV file to read.
    columns : list of str
        List of column names to validate in the CSV file.

    Returns
    -------
    pd.DataFrame or None
        DataFrame containing the CSV data, or None if file not found.

    Examples
    --------
    >>> df = read_csv("compounds.csv", ["smiles", "id", "affinity"])
    >>> df = read_csv("results.csv", ["name", "value"])
    """
    try:
        df = pd.read_csv(path)
        logger.info(f"Read {path} containing {len(df)} entries")

        for col in columns:
            if col and col in df.columns:
                logger.info(f"Column {col} in {path}")
            else:
                logger.warning(f"Column {col} not in {path}")

        return df

    except FileNotFoundError:
        logging.error(f"File not found: {path}")

def read_yaml(path):
    """Read and parse a YAML file.

    Loads a YAML file and logs its contents for debugging purposes.

    Parameters
    ----------
    path : str
        Path to the YAML file to read.

    Returns
    -------
    dict or None
        Dictionary containing the parsed YAML data, or None if an error occurs.

    Raises
    ------
    FileNotFoundError
        If the specified file does not exist.
    yaml.YAMLError
        If the YAML file is malformed or cannot be parsed.

    Examples
    --------
    >>> config = read_yaml("config.yaml")
    >>> system = read_yaml("system.yaml")
    """
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

def read_sdf(path):
    """Read molecules from an SDF file.

    Loads molecules from an SDF (Structure Data File) using RDKit,
    filtering out invalid molecules.

    Parameters
    ----------
    path : str
        Path to the SDF file to read.

    Returns
    -------
    list of rdkit.Chem.Mol or None
        List of valid RDKit molecule objects, or None if an error occurs.

    Examples
    --------
    >>> mols = read_sdf("ligands.sdf")
    >>> mols = read_sdf("compounds.sdf")

    Notes
    -----
    Invalid molecules (those that RDKit cannot parse) are automatically
    filtered out and not included in the returned list.
    """
    try:
        suppl = Chem.SDMolSupplier(path)
        mols = [mol for mol in suppl if mol is not None]

        logging.info(f"{path} loaded successfully. Total valid molecules: {len(mols)}")

        return mols

    except FileNotFoundError:
        logging.error(f"File not found: {path}")
    except Exception as e:
        logging.error(f"Error reading SDF: {e}")

def delete_last_line(file_path):
    """Delete the last line from a file in-place.

    Reads the entire file, removes the last line, and writes the
    modified content back to the same file.

    Parameters
    ----------
    file_path : str
        Path to the file to modify.

    Examples
    --------
    >>> delete_last_line("output.txt")
    >>> delete_last_line("results.log")

    Notes
    -----
    This operation modifies the file in-place. If the file is empty,
    no action is taken. For large files, consider using alternative
    methods to avoid loading the entire file into memory.
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()
    if lines:
        with open(file_path, 'w') as f:
            f.writelines(lines[:-1])

def generate_seeds(num_seeds: int, seed: Optional[int] = None) -> List[int]:
    """
    Generate a list of random integer seeds for reproducibility.

    Parameters
    ----------
    num_seeds : int
        Number of seeds to generate (must be >= 1).
    seed : int or None
        Optional global seed for reproducibility. If None, a random global
        seed is generated.

    Returns
    -------
    List[int]
        List of generated seeds.
    """
    if num_seeds < 1:
        raise ValueError(f"num_seeds must be >= 1, got {num_seeds}")

    # Generate or use the global seed
    global_seed = seed if seed is not None else random.randint(0, 2**32 - 1)
    rng = random.Random(global_seed)
    
    return [rng.randint(0, 2**32 - 1) for _ in range(num_seeds)]

def parse_censored_affinity(
    affinity_series: pd.Series,
    keep_sign: bool = True
) -> pd.DataFrame:
    """Parse affinity values with censoring signs into numeric and sign components.

    Separates censoring indicators (e.g., '>', '<', '>=', '<=') from the
    numeric affinity values, allowing for downstream analysis of censored data.

    Parameters
    ----------
    affinity_series : pd.Series
        Series of affinity values, possibly as strings with censoring signs
        (e.g., '>5.0', '<=3.2').
    keep_sign : bool, default=True
        Whether to retain the censoring sign in the output DataFrame.
        If False, the 'affinity_sign' column will be set to None.

    Returns
    -------
    pd.DataFrame
        DataFrame with two columns:
        - 'affinity_value' : float - The numeric part of the affinity
        - 'affinity_sign' : str or None - The censoring sign if present

    Examples
    --------
    >>> import pandas as pd
    >>> data = pd.Series(['>5.0', '3.2', '<=2.5'])
    >>> parse_censored_affinity(data)
       affinity_value affinity_sign
    0             5.0             >
    1             3.2          None
    2             2.5            <=
    """
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
    """Remove rows containing censoring signs in specified columns.

    Filters out rows where any of the specified columns contain censoring
    indicators (>, <, >=, <=), retaining only rows with purely numeric values.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to filter.
    cols : list of str
        List of column names to check for censoring signs.

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame containing only rows where all specified columns
        have numeric values without censoring signs.

    Examples
    --------
    >>> import pandas as pd
    >>> df = pd.DataFrame({
    ...     'exp': ['5.0', '>3.0', '2.5'],
    ...     'pred': ['4.8', '3.2', '2.3']
    ... })
    >>> remove_censored_affinity(df, ['exp'])
          exp pred
    0     5.0  4.8
    2     2.5  2.3
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
    """Remove censoring signs from values and convert to floats.

    Strips censoring indicators (>, <, >=, <=) from the beginning of values
    in specified columns and converts them to numeric floats. Non-numeric
    values after stripping are set to NaN.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to process.
    cols : list of str
        List of column names to strip censoring signs from.

    Returns
    -------
    pd.DataFrame
        DataFrame with censoring signs removed and values converted to float
        in the specified columns.

    Examples
    --------
    >>> import pandas as pd
    >>> df = pd.DataFrame({
    ...     'exp': ['>5.0', '3.0', '<=2.5'],
    ...     'pred': ['4.8', '3.2', '2.3']
    ... })
    >>> strip_censoring_signs(df, ['exp'])
          exp pred
    0     5.0  4.8
    1     3.0  3.2
    2     2.5  2.3
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
    censoring: str = 'remove'
) -> pd.DataFrame:
    """Prepare a DataFrame for affinity analysis by handling censoring signs.

    Processes a DataFrame to handle censored affinity data, either by
    removing censored rows entirely or by stripping the censoring signs
    and retaining the numeric values.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to process.
    cols : list of str
        List of column names to check and clean for censoring signs.
    censoring : {'remove', 'strip'}, default='remove'
        Strategy for handling censoring signs:
        - 'remove': Remove rows with censoring signs in any specified column
        - 'strip': Remove censoring signs and use the numeric part

    Returns
    -------
    pd.DataFrame
        Cleaned DataFrame ready for numeric analysis.

    Raises
    ------
    ValueError
        If censoring parameter is not 'remove' or 'strip'.

    Examples
    --------
    >>> import pandas as pd
    >>> df = pd.DataFrame({'exp': ['>5.0', '3.0'], 'pred': ['4.8', '3.2']})
    >>> prepare_affinity_dataframe(df, ['exp'], censoring='remove')
          exp pred
    1     3.0  3.2

    >>> prepare_affinity_dataframe(df, ['exp'], censoring='strip')
          exp pred
    0     5.0  4.8
    1     3.0  3.2
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
    """Convert Boltz affinity predictions to IC50 and pIC50 values.

    Transforms affinity predictions from log(IC50) in μM to IC50 (μM) and
    pIC50 (kcal/mol). Optionally saves the enriched DataFrame to a CSV file.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing affinity predictions.
    affinity_col : str, default='affinity_pred_value'
        Column name for affinity predictions. Must contain numeric values
        (no censoring signs).
    output_path : str, optional
        If provided, save the updated DataFrame to this path as a CSV file.

    Returns
    -------
    pd.DataFrame
        DataFrame with three additional columns:
        - 'affinity_value' : float - Copy of the input affinity values
        - 'IC50_uM' : float - IC50 in micromolar (10^affinity_value)
        - 'pIC50_kcal_per_mol' : float - pIC50 in kcal/mol

    Examples
    --------
    >>> import pandas as pd
    >>> df = pd.DataFrame({'affinity_pred_value': [6.0, 7.0, 5.5]})
    >>> result = convert_boltz_affinity_to_ic50(df)
    >>> result[['IC50_uM', 'pIC50_kcal_per_mol']]
       IC50_uM  pIC50_kcal_per_mol
    0      1.0                0.00
    1     10.0               -1.36
    2      3.16               0.68

    Notes
    -----
    The conversion formula for pIC50 is: (6 - affinity_value) * 1.364
    """
    df = df.copy()
    df['affinity_value'] = df[affinity_col].astype(float)
    df['IC50_uM'] = 10 ** df['affinity_value']
    df['pIC50_kcal_per_mol'] = (6 - df['affinity_value']) * 1.364
    if output_path is not None:
        df.to_csv(output_path, index=False)
    return df

def drop_and_log_nans(df: pd.DataFrame, cols: list, context: str = "") -> pd.DataFrame:
    """Drop rows with NaN values and log the number of rows removed.

    Removes rows containing NaN values in any of the specified columns
    and logs the count of removed rows for debugging and quality control.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to process.
    cols : list of str
        Column names to check for NaN values.
    context : str, default=""
        Optional context string to include in the log message for clarity.

    Returns
    -------
    pd.DataFrame
        DataFrame with rows containing NaN in specified columns removed.

    Examples
    --------
    >>> import pandas as pd
    >>> import numpy as np
    >>> df = pd.DataFrame({'a': [1, 2, np.nan], 'b': [4, np.nan, 6]})
    >>> drop_and_log_nans(df, ['a', 'b'], context="correlation analysis")
         a    b
    0  1.0  4.0
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
    censoring: str = 'strip'
) -> dict:
    """Calculate correlation metrics between predicted and experimental affinities.

    Computes multiple statistical metrics to assess the agreement between
    predicted and experimental affinity values, handling censored data
    according to the specified strategy.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing both predicted and experimental values.
    pred_col : str
        Column name for predicted affinity values.
    exp_col : str
        Column name for experimental (expected) affinity values.
    sample_size : int, optional
        If specified, randomly sample this many rows before calculating
        metrics (useful for large datasets).
    censoring : {'remove', 'strip'}, default='strip'
        Strategy for handling censoring signs:
        - 'remove': Remove rows with censoring signs in either column
        - 'strip': Remove censoring signs and use the numeric part

    Returns
    -------
    dict
        Dictionary containing the following correlation metrics:
        - 'r2' : float - Coefficient of determination (R²)
        - 'pearson' : float - Pearson correlation coefficient
        - 'spearman' : float - Spearman rank correlation coefficient
        - 'kendall' : float - Kendall tau correlation coefficient
        - 'rmse' : float - Root mean squared error
        - 'mae' : float - Mean absolute error

    Examples
    --------
    >>> import pandas as pd
    >>> df = pd.DataFrame({
    ...     'pred': [5.0, 6.0, 7.0],
    ...     'exp': [5.2, 5.8, 7.1]
    ... })
    >>> metrics = calculate_affinity_correlations(df, 'pred', 'exp')
    >>> print(f"R² = {metrics['r2']:.3f}")
    R² = 0.982
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
    """Create a correlation plot for predicted vs experimental affinities.

    Generates a seaborn jointplot showing the relationship between predicted
    and experimental affinity values, including marginal distributions and
    correlation metrics. Handles censored data according to the specified
    strategy.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing both predicted and experimental values.
    pred_col : str
        Column name for predicted affinity values.
    exp_col : str
        Column name for experimental (expected) affinity values.
    sample_size : int, optional
        If specified, randomly sample this many rows before plotting
        (useful for large datasets to improve performance).
    censoring : {'remove', 'strip'}, default='strip'
        Strategy for handling censoring signs:
        - 'remove': Remove rows with censoring signs in either column
        - 'strip': Remove censoring signs and use the numeric part
    output_path : str, optional
        If provided, save the plot to this path. If None, display the plot
        interactively using plt.show().

    Examples
    --------
    >>> import pandas as pd
    >>> df = pd.DataFrame({
    ...     'pred': [5.0, 6.0, 7.0, 5.5],
    ...     'exp': [5.2, 5.8, 7.1, 5.4]
    ... })
    >>> plot_affinity_correlation(df, 'pred', 'exp', output_path='correlation.png')

    >>> # For large datasets, use sampling
    >>> plot_affinity_correlation(df, 'pred', 'exp', sample_size=1000)

    Notes
    -----
    The plot includes:
    - Scatter plot with marginal distributions
    - Identity line (y=x) in red dashed
    - Legend with R², Pearson, Spearman, Kendall, RMSE, and MAE metrics
    - Axes labeled as 'Experimental Affinity' and 'Predicted Affinity'
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
