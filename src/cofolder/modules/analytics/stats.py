import logging

import pandas as pd
from scipy.stats import kendalltau, pearsonr, spearmanr
from sklearn.metrics import mean_absolute_error, r2_score, root_mean_squared_error

from cofolder.modules.analytics import dataset
from cofolder.modules.contracts import convert_metric_value
from cofolder.modules.utils import helpers, write

logger = logging.getLogger(__name__)


def convert_boltz_affinity_to_ic50(
    df: pd.DataFrame,
    affinity_col: str = "affinity_pred_value",
    output_path: str | None = None,
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
    df["affinity_value"] = df[affinity_col].astype(float)
    df["IC50_uM"] = df["affinity_value"].map(
        lambda value: convert_metric_value("affinity_pred_value", "IC50_uM", value)
    )
    df["pIC50_kcal_per_mol"] = df["affinity_value"].map(
        lambda value: convert_metric_value(
            "affinity_pred_value", "pIC50_kcal_per_mol", value
        )
    )
    if output_path is not None:
        write.write_csv(df, output_path, index=False)
    return df


def affinity_to_pic50_and_ic50(affinity_pred_value: float) -> tuple[float, float]:
    """
    Convert Boltz affinity_pred_value to pIC50 and IC50 (M).

    Definitions (Boltz):
    - pIC50 = 6 - affinity_pred_value
    - IC50 (M) = 10 ** (-pIC50)

    Parameters
    ----------
    affinity_pred_value : float
        Raw affinity prediction value.

    Returns
    -------
    tuple[float, float]
        (pIC50, IC50_M)
    """
    pIC50 = convert_metric_value("affinity_pred_value", "pIC50", affinity_pred_value)
    IC50_M = convert_metric_value("affinity_pred_value", "IC50_M", affinity_pred_value)
    return pIC50, IC50_M


def affinity_to_pic50_kcal_per_mol(affinity_pred_value: float) -> float:
    """
    Convert Boltz affinity_pred_value to pIC50-scaled kcal/mol.

    Definition (Boltz):
    - pIC50_kcal_per_mol = (6 - affinity_pred_value) * 1.364

    Parameters
    ----------
    affinity_pred_value : float
        Raw affinity prediction value.

    Returns
    -------
    float
        pIC50-scaled kcal/mol value.
    """
    return convert_metric_value(
        "affinity_pred_value", "pIC50_kcal_per_mol", affinity_pred_value
    )


def calculate_affinity_correlations(
    df: pd.DataFrame,
    pred_col: str,
    exp_col: str,
    sample_size: int | None = None,
    censoring: str = "strip",
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
    df = dataset.prepare_affinity_dataframe(
        df, [pred_col, exp_col], censoring=censoring
    )
    df = helpers.drop_and_log_nans(
        df, [pred_col, exp_col], context="correlation calculation"
    )
    if sample_size is not None and sample_size < len(df):
        df = df.sample(n=sample_size, random_state=42)
    x = df[pred_col]
    y = df[exp_col]
    metrics = {
        "r2": float(r2_score(y, x)),
        "pearson": float(pearsonr(x, y)[0]),
        "spearman": float(spearmanr(x, y)[0]),
        "kendall": float(kendalltau(x, y)[0]),
        "rmse": float(root_mean_squared_error(y, x)),
        "mae": float(mean_absolute_error(y, x)),
    }
    return metrics
