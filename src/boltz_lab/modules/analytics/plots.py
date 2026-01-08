import pandas as pd
from typing import Optional

import matplotlib.pyplot as plt
import seaborn as sns

from boltz_lab.modules.analytics import dataset, stats
from boltz_lab.modules.utils import helpers, write

import logging

logger = logging.getLogger(__name__)

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
    df = dataset.prepare_affinity_dataframe(df, [pred_col, exp_col], censoring=censoring)
    df = helpers.drop_and_log_nans(df, [pred_col, exp_col], context="correlation plotting")
    if sample_size is not None and sample_size < len(df):
        df = df.sample(n=sample_size, random_state=42)
    x = df[exp_col]
    y = df[pred_col]
    # Important to note here that we already sampled the DataFrame above,
    # so we don't need to sample again for metrics calculation.
    # Also censoring is handled in the prepare_affinity_dataframe function. so 'remove' or 'strip' is already applied.
    # And there is no need to pass
    metrics = stats.calculate_affinity_correlations(df, pred_col, exp_col, sample_size=None)
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
        write.save_plot(output_path)
    else:
        plt.show()
