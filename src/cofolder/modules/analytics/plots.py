from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pandas as pd
import seaborn as sns

from cofolder.modules.analytics import dataset, stats
from cofolder.modules.utils import helpers, write

logger = logging.getLogger(__name__)


def plot_affinity_correlation(
    df: pd.DataFrame,
    pred_col: str,
    exp_col: str,
    sample_size: Optional[int] = None,
    censoring: str = "strip",
    output_path: Optional[str] = None,
):
    """Create a correlation plot for predicted vs experimental affinities."""
    df = dataset.prepare_affinity_dataframe(df, [pred_col, exp_col], censoring=censoring)
    df = helpers.drop_and_log_nans(df, [pred_col, exp_col], context="correlation plotting")
    if sample_size is not None and sample_size < len(df):
        df = df.sample(n=sample_size, random_state=42)
    x = df[exp_col]
    y = df[pred_col]
    metrics = stats.calculate_affinity_correlations(df, pred_col, exp_col, sample_size=None)
    plt.figure(figsize=(7, 7))
    g = sns.jointplot(x=x, y=y, kind="scatter", marginal_kws=dict(bins=30, fill=True))
    g.ax_joint.plot([x.min(), x.max()], [x.min(), x.max()], "r--", alpha=0.5)
    legend = "\n".join(
        [
            f"R² = {metrics['r2']:.3f}",
            f"Pearson = {metrics['pearson']:.3f}",
            f"Spearman = {metrics['spearman']:.3f}",
            f"Kendall = {metrics['kendall']:.3f}",
            f"RMSE = {metrics['rmse']:.3f}",
            f"MAE = {metrics['mae']:.3f}",
        ]
    )
    g.ax_joint.legend([legend], loc="upper left", fontsize=9, frameon=True)
    g.set_axis_labels("Experimental Affinity", "Predicted Affinity")
    plt.tight_layout()
    if output_path:
        write.save_plot(output_path)
    else:
        plt.show()


def _save_figure(fig: plt.Figure, output_path: str | Path) -> Path:
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight")
    logger.info("Plot saved to %s", path)
    return path


def plot_reference_overlap_scatter(
    df: pd.DataFrame,
    *,
    output_dir: str | Path,
    file_stem: str,
    x_col: str,
    y_col: str,
    x_threshold: float,
    y_threshold: float,
    x_label: str,
    y_label: str,
    title: str,
    query_1_col: str,
    query_2_col: str,
) -> list[Path]:
    """Plot reference-overlap diagnostics for decision support."""
    plot_df = df.copy()
    if plot_df.empty:
        return []

    plot_df[x_col] = pd.to_numeric(plot_df.get(x_col), errors="coerce").clip(lower=0.0, upper=1.0)
    plot_df[y_col] = pd.to_numeric(plot_df.get(y_col), errors="coerce").clip(lower=0.0, upper=1.0)
    has_query_1 = (
        plot_df.get(query_1_col, pd.Series(dtype=object)).fillna("").astype(str).str.strip() != ""
    )
    has_query_2 = (
        plot_df.get(query_2_col, pd.Series(dtype=object)).fillna("").astype(str).str.strip() != ""
    )
    if not (has_query_1 & has_query_2).any():
        return []
    plot_df["_plot_x"] = plot_df[x_col]
    plot_df.loc[plot_df["_plot_x"].isna() & has_query_1, "_plot_x"] = 0.0
    plot_df["_plot_y"] = plot_df[y_col]
    plot_df.loc[plot_df["_plot_y"].isna() & has_query_2, "_plot_y"] = 0.0
    plot_df = plot_df[
        (plot_df["_plot_x"] >= x_threshold)
        | (plot_df["_plot_y"] >= y_threshold)
    ].copy()
    if plot_df.empty:
        return []

    plot_df["source"] = plot_df.get("source", pd.Series(dtype=object)).fillna("unknown").astype(str)

    fig, ax = plt.subplots(figsize=(7.2, 7.2))
    ax.add_patch(
        Rectangle(
            (0.0, 0.0),
            x_threshold,
            y_threshold,
            facecolor="#d9d9d9",
            edgecolor="none",
            alpha=0.35,
            zorder=0,
        )
    )

    palette = {
        "public": "#1f1f1f",
        "custom": "#d95f02",
        "mixed": "#1b9e77",
        "unknown": "#7570b3",
    }
    markers = {
        "public": "o",
        "custom": "^",
        "mixed": "s",
        "unknown": "D",
    }
    display_names = {
        "public": "Public references",
        "custom": "Custom references",
        "mixed": "Mixed provenance",
        "unknown": "Other provenance",
    }
    source_order = [source for source in ("public", "custom", "mixed", "unknown") if source in set(plot_df["source"])]
    source_order.extend(
        source for source in sorted(set(plot_df["source"])) if source not in set(source_order)
    )

    for source in source_order:
        source_group = plot_df[plot_df["source"] == source].copy()
        if source_group.empty:
            continue
        ax.scatter(
            source_group["_plot_x"],
            source_group["_plot_y"],
            marker=markers.get(source, "o"),
            linewidths=0.6,
            s=52,
            alpha=0.86,
            zorder=3,
            c=palette.get(source, "#555555"),
            edgecolors="white",
        )

    ax.axvline(x_threshold, color="#4d4d4d", linestyle="--", linewidth=1.2, zorder=2)
    ax.axhline(y_threshold, color="#4d4d4d", linestyle="--", linewidth=1.2, zorder=2)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, color="#ececec", linewidth=0.8, zorder=1)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)

    fig.tight_layout()
    output_root = Path(output_dir)
    saved_paths = [
        _save_figure(fig, output_root / f"{file_stem}.png"),
        _save_figure(fig, output_root / f"{file_stem}.pdf"),
    ]
    plt.close(fig)
    return saved_paths


def plot_bias_reference_overlap(
    df: pd.DataFrame,
    *,
    output_dir: str | Path,
    file_stem: str = "bias_reference_overlap_scatter",
    sequence_col: str = "plot_sequence_similarity",
    ligand_col: str = "plot_ecfp_similarity",
    sequence_threshold: float = 0.25,
    ligand_threshold: float = 0.35,
) -> list[Path]:
    return plot_reference_overlap_scatter(
        df,
        output_dir=output_dir,
        file_stem=file_stem,
        x_col=ligand_col,
        y_col=sequence_col,
        x_threshold=ligand_threshold,
        y_threshold=sequence_threshold,
        x_label="Ligand reference overlap (ECFP Tanimoto)",
        y_label="Protein reference overlap (sequence identity, normalized)",
        title="Bias reference-overlap diagnostic",
        query_1_col="query_ligand_chain_id",
        query_2_col="query_protein_chain_id",
    )
