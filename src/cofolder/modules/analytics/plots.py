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


def _annotation_label(row: pd.Series) -> str | None:
    for column in ("reference_label", "pdb_id", "protein_pdb_id", "ligand_pdb_id", "dataset_name"):
        value = row.get(column)
        if pd.notna(value) and str(value).strip():
            return str(value).strip()
    return None


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
    """Plot bias reference-overlap diagnostics for decision support."""
    plot_df = df.copy()
    if plot_df.empty:
        return []

    plot_df[sequence_col] = pd.to_numeric(plot_df.get(sequence_col), errors="coerce")
    plot_df[ligand_col] = pd.to_numeric(plot_df.get(ligand_col), errors="coerce")
    plot_df = plot_df.dropna(subset=[sequence_col, ligand_col]).copy()
    if plot_df.empty:
        return []

    plot_df[sequence_col] = plot_df[sequence_col].clip(lower=0.0, upper=1.0)
    plot_df[ligand_col] = plot_df[ligand_col].clip(lower=0.0, upper=1.0)
    plot_df["source"] = plot_df.get("source", pd.Series(dtype=object)).fillna("unknown").astype(str)

    fig, ax = plt.subplots(figsize=(7.2, 7.2))
    ax.add_patch(
        Rectangle(
            (0.0, 0.0),
            ligand_threshold,
            sequence_threshold,
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
        group = plot_df[plot_df["source"] == source].copy()
        if group.empty:
            continue
        ax.scatter(
            group[ligand_col],
            group[sequence_col],
            label=display_names.get(source, source),
            c=palette.get(source, "#555555"),
            marker=markers.get(source, "o"),
            edgecolors="white",
            linewidths=0.6,
            s=52,
            alpha=0.86,
            zorder=3,
        )

    ax.axvline(ligand_threshold, color="#4d4d4d", linestyle="--", linewidth=1.2, zorder=2)
    ax.axhline(sequence_threshold, color="#4d4d4d", linestyle="--", linewidth=1.2, zorder=2)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, color="#ececec", linewidth=0.8, zorder=1)
    ax.set_xlabel("Ligand reference overlap (ECFP Tanimoto)")
    ax.set_ylabel("Protein reference overlap (sequence identity, normalized)")
    ax.set_title("Bias reference-overlap diagnostic")
    ax.text(
        0.02,
        0.98,
        "Reference overlap only; not affinity, confidence, or cofolding-success prediction.",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        color="#333333",
        bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none", "pad": 3},
    )

    annotated_labels: set[str] = set()
    ranked = plot_df.sort_values(
        by=[sequence_col, ligand_col, "reference_label"],
        ascending=[False, False, True],
        na_position="last",
    )
    for candidate in (
        ranked.iloc[0] if not ranked.empty else None,
        ranked[ranked[sequence_col] >= 0.8].iloc[0] if (ranked[sequence_col] >= 0.8).any() else None,
    ):
        if candidate is None:
            continue
        label = _annotation_label(candidate)
        if not label or label in annotated_labels:
            continue
        annotated_labels.add(label)
        ax.annotate(
            label,
            (float(candidate[ligand_col]), float(candidate[sequence_col])),
            xytext=(6, 6),
            textcoords="offset points",
            fontsize=8,
            color="#1f1f1f",
            bbox={"facecolor": "white", "alpha": 0.75, "edgecolor": "none", "pad": 2},
        )

    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels, loc="lower right", frameon=True, framealpha=0.95)

    fig.tight_layout()
    output_root = Path(output_dir)
    saved_paths = [
        _save_figure(fig, output_root / f"{file_stem}.png"),
        _save_figure(fig, output_root / f"{file_stem}.pdf"),
    ]
    plt.close(fig)
    return saved_paths
