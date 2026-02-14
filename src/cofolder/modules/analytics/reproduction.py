"""Scaffold for model reproduction and bias metrics.

This module defines the output schema for model reproduction and bias
assessment. Metric computation is intentionally not implemented yet.
"""

from pathlib import Path
import logging

SYSTEM_REPRODUCTION_COLUMNS = (
    "ligand_rmsd_ref_mean",
    "protein_rmsd_ref_mean",
    "ifp_similarity_ref_mean",
    "pocket_coverage_ref",
    "ligand_pose_overlap_ref",
    "protein_sequence_identity_train_max",
    "ligand_fingerprint_similarity_train_max",
)

CHAIN_REPRODUCTION_COLUMNS = (
    "ligand_rmsd_ref",
    "protein_rmsd_ref",
    "ifp_similarity_ref",
    "pocket_coverage_ref",
    "ligand_pose_overlap_ref",
    "protein_sequence_identity_train",
    "ligand_fingerprint_similarity_train",
)


def scaffold_reproduction_metrics(system_df, chain_df, reference_path: Path | None, logger: logging.Logger | None = None):
    """Add placeholder columns for reproduction/bias metrics.

    Parameters
    ----------
    system_df : pandas.DataFrame
        System-level metrics table.
    chain_df : pandas.DataFrame
        Chain-level metrics table.
    reference_path : Path | None
        Optional path to reference structure for reproduction metrics.
    logger : logging.Logger | None
        Optional logger.

    Returns
    -------
    tuple[pandas.DataFrame, pandas.DataFrame]
        Updated system and chain tables with schema columns present.
    """
    for column in SYSTEM_REPRODUCTION_COLUMNS:
        if column not in system_df.columns:
            system_df[column] = None

    for column in CHAIN_REPRODUCTION_COLUMNS:
        if column not in chain_df.columns:
            chain_df[column] = None

    if logger is not None:
        if reference_path is None:
            logger.info(
                "Reproduction schema initialized without reference structure. "
                "Metrics are placeholders and not computed."
            )
        else:
            logger.info(
                "Reproduction schema initialized with reference structure at %s. "
                "Metric computation scaffold only.",
                reference_path
            )

    return system_df, chain_df
