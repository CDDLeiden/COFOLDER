"""Aggregate normalized runner outputs and workflow robustness metrics."""

import logging
import shutil
from collections.abc import Collection
from pathlib import Path

import numpy as np
import pandas as pd

from cofolder.modules.analytics import align
from cofolder.modules.contracts import AmbiguousIdentityError
from cofolder.modules.runners.contracts import (
    RUNNER_PROVENANCE_COLUMNS,
    RunnerChainIdentity,
)
from cofolder.modules.utils import read

# Preserve the established logger identity during the behavior-preserving move.
_DEFAULT_LOGGER = logging.getLogger(__name__)
logger = _DEFAULT_LOGGER


def gather_structures(
    base_dir: Path,
    system_name: str,
    repeats: int,
    logger: logging.Logger | None = None,
    repeat_ids: list[int] | None = None,
):
    """
    Gather all predicted structure files from repeats into a single folder.

    Parameters
    ----------
    base_dir : Path
        Top-level working directory (e.g., wrk_dir/raw).
    system_name : str
        Name of the system used in prediction outputs.
    repeats : int
        Number of repeats performed.
    logger : logging.Logger, optional
        Logger instance for reporting.

    Notes
    -----
    - Looks into: base_dir/repeat_{i}/normalized/structures/
    - Copies all .cif/.mmcif/.pdb files to: base_dir/results/structures/
    """
    if logger is None:
        logger = _DEFAULT_LOGGER

    target_dir = base_dir / "results" / "structures"
    target_dir.mkdir(parents=True, exist_ok=True)
    logger.info("Gathering structures into: %s", target_dir)

    for i in repeat_ids or range(1, repeats + 1):
        repeat_dir = base_dir / "raw" / f"repeat_{i}" / "normalized" / "structures"
        if not repeat_dir.exists():
            logger.warning("Normalized structures directory not found: %s", repeat_dir)
            continue

        for file_path in repeat_dir.iterdir():
            if file_path.suffix.lower() not in (".cif", ".mmcif", ".pdb"):
                continue

            new_name = file_path.name
            dest_path = target_dir / new_name
            try:
                shutil.copy2(file_path, dest_path)
                logger.debug("Copied %s -> %s", file_path, dest_path)
            except Exception as e:
                logger.error("Failed to copy %s: %s", file_path, e)

    logger.info(
        "Structure gathering complete. Total files in %s: %d",
        target_dir,
        len(list(target_dir.iterdir())),
    )


def merge_runner_results(
    raw_dir: Path,
    repeats: int,
    logger: logging.Logger | None = None,
    repeat_ids: list[int] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, list[dict]]:
    """Merge normalized per-repeat runner outputs into canonical DataFrames."""
    if logger is None:
        logger = _DEFAULT_LOGGER

    system_frames: list[pd.DataFrame] = []
    chain_frames: list[pd.DataFrame] = []
    manifests: list[dict] = []

    for repeat in repeat_ids or range(1, repeats + 1):
        normalized_dir = raw_dir / f"repeat_{repeat}" / "normalized"
        system_metrics_path = normalized_dir / "system_metrics.csv"
        chain_metrics_path = normalized_dir / "chain_metrics.csv"
        manifest_path = normalized_dir / "manifest.json"

        if manifest_path.exists():
            manifest = read.read_json(manifest_path)
            if manifest:
                manifests.append(manifest)

        if system_metrics_path.exists():
            system_frames.append(
                pd.read_csv(system_metrics_path, dtype={"effective_seed": "UInt64"})
            )
        else:
            logger.warning("Missing normalized system metrics: %s", system_metrics_path)

        if chain_metrics_path.exists():
            chain_frames.append(
                pd.read_csv(chain_metrics_path, dtype={"effective_seed": "UInt64"})
            )
        else:
            logger.warning("Missing normalized chain metrics: %s", chain_metrics_path)

    system_df = (
        pd.concat(system_frames, ignore_index=True) if system_frames else pd.DataFrame()
    )
    chain_df = (
        pd.concat(chain_frames, ignore_index=True) if chain_frames else pd.DataFrame()
    )

    if not system_df.empty:
        system_df = system_df.reset_index(drop=True)
        system_df.insert(0, "idx", range(len(system_df)))

    if not chain_df.empty:
        chain_df = chain_df.reset_index(drop=True)
        chain_df.insert(0, "idx", range(len(chain_df)))

    return system_df, chain_df, manifests


def add_chain_info(
    chain_df: pd.DataFrame,
    chain_identities: Collection[RunnerChainIdentity],
) -> pd.DataFrame:
    """
    Populate CHAIN_ID, ENTITY_TYPE, and ligand_molecule_id in chain_df
    from an explicit mapping produced during system normalization.

    Mapping is done by normalized entity position: ``conf_chain_id`` maps to the
    ordered system entity and yields both a chain-independent ``ENTITY_ID`` and the
    declared ``CHAIN_ID``.
    """

    identity_map = {item.conf_chain_id: item for item in chain_identities}
    if len(identity_map) != len(chain_identities):
        raise AmbiguousIdentityError(
            "Runner chain identity mapping contains duplicate conf_chain_id values."
        )
    if not identity_map and not chain_df.empty:
        raise AmbiguousIdentityError(
            "Runner chain records cannot be enriched without a chain identity mapping."
        )

    # --------------------------------------------------
    # Assign CHAIN_ID, ENTITY_TYPE, ligand_molecule_id
    # --------------------------------------------------
    chain_df["CHAIN_ID"] = None
    chain_df["ENTITY_ID"] = None
    chain_df["ENTITY_TYPE"] = None
    chain_df["ligand_molecule_id"] = None

    for idx, row in chain_df.iterrows():
        try:
            conf_id = int(row["conf_chain_id"])
        except (TypeError, ValueError) as exc:
            raise AmbiguousIdentityError(
                f"Runner chain record at row {idx} has invalid conf_chain_id "
                f"{row.get('conf_chain_id')!r}."
            ) from exc

        identity = identity_map.get(conf_id)
        if identity is None:
            raise AmbiguousIdentityError(
                f"Runner conf_chain_id {conf_id} at row {idx} cannot be mapped to "
                "a normalized system entity."
            )

        chain_df.at[idx, "CHAIN_ID"] = identity.chain_id
        chain_df.at[idx, "ENTITY_ID"] = identity.entity_id
        chain_df.at[idx, "ENTITY_TYPE"] = identity.entity_type
        chain_df.at[idx, "ligand_molecule_id"] = identity.ligand_molecule_id

    # --------------------------------------------------
    # Reorder columns
    # --------------------------------------------------
    cols = chain_df.columns.tolist()
    for col in ("CHAIN_ID", "ENTITY_ID", "ENTITY_TYPE", "ligand_molecule_id"):
        if col in cols:
            cols.remove(col)

    insert_at = cols.index("conf_chain_id")
    cols.insert(insert_at + 1, "CHAIN_ID")
    cols.insert(insert_at + 2, "ENTITY_ID")
    cols.insert(insert_at + 3, "ENTITY_TYPE")
    cols.insert(insert_at + 4, "ligand_molecule_id")

    return chain_df[cols]


def assess_numeric_variance(
    values: list[float],
    prefix: str,
) -> dict[str, float | None]:
    """
    Assess robustness of numeric values across repeats / diffusion samples.

    Parameters
    ----------
    values : list[float]
        Numeric values across runs (NaN/None allowed).
    prefix : str
        Prefix for output metric names.

    Returns
    -------
    dict
        Robustness statistics (mean, std).
    """
    s = pd.Series(values, dtype="float64").dropna()

    # No data
    if len(s) == 0:
        return {
            f"{prefix}_mean": None,
            f"{prefix}_std": None,
        }

    # Single value → no variance
    if len(s) == 1:
        val = float(s.iloc[0])
        return {
            f"{prefix}_mean": val,
            f"{prefix}_std": 0.0,
        }

    mean = float(s.mean())
    std = float(s.std(ddof=1))

    return {
        f"{prefix}_mean": mean,
        f"{prefix}_std": std,
    }


def assess_bitstring_similarity(
    bitstrings: list,
    prefix: str,
    id: str,
    wrk_dir: Path | str,
) -> dict[str, float | None]:
    """
    Assess robustness of binary fingerprints (e.g. IFPs) across runs
    using pairwise Tanimoto similarity.

    Saves a full NxN similarity matrix to
    Path(wrk_dir)/results/matrices/similarity_matrix_{prefix}.csv

    Parameters
    ----------
    bitstrings : list
        One fingerprint per run.
    prefix : str
        Prefix for output metric names.
    wrk_dir : Path or str
        Root working directory to save the matrix.

    Returns
    -------
    dict
        Mean and std of pairwise similarities (flattened off-diagonal).
    """

    def to_set(x):
        if x is None:
            return set()
        if isinstance(x, set):
            return x
        if isinstance(x, str):
            return {i for i, v in enumerate(x) if v == "1"}
        # assume iterable of bool/int
        return {i for i, v in enumerate(x) if bool(v)}

    sets = [to_set(x) for x in bitstrings]

    n = len(sets)
    if n < 2:
        return {
            f"{prefix}_mean": None,
            f"{prefix}_std": None,
        }

    # Compute NxN similarity matrix
    mat = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i, n):
            a, b = sets[i], sets[j]

            # both empty → identical
            if not a and not b:
                val = 1.0
            else:
                val = len(a & b) / len(a | b)
            mat[i, j] = val
            mat[j, i] = val

    # --- save full NxN similarity matrix ---
    matrices_dir = Path(wrk_dir) / "results" / "matrices"
    matrices_dir.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(
        mat,
        index=[f"run_{i}" for i in range(n)],
        columns=[f"run_{i}" for i in range(n)],
    )
    df.to_csv(matrices_dir / f"similarity_matrix_{prefix}_{id}.csv")

    # Flatten off-diagonal for summary statistics
    off_diag = mat[np.triu_indices(n, k=1)]
    pairwise_std = float(off_diag.std(ddof=1)) if len(off_diag) > 1 else None
    return {
        f"{prefix}_mean": float(off_diag.mean()),
        f"{prefix}_std": pairwise_std,
    }


def gather_robustness_results(
    system_df: pd.DataFrame,
    chain_df: pd.DataFrame,
    wrk_dir,
    reference_path: Path | str | None = None,
) -> pd.DataFrame:
    """
    Aggregate robustness metrics across repeats and diffusion samples.

    Adds structural robustness:
        - polymer chain RMSD (aligned on protein CA)
        - ligand pose RMSD (after global alignment)

    If reference_path is provided, all predicted structures are aligned to
    that reference structure before RMSD calculations.
    """

    rows: list[dict] = []

    logger.debug(
        "[ROBUSTNESS] Starting robustness aggregation "
        "(system_df columns=%d, chain_df columns=%d)",
        len(system_df.columns),
        len(chain_df.columns),
    )

    # ==============================================================
    # STRUCTURAL ALIGNMENT + RMSD PRECOMPUTE
    # ==============================================================

    cif_folder = Path(wrk_dir) / "results" / "structures"
    aligned_folder = Path(wrk_dir) / "results" / "structures_aligned"

    unique_cifs = chain_df["cif_file"].unique()
    cif_paths = [cif_folder / f for f in unique_cifs]
    if reference_path is not None:
        cif_paths = [Path(reference_path)] + cif_paths

    structures = align._load_structures(cif_paths)
    aligned_structs, ref_struct, ref_name = align._align_structures_on_protein_ca(
        structures, save_dir=aligned_folder
    )

    chain_rmsd_map = align._compute_chain_rmsd(chain_df, aligned_structs, wrk_dir)
    ligand_rmsd_map = align._compute_ligand_rmsd(chain_df, aligned_structs, wrk_dir)

    # ==============================================================
    # SYSTEM-LEVEL ROBUSTNESS
    # ==============================================================
    system_row = {
        "ENTITY_TYPE": "system",
        "ENTITY_ID": "system",
    }

    for col in system_df.columns:
        if is_metadata_column(col):
            continue

        series = system_df[col]

        if is_numeric_metric(series):
            stats = assess_numeric_variance(series.tolist(), prefix=col)
            system_row.update(stats)

    rows.append(system_row)

    # ==============================================================
    # CHAIN-LEVEL ROBUSTNESS (POLYMERS)
    # ==============================================================

    non_ligand_df = chain_df[chain_df["ENTITY_TYPE"] != "ligand"]

    for (chain_id, entity_type), group in non_ligand_df.groupby(
        ["CHAIN_ID", "ENTITY_TYPE"]
    ):
        chain_row = {
            "ENTITY_TYPE": entity_type,
            "ENTITY_ID": chain_id,
        }

        for col in group.columns:
            if is_metadata_column(col):
                continue

            series = group[col]

            if is_numeric_metric(series):
                stats = assess_numeric_variance(series.tolist(), prefix=col)
                chain_row.update(stats)

        # --- structural RMSD metric ---
        rmsd_vals = chain_rmsd_map.get((chain_id, entity_type))
        if rmsd_vals:
            stats = assess_numeric_variance(rmsd_vals, prefix="struct_rmsd")
            chain_row.update(stats)

        rows.append(chain_row)

    # --------------------------------------------------------------
    # LIGANDS — AGGREGATE BY CHAIN_ID TO PRESERVE PUBLIC IDENTITY
    # --------------------------------------------------------------

    ligand_df = chain_df[chain_df["ENTITY_TYPE"] == "ligand"]

    for chain_id, group in ligand_df.groupby("CHAIN_ID"):
        ligand_row = {
            "ENTITY_TYPE": "ligand",
            "ENTITY_ID": chain_id,
        }

        for col in group.columns:
            if is_metadata_column(col):
                continue

            series = group[col]

            if is_numeric_metric(series):
                stats = assess_numeric_variance(series.tolist(), prefix=col)
                ligand_row.update(stats)

            elif is_ifp_metric(col, series=series):
                stats = assess_bitstring_similarity(
                    series.tolist(), prefix=col, id=chain_id, wrk_dir=wrk_dir
                )
                ligand_row.update(stats)

        # --- structural RMSD metric ---
        molecule_ids = group["ligand_molecule_id"].dropna().astype(str).unique()
        rmsd_vals = (
            ligand_rmsd_map.get(molecule_ids[0]) if len(molecule_ids) == 1 else None
        )
        if rmsd_vals:
            stats = assess_numeric_variance(rmsd_vals, prefix="struct_rmsd")
            ligand_row.update(stats)

        rows.append(ligand_row)

    results_df = pd.DataFrame(rows)

    logger.debug(
        "[ROBUSTNESS] Finished robustness aggregation (rows=%d, columns=%d)",
        results_df.shape[0],
        results_df.shape[1],
    )

    return results_df


def is_metadata_column(column_name: str) -> bool:
    """
    Columns that should never be treated as metrics.
    """
    metadata_cols = {
        "idx",
        "model_name",
        "repeat",
        "diffusion_sample",
        "cif_file",
        "CHAIN_ID",
        "ENTITY_ID",
        "ENTITY_TYPE",
        "conf_chain_id",
        "ligand_molecule_id",
    } | set(RUNNER_PROVENANCE_COLUMNS)

    result = column_name in metadata_cols

    if result:
        logger.debug(
            "[ROBUSTNESS] Skipping metadata column '%s'",
            column_name,
        )

    return result


def is_numeric_metric(series: pd.Series) -> bool:
    """
    Determine whether a Series represents numeric values, even if stored as strings.

    A series is considered numeric if at least one non-null value can be
    safely converted to float, and all non-null values are convertible.
    """
    if series.empty:
        return False

    # Drop NaNs and empty strings
    s = series.dropna()
    s = s[s.astype(str).str.strip() != ""]

    if len(s) == 0:
        return False

    try:
        # Try vectorized conversion
        converted = pd.to_numeric(s, errors="raise")
    except Exception:
        return False

    # At least one real numeric value must exist
    return np.isfinite(converted).any()


def is_ifp_metric(column_name: str, series: pd.Series | None = None) -> bool:
    """
    Determine whether a column represents a binary interaction fingerprint (IFP).

    Returns True only if the column name suggests an IFP **and** the column has any non-empty values.
    """
    name = column_name.lower()
    result = name.startswith("ifp") or "fingerprint" in name

    # If a series is provided, check it is not all empty/NaN
    if series is not None and result:
        if series.dropna().empty:
            result = False

    logger.debug(
        "[ROBUSTNESS] Metric '%s' classified as IFP: %s",
        column_name,
        result,
    )

    return result
