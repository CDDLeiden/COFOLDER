import shutil
from pathlib import Path
import logging
import pandas as pd
import numpy as np

from cofolder.modules.utils import read, write
from cofolder.modules.analytics import stats, structure

logger = logging.getLogger(__name__)

def gather_structures(base_dir: Path, system_name: str, repeats: int, logger: logging.Logger = None):
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
    - Looks into: base_dir/repeat_{i}/predictions/{system_name}/
    - Copies all .cif and .pdb files to: base_dir/../results/structures/
    - Renames files as {index}_{system_name}_model_{idx}.cif
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    target_dir = base_dir / "results" / "structures"
    target_dir.mkdir(parents=True, exist_ok=True)
    logger.info("Gathering structures into: %s", target_dir)

    for i in range(1, repeats + 1):
        repeat_dir = base_dir / "raw" / f"repeat_{i}" \
            / f"boltz_results_{system_name}" \
            / "predictions" / system_name
        if not repeat_dir.exists():
            logger.warning("Predictions directory not found: %s", repeat_dir)
            continue

        for file_path in repeat_dir.iterdir():
            if file_path.suffix.lower() not in (".cif", ".pdb"):
                continue

            # Rename file as {i}_{system_name}_model_{idx}.cif
            # Assuming the original file name contains model info
            new_name = f"{i}_{file_path.name}"
            dest_path = target_dir / new_name
            try:
                shutil.copy2(file_path, dest_path)
                logger.debug("Copied %s -> %s", file_path, dest_path)
            except Exception as e:
                logger.error("Failed to copy %s: %s", file_path, e)

    logger.info("Structure gathering complete. Total files in %s: %d",
                target_dir, len(list(target_dir.iterdir())))

def initialize_results(
    raw_dir: Path,
    system_name: str,
    repeats: int,
    diffusion_samples: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Initialize empty system-level and chain-level result DataFrames.

    Notes
    -----
    - cif_id is intentionally NOT used.
    - conf_chain_id corresponds to the index used in confidence JSON files.
    - CHAIN_ID will be populated later from the system definition.
    """

    # -------------------------
    # System-level DataFrame
    # -------------------------
    system_rows = []
    idx = 0

    for repeat in range(1, repeats + 1):
        for sample in range(diffusion_samples):
            system_rows.append(
                {
                    "idx": idx,
                    "cif_file": f'{repeat}_{system_name}_model_{sample}.cif',
                    "model_name": system_name,
                    "repeat": repeat,
                    "diffusion_sample": sample,
                }
            )
            idx += 1

    system_df = pd.DataFrame(system_rows)

    # -------------------------
    # Chain-level DataFrame
    # -------------------------
    # Determine number of chains from confidence JSONs
    conf_path = (
        raw_dir
        / "repeat_1"
        / f"boltz_results_{system_name}"
        / "predictions"
        / system_name
        / f"confidence_{system_name}_model_0.json"
    )

    if not conf_path.exists():
        raise FileNotFoundError(f"Missing confidence file: {conf_path}")

    data = read.read_json(conf_path)
    if "chains_ptm" not in data:
        raise KeyError("chains_ptm missing from confidence JSON")

    chain_ids = sorted(data["chains_ptm"].keys(), key=int)

    chain_rows = []
    idx = 0

    for repeat in range(1, repeats + 1):
        for sample in range(diffusion_samples):
            for conf_chain_id in chain_ids:
                chain_rows.append(
                    {
                        "idx": idx,
                        "CHAIN_ID": None,               # filled later from system
                        "ENTITY_TYPE": None,            # filled later from system
                        "conf_chain_id": int(conf_chain_id),
                        "cif_file": f'{repeat}_{system_name}_model_{sample}.cif',
                        "model_name": system_name,
                        "repeat": repeat,
                        "diffusion_sample": sample,
                    }
                )
                idx += 1

    chain_df = pd.DataFrame(chain_rows)

    return system_df, chain_df

def add_chain_info(chain_df: pd.DataFrame, sys: "System") -> pd.DataFrame:
    """
    Populate CHAIN_ID, ENTITY_TYPE, and ligand_molecule_id in chain_df
    based on system sequences.

    Mapping is done by POSITION:
    conf_chain_id (0,1,2,...) → ordered system chain IDs.
    """

    sequences = sys.find_value(key="sequences") or []
    if not sequences:
        logger.warning("No sequences found in system to map chain information.")
        return chain_df

    # --------------------------------------------------
    # Build ordered chain list and metadata from system
    # --------------------------------------------------
    ordered_chain_ids: list[str] = []
    chain_to_entity: dict[str, str] = {}
    chain_to_molecule_id: dict[str, str] = {}

    for seq_entry in sequences:
        if not isinstance(seq_entry, dict):
            continue

        entity_type = next(iter(seq_entry))
        entity_data = seq_entry[entity_type]

        chain_ids = entity_data.get("id")
        if chain_ids is None:
            continue

        if not isinstance(chain_ids, list):
            chain_ids = [chain_ids]

        # ----------------------------------------------
        # Resolve molecule identifier
        # ----------------------------------------------
        if entity_type == "ligand":
            molecule_id = (
                entity_data.get("ccd")
                or entity_data.get("smiles")
                or "UNKNOWN_LIGAND"
            )
        else:
            molecule_id = None  # handled per-chain below

        for cid in chain_ids:
            cid = str(cid)
            ordered_chain_ids.append(cid)
            chain_to_entity[cid] = entity_type

            if entity_type == "ligand":
                chain_to_molecule_id[cid] = molecule_id
            else:
                chain_to_molecule_id[cid] = f"{entity_type}_{cid}"

    # --------------------------------------------------
    # Assign CHAIN_ID, ENTITY_TYPE, ligand_molecule_id
    # --------------------------------------------------
    chain_df["CHAIN_ID"] = None
    chain_df["ENTITY_TYPE"] = None
    chain_df["ligand_molecule_id"] = None

    for idx, row in chain_df.iterrows():
        conf_id = int(row["conf_chain_id"])

        if conf_id >= len(ordered_chain_ids):
            logger.warning(
                "conf_chain_id %d exceeds system chain count (%d)",
                conf_id,
                len(ordered_chain_ids),
            )
            continue

        chain_id = ordered_chain_ids[conf_id]
        entity_type = chain_to_entity[chain_id]

        chain_df.at[idx, "CHAIN_ID"] = chain_id
        chain_df.at[idx, "ENTITY_TYPE"] = entity_type
        chain_df.at[idx, "ligand_molecule_id"] = chain_to_molecule_id[chain_id]

    # --------------------------------------------------
    # Reorder columns
    # --------------------------------------------------
    cols = chain_df.columns.tolist()
    for col in ("CHAIN_ID", "ENTITY_TYPE", "ligand_molecule_id"):
        if col in cols:
            cols.remove(col)

    insert_at = cols.index("conf_chain_id")
    cols.insert(insert_at + 1, "CHAIN_ID")
    cols.insert(insert_at + 2, "ENTITY_TYPE")
    cols.insert(insert_at + 3, "ligand_molecule_id")

    return chain_df[cols]

def gather_confidence_metrics(
    raw_dir: Path,
    system_df: pd.DataFrame,
    chain_df: pd.DataFrame,
    system_name: str,
    repeats: int,
    diffusion_samples: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Gather confidence metrics and append them as new columns to existing DataFrames.

    Parameters
    ----------
    raw_dir : Path
        Path to wrk_dir/raw.
    system_df : pd.DataFrame
        DataFrame with columns ['idx', 'model_name', 'repeat', 'diffusion_sample']
    chain_df : pd.DataFrame
        DataFrame with columns ['idx', 'CHAIN_ID', 'ENTITY_TYPE', 'conf_chain_id', 'model_name', 'repeat', 'diffusion_sample']
    system_name : str
        Name of the system.
    repeats : int
        Number of repeats.
    diffusion_samples : int
        Number of diffusion samples.

    Returns
    -------
    system_df : pd.DataFrame
        Original DataFrame with new model-level metric columns appended.
    chain_df : pd.DataFrame
        Original DataFrame with new chain-level metric columns appended.
    """

    for repeat in range(1, repeats + 1):
        for sample_idx in range(diffusion_samples):
            conf_path = (
                raw_dir
                / f"repeat_{repeat}"
                / f"boltz_results_{system_name}"
                / "predictions"
                / system_name
                / f"confidence_{system_name}_model_{sample_idx}.json"
            )

            if not conf_path.exists():
                logger.warning("Missing confidence file: %s", conf_path)
                continue

            data = read.read_json(conf_path)
            if not data:
                continue

            # --- Append model-level metrics ---
            model_mask = (
                (system_df["repeat"] == repeat) &
                (system_df["diffusion_sample"] == sample_idx)
            )

            for key, value in data.items():
                if key in {"chains_ptm", "pair_chains_iptm"}:
                    continue
                elif not isinstance(value, (dict, list)):
                    system_df.loc[model_mask, key] = value

            # --- Append chain-level metrics ---
            chains_ptm = data.get("chains_ptm", {})
            pair_chains_iptm = data.get("pair_chains_iptm", {})

            for idx, row in chain_df.iterrows():
                # Only update rows for the current repeat and diffusion_sample
                if row["repeat"] != repeat or row["diffusion_sample"] != sample_idx:
                    continue

                cid = str(row["conf_chain_id"])

                # Add chains_ptm if present
                if cid in chains_ptm:
                    chain_df.at[idx, "chains_ptm"] = chains_ptm[cid]

                # Flatten all inner values of pair_chains_iptm for this chain
                if cid in pair_chains_iptm:
                    inner_dict = pair_chains_iptm[cid]
                    for other_chain_id, val in inner_dict.items():
                        col_name = f"pair_chains_iptm_{other_chain_id}"
                        chain_df.at[idx, col_name] = val

    logger.info(
        "Confidence metrics appended — system_df: %d rows, chain_df: %d rows",
        len(system_df),
        len(chain_df),
    )

    return system_df, chain_df

def gather_affinity_metrics(
    raw_dir: Path,
    chain_df: pd.DataFrame,
    system_name: str,
    sys: "System",
    repeats: int,
    extended: bool = False,
) -> pd.DataFrame:
    """
    Gather affinity metrics and append them to chain_df.

    Rules:
    - Affinity is enabled only if specified in system properties
    - Affinity is repeat-specific
    - All diffusion samples within a repeat share the same affinity values
    - Only the binder ligand chain receives affinity values
    """

    # --------------------------------------------------
    # Check if affinity prediction is enabled
    # --------------------------------------------------
    properties = sys.find_value(key="properties") or []
    affinity_props = None

    for prop in properties:
        if isinstance(prop, dict) and "affinity" in prop:
            affinity_props = prop["affinity"]
            break

    if affinity_props is None:
        return chain_df

    binder_chain_id = affinity_props.get("binder")
    if binder_chain_id is None:
        logger.warning("Affinity specified but no binder chain ID found.")
        return chain_df

    binder_chain_id = str(binder_chain_id)

    # --------------------------------------------------
    # Ensure output columns exist
    # --------------------------------------------------
    base_columns = (
        "affinity_pred_value",
        "affinity_probability_binary",
    )

    extended_columns = (
        "pIC50",
        "IC50_M",
        "pIC50_kcal_per_mol",
    )

    for col in base_columns:
        if col not in chain_df.columns:
            chain_df[col] = None

    if extended:
        for col in extended_columns:
            if col not in chain_df.columns:
                chain_df[col] = None

   # --------------------------------------------------
    # Repeat loop
    # --------------------------------------------------
    for repeat in range(1, repeats + 1):

        affinity_path = (
            raw_dir
            / f"repeat_{repeat}"
            / f"boltz_results_{system_name}"
            / "predictions"
            / system_name
            / f"affinity_{system_name}.json"
        )

        if not affinity_path.exists():
            logger.warning("Missing affinity file: %s", affinity_path)
            continue

        affinity_data = read.read_json(affinity_path)
        if not affinity_data:
            continue

        affinity_pred_value = affinity_data.get("affinity_pred_value")
        affinity_probability_binary = affinity_data.get(
            "affinity_probability_binary"
        )

        if affinity_pred_value is None or affinity_probability_binary is None:
            logger.warning(
                "Affinity values missing in file: %s", affinity_path
            )
            continue

        # --------------------------------------------------
        # Derived conversions (extended metrics only)
        # --------------------------------------------------
        if extended:
            pIC50, IC50_M = stats.affinity_to_pic50_and_ic50(
                affinity_pred_value
            )
            pIC50_kcal_per_mol = stats.affinity_to_pic50_kcal_per_mol(
                affinity_pred_value
            )

        # --------------------------------------------------
        # Assign to chain_df
        # --------------------------------------------------
        for idx, row in chain_df.iterrows():
            if row["repeat"] != repeat:
                continue
            if row["CHAIN_ID"] != binder_chain_id:
                continue

            chain_df.at[idx, "affinity_pred_value"] = affinity_pred_value
            chain_df.at[idx, "affinity_probability_binary"] = (
                affinity_probability_binary
            )

            if extended:
                chain_df.at[idx, "pIC50"] = pIC50
                chain_df.at[idx, "IC50_M"] = IC50_M
                chain_df.at[idx, "pIC50_kcal_per_mol"] = pIC50_kcal_per_mol

            logger.info(
                "Affinity metrics appended (extended=%s) — chain_df: %d rows",
                extended,
                len(chain_df),
            )

    return chain_df

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
    var = std ** 2

    return {
        f"{prefix}_mean": mean,
        f"{prefix}_std": std,
    }

def assess_bitstring_similarity(
    bitstrings: list,
    prefix: str,
) -> dict[str, float | None]:
    """
    Assess robustness of binary fingerprints (e.g. IFPs) across runs
    using pairwise Tanimoto similarity.

    Parameters
    ----------
    bitstrings : list
        One fingerprint per run.
    prefix : str
        Prefix for output metric names.

    Returns
    -------
    dict
        Mean, std, min, max pairwise similarity.
    """

    def to_set(x):
        if x is None:
            return None
        if isinstance(x, set):
            return x
        if isinstance(x, str):
            return {i for i, v in enumerate(x) if v == "1"}
        # assume iterable of bool/int
        return {i for i, v in enumerate(x) if bool(v)}

    sets = [to_set(x) for x in bitstrings if x is not None]

    if len(sets) < 2:
        return {
            f"{prefix}_mean": None,
            f"{prefix}_std": None,
        }

    sims = []

    for i in range(len(sets)):
        for j in range(i + 1, len(sets)):
            a, b = sets[i], sets[j]

            # both empty → identical
            if not a and not b:
                sims.append(1.0)
            else:
                sims.append(len(a & b) / len(a | b))

    s = pd.Series(sims, dtype="float64")

    return {
        f"{prefix}_mean": float(s.mean()),
        f"{prefix}_std": float(s.std(ddof=1)),
    }

def assess_rmsd_robustness_placeholder(
    wrk_dir,
    entity_type: str,
    entity_id: str,
) -> dict[str, float | None]:
    """
    Placeholder for RMSD robustness assessment across runs.

    Parameters
    ----------
    wrk_dir : Path
        Working directory containing CIF files.
    entity_type : str
        'system', 'protein', or 'ligand'
    entity_id : str
        Chain ID or 'system'

    Returns
    -------
    dict
        RMSD robustness metrics (currently empty).
    """
    # TODO:
    # - collect CIFs across repeats / diffusion samples
    # - align structures
    # - compute pairwise RMSDs
    # - feed RMSDs into assess_numeric_variance
    return {
        "rmsd_mean": None,
        "rmsd_std": None,
    }

def gather_robustness_results(
    system_df: pd.DataFrame,
    chain_df: pd.DataFrame,
    wrk_dir,
) -> pd.DataFrame:
    """
    Aggregate robustness metrics across repeats and diffusion samples.

    Robustness is computed for *all calculated metrics* present in the
    system_df and chain_df, excluding metadata columns.

    Ligands are aggregated by ligand_molecule_id.
    Polymers (protein / RNA / DNA) are aggregated per chain.
    """

    rows: list[dict] = []

    logger.debug(
        "[ROBUSTNESS] Starting robustness aggregation "
        "(system_df columns=%d, chain_df columns=%d)",
        len(system_df.columns),
        len(chain_df.columns),
    )

    # ==============================================================
    # SYSTEM-LEVEL ROBUSTNESS
    # ==============================================================
    system_row = {
        "ENTITY_TYPE": "system",
        "ENTITY_ID": "system",
    }

    for col in system_df.columns:
        logger.debug("[ROBUSTNESS][SYSTEM] Evaluating column '%s'", col)

        if is_metadata_column(col):
            continue

        series = system_df[col]

        if is_numeric_metric(series):
            logger.debug(
                "[ROBUSTNESS][SYSTEM] Numeric variance assessment for '%s'",
                col,
            )
            stats = assess_numeric_variance(
                series.tolist(),
                prefix=col,
            )
            system_row.update(stats)
        else:
            logger.debug(
                "[ROBUSTNESS][SYSTEM] Column '%s' skipped (non-numeric)",
                col,
            )

    # RMSD placeholder (system)
    system_row.update(
        assess_rmsd_robustness_placeholder(
            wrk_dir=wrk_dir,
            entity_type="system",
            entity_id="system",
        )
    )

    rows.append(system_row)

    # ==============================================================
    # CHAIN-LEVEL ROBUSTNESS
    # ==============================================================
    ligand_df = chain_df[chain_df["ENTITY_TYPE"] == "ligand"]
    non_ligand_df = chain_df[chain_df["ENTITY_TYPE"] != "ligand"]

    # --------------------------------------------------------------
    # NON-LIGAND CHAINS (protein / RNA / DNA)
    # --------------------------------------------------------------
    for (chain_id, entity_type), group in non_ligand_df.groupby(
        ["CHAIN_ID", "ENTITY_TYPE"]
    ):
        logger.debug(
            "[ROBUSTNESS][CHAIN] Processing chain '%s' (%s)",
            chain_id,
            entity_type,
        )

        chain_row = {
            "ENTITY_TYPE": entity_type,
            "ENTITY_ID": chain_id,
        }

        for col in group.columns:
            logger.debug(
                "[ROBUSTNESS][CHAIN %s] Evaluating column '%s'",
                chain_id,
                col,
            )

            if is_metadata_column(col):
                continue

            series = group[col]

            if is_numeric_metric(series):
                stats = assess_numeric_variance(
                    series.tolist(),
                    prefix=col,
                )
                chain_row.update(stats)

            elif is_ifp_metric(col):
                stats = assess_bitstring_similarity(
                    series.tolist(),
                    prefix=col,
                )
                chain_row.update(stats)

        # RMSD placeholder (chain)
        chain_row.update(
            assess_rmsd_robustness_placeholder(
                wrk_dir=wrk_dir,
                entity_type=entity_type,
                entity_id=chain_id,
            )
        )

        rows.append(chain_row)

    # --------------------------------------------------------------
    # LIGANDS — AGGREGATE BY ligand_molecule_id
    # --------------------------------------------------------------
    for ligand_molecule_id, group in ligand_df.groupby("ligand_molecule_id"):
        logger.debug(
            "[ROBUSTNESS][LIGAND] Processing molecule '%s'",
            ligand_molecule_id,
        )

        ligand_row = {
            "ENTITY_TYPE": "ligand",
            "ENTITY_ID": ligand_molecule_id,
        }

        for col in group.columns:
            logger.debug(
                "[ROBUSTNESS][LIGAND %s] Evaluating column '%s'",
                ligand_molecule_id,
                col,
            )

            if is_metadata_column(col):
                continue

            series = group[col]

            if is_numeric_metric(series):
                stats = assess_numeric_variance(
                    series.tolist(),
                    prefix=col,
                )
                ligand_row.update(stats)

            elif is_ifp_metric(col):
                stats = assess_bitstring_similarity(
                    series.tolist(),
                    prefix=col,
                )
                ligand_row.update(stats)

        # RMSD placeholder (ligand molecule)
        ligand_row.update(
            assess_rmsd_robustness_placeholder(
                wrk_dir=wrk_dir,
                entity_type="ligand",
                entity_id=ligand_molecule_id,
            )
        )

        rows.append(ligand_row)

    results_df = pd.DataFrame(rows)

    logger.debug(
        "[ROBUSTNESS] Finished robustness aggregation "
        "(rows=%d, columns=%d)",
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
        "ENTITY_TYPE",
        "conf_chain_id",
    }

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

def is_ifp_metric(column_name: str) -> bool:
    """
    Determine whether a column represents a binary interaction fingerprint.
    """
    name = column_name.lower()
    result = name.startswith("ifp") or "fingerprint" in name

    logger.debug(
        "[ROBUSTNESS] Metric '%s' classified as IFP: %s",
        column_name,
        result,
    )

    return result