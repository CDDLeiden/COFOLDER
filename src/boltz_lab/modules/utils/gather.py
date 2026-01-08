import shutil
from pathlib import Path
from typing import Tuple
import logging
import pandas as pd

from boltz_lab.modules.utils import read, write

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
    Populate CHAIN_ID and ENTITY_TYPE in chain_df based on system sequences.

    Mapping is done by POSITION:
    conf_chain_id (0,1,2,...) → ordered system chain IDs.
    """

    sequences = sys.find_value(key="sequences") or []
    if not sequences:
        logger.warning("No sequences found in system to map chain information.")
        return chain_df

    # --------------------------------------------------
    # Build ordered chain list from system
    # --------------------------------------------------
    ordered_chain_ids: list[str] = []
    chain_to_entity: dict[str, str] = {}

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

        for cid in chain_ids:
            cid = str(cid)
            ordered_chain_ids.append(cid)
            chain_to_entity[cid] = entity_type

    # --------------------------------------------------
    # Assign CHAIN_ID + ENTITY_TYPE by conf_chain_id
    # --------------------------------------------------
    chain_df["CHAIN_ID"] = None
    chain_df["ENTITY_TYPE"] = None

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
        chain_df.at[idx, "CHAIN_ID"] = chain_id
        chain_df.at[idx, "ENTITY_TYPE"] = chain_to_entity[chain_id]

    # --------------------------------------------------
    # Reorder columns: CHAIN_ID, ENTITY_TYPE together
    # --------------------------------------------------
    cols = chain_df.columns.tolist()
    for col in ("CHAIN_ID", "ENTITY_TYPE"):
        if col in cols:
            cols.remove(col)

    insert_at = cols.index("conf_chain_id")
    cols.insert(insert_at, "CHAIN_ID")
    cols.insert(insert_at + 1, "ENTITY_TYPE")

    chain_df = chain_df[cols]

    return chain_df

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