"""Shared helpers for reading and qualifying Validate metric outputs."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import pandas as pd

SYSTEM_METADATA_COLUMNS = frozenset(
    {"idx", "cif_file", "model_name", "repeat", "diffusion_sample"}
)
CHAIN_METADATA_COLUMNS = frozenset(
    {
        "idx",
        "conf_chain_id",
        "CHAIN_ID",
        "ENTITY_TYPE",
        "ligand_molecule_id",
        "cif_file",
        "model_name",
        "repeat",
        "diffusion_sample",
    }
)


def read_metric_frames(run_dir: str | Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read Validate's system and chain CSVs, returning empty frames if absent."""

    results_dir = Path(run_dir) / "results"
    system_path = results_dir / "system_metrics.csv"
    chain_path = results_dir / "chain_metrics.csv"
    system_df = pd.read_csv(system_path) if system_path.exists() else pd.DataFrame()
    chain_df = pd.read_csv(chain_path) if chain_path.exists() else pd.DataFrame()
    return system_df, chain_df


def collect_qualified_metric_values(
    system_df: pd.DataFrame,
    chain_df: pd.DataFrame,
) -> dict[str, tuple[Any, ...]]:
    """Return all metric values under Screen-compatible qualified names."""

    values: dict[str, list[Any]] = {}
    for column in system_df.columns:
        if column in SYSTEM_METADATA_COLUMNS:
            continue
        values[f"system__{column}"] = system_df[column].tolist()

    if not chain_df.empty and not {"CHAIN_ID", "ENTITY_TYPE"}.issubset(
        chain_df.columns
    ):
        # Retain compatibility with older/minimal Validate outputs that did not
        # carry chain identity columns. Bare selectors can still resolve these;
        # callers need modern metadata to use entity-qualified selectors.
        for column in chain_df.columns:
            if column in CHAIN_METADATA_COLUMNS:
                continue
            values[f"chain__{column}"] = chain_df[column].tolist()

    if not chain_df.empty and {"CHAIN_ID", "ENTITY_TYPE"}.issubset(chain_df.columns):
        conf_chain_ids: dict[str, str] = {}
        if "conf_chain_id" in chain_df.columns:
            for _, row in chain_df.iterrows():
                if pd.isna(row.get("conf_chain_id")):
                    continue
                chain_id = str(row.get("CHAIN_ID", "")).strip()
                if not chain_id:
                    continue
                try:
                    conf_id = str(int(row["conf_chain_id"]))
                except (TypeError, ValueError):
                    continue
                conf_chain_ids.setdefault(conf_id, chain_id)

        for _, row in chain_df.iterrows():
            chain_id = str(row.get("CHAIN_ID", "")).strip() or "NA"
            entity_type = str(row.get("ENTITY_TYPE", "")).strip() or "unknown"
            prefix = f"{entity_type}_{chain_id}"
            for column, value in row.items():
                if column in CHAIN_METADATA_COLUMNS:
                    continue
                values.setdefault(f"{prefix}__{column}", []).append(value)
                match = re.fullmatch(r"pair_chains_iptm_(\d+)", str(column))
                if match and match.group(1) in conf_chain_ids:
                    semantic = f"pair_chains_iptm_{conf_chain_ids[match.group(1)]}"
                    values.setdefault(f"{prefix}__{semantic}", []).append(value)

    return {name: tuple(items) for name, items in values.items()}


def primary_metric_values(
    system_df: pd.DataFrame,
    chain_df: pd.DataFrame,
) -> dict[str, Any]:
    """Return the preferred repeat-1/sample-0 value for each qualified metric."""

    primary_system = _primary_rows(system_df, per_chain=False)
    primary_chain = _primary_rows(chain_df, per_chain=True)
    values = collect_qualified_metric_values(primary_system, primary_chain)
    return {name: items[0] for name, items in values.items() if items}


def _primary_rows(df: pd.DataFrame, *, per_chain: bool) -> pd.DataFrame:
    if df.empty:
        return df
    if {"repeat", "diffusion_sample"}.issubset(df.columns):
        selected = df[(df["repeat"] == 1) & (df["diffusion_sample"] == 0)]
        if not selected.empty:
            return selected
    if per_chain and "CHAIN_ID" in df.columns:
        return df.drop_duplicates(subset=["CHAIN_ID"], keep="first")
    return df.iloc[:1]
