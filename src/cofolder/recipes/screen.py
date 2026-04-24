"""Screening workflow implemented as a Validate wrapper over CSV inputs."""

from __future__ import annotations

import copy
import logging
import re
from pathlib import Path
from typing import Any

import pandas as pd

from cofolder.modules.input import system
from cofolder.modules.utils import read, write
from cofolder.recipes.validate import Validate

logger = logging.getLogger(__name__)


class Screen:
    """Run Validate for each molecule row in a CSV by adapting system.yaml."""

    def __init__(
        self,
        wrk_dir: str,
        system_path: str,
        options_path: str,
        runner: str = "boltz",
        variable: list[str] | None = None,
        variable_csv: str | None = None,
        col_variable: list[str] | None = None,
        col_id: str | None = None,
        merge_data: str | None = None,
        repeats: int = 1,
        seed: int | None = None,
        scoring_functions: list[str] | None = None,
        assess_robustness: bool = True,
        assess_bias: bool = False,
        protein_training_data_path: str | None = None,
        ligand_training_data_path: str | None = None,
        bias_release_cutoff: str = "2023-06-01",
        bias_ligand_similarity_threshold: float = 0.35,
        bias_chains: list[str] | None = None,
        build_bias_training_data: bool = False,
        bias_training_components_cif: str | None = None,
        conformers: str | None = None,
        sdf_file: str | None = None,
        reference_path: str | None = None,
        pocket_coverage_reference: str | None = None,
        reproduction_metrics: list[str] | None = None,
    ):
        self.wrk_dir = Path(wrk_dir)
        self.system_path = Path(system_path)
        self.options_path = Path(options_path)
        self.runner = str(runner)

        self.variable_raw = variable or []
        self.col_variable = col_variable or []
        self.col_id = col_id
        self.variable_csv = Path(variable_csv) if variable_csv else None
        self.merge_data = self._parse_list(merge_data)

        self.variable_paths = [self._parse_path(v) for v in self.variable_raw]

        self.validate_kwargs: dict[str, Any] = {
            "repeats": repeats,
            "seed": seed,
            "scoring_functions": scoring_functions,
            "assess_robustness": assess_robustness,
            "assess_bias": assess_bias,
            "protein_training_data_path": protein_training_data_path,
            "ligand_training_data_path": ligand_training_data_path,
            "bias_release_cutoff": bias_release_cutoff,
            "bias_ligand_similarity_threshold": bias_ligand_similarity_threshold,
            "bias_chains": bias_chains,
            "build_bias_training_data": build_bias_training_data,
            "bias_training_components_cif": bias_training_components_cif,
            "conformers": conformers,
            "sdf_file": sdf_file,
            "reference_path": reference_path,
            "pocket_coverage_reference": pocket_coverage_reference,
            "reproduction_metrics": reproduction_metrics,
        }

        self.base_system = read.read_yaml(path=self.system_path)
        self.logger = logging.getLogger("cofolder.screen")

        self._validate_config()

    def _validate_config(self) -> None:
        if self.variable_csv is None:
            raise ValueError("--variable_csv is required.")
        if not self.variable_csv.exists():
            raise ValueError(f"--variable_csv does not exist: {self.variable_csv}")
        if not self.variable_csv.is_file():
            raise ValueError(f"--variable_csv is not a file: {self.variable_csv}")
        if self.col_id is None or not str(self.col_id).strip():
            raise ValueError("--col_id is required.")
        if not self.variable_raw or not self.col_variable:
            raise ValueError("At least one --variable/--col_variable pair is required.")
        if len(self.variable_raw) != len(self.col_variable):
            raise ValueError(
                "Number of --variable entries must match number of --col_variable entries. "
                f"Got variable={len(self.variable_raw)} col_variable={len(self.col_variable)}."
            )

    def run(self) -> None:
        self.wrk_dir.mkdir(parents=True, exist_ok=True)
        df = pd.read_csv(self.variable_csv)

        required_cols = [self.col_id, *self.col_variable, *self.merge_data]
        missing = [c for c in required_cols if c not in df.columns]
        if missing:
            raise ValueError(f"CSV missing required columns: {', '.join(missing)}")

        records: list[dict[str, Any]] = []
        records_with_scores: list[dict[str, Any]] = []
        total = len(df)
        for i, (_, row) in enumerate(df.iterrows(), 1):
            compound_id = str(row[self.col_id])
            safe_id = self._safe_name(compound_id)
            run_dir = self.wrk_dir / f"{i}_{safe_id}"
            run_dir.mkdir(parents=True, exist_ok=True)
            row_system_path = run_dir / "screen_system.yaml"

            summary = {
                "index": i,
                self.col_id: compound_id,
                "status": "success",
                "error_message": "",
                "run_dir": str(run_dir),
            }
            detailed: dict[str, Any] = {k: row.get(k) for k in df.columns}
            detailed.update(summary)
            for col in self.col_variable:
                summary[col] = row.get(col)
            for col in self.merge_data:
                summary[col] = row.get(col)

            self.logger.info("(%d/%d) screening %s", i, total, compound_id)

            try:
                sys_obj = system.System(system=copy.deepcopy(self.base_system))
                for path, col in zip(self.variable_paths, self.col_variable):
                    value = row[col]
                    if pd.isna(value) or (isinstance(value, str) and value.strip() == ""):
                        raise ValueError(
                            f"Empty value for mapped column '{col}' in row {i} ({compound_id})."
                        )
                    if self._path_targets_smiles(path) and not self._is_valid_smiles(str(value)):
                        raise ValueError(
                            f"Invalid SMILES in column '{col}' for row {i} ({compound_id}): {value}"
                        )
                    sys_obj.update_system(value=value, path=path)

                write.write_yaml(sys_obj, path=row_system_path)

                row_validate_kwargs = dict(self.validate_kwargs)
                if (
                    row_validate_kwargs.get("assess_bias")
                    and row_validate_kwargs.get("build_bias_training_data")
                ):
                    bias_train_dir = run_dir / "results" / "bias_train"
                    bias_train_dir.mkdir(parents=True, exist_ok=True)
                    # Build ligand references per-screened system after row system YAML exists.
                    row_validate_kwargs["ligand_training_data_path"] = str(
                        bias_train_dir / "ligand_training_data.csv"
                    )

                validator = Validate(
                    wrk_dir=str(run_dir),
                    system_path=str(row_system_path),
                    options_path=str(self.options_path),
                    runner=self.runner,
                    **row_validate_kwargs,
                )
                validator.run()
                detailed.update(self._collect_score_columns(run_dir=run_dir))
            except Exception as exc:
                summary["status"] = "failed"
                summary["error_message"] = str(exc)
                detailed["status"] = "failed"
                detailed["error_message"] = str(exc)
                self.logger.exception("Screen row failed (%s): %s", compound_id, exc)

            records.append(summary)
            records_with_scores.append(detailed)

        out_csv = self.wrk_dir / "screen_results.csv"
        pd.DataFrame(records).to_csv(out_csv, index=False)
        out_scores_csv = self.wrk_dir / "screen_results_with_scores.csv"
        pd.DataFrame(records_with_scores).to_csv(out_scores_csv, index=False)

        failures = sum(1 for r in records if r["status"] == "failed")
        successes = len(records) - failures
        self.logger.info(
            "Screen complete: total=%d success=%d failed=%d summary=%s merged=%s",
            len(records),
            successes,
            failures,
            out_csv,
            out_scores_csv,
        )

    def _collect_score_columns(self, run_dir: Path) -> dict[str, Any]:
        """Collect computed scores from per-row validate outputs."""
        out: dict[str, Any] = {}
        system_csv = run_dir / "results" / "system_metrics.csv"
        chain_csv = run_dir / "results" / "chain_metrics.csv"

        if system_csv.exists():
            try:
                sdf = pd.read_csv(system_csv)
                if not sdf.empty:
                    row = self._select_metrics_row(sdf)
                    for col, value in row.items():
                        if col in {"idx", "cif_file", "model_name", "repeat", "diffusion_sample"}:
                            continue
                        out[f"system__{col}"] = value
            except Exception as exc:
                self.logger.warning("Failed reading system metrics from %s: %s", system_csv, exc)

        if chain_csv.exists():
            try:
                cdf = pd.read_csv(chain_csv)
                if not cdf.empty and {"CHAIN_ID", "ENTITY_TYPE"}.issubset(cdf.columns):
                    cdf = self._select_chain_metrics_rows(cdf)
                    for _, crow in cdf.iterrows():
                        chain_id = str(crow.get("CHAIN_ID", "")).strip() or "NA"
                        entity_type = str(crow.get("ENTITY_TYPE", "")).strip() or "unknown"
                        prefix = f"{entity_type}_{chain_id}"
                        for col, value in crow.items():
                            if col in {
                                "idx",
                                "conf_chain_id",
                                "CHAIN_ID",
                                "ENTITY_TYPE",
                                "ligand_molecule_id",
                                "cif_file",
                                "model_name",
                                "repeat",
                                "diffusion_sample",
                            }:
                                continue
                            out[f"{prefix}__{col}"] = value
            except Exception as exc:
                self.logger.warning("Failed reading chain metrics from %s: %s", chain_csv, exc)

        return out

    @staticmethod
    def _select_metrics_row(df: pd.DataFrame) -> pd.Series:
        """Prefer repeat=1/sample=0 row if available; otherwise first row."""
        if {"repeat", "diffusion_sample"}.issubset(df.columns):
            sub = df[(df["repeat"] == 1) & (df["diffusion_sample"] == 0)]
            if not sub.empty:
                return sub.iloc[0]
        return df.iloc[0]

    @staticmethod
    def _select_chain_metrics_rows(df: pd.DataFrame) -> pd.DataFrame:
        """Prefer repeat=1/sample=0 chain rows when available; otherwise de-duplicate by chain."""
        if {"repeat", "diffusion_sample"}.issubset(df.columns):
            sub = df[(df["repeat"] == 1) & (df["diffusion_sample"] == 0)]
            if not sub.empty:
                return sub
        if "CHAIN_ID" in df.columns:
            return df.drop_duplicates(subset=["CHAIN_ID"], keep="first")
        return df

    @staticmethod
    def _parse_path(path_str: str) -> list[Any]:
        if not path_str or not str(path_str).strip():
            raise ValueError("--variable path cannot be empty.")
        return [int(v.strip()) if v.strip().isdigit() else v.strip() for v in str(path_str).split(",")]

    @staticmethod
    def _parse_list(input_str: str | None) -> list[str]:
        if not input_str:
            return []
        return [v.strip() for v in str(input_str).split(",") if v.strip()]

    @staticmethod
    def _safe_name(value: str) -> str:
        token = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())
        return token or "item"

    @staticmethod
    def _path_targets_smiles(path: list[Any]) -> bool:
        return bool(path) and str(path[-1]).strip().lower() == "smiles"

    @staticmethod
    def _is_valid_smiles(smiles: str) -> bool:
        try:
            from rdkit import Chem
        except Exception:
            # If RDKit is unavailable, do not block screen row execution here.
            return True
        s = str(smiles).strip()
        if not s:
            return False
        return Chem.MolFromSmiles(s) is not None
