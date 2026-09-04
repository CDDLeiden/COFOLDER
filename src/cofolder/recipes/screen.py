"""Screening workflow implemented as a Validate wrapper over CSV inputs."""

from __future__ import annotations

import copy
import json
import logging
import math
import re
from pathlib import Path
from typing import Any

import pandas as pd

from cofolder.modules.analytics.ifp_clustering import (
    SUMMARY_COLUMNS as IFP_CLUSTER_SUMMARY_COLUMNS,
)
from cofolder.modules.analytics.ifp_clustering import cluster_binary_ifps
from cofolder.modules.analytics.reproduction import (
    _build_reference_ifp_from_custom,
    _parse_custom_pocket_reference,
)
from cofolder.modules.input import system
from cofolder.modules.input.system import iter_system_chains
from cofolder.modules.runners import get_runner
from cofolder.modules.runners.msa import resolve_declared_msa_paths
from cofolder.modules.utils import read, write
from cofolder.recipes._metrics import primary_metric_values, read_metric_frames
from cofolder.recipes.validate import Validate

logger = logging.getLogger(__name__)


SYSTEM_SCREEN_METRICS = (
    "confidence_score",
    "ptm",
    "iptm",
    "bias_prot_sim_train_max",
    "bias_prot_sim_train_pairwise_max",
    "bias_lig_sim_train_max",
    "pocket_coverage_ref",
    "pocket_coverage_ref_mean",
    "pocket_coverage_custom",
    "pocket_coverage_custom_mean",
)
PROTEIN_SCREEN_METRICS = (
    "chains_ptm",
    "bias_prot_sim_train",
    "bias_prot_sim_train_pairwise",
)
LIGAND_SCREEN_METRICS = (
    "chains_ptm",
    "affinity_pred_value",
    "affinity_probability_binary",
    "pIC50",
    "IC50_M",
    "pIC50_kcal_per_mol",
    "sasa",
    "sasa_norm_heavy",
    "ifp_distance",
    "ifp_prolif",
    "bias_lig_sim_train",
    "pocket_coverage_ref",
    "pocket_coverage_custom",
    "ligand_pose_overlap_ref",
)


class Screen:
    """Run Validate for each CSV row and consolidate optional IFP analyses.

    ``cluster_ifps`` applies deterministic average-linkage clustering to compatible
    binary distance IFPs after prediction. ``ifp_cluster_similarity_threshold`` is
    the inclusive Jaccard-similarity cut and defaults to ``0.5``. Both consolidated
    CSVs always contain cluster ID/status columns; an enabled run additionally writes
    ``ifp_cluster_summary.csv``.
    """

    def __init__(
        self,
        wrk_dir: str,
        system_path: str,
        options_path: str,
        runner: str = "boltz2",
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
        ifp_filter_threshold: float | None = None,
        ifp_ligand_chain: str | None = None,
        cluster_ifps: bool = False,
        ifp_cluster_similarity_threshold: float = 0.5,
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
        self.ifp_filter_threshold = ifp_filter_threshold
        self.ifp_ligand_chain = (
            str(ifp_ligand_chain).strip() if ifp_ligand_chain is not None else None
        )
        self.pocket_coverage_reference = pocket_coverage_reference
        self.cluster_ifps = bool(cluster_ifps)
        self.ifp_cluster_similarity_threshold = ifp_cluster_similarity_threshold
        self._ifp_filter_reference_spec: dict[str, Any] | None = None

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
        self.base_system_obj = system.System(system=self.base_system)
        resolve_declared_msa_paths(
            self.base_system_obj,
            base_dir=self.system_path.parent,
        )
        self.base_system = self.base_system_obj.system
        self.runner_impl = get_runner(self.runner)
        self.reusable_msa_dir = (
            self.wrk_dir / "shared" / "msa" / self.runner
            if getattr(self.runner_impl, "supports_msa_reuse", False)
            else None
        )
        if self.reusable_msa_dir is not None:
            options_obj = self.runner_impl.load_options(self.options_path)
            resolve_msa_reuse_settings = getattr(
                self.runner_impl, "msa_reuse_settings", None
            )
            self.msa_reuse_settings = (
                resolve_msa_reuse_settings(options_obj)
                if callable(resolve_msa_reuse_settings)
                else {}
            )
        else:
            self.msa_reuse_settings = {}
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

        try:
            self.ifp_cluster_similarity_threshold = float(
                self.ifp_cluster_similarity_threshold
            )
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "--ifp_cluster_similarity_threshold must be a number in [0, 1]."
            ) from exc
        if (
            not math.isfinite(self.ifp_cluster_similarity_threshold)
            or not 0 <= self.ifp_cluster_similarity_threshold <= 1
        ):
            raise ValueError("--ifp_cluster_similarity_threshold must be in [0, 1].")

        scoring_functions = self.validate_kwargs.get("scoring_functions")
        if (
            self.cluster_ifps
            and scoring_functions is not None
            and "ifp_distance" not in scoring_functions
        ):
            raise ValueError(
                "--cluster_ifps requires distance IFP scoring; include "
                "'ifp_distance' in --scoring_functions."
            )

        if self.ifp_filter_threshold is None:
            if self.cluster_ifps:
                self._validate_ifp_ligand_selection("--cluster_ifps")
            return

        try:
            self.ifp_filter_threshold = float(self.ifp_filter_threshold)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "--ifp_filter_threshold must be a number in [0, 1]."
            ) from exc
        if (
            not math.isfinite(self.ifp_filter_threshold)
            or not 0 <= self.ifp_filter_threshold <= 1
        ):
            raise ValueError("--ifp_filter_threshold must be in [0, 1].")

        if scoring_functions is not None and "ifp_distance" not in scoring_functions:
            raise ValueError(
                "--ifp_filter_threshold requires distance IFP scoring; include "
                "'ifp_distance' in --scoring_functions."
            )

        self._validate_filter_reference()
        self._validate_ifp_ligand_selection("--ifp_filter_threshold")

    def _validate_ifp_ligand_selection(self, option: str) -> None:
        ligand_count = self._configured_ligand_count()
        if ligand_count == 0:
            raise ValueError(f"{option} requires a system containing a ligand.")
        if ligand_count > 1 and not self.ifp_ligand_chain:
            raise ValueError(
                f"--ifp_ligand_chain is required with {option} when the system "
                "contains multiple ligand chains."
            )

    def _validate_filter_reference(self) -> None:
        """Parse the filter reference before any prediction work starts."""
        value = self.pocket_coverage_reference
        if value is None or not str(value).strip():
            raise ValueError(
                "--ifp_filter_threshold requires --pocket_coverage_reference."
            )

        raw = str(value).strip()
        path = Path(raw).expanduser()
        looks_like_path = bool(path.suffix) or "/" in raw or "\\" in raw
        if looks_like_path and not path.exists():
            raise ValueError(f"--pocket_coverage_reference file does not exist: {path}")
        if path.exists() and not path.is_file():
            raise ValueError(f"--pocket_coverage_reference is not a file: {path}")

        try:
            reference_value = str(path) if path.exists() else raw
            self._ifp_filter_reference_spec = _parse_custom_pocket_reference(
                reference_value
            )
        except (OSError, ValueError) as exc:
            raise ValueError(f"Invalid --pocket_coverage_reference: {exc}") from exc
        if self._ifp_filter_reference_spec is None:
            raise ValueError("--pocket_coverage_reference must not be empty.")

    def _configured_ligand_count(self) -> int:
        count = 0
        for entry in self.base_system.get("sequences", []):
            if not isinstance(entry, dict) or "ligand" not in entry:
                continue
            ligand = entry.get("ligand")
            identifier = ligand.get("id") if isinstance(ligand, dict) else None
            count += len(identifier) if isinstance(identifier, list) else 1
        return count

    def run(self) -> pd.DataFrame:
        """Run the screen, write both CSV summaries, and return the merged results."""

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
            filter_result = self._default_filter_result()
            summary.update(filter_result)
            detailed: dict[str, Any] = {k: row.get(k) for k in df.columns}
            detailed.update(summary)
            self._ensure_screen_metric_schema(detailed)
            detailed.update(self._default_cluster_result())
            summary.update(self._default_cluster_result())
            for col in self.col_variable:
                summary[col] = row.get(col)
            for col in self.merge_data:
                summary[col] = row.get(col)

            self.logger.info("(%d/%d) screening %s", i, total, compound_id)

            sys_obj: system.System | None = None
            try:
                sys_obj = system.System(system=copy.deepcopy(self.base_system))
                for path, col in zip(self.variable_paths, self.col_variable):
                    value = row[col]
                    if pd.isna(value) or (
                        isinstance(value, str) and value.strip() == ""
                    ):
                        raise ValueError(
                            f"Empty value for mapped column '{col}' in row {i} ({compound_id})."
                        )
                    if self._path_targets_smiles(path) and not self._is_valid_smiles(
                        str(value)
                    ):
                        raise ValueError(
                            f"Invalid SMILES in column '{col}' for row {i} ({compound_id}): {value}"
                        )
                    sys_obj.update_system(value=value, path=path)

                if self.reusable_msa_dir is not None:
                    injected = self.runner_impl.inject_reusable_msas(
                        sys_obj,
                        self.reusable_msa_dir,
                        settings=self.msa_reuse_settings,
                    )
                    if injected:
                        self.logger.info(
                            "Reusing %d shared protein MSA(s) for %s.",
                            injected,
                            compound_id,
                        )

                write.write_yaml(sys_obj, path=row_system_path)

                row_validate_kwargs = dict(self.validate_kwargs)
                if row_validate_kwargs.get("assess_bias") and row_validate_kwargs.get(
                    "build_bias_training_data"
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
                    reusable_msa_dir=(
                        str(self.reusable_msa_dir)
                        if self.reusable_msa_dir is not None
                        else None
                    ),
                    **row_validate_kwargs,
                )
                validator.run()
                if self.reusable_msa_dir is not None:
                    self.runner_impl.inject_reusable_msas(
                        sys_obj,
                        self.reusable_msa_dir,
                        settings=self.msa_reuse_settings,
                    )
                    write.write_yaml(sys_obj, path=row_system_path)
                detailed.update(self._collect_score_columns(run_dir=run_dir))
                filter_result = self._evaluate_ifp_filter(run_dir=run_dir)
                summary.update(filter_result)
                detailed.update(filter_result)
            except Exception as exc:
                if self.reusable_msa_dir is not None and sys_obj is not None:
                    try:
                        injected = self.runner_impl.inject_reusable_msas(
                            sys_obj,
                            self.reusable_msa_dir,
                            settings=self.msa_reuse_settings,
                        )
                        if injected and row_system_path.exists():
                            write.write_yaml(sys_obj, path=row_system_path)
                    except Exception as cache_exc:
                        self.logger.warning(
                            "Could not update failed row YAML from reusable MSA cache: %s",
                            cache_exc,
                        )
                summary["status"] = "failed"
                summary["error_message"] = str(exc)
                detailed["status"] = "failed"
                detailed["error_message"] = str(exc)
                if self.ifp_filter_threshold is not None:
                    filter_result = self._not_evaluable_filter_result("row_failed")
                    summary.update(filter_result)
                    detailed.update(filter_result)
                self.logger.exception("Screen row failed (%s): %s", compound_id, exc)

            records.append(summary)
            records_with_scores.append(detailed)

        summary_df = pd.DataFrame(records)
        results_df = pd.DataFrame(records_with_scores)
        self._apply_ifp_clustering(summary_df, results_df)

        out_csv = self.wrk_dir / "screen_results.csv"
        self._set_filter_dtypes(summary_df)
        summary_df.to_csv(out_csv, index=False)
        out_scores_csv = self.wrk_dir / "screen_results_with_scores.csv"
        self._set_filter_dtypes(results_df)
        results_df.to_csv(out_scores_csv, index=False)

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
        return results_df

    def _default_cluster_result(self) -> dict[str, Any]:
        return {
            "ifp_cluster_id": None,
            "ifp_cluster_status": (
                "not_evaluable" if self.cluster_ifps else "not_applied"
            ),
        }

    def _apply_ifp_clustering(
        self,
        summary_df: pd.DataFrame,
        results_df: pd.DataFrame,
    ) -> None:
        """Annotate both outputs and write a deterministic cluster summary."""

        if not self.cluster_ifps:
            return

        parsed_rows: list[tuple[int, list[int], str]] = []
        expected_width: int | None = None
        for position, row in summary_df.iterrows():
            if row.get("status") != "success":
                continue
            fingerprint, _ = self._load_selected_ligand_ifp(Path(row["run_dir"]))
            if fingerprint is None:
                continue
            if expected_width is None:
                expected_width = len(fingerprint)
            if len(fingerprint) != expected_width:
                continue
            parsed_rows.append((position, fingerprint, str(row[self.col_id])))

        if parsed_rows:
            clustered = cluster_binary_ifps(
                [row[1] for row in parsed_rows],
                [row[2] for row in parsed_rows],
                similarity_threshold=self.ifp_cluster_similarity_threshold,
            )
            for (position, _, _), cluster_id in zip(parsed_rows, clustered.cluster_ids):
                for frame in (summary_df, results_df):
                    frame.at[position, "ifp_cluster_id"] = cluster_id
                    frame.at[position, "ifp_cluster_status"] = "clustered"
            cluster_summary = clustered.summary
        else:
            cluster_summary = pd.DataFrame(columns=IFP_CLUSTER_SUMMARY_COLUMNS)

        cluster_summary.to_csv(self.wrk_dir / "ifp_cluster_summary.csv", index=False)

    def _load_selected_ligand_ifp(self, run_dir: Path) -> tuple[list[int] | None, str]:
        chain_csv = run_dir / "results" / "chain_metrics.csv"
        if not chain_csv.exists():
            return None, "missing_chain_metrics"
        try:
            chain_df = pd.read_csv(chain_csv)
        except Exception:
            return None, "malformed_chain_metrics"
        required = {"CHAIN_ID", "ENTITY_TYPE", "ifp_distance"}
        if chain_df.empty or not required.issubset(chain_df.columns):
            return None, "missing_ifp"
        chain_df = self._select_chain_metrics_rows(chain_df)
        ligand_rows = chain_df[
            chain_df["ENTITY_TYPE"].astype(str).str.lower() == "ligand"
        ]
        if self.ifp_ligand_chain:
            ligand_rows = ligand_rows[
                ligand_rows["CHAIN_ID"].astype(str) == self.ifp_ligand_chain
            ]
        elif ligand_rows["CHAIN_ID"].astype(str).nunique() != 1:
            return None, "ambiguous_ligand_chain"
        if ligand_rows.empty:
            return None, "ligand_chain_not_found"
        return self._parse_binary_ifp(ligand_rows.iloc[0].get("ifp_distance"))

    def _default_filter_result(self) -> dict[str, Any]:
        if self.ifp_filter_threshold is not None:
            return self._not_evaluable_filter_result("missing_ifp")
        return {
            "ifp_filter_pass": pd.NA,
            "ifp_filter_status": "not_applied",
            "ifp_filter_reason": "filtering_disabled",
            "ifp_filter_overlap": None,
            "ifp_filter_threshold": None,
            "ifp_filter_reference": None,
        }

    def _not_evaluable_filter_result(self, reason: str) -> dict[str, Any]:
        return {
            "ifp_filter_pass": pd.NA,
            "ifp_filter_status": "not_evaluable",
            "ifp_filter_reason": reason,
            "ifp_filter_overlap": None,
            "ifp_filter_threshold": self.ifp_filter_threshold,
            "ifp_filter_reference": self.pocket_coverage_reference,
        }

    def _evaluate_ifp_filter(self, run_dir: Path) -> dict[str, Any]:
        """Evaluate strict reference overlap for one completed screen row."""
        if self.ifp_filter_threshold is None:
            return self._default_filter_result()

        chain_csv = run_dir / "results" / "chain_metrics.csv"
        if not chain_csv.exists():
            return self._not_evaluable_filter_result("missing_chain_metrics")
        try:
            chain_df = pd.read_csv(chain_csv)
        except Exception:
            return self._not_evaluable_filter_result("malformed_chain_metrics")
        required = {"CHAIN_ID", "ENTITY_TYPE", "ifp_distance"}
        if chain_df.empty or not required.issubset(chain_df.columns):
            return self._not_evaluable_filter_result("missing_ifp")

        chain_df = self._select_chain_metrics_rows(chain_df)
        ligand_rows = chain_df[
            chain_df["ENTITY_TYPE"].astype(str).str.lower() == "ligand"
        ]
        if self.ifp_ligand_chain:
            ligand_rows = ligand_rows[
                ligand_rows["CHAIN_ID"].astype(str) == self.ifp_ligand_chain
            ]
        elif ligand_rows["CHAIN_ID"].astype(str).nunique() != 1:
            return self._not_evaluable_filter_result("ambiguous_ligand_chain")
        if ligand_rows.empty:
            return self._not_evaluable_filter_result("ligand_chain_not_found")

        ligand_row = ligand_rows.iloc[0]
        pred_ifp, parse_reason = self._parse_binary_ifp(ligand_row.get("ifp_distance"))
        if pred_ifp is None:
            return self._not_evaluable_filter_result(parse_reason)

        ref_ifp = self._resolve_filter_reference(
            run_dir=run_dir,
            chain_df=chain_df,
            ligand_row=ligand_row,
        )
        if ref_ifp is None:
            return self._not_evaluable_filter_result("reference_resolution_failed")
        if len(pred_ifp) != len(ref_ifp):
            return self._not_evaluable_filter_result("incompatible_vector_lengths")
        if not any(ref_ifp):
            return self._not_evaluable_filter_result("empty_reference")

        overlap = sum(p and r for p, r in zip(pred_ifp, ref_ifp)) / sum(ref_ifp)
        passed = overlap >= self.ifp_filter_threshold
        return {
            "ifp_filter_pass": bool(passed),
            "ifp_filter_status": "accepted" if passed else "rejected",
            "ifp_filter_reason": "threshold_met" if passed else "below_threshold",
            "ifp_filter_overlap": float(overlap),
            "ifp_filter_threshold": self.ifp_filter_threshold,
            "ifp_filter_reference": self.pocket_coverage_reference,
        }

    def _resolve_filter_reference(
        self,
        run_dir: Path,
        chain_df: pd.DataFrame,
        ligand_row: pd.Series,
    ) -> list[int] | None:
        spec = self._ifp_filter_reference_spec
        if spec is None:
            return None
        if spec["kind"] == "bits":
            return list(spec["bits"])

        cif_name = ligand_row.get("cif_file")
        if pd.isna(cif_name) or "cif_file" not in chain_df.columns:
            return None
        model_rows = chain_df[chain_df["cif_file"] == cif_name]
        receptor_rows = model_rows[
            model_rows["ENTITY_TYPE"].astype(str).str.lower() == "protein"
        ]
        receptor_ids = receptor_rows["CHAIN_ID"].dropna().astype(str).unique()
        if len(receptor_ids) != 1:
            return None

        structure_path = run_dir / "results" / "structures" / str(cif_name)
        if not structure_path.exists():
            return None
        try:
            import gemmi

            model = gemmi.read_structure(str(structure_path))[0]
            receptor_chain = model[receptor_ids[0]]
            residues = list(receptor_chain)
            residue_order = [(" ", int(res.seqid.num), " ") for res in residues]
            residue_names = {
                residue_id: str(res.name)
                for residue_id, res in zip(residue_order, residues)
            }
        except Exception as exc:
            self.logger.warning(
                "Unable to resolve filter reference from %s: %s", structure_path, exc
            )
            return None

        return _build_reference_ifp_from_custom(
            spec,
            residue_order=residue_order,
            residue_names=residue_names,
            logger=self.logger,
            warn_key=f"cif={cif_name},chain={ligand_row.get('CHAIN_ID')}",
        )

    @staticmethod
    def _parse_binary_ifp(value: Any) -> tuple[list[int] | None, str]:
        if value is None or (not isinstance(value, (list, tuple)) and pd.isna(value)):
            return None, "missing_ifp"
        parsed = value
        if isinstance(value, str):
            if not value.strip():
                return None, "missing_ifp"
            try:
                parsed = json.loads(value)
            except (TypeError, json.JSONDecodeError):
                return None, "malformed_ifp"
        if not isinstance(parsed, (list, tuple)):
            return None, "malformed_ifp"
        if not parsed:
            return None, "missing_ifp"
        if any(item not in (0, 1, False, True) for item in parsed):
            return None, "malformed_ifp"
        return [int(item) for item in parsed], ""

    @staticmethod
    def _set_filter_dtypes(df: pd.DataFrame) -> None:
        if "ifp_filter_pass" in df.columns:
            df["ifp_filter_pass"] = df["ifp_filter_pass"].astype("boolean")

    def _collect_score_columns(self, run_dir: Path) -> dict[str, Any]:
        """Collect computed scores from per-row validate outputs."""
        out: dict[str, Any] = {}
        try:
            system_df, chain_df = read_metric_frames(run_dir)
            out.update(primary_metric_values(system_df, chain_df))
        except (OSError, UnicodeError, pd.errors.ParserError) as exc:
            self.logger.warning("Failed reading metrics from %s: %s", run_dir, exc)

        self._ensure_screen_metric_schema(out)
        return out

    def _ensure_screen_metric_schema(self, output: dict[str, Any]) -> None:
        """Populate stable manuscript-facing score columns, using nulls when unavailable."""

        for column in SYSTEM_SCREEN_METRICS:
            output.setdefault(f"system__{column}", None)

        chains = list(iter_system_chains(self.base_system_obj))
        protein_ids = [
            chain.chain_id for chain in chains if chain.entity_type == "protein"
        ]
        for chain in chains:
            prefix = f"{chain.entity_type}_{chain.chain_id}"
            metrics: tuple[str, ...] = ()
            if chain.entity_type == "protein":
                metrics = PROTEIN_SCREEN_METRICS
            elif chain.entity_type == "ligand":
                metrics = LIGAND_SCREEN_METRICS
            for column in metrics:
                output.setdefault(f"{prefix}__{column}", None)
            if chain.entity_type == "ligand":
                for protein_id in protein_ids:
                    output.setdefault(
                        f"{prefix}__pair_chains_iptm_{protein_id}",
                        None,
                    )

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
        return [
            int(v.strip()) if v.strip().isdigit() else v.strip()
            for v in str(path_str).split(",")
        ]

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
