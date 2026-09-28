"""Private metric shaping, filtering, and clustering for the Screen recipe."""

from __future__ import annotations

import json
import logging
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Any

import pandas as pd

from cofolder.modules.analytics.ifp_clustering import (
    LEAF_ORDER_COLUMNS as IFP_CLUSTER_LEAF_ORDER_COLUMNS,
)
from cofolder.modules.analytics.ifp_clustering import (
    LINKAGE_COLUMNS as IFP_CLUSTER_LINKAGE_COLUMNS,
)
from cofolder.modules.analytics.ifp_clustering import (
    SUMMARY_COLUMNS as IFP_CLUSTER_SUMMARY_COLUMNS,
)
from cofolder.modules.analytics.ifp_clustering import cluster_interaction_fingerprints
from cofolder.modules.analytics.ifp_filtering import (
    IFPFilterMode,
    ReferenceIFPFilterPolicy,
    evaluate_reference_ifp_filter,
)
from cofolder.modules.analytics.reference_ifp import (
    IFPExtractionConfig,
    IFPSimilarityMetric,
    IFPTaxonomy,
    InteractionFingerprint,
    LigandSelector,
    ReferenceEntitySelectionError,
    ReferenceIFPError,
    compare_interaction_fingerprints,
    extract_interaction_fingerprint,
    map_reference_identities,
)
from cofolder.modules.analytics.reproduction import _build_reference_ifp_from_custom
from cofolder.modules.contracts import (
    SCREEN_METRIC_PROFILES,
    ScreenExecutionCardinalityError,
)
from cofolder.modules.input.system import iter_system_chains
from cofolder.modules.runners import PlannedExecution, RunnerExecutionPlan
from cofolder.recipes._metrics import primary_metric_values, read_metric_frames

SYSTEM_SCREEN_METRICS = SCREEN_METRIC_PROFILES["system"]
PROTEIN_SCREEN_METRICS = SCREEN_METRIC_PROFILES["protein"]
LIGAND_SCREEN_METRICS = SCREEN_METRIC_PROFILES["ligand"]


@dataclass(frozen=True, slots=True)
class ScreenPostprocessConfig:
    """Resolved, immutable settings consumed after Screen input validation."""

    wrk_dir: Path
    ligand_chain: str
    base_system_obj: Any
    col_id: str
    cluster_ifps: bool
    ifp_cluster_similarity_threshold: float
    ifp_filter_threshold: float | None
    ifp_taxonomy: IFPTaxonomy
    pocket_coverage_reference: str | None
    resolved_ifp_filter_source: str | None
    ifp_filter_reference_spec: Mapping[str, Any] | None
    ifp_filter_policy_config: ReferenceIFPFilterPolicy | None
    reference_ligand_selector: LigandSelector | None
    ifp_reference_receptor_chains: tuple[str, ...]
    validate_kwargs: Mapping[str, Any]
    execution_plan: RunnerExecutionPlan | None

    @classmethod
    def from_screen(cls, screen: Any) -> "ScreenPostprocessConfig":
        reference_spec = (
            MappingProxyType(dict(screen._ifp_filter_reference_spec))
            if screen._ifp_filter_reference_spec is not None
            else None
        )
        return cls(
            wrk_dir=screen.wrk_dir,
            ligand_chain=screen.ligand_chain,
            base_system_obj=screen.base_system_obj,
            col_id=screen.col_id,
            cluster_ifps=screen.cluster_ifps,
            ifp_cluster_similarity_threshold=screen.ifp_cluster_similarity_threshold,
            ifp_filter_threshold=screen.ifp_filter_threshold,
            ifp_taxonomy=screen.ifp_taxonomy,
            pocket_coverage_reference=screen.pocket_coverage_reference,
            resolved_ifp_filter_source=screen._resolved_ifp_filter_source,
            ifp_filter_reference_spec=reference_spec,
            ifp_filter_policy_config=screen._ifp_filter_policy_config,
            reference_ligand_selector=screen._reference_ligand_selector,
            ifp_reference_receptor_chains=screen.ifp_reference_receptor_chains,
            validate_kwargs=MappingProxyType(dict(screen.validate_kwargs)),
            execution_plan=screen.execution_plan,
        )


class ScreenPostprocessor:
    """Apply Screen-owned metric shaping and optional IFP decisions."""

    _CONFIG_ALIASES = {
        "_resolved_ifp_filter_source": "resolved_ifp_filter_source",
        "_ifp_filter_reference_spec": "ifp_filter_reference_spec",
        "_ifp_filter_policy_config": "ifp_filter_policy_config",
        "_reference_ligand_selector": "reference_ligand_selector",
    }

    def __init__(
        self,
        config: ScreenPostprocessConfig,
        *,
        logger: logging.Logger,
    ) -> None:
        self.config = config
        self.logger = logger
        self._reference_ifp: InteractionFingerprint | None = None
        self._reference_ifp_failure: str | None = None
        self._prediction_ifp_cache: dict[
            tuple[Path, int | None, int | None],
            tuple[InteractionFingerprint | None, str],
        ] = {}
        self._row_ligand_chains: dict[Path, str] = {}

    def register_row_ligand_chain(self, run_dir: Path, chain_id: str) -> None:
        """Record the analysis ligand selected for one mapped system row."""

        self._row_ligand_chains[Path(run_dir)] = str(chain_id)

    def _ligand_chain_for(self, run_dir: Path) -> str:
        return self._row_ligand_chains.get(Path(run_dir), self.ligand_chain)

    def __getattr__(self, name: str) -> Any:
        target = self._CONFIG_ALIASES.get(name, name)
        fields = ScreenPostprocessConfig.__dataclass_fields__
        if target in fields:
            return getattr(self.config, target)
        raise AttributeError(name)

    def _ifp_filter_enabled(self) -> bool:
        return self._ifp_filter_policy_config is not None

    def _prepare_reference_ifp(self) -> None:
        if (
            not self._ifp_filter_enabled()
            or self._resolved_ifp_filter_source != "reference_complex"
        ):
            return
        reference_path = Path(self.validate_kwargs["reference_path"])
        try:
            self._reference_ifp = extract_interaction_fingerprint(
                reference_path,
                ligand=self._reference_ligand_selector,
                receptor_chains=self.ifp_reference_receptor_chains or None,
                config=IFPExtractionConfig(taxonomy=self.ifp_taxonomy),
            )
        except ReferenceEntitySelectionError as exc:
            self._reference_ifp_failure = str(exc)
            return
        except ReferenceIFPError as exc:
            if self.ifp_taxonomy is not IFPTaxonomy.PROLIF:
                raise
            message = str(exc)
            self._reference_ifp_failure = (
                "prolif_worker_timeout"
                if "timeout" in message
                else "prolif_worker_crashed"
            )
            return
        try:
            policy = self._ifp_filter_policy_config
            if policy is not None and policy.required_interactions:
                unknown = (
                    policy.required_interactions - self._reference_ifp.interactions
                )
                if unknown:
                    raise ValueError(
                        "Required interactions are absent from the reference fingerprint: "
                        + ", ".join(str(item) for item in sorted(unknown))
                    )
        except AttributeError as exc:
            raise ValueError("Reference IFP was not prepared.") from exc

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
            for name in (
                "ifp_cluster_summary.csv",
                "ifp_cluster_linkage.csv",
                "ifp_cluster_leaf_order.csv",
            ):
                (self.wrk_dir / "results" / name).unlink(missing_ok=True)
            return

        parsed_rows: list[tuple[int, InteractionFingerprint, str]] = []
        for position, row in summary_df.iterrows():
            if row.get("status") != "success":
                continue
            fingerprint, _ = self._load_selected_interaction_fingerprint(
                Path(row["run_dir"]),
                repeat_id=int(row["repeat_id"]),
                sample_id=(
                    int(row["sample_id"]) if pd.notna(row.get("sample_id")) else None
                ),
            )
            if fingerprint is None:
                continue
            parsed_rows.append((position, fingerprint, str(row["execution_key"])))

        if parsed_rows:
            clustered = cluster_interaction_fingerprints(
                [row[1] for row in parsed_rows],
                [row[2] for row in parsed_rows],
                similarity_threshold=self.ifp_cluster_similarity_threshold,
            )
            for (position, _, _), cluster_id in zip(parsed_rows, clustered.cluster_ids):
                for frame in (summary_df, results_df):
                    frame.at[position, "ifp_cluster_id"] = cluster_id
                    frame.at[position, "ifp_cluster_status"] = "clustered"
            cluster_summary = clustered.summary.copy()
            cluster_linkage = pd.DataFrame(
                [
                    {
                        "merge_index": len(clustered.member_ids) + offset,
                        "left_child": int(row[0]),
                        "right_child": int(row[1]),
                        "jaccard_distance": row[2],
                        "member_count": int(row[3]),
                    }
                    for offset, row in enumerate(clustered.linkage_matrix)
                ],
                columns=IFP_CLUSTER_LINKAGE_COLUMNS,
            )
            cluster_leaf_order = pd.DataFrame(
                [
                    {
                        "leaf_position": position,
                        "input_index": input_index,
                        "member_id": clustered.member_ids[input_index],
                        "ifp_cluster_id": clustered.cluster_ids[input_index],
                    }
                    for position, input_index in enumerate(clustered.leaf_indices)
                ],
                columns=IFP_CLUSTER_LEAF_ORDER_COLUMNS,
            )
            cluster_summary["ifp_taxonomy"] = self.ifp_taxonomy.value
            if self._ifp_filter_enabled():
                annotations = []
                for cluster_id in cluster_summary["ifp_cluster_id"]:
                    members = summary_df[summary_df["ifp_cluster_id"] == cluster_id]
                    similarities = pd.to_numeric(
                        members.get("ifp_filter_similarity", pd.Series(dtype=float)),
                        errors="coerce",
                    ).dropna()
                    accepted = members[
                        members.get(
                            "ifp_filter_status",
                            pd.Series(index=members.index, dtype=object),
                        )
                        == "accepted"
                    ]
                    annotations.append(
                        {
                            "ifp_cluster_id": cluster_id,
                            "reference_evaluable_count": int(similarities.size),
                            "reference_accepted_count": len(accepted),
                            "reference_accepted_member_ids": json.dumps(
                                accepted[self.col_id].astype(str).tolist()
                            ),
                            "mean_reference_similarity": (
                                float(similarities.mean())
                                if not similarities.empty
                                else None
                            ),
                            "max_reference_similarity": (
                                float(similarities.max())
                                if not similarities.empty
                                else None
                            ),
                        }
                    )
                cluster_summary = cluster_summary.merge(
                    pd.DataFrame(annotations), on="ifp_cluster_id", how="left"
                )
        else:
            cluster_summary = pd.DataFrame(columns=IFP_CLUSTER_SUMMARY_COLUMNS)
            cluster_linkage = pd.DataFrame(columns=IFP_CLUSTER_LINKAGE_COLUMNS)
            cluster_leaf_order = pd.DataFrame(columns=IFP_CLUSTER_LEAF_ORDER_COLUMNS)

        public_results_dir = self.wrk_dir / "results"
        public_results_dir.mkdir(parents=True, exist_ok=True)
        cluster_summary.to_csv(
            public_results_dir / "ifp_cluster_summary.csv", index=False
        )
        cluster_linkage.to_csv(
            public_results_dir / "ifp_cluster_linkage.csv", index=False
        )
        cluster_leaf_order.to_csv(
            public_results_dir / "ifp_cluster_leaf_order.csv", index=False
        )

    def _publish_prolif_events(self, summary_df: pd.DataFrame) -> None:
        """Publish occurrence-level ProLIF events with stable Screen identities."""

        path = self.wrk_dir / "results" / "ifp_interaction_events.jsonl"
        if self.ifp_taxonomy is not IFPTaxonomy.PROLIF or not (
            self.cluster_ifps or self._ifp_filter_enabled()
        ):
            path.unlink(missing_ok=True)
            return

        records: list[dict[str, object]] = []
        for _, row in summary_df.iterrows():
            if row.get("status") != "success":
                continue
            fingerprint, _ = self._load_selected_interaction_fingerprint(
                Path(row["run_dir"]),
                repeat_id=int(row["repeat_id"]),
                sample_id=(
                    int(row["sample_id"]) if pd.notna(row.get("sample_id")) else None
                ),
            )
            if fingerprint is None:
                continue
            identity = {
                "compound_id": str(row[self.col_id]),
                "execution_key": str(row["execution_key"]),
                "repeat_id": int(row["repeat_id"]),
                "model_id": str(row["model_id"]),
                "sample_id": (
                    int(row["sample_id"]) if pd.notna(row.get("sample_id")) else None
                ),
            }
            records.extend(
                {**identity, **event} for event in fingerprint.serialized_events()
            )

        path.parent.mkdir(parents=True, exist_ok=True)
        temporary_path = path.with_name(f".{path.name}.tmp")
        try:
            temporary_path.write_text(
                "".join(
                    json.dumps(record, sort_keys=True) + "\n" for record in records
                ),
                encoding="utf-8",
            )
            temporary_path.replace(path)
        finally:
            temporary_path.unlink(missing_ok=True)

    def _load_selected_interaction_fingerprint(
        self,
        run_dir: Path,
        *,
        repeat_id: int | None = None,
        sample_id: int | None = None,
    ) -> tuple[InteractionFingerprint | None, str]:
        run_dir = Path(run_dir)
        ligand_chain = self._ligand_chain_for(run_dir)
        cache_key = (run_dir, repeat_id, sample_id)
        if cache_key in self._prediction_ifp_cache:
            return self._prediction_ifp_cache[cache_key]
        _, chain_df = read_metric_frames(run_dir)
        if chain_df.empty or not {"CHAIN_ID", "ENTITY_TYPE", "cif_file"}.issubset(
            chain_df.columns
        ):
            result = (None, "missing_chain_metrics")
            self._prediction_ifp_cache[cache_key] = result
            return result
        selected = self._select_chain_metrics_rows(
            chain_df, repeat_id=repeat_id, sample_id=sample_id
        )
        ligand_rows = selected[
            (selected["ENTITY_TYPE"].astype(str).str.lower() == "ligand")
            & (selected["CHAIN_ID"].astype(str) == ligand_chain)
        ]
        if ligand_rows.empty:
            result = (None, "ligand_chain_not_found")
            self._prediction_ifp_cache[cache_key] = result
            return result
        cif_name = ligand_rows.iloc[0].get("cif_file")
        path = run_dir / "results" / "structures" / str(cif_name)
        receptor_chains = tuple(
            selected.loc[
                selected["ENTITY_TYPE"].astype(str).str.lower() == "protein", "CHAIN_ID"
            ]
            .dropna()
            .astype(str)
            .unique()
        )
        try:
            result = (
                extract_interaction_fingerprint(
                    path,
                    ligand=LigandSelector(chain_id=ligand_chain),
                    receptor_chains=receptor_chains or None,
                    config=IFPExtractionConfig(taxonomy=self.ifp_taxonomy),
                ),
                "",
            )
        except ReferenceIFPError as exc:
            result = (None, str(exc))
        self._prediction_ifp_cache[cache_key] = result
        return result

    def _default_filter_result(self) -> dict[str, Any]:
        if self._ifp_filter_enabled():
            return self._not_evaluable_filter_result("missing_ifp")
        return {
            "ifp_filter_pass": pd.NA,
            "ifp_filter_status": "not_applied",
            "ifp_filter_reason": "filtering_disabled",
            "ifp_filter_overlap": None,
            "ifp_filter_similarity": None,
            "ifp_filter_threshold": None,
            "ifp_filter_reference": None,
            "ifp_filter_similarity_metric": None,
            "ifp_filter_policy": None,
            "ifp_filter_taxonomy": None,
            "ifp_filter_required_interactions": None,
            "ifp_filter_missing_interactions": None,
            "ifp_filter_mapping_status": None,
            "ifp_filter_mapping_failures": None,
        }

    def _not_evaluable_filter_result(self, reason: str) -> dict[str, Any]:
        return {
            "ifp_filter_pass": pd.NA,
            "ifp_filter_status": "not_evaluable",
            "ifp_filter_reason": reason,
            "ifp_filter_overlap": None,
            "ifp_filter_similarity": None,
            "ifp_filter_threshold": self.ifp_filter_threshold,
            "ifp_filter_reference": (
                str(self.validate_kwargs.get("reference_path"))
                if self._resolved_ifp_filter_source == "reference_complex"
                else self.pocket_coverage_reference
            ),
            "ifp_filter_similarity_metric": (
                self._ifp_filter_policy_config.similarity_metric.value
                if self._ifp_filter_policy_config
                else None
            ),
            "ifp_filter_policy": (
                self._ifp_filter_policy_config.mode.value
                if self._ifp_filter_policy_config
                else None
            ),
            "ifp_filter_taxonomy": self.ifp_taxonomy.value,
            "ifp_filter_required_interactions": (
                json.dumps(
                    [
                        str(item)
                        for item in sorted(
                            self._ifp_filter_policy_config.required_interactions or ()
                        )
                    ]
                )
                if self._ifp_filter_policy_config
                else None
            ),
            "ifp_filter_missing_interactions": json.dumps([]),
            "ifp_filter_mapping_status": "unmappable",
            "ifp_filter_mapping_failures": json.dumps([reason]),
        }

    def _evaluate_ifp_filter(
        self,
        run_dir: Path,
        *,
        repeat_id: int | None = None,
        sample_id: int | None = None,
    ) -> dict[str, Any]:
        """Evaluate strict reference overlap for one completed execution."""
        if not self._ifp_filter_enabled():
            return self._default_filter_result()
        if self._resolved_ifp_filter_source == "reference_complex":
            return self._evaluate_reference_complex_filter(
                run_dir, repeat_id=repeat_id, sample_id=sample_id
            )

        _, chain_df = read_metric_frames(run_dir)
        ligand_chain = self._ligand_chain_for(run_dir)
        if chain_df.empty:
            return self._not_evaluable_filter_result("missing_chain_metrics")
        required = {"CHAIN_ID", "ENTITY_TYPE", "ifp_distance"}
        if chain_df.empty or not required.issubset(chain_df.columns):
            return self._not_evaluable_filter_result("missing_ifp")

        chain_df = self._select_chain_metrics_rows(
            chain_df, repeat_id=repeat_id, sample_id=sample_id
        )
        ligand_rows = chain_df[
            chain_df["ENTITY_TYPE"].astype(str).str.lower() == "ligand"
        ]
        if ligand_chain:
            ligand_rows = ligand_rows[
                ligand_rows["CHAIN_ID"].astype(str) == ligand_chain
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

        intersection = sum(p and r for p, r in zip(pred_ifp, ref_ifp))
        overlap = intersection / sum(ref_ifp)
        union = sum(bool(p or r) for p, r in zip(pred_ifp, ref_ifp))
        jaccard = intersection / union if union else 0.0
        assert self._ifp_filter_policy_config is not None
        score = (
            jaccard
            if self._ifp_filter_policy_config.similarity_metric
            is IFPSimilarityMetric.JACCARD
            else overlap
        )
        passed = score >= float(self.ifp_filter_threshold)
        return {
            "ifp_filter_pass": bool(passed),
            "ifp_filter_status": "accepted" if passed else "rejected",
            "ifp_filter_reason": "threshold_met" if passed else "below_threshold",
            "ifp_filter_overlap": float(overlap),
            "ifp_filter_similarity": float(score),
            "ifp_filter_threshold": self.ifp_filter_threshold,
            "ifp_filter_reference": self.pocket_coverage_reference,
            "ifp_filter_similarity_metric": self._ifp_filter_policy_config.similarity_metric.value,
            "ifp_filter_policy": self._ifp_filter_policy_config.mode.value,
            "ifp_filter_taxonomy": self.ifp_taxonomy.value,
            "ifp_filter_required_interactions": json.dumps([]),
            "ifp_filter_missing_interactions": json.dumps([]),
            "ifp_filter_mapping_status": "not_applicable",
            "ifp_filter_mapping_failures": json.dumps([]),
        }

    def _evaluate_reference_complex_filter(
        self,
        run_dir: Path,
        *,
        repeat_id: int | None = None,
        sample_id: int | None = None,
    ) -> dict[str, Any]:
        if self._reference_ifp_failure:
            return self._not_evaluable_filter_result(self._reference_ifp_failure)
        if self._reference_ifp is None or self._ifp_filter_policy_config is None:
            return self._not_evaluable_filter_result("reference_ifp_extraction_failed")
        try:
            prediction, extraction_reason = self._load_selected_interaction_fingerprint(
                run_dir, repeat_id=repeat_id, sample_id=sample_id
            )
            if prediction is None:
                reason = (
                    "prolif_worker_timeout"
                    if "timeout" in extraction_reason
                    else (
                        "prolif_worker_crashed"
                        if "crashed" in extraction_reason
                        else extraction_reason or "prediction_ifp_extraction_failed"
                    )
                )
                return self._not_evaluable_filter_result(reason)
            predicted_path = prediction.source_path
            mapping = map_reference_identities(
                self._reference_ifp,
                prediction,
                reference_structure_path=Path(self.validate_kwargs["reference_path"]),
                predicted_structure_path=predicted_path,
            )
            comparison = compare_interaction_fingerprints(
                self._reference_ifp, prediction, mapping
            )
            outcome = evaluate_reference_ifp_filter(
                self._reference_ifp, comparison, self._ifp_filter_policy_config
            )
        except ReferenceIFPError as exc:
            message = str(exc)
            reason = (
                "prolif_worker_timeout"
                if "prolif_worker_timeout" in message
                else (
                    "prolif_worker_crashed"
                    if "prolif_worker_crashed" in message
                    else "prediction_ifp_extraction_failed"
                )
            )
            return self._not_evaluable_filter_result(reason)

        selected_score = comparison.similarities.get(
            self._ifp_filter_policy_config.similarity_metric
        )
        return {
            "ifp_filter_pass": outcome.passed if outcome.passed is not None else pd.NA,
            "ifp_filter_status": outcome.status,
            "ifp_filter_reason": outcome.reason,
            "ifp_filter_overlap": comparison.similarities.get(
                IFPSimilarityMetric.REFERENCE_COVERAGE
            ),
            "ifp_filter_similarity": selected_score,
            "ifp_filter_threshold": self._ifp_filter_policy_config.threshold,
            "ifp_filter_reference": str(self.validate_kwargs["reference_path"]),
            "ifp_filter_similarity_metric": self._ifp_filter_policy_config.similarity_metric.value,
            "ifp_filter_policy": self._ifp_filter_policy_config.mode.value,
            "ifp_filter_taxonomy": self.ifp_taxonomy.value,
            "ifp_filter_required_interactions": json.dumps(
                [
                    str(item)
                    for item in sorted(
                        self._ifp_filter_policy_config.required_interactions
                        or self._reference_ifp.interactions
                        if self._ifp_filter_policy_config.mode is IFPFilterMode.REQUIRED
                        else ()
                    )
                ]
            ),
            "ifp_filter_missing_interactions": json.dumps(
                [str(item) for item in comparison.missing_interactions]
            ),
            "ifp_filter_mapping_status": comparison.mapping.status.value,
            "ifp_filter_mapping_failures": json.dumps(
                list(comparison.mapping.failures)
            ),
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

    def _collect_execution_score_columns(
        self,
        system_df: pd.DataFrame,
        chain_df: pd.DataFrame,
        slot: PlannedExecution,
    ) -> tuple[dict[str, Any], bool]:
        """Collect metrics for exactly one planned repeat/model/sample slot."""

        def select(frame: pd.DataFrame) -> pd.DataFrame:
            if frame.empty:
                return frame
            mask = pd.Series(True, index=frame.index)
            if "repeat" in frame:
                mask &= pd.to_numeric(frame["repeat"], errors="coerce").eq(
                    slot.repeat_id
                )
            if "diffusion_sample" in frame and slot.sample_id is not None:
                mask &= pd.to_numeric(frame["diffusion_sample"], errors="coerce").eq(
                    slot.sample_id
                )
            if "model_name" in frame:
                observed_models = set(frame["model_name"].dropna().astype(str))
                if slot.model_id in observed_models:
                    mask &= frame["model_name"].astype(str).eq(slot.model_id)
            return frame.loc[mask]

        selected_system = select(system_df)
        selected_chain = select(chain_df)
        if len(selected_system) > 1:
            raise ScreenExecutionCardinalityError(
                "Validate produced duplicate system rows for execution "
                f"repeat={slot.repeat_id}, model={slot.model_id!r}, "
                f"sample={slot.sample_id!r}."
            )
        observed = not selected_system.empty or not selected_chain.empty
        values = (
            primary_metric_values(selected_system, selected_chain) if observed else {}
        )
        self._ensure_screen_metric_schema(values)
        return values, observed

    def _validate_observed_execution_keys(self, system_df: pd.DataFrame) -> None:
        """Reject normalized system outputs outside the planned execution matrix."""
        if system_df.empty or self.execution_plan is None:
            return
        required = {"repeat", "diffusion_sample"}
        if not required.issubset(system_df.columns):
            raise ScreenExecutionCardinalityError(
                "Validate system output lacks repeat/diffusion_sample identity."
            )
        planned = {
            (slot.repeat_id, slot.sample_id) for slot in self.execution_plan.executions
        }
        execution_rows = system_df.dropna(subset=["repeat", "diffusion_sample"])
        observed = {
            (int(row["repeat"]), int(row["diffusion_sample"]))
            for _, row in execution_rows.iterrows()
        }
        extras = sorted(observed - planned)
        if extras:
            raise ScreenExecutionCardinalityError(
                f"Validate produced executions outside the planned matrix: {extras}."
            )

    def _ensure_screen_metric_schema(
        self, output: dict[str, Any], system_obj: Any | None = None
    ) -> None:
        """Populate stable manuscript-facing score columns, using nulls when unavailable."""

        for column in SYSTEM_SCREEN_METRICS:
            output.setdefault(f"system__{column}", None)

        chains = list(iter_system_chains(system_obj or self.base_system_obj))
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
    def _select_chain_metrics_rows(
        df: pd.DataFrame,
        *,
        repeat_id: int | None = None,
        sample_id: int | None = None,
    ) -> pd.DataFrame:
        """Select one execution's chain rows, retaining legacy primary selection."""
        if repeat_id is not None and "repeat" in df.columns:
            selected = df[pd.to_numeric(df["repeat"], errors="coerce").eq(repeat_id)]
            if sample_id is not None and "diffusion_sample" in selected.columns:
                selected = selected[
                    pd.to_numeric(selected["diffusion_sample"], errors="coerce").eq(
                        sample_id
                    )
                ]
            return selected
        if {"repeat", "diffusion_sample"}.issubset(df.columns):
            sub = df[(df["repeat"] == 1) & (df["diffusion_sample"] == 0)]
            if not sub.empty:
                return sub
        if "CHAIN_ID" in df.columns:
            return df.drop_duplicates(subset=["CHAIN_ID"], keep="first")
        return df
