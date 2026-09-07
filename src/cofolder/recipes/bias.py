from __future__ import annotations

import logging
from contextlib import nullcontext
from datetime import date
from pathlib import Path
from uuid import uuid4

import pandas as pd

from cofolder.modules.analytics.bias import apply_bias_metrics
from cofolder.modules.analytics.bias_training import run_build_bias_training_data
from cofolder.modules.contracts import (
    PUBLIC_SCHEMA_VERSION,
    ArtifactReference,
    EvidenceSource,
    FailureStage,
    OutputIdentity,
    PublicManifest,
    PublicOutputBundle,
    WorkflowExecutionError,
    WorkflowKind,
    bundle_from_frames,
    failure_from_exception,
    write_public_bundle,
)
from cofolder.modules.input import system
from cofolder.modules.utils import read
from cofolder.modules.utils.timing import DebugTimingCollector

logger = logging.getLogger(__name__)


def normalize_bias_chains(
    bias_chains: set[str] | list[str] | tuple[str, ...] | None,
) -> set[str] | None:
    normalized: set[str] = set()
    for value in bias_chains or []:
        for token in str(value).split(","):
            token = token.strip()
            if token:
                normalized.add(token.upper())
    return normalized or None


def _iter_chain_ids(raw_ids) -> list[str]:
    if raw_ids is None:
        return []
    if not isinstance(raw_ids, list):
        raw_ids = [raw_ids]
    return [str(chain_id).strip() for chain_id in raw_ids if str(chain_id).strip()]


def _resolve_ligand_molecule_id(ligand_data: dict) -> str:
    ccd_value = ligand_data.get("ccd")
    if ccd_value is None:
        ccd_ids = []
    else:
        ccd_ids = ccd_value if isinstance(ccd_value, list) else [ccd_value]
    for ccd_id in ccd_ids:
        ccd_text = str(ccd_id).strip()
        if ccd_text:
            return ccd_text.upper()

    smiles = str(ligand_data.get("smiles", "")).strip()
    if smiles:
        return smiles
    return "UNKNOWN_LIGAND"


def build_bias_dataframes(
    sys_obj: system.System,
    *,
    system_name: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    system_df = pd.DataFrame(
        [
            {
                "model_name": system_name,
                "repeat": 1,
                "diffusion_sample": 0,
            }
        ]
    )

    chain_rows: list[dict[str, object]] = []
    sequences = sys_obj.find_value(key="sequences") or []
    conf_chain_id = 0

    for entity_position, seq_entry in enumerate(sequences):
        if not isinstance(seq_entry, dict):
            continue
        if "protein" in seq_entry:
            protein_data = seq_entry.get("protein") or {}
            for chain_id in _iter_chain_ids(protein_data.get("id")):
                chain_rows.append(
                    {
                        "conf_chain_id": conf_chain_id,
                        "CHAIN_ID": chain_id,
                        "ENTITY_ID": f"entity:{entity_position}",
                        "ENTITY_TYPE": "protein",
                        "ligand_molecule_id": f"protein_{chain_id}",
                        "model_name": system_name,
                        "repeat": 1,
                        "diffusion_sample": 0,
                    }
                )
                conf_chain_id += 1

        if "ligand" in seq_entry:
            ligand_data = seq_entry.get("ligand") or {}
            molecule_id = _resolve_ligand_molecule_id(ligand_data)
            for chain_id in _iter_chain_ids(ligand_data.get("id")):
                chain_rows.append(
                    {
                        "conf_chain_id": conf_chain_id,
                        "CHAIN_ID": chain_id,
                        "ENTITY_ID": f"entity:{entity_position}",
                        "ENTITY_TYPE": "ligand",
                        "ligand_molecule_id": molecule_id,
                        "model_name": system_name,
                        "repeat": 1,
                        "diffusion_sample": 0,
                    }
                )
                conf_chain_id += 1

    chain_df = pd.DataFrame(
        chain_rows,
        columns=[
            "conf_chain_id",
            "CHAIN_ID",
            "ENTITY_ID",
            "ENTITY_TYPE",
            "ligand_molecule_id",
            "model_name",
            "repeat",
            "diffusion_sample",
        ],
    )
    return system_df, chain_df


def find_invalid_ligand_smiles_chain_ids(sys_obj: system.System) -> set[str]:
    try:
        from rdkit import Chem
    except Exception:
        return set()

    invalid: set[str] = set()
    sequences = sys_obj.find_value(key="sequences") or []
    for entry in sequences:
        if not isinstance(entry, dict) or "ligand" not in entry:
            continue
        ligand_data = entry.get("ligand") or {}
        smiles = str(ligand_data.get("smiles", "")).strip()
        if not smiles:
            continue
        if Chem.MolFromSmiles(smiles) is not None:
            continue
        invalid.update(
            chain_id.upper() for chain_id in _iter_chain_ids(ligand_data.get("id"))
        )
    return invalid


class BiasAssessmentWorkflow:
    def __init__(
        self,
        *,
        wrk_dir: Path,
        system_path: Path,
        sys_obj: system.System,
        protein_training_data_path: Path | None,
        ligand_training_data_path: Path | None,
        bias_release_cutoff: str,
        bias_ligand_similarity_threshold: float,
        bias_chains: set[str] | list[str] | None,
        build_bias_training_data: bool,
        bias_training_components_cif: Path | None,
        custom_protein_reference_path: Path | None = None,
        custom_ligand_reference_path: Path | None = None,
        logger: logging.Logger | None = None,
        timings: DebugTimingCollector | None = None,
        strict_training_data: bool = False,
    ):
        self.wrk_dir = Path(wrk_dir)
        self.system_path = Path(system_path)
        self.run_id = str(uuid4())
        self.sys = sys_obj
        self.protein_training_data_path = protein_training_data_path
        self.ligand_training_data_path = ligand_training_data_path
        self.custom_protein_reference_path = custom_protein_reference_path
        self.custom_ligand_reference_path = custom_ligand_reference_path
        self.bias_release_cutoff = str(bias_release_cutoff)
        self.bias_ligand_similarity_threshold = float(bias_ligand_similarity_threshold)
        self.bias_chains = normalize_bias_chains(bias_chains)
        self.build_bias_training_data = bool(build_bias_training_data)
        self.bias_training_components_cif = bias_training_components_cif
        self.logger = logger or logging.getLogger(__name__)
        self.timings = timings
        self.strict_training_data = strict_training_data

    @property
    def output_dir(self) -> Path:
        return self.wrk_dir / "results" / "bias_train"

    def _debug_timer(self, label: str):
        if self.timings is None:
            return nullcontext()
        return self.timings.measure(label, logger=self.logger)

    def _validate_release_cutoff(self) -> None:
        try:
            date.fromisoformat(self.bias_release_cutoff)
        except ValueError as exc:
            raise ValueError(
                "--bias_release_cutoff must be a valid ISO date (YYYY-MM-DD)."
            ) from exc

    def _resolve_components_cif_path(self) -> Path | None:
        if self.bias_training_components_cif is not None:
            return self.bias_training_components_cif
        if self.protein_training_data_path is None:
            return None
        default_path = self.protein_training_data_path.parent / "ccd" / "components.cif"
        if default_path.exists() and default_path.is_file():
            return default_path
        return None

    def apply(
        self,
        *,
        system_df: pd.DataFrame,
        chain_df: pd.DataFrame,
        boltz_cache_path: Path | str | None,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        self._validate_release_cutoff()

        protein_metrics_path = self.protein_training_data_path
        ligand_metrics_path = self.ligand_training_data_path

        if self.build_bias_training_data:
            protein_metrics_path, ligand_metrics_path = self._build_training_data(
                chain_df=chain_df,
                ligand_metrics_path=ligand_metrics_path,
            )

        scoped_chain_df = chain_df
        if self.bias_chains:
            scoped_chain_df = chain_df[
                chain_df["CHAIN_ID"]
                .astype(str)
                .str.strip()
                .str.upper()
                .isin(self.bias_chains)
            ].copy()

        needs_protein = bool(
            (
                scoped_chain_df.get("ENTITY_TYPE", pd.Series(dtype=object)).astype(str)
                == "protein"
            ).any()
        )
        needs_ligand = bool(
            (
                scoped_chain_df.get("ENTITY_TYPE", pd.Series(dtype=object)).astype(str)
                == "ligand"
            ).any()
        )

        has_protein_file = (
            protein_metrics_path is not None
            and protein_metrics_path.exists()
            and protein_metrics_path.is_file()
        )
        has_ligand_file = (
            ligand_metrics_path is not None
            and ligand_metrics_path.exists()
            and ligand_metrics_path.is_file()
        )
        has_protein_source = (
            has_protein_file or self.custom_protein_reference_path is not None
        )
        has_ligand_source = (
            has_ligand_file or self.custom_ligand_reference_path is not None
        )

        missing_sources: list[str] = []
        if needs_protein and not has_protein_source:
            missing_sources.append("protein")
        if needs_ligand and not has_ligand_source:
            missing_sources.append("ligand")

        if not missing_sources:
            with self._debug_timer("scores.bias_metrics.total"):
                return apply_bias_metrics(
                    system_df=system_df,
                    chain_df=chain_df,
                    sys_obj=self.sys,
                    protein_training_data_path=protein_metrics_path,
                    ligand_training_data_path=ligand_metrics_path,
                    custom_protein_reference_path=self.custom_protein_reference_path,
                    custom_ligand_reference_path=self.custom_ligand_reference_path,
                    release_cutoff=self.bias_release_cutoff,
                    bias_chains=self.bias_chains,
                    protein_top_n=100,
                    boltz_cache_path=boltz_cache_path or "~/.boltz",
                    components_cif_path=self._resolve_components_cif_path(),
                    output_dir=self.output_dir,
                    logger=self.logger,
                    timings=self.timings,
                    strict_query_resolution=self.strict_training_data,
                )

        message = (
            "Bias assessment is missing required reference sources for the selected query entities: "
            f"{', '.join(missing_sources)}. "
            f"Resolved public protein={protein_metrics_path} public ligand={ligand_metrics_path} "
            f"custom protein={self.custom_protein_reference_path} custom ligand={self.custom_ligand_reference_path}."
        )
        if self.strict_training_data:
            raise ValueError(message)
        self.logger.warning("%s Skipping bias metrics.", message)
        return system_df, chain_df

    def _build_training_data(
        self,
        *,
        chain_df: pd.DataFrame,
        ligand_metrics_path: Path | None,
    ) -> tuple[Path, Path]:
        if self.protein_training_data_path is None:
            raise ValueError(
                "--build_bias_training_data requires --protein_training_data_path."
            )

        if ligand_metrics_path is None:
            ligand_metrics_path = self.output_dir / "ligand_training_data.csv"
            ligand_metrics_path.parent.mkdir(parents=True, exist_ok=True)

        components_cif = self.bias_training_components_cif
        if components_cif is None:
            components_cif = (
                self.protein_training_data_path.parent / "ccd" / "components.cif"
            )
        if not components_cif.exists():
            raise ValueError(
                f"components.cif not found for bias-training build: {components_cif}. "
                "Provide --bias_training_components_cif or prepare the training-data root."
            )

        available_protein_chains = {
            str(row["CHAIN_ID"]).strip().upper()
            for _, row in chain_df.iterrows()
            if str(row.get("ENTITY_TYPE")) == "protein"
            and str(row.get("CHAIN_ID", "")).strip()
        }
        available_ligand_chains = {
            str(row["CHAIN_ID"]).strip().upper()
            for _, row in chain_df.iterrows()
            if str(row.get("ENTITY_TYPE")) == "ligand"
            and str(row.get("CHAIN_ID", "")).strip()
        }
        selected_chains = self.bias_chains
        selected_protein_chains = (
            available_protein_chains & selected_chains
            if selected_chains is not None
            else available_protein_chains
        )
        selected_ligand_chains = (
            available_ligand_chains & selected_chains
            if selected_chains is not None
            else available_ligand_chains
        )

        invalid_smiles_chains = find_invalid_ligand_smiles_chain_ids(self.sys)
        invalid_selected_ligand_chains = selected_ligand_chains & invalid_smiles_chains
        if invalid_selected_ligand_chains:
            self.logger.warning(
                "Invalid ligand SMILES detected for chains=%s. "
                "Skipping ligand ECFP protocol for this bias-build run.",
                sorted(invalid_selected_ligand_chains),
            )
            selected_ligand_chains = (
                selected_ligand_chains - invalid_selected_ligand_chains
            )

        run_protein_protocol = bool(selected_protein_chains)
        run_ligand_protocol = bool(selected_ligand_chains)
        build_protein_training_path = self.protein_training_data_path
        build_protein_training_path.parent.mkdir(parents=True, exist_ok=True)
        if not run_protein_protocol:
            build_protein_training_path = (
                self.output_dir / "_protein_training_data_build_tmp.csv"
            )
            build_protein_training_path.parent.mkdir(parents=True, exist_ok=True)

        self.logger.info(
            "Bias build protocol selection: run_protein=%s chains=%s | "
            "run_ligand=%s chains=%s | protein_build_csv=%s | ligand_build_csv=%s",
            run_protein_protocol,
            sorted(selected_protein_chains),
            run_ligand_protocol,
            sorted(selected_ligand_chains),
            build_protein_training_path,
            ligand_metrics_path,
        )
        with self._debug_timer("bias.training_data.build"):
            run_build_bias_training_data(
                system_path=self.system_path,
                components_cif=components_cif,
                output_protein_csv=build_protein_training_path,
                output_ligand_csv=ligand_metrics_path,
                release_cutoff=self.bias_release_cutoff,
                ligand_similarity_threshold=self.bias_ligand_similarity_threshold,
                overwrite=True,
                skip_bias_csv=True,
                skip_protein_mmseqs=(not run_protein_protocol),
                skip_ligand_ecfp=(not run_ligand_protocol),
                ligand_chains=selected_ligand_chains if run_ligand_protocol else None,
                timings=self.timings,
            )

        protein_metrics_path = (
            build_protein_training_path
            if run_protein_protocol
            else self.protein_training_data_path
        )
        return protein_metrics_path, ligand_metrics_path


class Bias:
    """High-level orchestrator for standalone bias workflow."""

    def __init__(
        self,
        wrk_dir: str,
        system_path: str,
        protein_training_data_path: str | None = None,
        ligand_training_data_path: str | None = None,
        custom_protein_reference_path: str | None = None,
        custom_ligand_reference_path: str | None = None,
        bias_release_cutoff: str = "2023-06-01",
        bias_ligand_similarity_threshold: float = 0.35,
        bias_chains: list[str] | None = None,
        build_bias_training_data: bool = False,
        bias_training_components_cif: str | None = None,
    ):
        self.wrk_dir = Path(wrk_dir)
        self.system_path = Path(system_path)
        self.run_id = str(uuid4())
        self.protein_training_data_path = (
            Path(protein_training_data_path) if protein_training_data_path else None
        )
        self.ligand_training_data_path = (
            Path(ligand_training_data_path) if ligand_training_data_path else None
        )
        self.custom_protein_reference_path = (
            Path(custom_protein_reference_path)
            if custom_protein_reference_path
            else None
        )
        self.custom_ligand_reference_path = (
            Path(custom_ligand_reference_path) if custom_ligand_reference_path else None
        )
        self.bias_release_cutoff = str(bias_release_cutoff)
        self.bias_ligand_similarity_threshold = float(bias_ligand_similarity_threshold)
        self.bias_chains = normalize_bias_chains(bias_chains)
        self.build_bias_training_data = bool(build_bias_training_data)
        self.bias_training_components_cif = (
            Path(bias_training_components_cif) if bias_training_components_cif else None
        )

        self.logger = logging.getLogger(__name__)
        self._system = read.read_yaml(path=self.system_path)
        self.sys = system.System(system=self._system)
        self.timings = DebugTimingCollector(logger=self.logger)

    def _debug_timer(self, label: str):
        return self.timings.measure(label, logger=self.logger)

    def _log_timing_summary(self) -> None:
        self.timings.log_summary(logger=self.logger)

    def run(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        try:
            return self._run_impl()
        except WorkflowExecutionError:
            raise
        except Exception as exc:
            identity = OutputIdentity(
                workflow=WorkflowKind.BIAS,
                run_id=self.run_id,
                system_id=self.system_path.stem,
            )
            failure = failure_from_exception(
                exc,
                identity=identity,
                stage=FailureStage.ANALYTICS,
                error_code="bias_analytics_failed",
            )
            output_dir = self.wrk_dir / "results"
            write_public_bundle(
                PublicOutputBundle(
                    manifest=PublicManifest(
                        schema_version=PUBLIC_SCHEMA_VERSION,
                        identity=identity,
                        status="failed",
                    ),
                    records=(failure,),
                ),
                output_dir,
            )
            raise WorkflowExecutionError(
                str(exc), failures=(failure,), output_dir=output_dir
            ) from exc

    def _run_impl(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        with self._debug_timer("bias.total"):
            self.wrk_dir.mkdir(parents=True, exist_ok=True)
            output_dir = self.wrk_dir / "results" / "bias_train"
            output_dir.mkdir(parents=True, exist_ok=True)
            for legacy_name in ("system_metrics.csv", "chain_metrics.csv"):
                (output_dir / legacy_name).unlink(missing_ok=True)

            system_df, chain_df = build_bias_dataframes(
                self.sys,
                system_name=self.system_path.stem,
            )
            if chain_df.empty:
                raise ValueError(
                    "The system definition must include at least one protein or ligand "
                    "component with an explicit chain ID for bias assessment."
                )

            workflow = BiasAssessmentWorkflow(
                wrk_dir=self.wrk_dir,
                system_path=self.system_path,
                sys_obj=self.sys,
                protein_training_data_path=self.protein_training_data_path,
                ligand_training_data_path=self.ligand_training_data_path,
                custom_protein_reference_path=self.custom_protein_reference_path,
                custom_ligand_reference_path=self.custom_ligand_reference_path,
                bias_release_cutoff=self.bias_release_cutoff,
                bias_ligand_similarity_threshold=self.bias_ligand_similarity_threshold,
                bias_chains=self.bias_chains,
                build_bias_training_data=self.build_bias_training_data,
                bias_training_components_cif=self.bias_training_components_cif,
                logger=self.logger,
                timings=self.timings,
                strict_training_data=True,
            )
            system_df, chain_df = workflow.apply(
                system_df=system_df,
                chain_df=chain_df,
                boltz_cache_path="~/.boltz",
            )

            evidence = [
                EvidenceSource(
                    kind="training_dataset", identifier=path.name, path=str(path)
                )
                for path in (
                    self.protein_training_data_path,
                    self.ligand_training_data_path,
                    self.custom_protein_reference_path,
                    self.custom_ligand_reference_path,
                )
                if path is not None
            ]
            artifacts = (
                ArtifactReference("bias_supporting_outputs", "bias_train", "directory"),
            )
            bundle = bundle_from_frames(
                system_df,
                chain_df,
                identity=OutputIdentity(
                    workflow=WorkflowKind.BIAS,
                    run_id=self.run_id,
                    system_id=self.system_path.stem,
                ),
                evidence=evidence,
                requested_metrics={"bias_metrics"},
                artifacts=artifacts,
            )
            write_public_bundle(bundle, self.wrk_dir / "results")
            self.logger.info(
                "Standalone bias public outputs written to %s", self.wrk_dir / "results"
            )

        self._log_timing_summary()
        return system_df, chain_df
