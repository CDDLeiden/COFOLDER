from __future__ import annotations

import logging
from contextlib import nullcontext
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
from cofolder.modules.input import InputValidationError, system
from cofolder.modules.utils import read
from cofolder.modules.utils.timing import DebugTimingCollector
from cofolder.recipes._completion import report_completion

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
        bias_protein_similarity_threshold: float,
        bias_ligand_similarity_threshold: float,
        bias_chains: set[str] | list[str] | None,
        build_bias_training_data: bool,
        bias_training_components_cif: Path | None,
        bias_training_data_protein_path: Path | None = None,
        bias_training_data_ligand_path: Path | None = None,
        bias_query_cache_path: Path | None = None,
        custom_bias_reference_path: Path | None = None,
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
        self.bias_training_data_protein_path = bias_training_data_protein_path
        self.bias_training_data_ligand_path = bias_training_data_ligand_path
        self.bias_query_cache_path = bias_query_cache_path
        self.custom_bias_reference_path = custom_bias_reference_path
        self.custom_protein_reference_path = custom_protein_reference_path
        self.custom_ligand_reference_path = custom_ligand_reference_path
        self.bias_release_cutoff = str(bias_release_cutoff)
        self.bias_protein_similarity_threshold = float(
            bias_protein_similarity_threshold
        )
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
        from cofolder.modules.analytics.bias_database import parse_bias_release_policy

        parse_bias_release_policy(self.bias_release_cutoff)
        for name, value in (
            ("protein", self.bias_protein_similarity_threshold),
            ("ligand", self.bias_ligand_similarity_threshold),
        ):
            if not 0.0 <= value <= 1.0:
                raise ValueError(
                    f"Bias {name} similarity threshold must be between 0 and 1."
                )

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
        custom_protein_reference_path = self.custom_protein_reference_path
        custom_ligand_reference_path = self.custom_ligand_reference_path
        custom_fingerprint = None
        if self.custom_bias_reference_path is not None:
            from cofolder.modules.analytics.bias_database import validate_custom_bias_reference_bundle

            if self.custom_protein_reference_path or self.custom_ligand_reference_path:
                raise ValueError(
                    "--custom_bias_reference_path cannot be combined with paired custom reference paths."
                )
            custom_bundle = validate_custom_bias_reference_bundle(
                self.custom_bias_reference_path
            )
            custom_protein_reference_path = custom_bundle.protein_path
            custom_ligand_reference_path = custom_bundle.ligand_path
            custom_fingerprint = custom_bundle.fingerprint

        database_mode = bool(
            self.bias_training_data_protein_path
            or self.bias_training_data_ligand_path
            or (
                self.protein_training_data_path is None
                and self.ligand_training_data_path is None
                and custom_protein_reference_path is None
                and custom_ligand_reference_path is None
                and not self.build_bias_training_data
            )
        )
        if database_mode:
            if self.protein_training_data_path or self.ligand_training_data_path or self.build_bias_training_data:
                raise ValueError(
                    "Database-backed bias paths cannot be combined with legacy "
                    "--protein_training_data_path, --ligand_training_data_path, or "
                    "--build_bias_training_data."
                )
            from cofolder.modules.analytics.bias_database import (
                DEFAULT_LIGAND_DATABASE,
                DEFAULT_PROTEIN_DATABASE,
                materialize_bias_references,
                validate_bias_database_bundle,
            )

            entity_types = set(chain_df.get("ENTITY_TYPE", pd.Series(dtype=str)).astype(str))
            protein_bundle = (
                validate_bias_database_bundle(
                    self.bias_training_data_protein_path or DEFAULT_PROTEIN_DATABASE,
                    "protein",
                )
                if "protein" in entity_types
                else None
            )
            ligand_bundle = (
                validate_bias_database_bundle(
                    self.bias_training_data_ligand_path or DEFAULT_LIGAND_DATABASE,
                    "ligand",
                )
                if "ligand" in entity_types
                else None
            )
            references = materialize_bias_references(
                system_obj=self.sys,
                protein_bundle=protein_bundle,
                ligand_bundle=ligand_bundle,
                output_dir=self.output_dir,
                release_cutoff=self.bias_release_cutoff,
                protein_similarity_threshold=self.bias_protein_similarity_threshold,
                ligand_similarity_threshold=self.bias_ligand_similarity_threshold,
                selected_chains=self.bias_chains,
                query_cache_path=self.bias_query_cache_path,
                custom_fingerprint=custom_fingerprint,
            )
            protein_metrics_path = references.protein_path
            ligand_metrics_path = references.ligand_path

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
            has_protein_file or custom_protein_reference_path is not None
        )
        has_ligand_source = (
            has_ligand_file or custom_ligand_reference_path is not None
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
                    custom_protein_reference_path=custom_protein_reference_path,
                    custom_ligand_reference_path=custom_ligand_reference_path,
                    release_cutoff=self.bias_release_cutoff,
                    bias_chains=self.bias_chains,
                    protein_top_n=100,
                    protein_similarity_threshold=self.bias_protein_similarity_threshold,
                    ligand_similarity_threshold=self.bias_ligand_similarity_threshold,
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
            f"custom protein={custom_protein_reference_path} custom ligand={custom_ligand_reference_path}."
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
                protein_similarity_threshold=(
                    self.bias_protein_similarity_threshold * 100.0
                ),
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
        bias_training_data_protein_path: str | None = None,
        bias_training_data_ligand_path: str | None = None,
        bias_query_cache_path: str | None = None,
        custom_bias_reference_path: str | None = None,
        custom_protein_reference_path: str | None = None,
        custom_ligand_reference_path: str | None = None,
        bias_release_cutoff: str = "2023-06-01",
        bias_protein_similarity_threshold: float = 0.25,
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
        self.bias_training_data_protein_path = (
            Path(bias_training_data_protein_path)
            if bias_training_data_protein_path
            else None
        )
        self.bias_training_data_ligand_path = (
            Path(bias_training_data_ligand_path)
            if bias_training_data_ligand_path
            else None
        )
        self.bias_query_cache_path = (
            Path(bias_query_cache_path) if bias_query_cache_path else None
        )
        self.custom_bias_reference_path = (
            Path(custom_bias_reference_path) if custom_bias_reference_path else None
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
        self.bias_protein_similarity_threshold = float(
            bias_protein_similarity_threshold
        )
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

    def preflight(self):
        """Validate the standalone bias request without searches or writes."""
        from cofolder.modules.analytics.bias_database import (
            bias_query_cache_readiness,
            DEFAULT_LIGAND_DATABASE,
            DEFAULT_PROTEIN_DATABASE,
            _selected_queries,
            parse_bias_release_policy,
            validate_custom_bias_reference_bundle,
            validate_bias_database_bundle,
        )
        from cofolder.modules.input.system import iter_system_chains
        from cofolder.recipes.preflight import PreflightReport

        try:
            self._validate_source_mode()
            policy = parse_bias_release_policy(self.bias_release_cutoff)
            for name, value in (
                ("protein", self.bias_protein_similarity_threshold),
                ("ligand", self.bias_ligand_similarity_threshold),
            ):
                if not 0.0 <= value <= 1.0:
                    raise ValueError(
                        f"Bias {name} similarity threshold must be between 0 and 1."
                    )
            sources: list[str] = []
            protein_bundle = None
            ligand_bundle = None
            database_mode = not (
                self.protein_training_data_path
                or self.ligand_training_data_path
                or self.custom_protein_reference_path
                or self.custom_ligand_reference_path
                or self.custom_bias_reference_path
                or self.build_bias_training_data
            ) or bool(
                self.bias_training_data_protein_path
                or self.bias_training_data_ligand_path
            )
            entity_types = {
                chain.entity_type
                for chain in iter_system_chains(self.sys)
                if not self.bias_chains
                or chain.chain_id.strip().upper() in self.bias_chains
            }
            if database_mode and "protein" in entity_types:
                protein_bundle = validate_bias_database_bundle(
                    self.bias_training_data_protein_path or DEFAULT_PROTEIN_DATABASE,
                    "protein",
                )
                sources.append(
                    f"protein:{protein_bundle.root}:{protein_bundle.fingerprint}"
                )
            if database_mode and "ligand" in entity_types:
                ligand_bundle = validate_bias_database_bundle(
                    self.bias_training_data_ligand_path or DEFAULT_LIGAND_DATABASE,
                    "ligand",
                )
                sources.append(
                    f"ligand:{ligand_bundle.root}:{ligand_bundle.fingerprint}"
                )
            if self.custom_bias_reference_path:
                custom = validate_custom_bias_reference_bundle(self.custom_bias_reference_path)
                sources.append(f"custom:{custom.root}:{custom.fingerprint}")
            cache_root, cache_ready, cache_state = bias_query_cache_readiness(
                self.bias_query_cache_path
            )
            ligand_smiles_by_id: dict[str, str] = {}
            if ligand_bundle is not None:
                source = pd.read_csv(ligand_bundle.data_path)
                if {"ligand_id", "smiles"} <= set(source.columns):
                    ligand_smiles_by_id = {
                        str(row["ligand_id"]).upper(): str(row["smiles"])
                        for _, row in source.dropna(
                            subset=["ligand_id", "smiles"]
                        ).iterrows()
                    }
            proteins, ligands = _selected_queries(
                self.sys, self.bias_chains, ligand_smiles_by_id
            )
            import hashlib
            from rdkit import Chem

            query_hashes = [
                f"protein:{chain}:{hashlib.sha256(''.join(sequence.split()).upper().encode()).hexdigest()}"
                for chain, sequence in sorted(proteins.items())
            ]
            for chain, smiles in sorted(ligands.items()):
                molecule = Chem.MolFromSmiles(smiles)
                canonical = (
                    Chem.MolToSmiles(molecule, isomericSmiles=True)
                    if molecule
                    else smiles
                )
                query_hashes.append(
                    f"ligand:{chain}:{hashlib.sha256(canonical.encode()).hexdigest()}"
                )
            return PreflightReport(
                workflow="bias",
                ready=True,
                selected_chains=tuple(
                    sorted(chain.chain_id for chain in iter_system_chains(self.sys))
                ),
                database_sources=tuple(sources),
                release_cutoff=self.bias_release_cutoff,
                release_policy=policy.mode,
                query_hashes=tuple(query_hashes),
                bias_query_cache=f"{cache_root} ({cache_state})",
                bias_cache_ready=cache_ready,
                expected_bias_backfill=bool(
                    protein_bundle is not None
                    and protein_bundle.sequence_index_path is not None
                    and proteins
                    and ligands
                ),
                protein_similarity_threshold=self.bias_protein_similarity_threshold,
                ligand_similarity_threshold=self.bias_ligand_similarity_threshold,
                output_dir=str(self.wrk_dir / "results"),
            )
        except Exception as exc:
            return PreflightReport(
                workflow="bias",
                ready=False,
                release_cutoff=self.bias_release_cutoff,
                output_dir=str(self.wrk_dir / "results"),
                messages=(str(exc),),
            )

    def _validate_source_mode(self) -> None:
        if self.custom_bias_reference_path and (
            self.custom_protein_reference_path or self.custom_ligand_reference_path
        ):
            raise ValueError(
                "--custom_bias_reference_path cannot be combined with paired custom reference paths."
            )
        if (
            self.bias_training_data_protein_path
            or self.bias_training_data_ligand_path
        ) and (
            self.protein_training_data_path
            or self.ligand_training_data_path
            or self.build_bias_training_data
        ):
            raise ValueError(
                "Database-backed bias paths cannot be combined with legacy public-reference paths/build mode."
            )

    @report_completion(WorkflowKind.BIAS)
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
                stage=(
                    FailureStage.INPUT_VALIDATION
                    if isinstance(exc, InputValidationError)
                    or "database" in str(exc).lower()
                    else FailureStage.ANALYTICS
                ),
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
            self._validate_source_mode()
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
                bias_training_data_protein_path=self.bias_training_data_protein_path,
                bias_training_data_ligand_path=self.bias_training_data_ligand_path,
                bias_query_cache_path=self.bias_query_cache_path,
                custom_bias_reference_path=self.custom_bias_reference_path,
                custom_protein_reference_path=self.custom_protein_reference_path,
                custom_ligand_reference_path=self.custom_ligand_reference_path,
                bias_release_cutoff=self.bias_release_cutoff,
                bias_protein_similarity_threshold=self.bias_protein_similarity_threshold,
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
                    self.bias_training_data_protein_path,
                    self.bias_training_data_ligand_path,
                    self.custom_bias_reference_path,
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

        self._log_timing_summary()
        return system_df, chain_df
