from __future__ import annotations

import copy
import logging
import os
from pathlib import Path

from cofolder.modules.analytics.bias import apply_bias_metrics
from cofolder.modules.analytics.bias_training import run_build_bias_training_data
from cofolder.modules.analytics.reproduction import scaffold_reproduction_metrics
from cofolder.modules.analytics.structure import Structure
from cofolder.modules.input import system
from cofolder.modules.runners import RunnerRequest, get_runner
from cofolder.modules.utils import gather, helpers, read, write
from cofolder.modules.utils.timing import DebugTimingCollector

logger = logging.getLogger(__name__)

DEFAULT_SCORING_FUNCTIONS = {
    "confidence_metrics",
    "affinity_metrics",
    "affinity_metrics_ext",
    "ifp_distance",
    "sasa",
    "sasa_normalized",
}
DEFAULT_REPRODUCTION_METRICS = {
    "protein_rmsd",
    "ligand_rmsd",
    "sucos",
    "pocket_coverage",
}
RUNNER_METRIC_GROUPS = {
    "confidence_metrics",
    "affinity_metrics",
    "affinity_metrics_ext",
}


class Validate(object):
    """High-level orchestrator for validation workflow."""

    def __init__(
        self,
        wrk_dir: str,
        system_path: str,
        options_path: str,
        runner: str = "boltz",
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
        self.runner_name = str(runner)
        self.repeats = repeats
        self.seed = seed
        self.assess_robustness = assess_robustness
        self.assess_bias = assess_bias
        self.protein_training_data_path = (
            Path(protein_training_data_path) if protein_training_data_path else None
        )
        self.ligand_training_data_path = (
            Path(ligand_training_data_path) if ligand_training_data_path else None
        )
        self.bias_release_cutoff = bias_release_cutoff
        self.bias_ligand_similarity_threshold = float(bias_ligand_similarity_threshold)
        self.build_bias_training_data = build_bias_training_data
        self.bias_training_components_cif = (
            Path(bias_training_components_cif) if bias_training_components_cif else None
        )
        normalized_bias_chains: set[str] = set()
        for value in (bias_chains or []):
            for token in str(value).split(","):
                token = token.strip()
                if token:
                    normalized_bias_chains.add(token.upper())
        self.bias_chains = normalized_bias_chains or None
        self.conformers = conformers
        self.sdf_file = Path(sdf_file) if sdf_file else None
        self.reference_path = Path(reference_path) if reference_path else None
        self.pocket_coverage_reference = pocket_coverage_reference
        self.reproduction_metrics = set(reproduction_metrics or DEFAULT_REPRODUCTION_METRICS)

        if self.reference_path is not None:
            if not self.reference_path.exists():
                raise ValueError(f"reference_path does not exist: {self.reference_path}")
            if not self.reference_path.is_file():
                raise ValueError(f"reference_path is not a file: {self.reference_path}")

        self.logger = logging.getLogger(__name__)
        self.scoring_functions = (
            set(scoring_functions)
            if scoring_functions is not None
            else set(DEFAULT_SCORING_FUNCTIONS)
        )

        self.runner = get_runner(self.runner_name)
        self.runner.ensure_available()

        self._system = read.read_yaml(path=self.system_path)
        self.base_system = system.System(system=self._system)

        self.logger.debug(
            "Initializing Validate with parameters: %s",
            {
                "wrk_dir": self.wrk_dir,
                "system_path": self.system_path,
                "options_path": self.options_path,
                "runner": self.runner_name,
                "repeats": self.repeats,
                "seed": self.seed,
                "conformers": self.conformers,
                "sdf_file": self.sdf_file,
                "assess_bias": self.assess_bias,
                "protein_training_data_path": self.protein_training_data_path,
                "ligand_training_data_path": self.ligand_training_data_path,
                "bias_release_cutoff": self.bias_release_cutoff,
                "bias_ligand_similarity_threshold": self.bias_ligand_similarity_threshold,
                "bias_chains": sorted(self.bias_chains) if self.bias_chains else None,
                "build_bias_training_data": self.build_bias_training_data,
                "bias_training_components_cif": self.bias_training_components_cif,
                "reference_path": self.reference_path,
                "pocket_coverage_reference": self.pocket_coverage_reference,
                "reproduction_metrics": sorted(self.reproduction_metrics),
                "scoring_functions": sorted(self.scoring_functions),
            },
        )

        self.timings = DebugTimingCollector(logger=self.logger)

    def _log_timing_summary(self) -> None:
        self.timings.log_summary(logger=self.logger)

    def _debug_timer(self, label: str):
        return self.timings.measure(label, logger=self.logger)

    def run(self):
        with self._debug_timer("validate.total"):
            self.wrk_dir.mkdir(parents=True, exist_ok=True)
            self.raw_dir = self.wrk_dir / "raw"
            self.raw_dir.mkdir(parents=True, exist_ok=True)
            self.logger.debug("Using raw directory for outputs: %s", self.raw_dir)

            unsupported_metric_groups = sorted(
                (self.scoring_functions & RUNNER_METRIC_GROUPS) - set(self.runner.capabilities)
            )
            if unsupported_metric_groups:
                self.logger.warning(
                    "Runner '%s' does not support scoring groups %s. "
                    "The run will continue and the related output columns will be empty.",
                    self.runner_name,
                    unsupported_metric_groups,
                )

            with self._debug_timer("seeds.resolve"):
                self.seed, self.run_seeds = helpers.get_seeds(
                    repeats=self.repeats,
                    seed=self.seed,
                    logger=self.logger,
                )

            with self._debug_timer("runner.options.load"):
                runner_options = self.runner.load_options(self.options_path)

            self.sys = system.System(system=copy.deepcopy(self.base_system.system))
            with self._debug_timer("runner.prepare_system"):
                preparation = self.runner.prepare_system(
                    # keep backend-specific prep behind the runner boundary
                    system_obj=self.sys,
                    options_obj=runner_options,
                    wrk_dir=self.raw_dir,
                    conformers=self.conformers,
                    sdf_file=self.sdf_file,
                    logger=self.logger,
                )
            self.sys = preparation.system_obj
            runner_options = preparation.options_obj
            for warning in preparation.warnings:
                self.logger.warning("%s", warning)

            yaml_path = os.path.join(self.raw_dir, str(self.system_path.name))
            with self._debug_timer("system_yaml.write"):
                write.write_yaml(self.sys, path=yaml_path)
            self.logger.info("Updated system saved to YAML: %s", yaml_path)

            runner_results = []
            for i, seed in enumerate(self.run_seeds, 1):
                self.logger.info(
                    "Running repeat %d/%d with seed %d using runner '%s'",
                    i,
                    self.repeats,
                    seed,
                    self.runner_name,
                )
                request = RunnerRequest(
                    runner_name=self.runner_name,
                    system_name=self.system_path.stem,
                    system_path=Path(yaml_path),
                    system_obj=self.sys,
                    options_path=self.options_path,
                    options_obj=runner_options,
                    repeat=i,
                    seed=seed,
                    repeat_dir=self.raw_dir / f"repeat_{i}",
                    raw_dir=self.raw_dir,
                    logger=self.logger,
                    timings=self.timings,
                    label_prefix=f"repeat_{i}",
                )
                runner_results.append(self.runner.run(request))

            with self._debug_timer("structures.gather"):
                gather.gather_structures(
                    base_dir=self.wrk_dir,
                    system_name=self.system_path.stem,
                    repeats=self.repeats,
                    logger=self.logger,
                )

            with self._debug_timer("results.merge_runner_outputs"):
                system_df, chain_df, manifests = gather.merge_runner_results(
                    raw_dir=self.raw_dir,
                    repeats=self.repeats,
                    logger=self.logger,
                )

            if not chain_df.empty:
                with self._debug_timer("results.add_chain_info"):
                    chain_df = gather.add_chain_info(chain_df, self.sys)

            system_df, chain_df = self._apply_unsupported_metric_columns(
                system_df=system_df,
                chain_df=chain_df,
                unsupported_metric_groups=unsupported_metric_groups,
            )

            structure_metrics = {
                "ifp_distance",
                "ifp_prolif",
                "sasa",
                "sasa_normalized",
            }
            if self.scoring_functions & structure_metrics and not chain_df.empty:
                with self._debug_timer("structure.init"):
                    structure = Structure(
                        wrk_dir=self.wrk_dir,
                        chain_df=chain_df,
                        cif_folder=Path(self.wrk_dir / "results" / "structures"),
                    )

                if "ifp_distance" in self.scoring_functions:
                    with self._debug_timer("scores.ifp_distance"):
                        chain_df = structure.add_ifp_distance()

                if "ifp_prolif" in self.scoring_functions:
                    with self._debug_timer("scores.ifp_prolif"):
                        chain_df = structure.add_ifp_prolif()

                if {"sasa", "sasa_normalized"} & self.scoring_functions:
                    with self._debug_timer("scores.sasa"):
                        chain_df = structure.add_sasa(
                            absolute="sasa" in self.scoring_functions,
                            normalized="sasa_normalized" in self.scoring_functions,
                        )

            with self._debug_timer("scores.reproduction_metrics"):
                system_df, chain_df = scaffold_reproduction_metrics(
                    system_df=system_df,
                    chain_df=chain_df,
                    reference_path=self.reference_path,
                    wrk_dir=self.wrk_dir,
                    pocket_coverage_reference=self.pocket_coverage_reference,
                    reproduction_metrics=sorted(self.reproduction_metrics),
                    logger=self.logger,
                )

            runtime_context = dict(preparation.runtime_context)
            for result in runner_results:
                runtime_context.update(result.runtime_context)
            for manifest in manifests:
                runtime_context.update((manifest.get("runtime_context") or {}))

            if self.assess_bias:
                system_df, chain_df = self._apply_bias_metrics(
                    system_df=system_df,
                    chain_df=chain_df,
                    runtime_context=runtime_context,
                )

            with self._debug_timer("results.write.system_chain"):
                write.write_csv(system_df, output_path=self.wrk_dir / "results" / "system_metrics.csv")
                write.write_csv(chain_df, output_path=self.wrk_dir / "results" / "chain_metrics.csv")

            diffusion_samples = self._resolve_diffusion_samples(system_df, runtime_context)
            if self.assess_robustness and (self.repeats > 1 or diffusion_samples > 1):
                with self._debug_timer("scores.robustness_metrics"):
                    results_df = gather.gather_robustness_results(
                        system_df=system_df,
                        chain_df=chain_df,
                        wrk_dir=self.wrk_dir,
                        reference_path=self.reference_path,
                    )
                    write.write_csv(
                        results_df,
                        output_path=self.wrk_dir / "results" / "robustness_metrics.csv",
                    )

        self._log_timing_summary()

    def _apply_unsupported_metric_columns(
        self,
        system_df,
        chain_df,
        unsupported_metric_groups: list[str],
    ):
        if "confidence_metrics" in unsupported_metric_groups:
            for column in ("ptm", "iptm", "confidence_score"):
                if column not in system_df.columns:
                    system_df[column] = None
            if "chains_ptm" not in chain_df.columns:
                chain_df["chains_ptm"] = None

        if "affinity_metrics" in unsupported_metric_groups:
            for column in ("affinity_pred_value", "affinity_probability_binary"):
                if column not in chain_df.columns:
                    chain_df[column] = None

        if "affinity_metrics_ext" in unsupported_metric_groups:
            for column in ("pIC50", "IC50_M", "pIC50_kcal_per_mol"):
                if column not in chain_df.columns:
                    chain_df[column] = None

        return system_df, chain_df

    def _apply_bias_metrics(self, system_df, chain_df, runtime_context: dict[str, object]):
        protein_metrics_path = self.protein_training_data_path
        ligand_metrics_path = self.ligand_training_data_path
        if self.build_bias_training_data:
            if self.protein_training_data_path is None:
                raise ValueError(
                    "--build_bias_training_data requires --protein_training_data_path."
                )
            if ligand_metrics_path is None:
                ligand_metrics_path = (
                    self.wrk_dir / "results" / "bias_train" / "ligand_training_data.csv"
                )
                ligand_metrics_path.parent.mkdir(parents=True, exist_ok=True)
            components_cif = self.bias_training_components_cif
            if components_cif is None:
                components_cif = self.protein_training_data_path.parent / "ccd" / "components.cif"
            if not components_cif.exists():
                raise ValueError(
                    f"components.cif not found for in-validate bias build: {components_cif}. "
                    "Provide --bias_training_components_cif or prepare training data root."
                )
            self.logger.info(
                "Building bias training data in validate(): components_cif=%s protein_csv=%s ligand_csv=%s",
                components_cif,
                self.protein_training_data_path,
                ligand_metrics_path,
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
            selected_chains = (
                {str(x).strip().upper() for x in self.bias_chains}
                if self.bias_chains
                else None
            )
            selected_protein_chains = (
                (available_protein_chains & selected_chains)
                if selected_chains is not None
                else available_protein_chains
            )
            selected_ligand_chains = (
                (available_ligand_chains & selected_chains)
                if selected_chains is not None
                else available_ligand_chains
            )
            invalid_smiles_chains = self._find_invalid_ligand_smiles_chain_ids()
            invalid_selected_ligand_chains = selected_ligand_chains & invalid_smiles_chains
            if invalid_selected_ligand_chains:
                self.logger.warning(
                    "Invalid ligand SMILES detected for chains=%s. "
                    "Skipping ligand ECFP protocol for this bias-build run.",
                    sorted(invalid_selected_ligand_chains),
                )
                selected_ligand_chains = selected_ligand_chains - invalid_selected_ligand_chains
            run_protein_protocol = bool(selected_protein_chains)
            run_ligand_protocol = bool(selected_ligand_chains)
            build_protein_training_path = self.protein_training_data_path
            if not run_protein_protocol:
                build_protein_training_path = (
                    self.wrk_dir / "results" / "bias_train" / "_protein_training_data_build_tmp.csv"
                )
                build_protein_training_path.parent.mkdir(parents=True, exist_ok=True)
            self.logger.info(
                "Bias build protocol selection: run_protein=%s chains=%s | run_ligand=%s chains=%s | "
                "protein_build_csv=%s | ligand_build_csv=%s",
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
            if run_protein_protocol:
                protein_metrics_path = build_protein_training_path

        protein_ok = (
            protein_metrics_path is not None
            and protein_metrics_path.exists()
            and protein_metrics_path.is_file()
        )
        ligand_ok = (
            ligand_metrics_path is not None
            and ligand_metrics_path.exists()
            and ligand_metrics_path.is_file()
        )

        if protein_ok and ligand_ok:
            with self._debug_timer("scores.bias_metrics.total"):
                system_df, chain_df = apply_bias_metrics(
                    system_df=system_df,
                    chain_df=chain_df,
                    sys_obj=self.sys,
                    protein_training_data_path=protein_metrics_path,
                    ligand_training_data_path=ligand_metrics_path,
                    release_cutoff=self.bias_release_cutoff,
                    bias_chains=self.bias_chains,
                    protein_top_n=100,
                    boltz_cache_path=runtime_context.get("cache_path") or "~/.boltz",
                    output_dir=self.wrk_dir / "results" / "bias_train",
                    logger=self.logger,
                    timings=self.timings,
                )
        else:
            self.logger.warning(
                "Bias assessment requested but training data files are missing/unavailable "
                "(protein=%s, ligand=%s). Skipping bias metrics. "
                "To bootstrap bias references, run: "
                "`python scripts/fetch_bias_training_data.py "
                "--output_root <bias_data_dir>` then "
                "`python scripts/build_bias_training_data.py "
                "--system_path <system.yaml> "
                "--components_cif <bias_data_dir>/ccd/components.cif "
                "--output_protein_csv <protein_training.csv> "
                "--output_ligand_csv <ligand_training.csv> "
                "--release_cutoff %s`",
                protein_metrics_path,
                ligand_metrics_path,
                self.bias_release_cutoff,
            )
        return system_df, chain_df

    @staticmethod
    def _resolve_diffusion_samples(system_df, runtime_context: dict[str, object]) -> int:
        if not system_df.empty and "diffusion_sample" in system_df.columns:
            try:
                return int(system_df["diffusion_sample"].max()) + 1
            except Exception:
                pass
        try:
            return int(runtime_context.get("diffusion_samples") or 1)
        except Exception:
            return 1

    def _find_invalid_ligand_smiles_chain_ids(self) -> set[str]:
        try:
            from rdkit import Chem
        except Exception:
            return set()

        invalid: set[str] = set()
        sequences = self.sys.find_value(key="sequences") or []
        for entry in sequences:
            if not isinstance(entry, dict) or "ligand" not in entry:
                continue
            ligand_data = entry.get("ligand") or {}
            smiles = str(ligand_data.get("smiles", "")).strip()
            if not smiles:
                continue
            if Chem.MolFromSmiles(smiles) is not None:
                continue
            chain_ids = ligand_data.get("id")
            if chain_ids is None:
                continue
            if isinstance(chain_ids, list):
                for cid in chain_ids:
                    invalid.add(str(cid).strip().upper())
            else:
                invalid.add(str(chain_ids).strip().upper())
        return {c for c in invalid if c}
