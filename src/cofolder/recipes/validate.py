import logging
import os
from pathlib import Path

from cofolder.modules.input import command, system
from cofolder.modules.entities import ligand
from cofolder.modules.runners.boltz_runner import run_boltz
from cofolder.modules.utils import helpers, gather, read, write
from cofolder.modules.analytics.structure import Structure
from cofolder.modules.analytics.reproduction import scaffold_reproduction_metrics
from cofolder.modules.analytics.bias import apply_bias_metrics
from cofolder.modules.analytics.bias_training import run_build_bias_training_data

logger = logging.getLogger(__name__)

DEFAULT_SCORING_FUNCTIONS = {
    "boltz_confidence_metrics",
    "boltz_affinity_metrics",
    "boltz_affinity_metrics_ext",
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

class Validate(object):
    """High-level orchestrator for validation workflow.

    This class coordinates the validation of protein-ligand co-folding
    for a single system using Boltz.

    Parameters
    ----------
    wrk_dir : str
        Working directory for output files.
    system_path : str
        Path to system YAML file defining the molecular system.
    options_path : str
        Path to Boltz options YAML file.
    repeats : int
        Number of repeats.
    seeds : list[int] or None
        Optional list of seeds for reproducibility.
    conformers : str or None
        '2D', '3D', or 'sdf' conformer generation.
    sdf_file : str or None
        Path to existing SDF file if conformers='sdf' or save location if conformers='2D' or '3D'.
    reference_path : str or None
        Optional path to reference PDB/CIF for model reproduction metrics.
    pocket_coverage_reference : str or None
        Optional custom pocket-coverage reference ('2 8 10', 'A2 S8 T10', or '0100000101').

    Examples
    --------
    >>> validator = Validate(
    ...     wrk_dir="./output",
    ...     system_path="system.yaml",
    ...     options_path="options.yaml"
    ... )
    >>> validator.run()
    """
    def __init__(
        self,
        wrk_dir: str,
        system_path: str,
        options_path: str,
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
        self.repeats = repeats
        self.seed = seed
        self.scoring_functions = scoring_functions
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
        # Setup logger
        self.logger = logging.getLogger(__name__)
        self.logger.debug("Initializing Validate with parameters: %s", {
            "wrk_dir": self.wrk_dir,
            "system_path": self.system_path,
            "options_path": self.options_path,
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
        })

        # Allow direct class usage without explicitly passing scoring functions.
        self.scoring_functions = (
            set(scoring_functions)
            if scoring_functions is not None
            else set(DEFAULT_SCORING_FUNCTIONS)
        )

        self.logger.debug(
            "Enabled scoring functions: %s",
            sorted(self.scoring_functions)
        )

        # Load YAML options and system
        self._options = read.read_yaml(path=self.options_path)
        self.opt = command.Command(options=self._options)

        self._system = read.read_yaml(path=self.system_path)
        self.sys = system.System(system=self._system)

        self.logger.debug("Validate initialization complete.")

    def run(self):
        """Execute the validation workflow.

        Runs Boltz on the configured system and saves
        results to the working directory.
        """
        # Ensure working directory exists
        self.wrk_dir.mkdir(parents=True, exist_ok=True)
        self.raw_dir = self.wrk_dir / "raw"
        self.raw_dir.mkdir(parents=True, exist_ok=True)
        self.logger.debug("Using raw directory for outputs: %s", self.raw_dir)

        # Get run seeds
        self.seed, self.run_seeds = helpers.get_seeds(
            repeats=self.repeats,
            seed=self.seed,
            logger=self.logger
        )

        # Generate conformers and store in cache
        if not self.conformers:
            self.resname = None
            logger.debug("No conformer generation requested. Using SMILES.")
        else:
            self.resname = ligand.handle_conformers(
                sys_obj=self.sys,
                opt_obj=self.opt,
                wrk_dir=self.raw_dir,
                conformers=self.conformers,
                sdf_file=self.sdf_file,
                logger=self.logger,
            )

        # --- Save updated system to YAML ---
        yaml_path = os.path.join(self.raw_dir, str(self.system_path.name))
        try:
            write.write_yaml(self.sys, path=yaml_path)
            self.logger.info("Updated system saved to YAML: %s", yaml_path)
        except Exception as e:
            self.logger.error("Failed to save updated system YAML: %s", e)

        # Set and run command
        for i, seed in enumerate(self.run_seeds, 1):
            logger.info("Running repeat %d/%d with seed %d", i, self.repeats, seed)
            
            # Update the options
            self.opt.seed = seed
            self.opt.out_dir = self.raw_dir / f"repeat_{i}"
            self.opt.system_path = yaml_path

            # Build the command for this repeat
            cmd = self.opt.set_command(system=self.sys)
            logger.info("Running: %s", " ".join(map(str, cmd)))

            # Execute the command
            run_boltz(cmd)

        # Gather structures from all repeats
        gather.gather_structures(
            base_dir=self.wrk_dir,
            system_name=self.system_path.stem,
            repeats=self.repeats,
            logger=self.logger
        )

        # initialize results dataframes
        system_df, chain_df = gather.initialize_results(
            raw_dir=self.raw_dir,
            system_name=self.system_path.stem,
            repeats=self.repeats,
            diffusion_samples=self.opt.find_value(key="diffusion_samples") if self.opt.find_value(key="diffusion_samples") else 1
        )
        chain_df = gather.add_chain_info(chain_df, self.sys)

        # gather confidence metrics
        if "boltz_confidence_metrics" in self.scoring_functions:
            system_df, chain_df = gather.gather_confidence_metrics(
                raw_dir=self.raw_dir,
                system_df=system_df,
                chain_df=chain_df,
                system_name=self.system_path.stem,
                repeats=self.repeats,
                diffusion_samples=(
                    self.opt.find_value(key="diffusion_samples")
                    if self.opt.find_value(key="diffusion_samples")
                    else 1
                )
            )

        # gather affinity metrics
        if {
            "boltz_affinity_metrics",
            "boltz_affinity_metrics_ext",
        } & self.scoring_functions:
            chain_df = gather.gather_affinity_metrics(
                raw_dir=self.raw_dir,
                chain_df=chain_df,
                system_name=self.system_path.stem,
                sys=self.sys,
                repeats=self.repeats,
                extended="boltz_affinity_metrics_ext" in self.scoring_functions,
            )

        # gather structure-based metrics
        structure_metrics = {
            "ifp_distance",
            "ifp_prolif",
            "sasa",
            "sasa_normalized",
        }

        if self.scoring_functions & structure_metrics:
            structure = Structure(
                wrk_dir=self.wrk_dir,
                chain_df=chain_df,
                cif_folder=Path(self.wrk_dir / "results" / "structures"),
            )

            if "ifp_distance" in self.scoring_functions:
                chain_df = structure.add_ifp_distance()

            if "ifp_prolif" in self.scoring_functions:
                chain_df = structure.add_ifp_prolif()

        if {
            "sasa",
            "sasa_normalized",
        } & self.scoring_functions:
            chain_df = structure.add_sasa(
                absolute="sasa" in self.scoring_functions,
                normalized="sasa_normalized" in self.scoring_functions,
            )

        # Scaffold model reproduction and bias schema.
        system_df, chain_df = scaffold_reproduction_metrics(
            system_df=system_df,
            chain_df=chain_df,
            reference_path=self.reference_path,
            wrk_dir=self.wrk_dir,
            pocket_coverage_reference=self.pocket_coverage_reference,
            reproduction_metrics=sorted(self.reproduction_metrics),
            logger=self.logger,
        )

        if self.assess_bias:
            if self.build_bias_training_data:
                if self.protein_training_data_path is None or self.ligand_training_data_path is None:
                    raise ValueError(
                        "--build_bias_training_data requires both "
                        "--protein_training_data_path and --ligand_training_data_path."
                    )
                components_cif = self.bias_training_components_cif
                if components_cif is None:
                    components_cif = self.ligand_training_data_path.parent / "ccd" / "components.cif"
                if not components_cif.exists():
                    raise ValueError(
                        f"components.cif not found for in-validate bias build: {components_cif}. "
                        "Provide --bias_training_components_cif or prepare training data root."
                    )
                self.logger.info(
                    "Building bias training data in validate(): components_cif=%s protein_csv=%s ligand_csv=%s",
                    components_cif,
                    self.protein_training_data_path,
                    self.ligand_training_data_path,
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
                self.logger.info(
                    "Bias build protocol selection: run_protein=%s chains=%s | run_ligand=%s chains=%s",
                    run_protein_protocol,
                    sorted(selected_protein_chains),
                    run_ligand_protocol,
                    sorted(selected_ligand_chains),
                )
                run_build_bias_training_data(
                    system_path=self.system_path,
                    components_cif=components_cif,
                    output_protein_csv=self.protein_training_data_path,
                    output_ligand_csv=self.ligand_training_data_path,
                    release_cutoff=self.bias_release_cutoff,
                    ligand_similarity_threshold=self.bias_ligand_similarity_threshold,
                    overwrite=True,
                    skip_bias_csv=True,
                    skip_protein_mmseqs=(not run_protein_protocol),
                    skip_ligand_ecfp=(not run_ligand_protocol),
                    ligand_chains=selected_ligand_chains if run_ligand_protocol else None,
                )

            protein_ok = (
                self.protein_training_data_path is not None
                and self.protein_training_data_path.exists()
                and self.protein_training_data_path.is_file()
            )
            ligand_ok = (
                self.ligand_training_data_path is not None
                and self.ligand_training_data_path.exists()
                and self.ligand_training_data_path.is_file()
            )

            if protein_ok and ligand_ok:
                system_df, chain_df = apply_bias_metrics(
                    system_df=system_df,
                    chain_df=chain_df,
                    sys_obj=self.sys,
                    protein_training_data_path=self.protein_training_data_path,
                    ligand_training_data_path=self.ligand_training_data_path,
                    release_cutoff=self.bias_release_cutoff,
                    bias_chains=self.bias_chains,
                    protein_top_n=100,
                    boltz_cache_path=self.opt.find_value(key="cache") or "~/.boltz",
                    output_dir=self.wrk_dir / "results" / "bias_train",
                    logger=self.logger,
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
                    self.protein_training_data_path,
                    self.ligand_training_data_path,
                    self.bias_release_cutoff,
                )

        write.write_csv(system_df, output_path=self.wrk_dir / "results" / "system_metrics.csv")
        write.write_csv(chain_df, output_path=self.wrk_dir / "results" / "chain_metrics.csv")

        if self.assess_robustness and \
            (self.repeats > 1 or \
             (self.opt.find_value(key="diffusion_samples") or 1) > 1):
            results_df = gather.gather_robustness_results(
                system_df=system_df,
                chain_df=chain_df,
                wrk_dir=self.wrk_dir,
                reference_path=self.reference_path,
            )

            write.write_csv(results_df, output_path=self.wrk_dir / "results" / "robustness_metrics.csv")

    def _find_invalid_ligand_smiles_chain_ids(self) -> set[str]:
        """Return ligand chain IDs with invalid/empty SMILES in current system."""
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
