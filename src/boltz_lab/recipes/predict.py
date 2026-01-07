import logging
import os
from pathlib import Path

from boltz_lab.modules.input import command, system
from boltz_lab.modules.entities import ligand
from boltz_lab.modules.runners.boltz_runner import run_boltz
from boltz_lab.modules.utils import helpers  

logger = logging.getLogger(__name__)

class Predict(object):
    """High-level orchestrator for prediction workflow.

    This class coordinates the prediction of protein-ligand co-folding
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

    Examples
    --------
    >>> predictor = Predict(
    ...     wrk_dir="./output",
    ...     system_path="system.yaml",
    ...     options_path="options.yaml"
    ... )
    >>> predictor.run()
    """
    def __init__(
        self,
        wrk_dir: str,
        system_path: str,
        options_path: str,
        repeats: int = 1,
        seed: int | None = None,
        conformers: str | None = None,
        sdf_file: str | None = None
    ):
        self.wrk_dir = Path(wrk_dir)
        self.system_path = Path(system_path)
        self.options_path = Path(options_path)
        self.repeats = repeats
        self.seed = seed
        self.conformers = conformers
        self.sdf_file = Path(sdf_file) if sdf_file else None

        # Setup logger
        self.logger = logging.getLogger('boltz-lab.recipies.predict')
        self.logger.debug("Initializing Predict with parameters: %s", {
            "wrk_dir": self.wrk_dir,
            "system_path": self.system_path,
            "options_path": self.options_path,
            "repeats": self.repeats,
            "seed": self.seed,
            "conformers": self.conformers,
            "sdf_file": self.sdf_file
        })

        # Load YAML options and system
        self._options = helpers.read_yaml(path=self.options_path)
        self.opt = command.Command(options=self._options)

        self._system = helpers.read_yaml(path=self.system_path)
        self.sys = system.System(system=self._system)

        self.logger.debug("Predict initialization complete.")

    def run(self):
        """Execute the prediction workflow.

        Runs Boltz prediction on the configured system and saves
        results to the working directory.
        """
        # Ensure working directory exists
        self.wrk_dir.mkdir(parents=True, exist_ok=True)
        
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
                wrk_dir=self.wrk_dir,
                conformers=self.conformers,
                sdf_file=self.sdf_file,
                global_seed=self.seed,
                logger=self.logger,
            )

        # Set and run command
        for i, seed in enumerate(self.run_seeds, 1):
            logger.info("Running repeat %d/%d with seed %d", i, self.repeats, seed)
            
            # Update the options
            self.opt.seed = seed
            self.opt.out_dir = self.wrk_dir / f"repeat_{i}"
            self.opt.system_path = self.system_path

            # Build the command for this repeat
            cmd = self.opt.set_command(system=self.sys)

            # Log the command safely (convert all Path objects to strings)
            logger.info("Running: %s", " ".join(map(str, cmd)))

            # Execute the command
            run_boltz(cmd)