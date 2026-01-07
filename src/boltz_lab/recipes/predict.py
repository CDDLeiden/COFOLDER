import logging
from pathlib import Path

from boltz_lab.modules.input import command, system
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
    generate_conformers : str or None
        '2D', '3D', or 'sdf' conformer generation.
    sdf_file : str or None
        Path to SDF file if generate_conformers='sdf'.

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
        generate_conformers: str | None = None,
        sdf_file: str | None = None
    ):
        self.wrk_dir = Path(wrk_dir)
        self.system_path = Path(system_path)
        self.options_path = Path(options_path)
        self.repeats = repeats
        self.seed = seed
        self.generate_conformers = generate_conformers
        self.sdf_file = Path(sdf_file) if sdf_file else None

        # Setup logger
        self.logger = logging.getLogger('boltz-lab.recipies.predict')
        self.logger.debug("Initializing Predict with parameters: %s", {
            "wrk_dir": self.wrk_dir,
            "system_path": self.system_path,
            "options_path": self.options_path,
            "repeats": self.repeats,
            "seed": self.seed,
            "generate_conformers": self.generate_conformers,
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
        
        # Get run seeds seeds
        if self.seed is None:
            # generate a random global seed if none provided
            self.seed = helpers.generate_seeds(num_seeds=1, seed=None)[0]
            self.logger.info("No global seed provided. Generated random global seed: %d", self.seed)
        else:
            self.logger.info("Using provided global seed: %d", self.seed)

        if self.repeats == 1:
            self.run_seeds = [self.seed]
            self.logger.debug("Single repeat: using global seed as run seed: %s", self.run_seeds)
        else:
            self.run_seeds = helpers.generate_seeds(
                num_seeds=self.repeats,
                seed=self.seed
            )
            self.logger.debug(
                "Multiple repeats: %d run seeds generated from global seed %d: %s",
                self.repeats,
                self.seed,
                self.run_seeds
            )

        self.logger.info("Run seeds to be used for this workflow: %s", self.run_seeds)

        # Handle conformer generation info
        #if self.generate_conformers in {'2D', '3D'}:
        #    self.logger.info("Generating %s conformers.", self.generate_conformers)

        # Set options
        self.opt.out_dir = self.wrk_dir
        self.opt.system_path = self.system_path

        # Set and run command

        # add repeats #TODO
        cmd = self.opt.set_command(system=self.sys)
        run_boltz(cmd)