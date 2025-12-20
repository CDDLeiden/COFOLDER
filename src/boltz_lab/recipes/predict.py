import subprocess
import time
import logging

# legacy imports
from boltz_lab.modules.input import command, system
from boltz_lab.modules.utils import helpers  

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
    debug : bool, optional
        Enable debug logging (default: False).

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
        debug: bool = False,
        # generate_conformers: str | None = None  # Uncomment when implemented
    ):
        self.wrk_dir = wrk_dir
        self.system_path = system_path
        self.options_path = options_path
        # self.generate_conformers = generate_conformers  # TODO

        # Setup logger
        self.logger = logging.getLogger('boltz-lab.prediction.Predict')
        self.logger.setLevel(logging.DEBUG if debug else logging.INFO)
        self.logger.debug("Initializing Predict with parameters: %s", {
            "wrk_dir": wrk_dir,
            "system_path": system_path,
            "options_path": options_path,
            # "generate_conformers": generate_conformers
        })

        # Load options and system
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
        start_time = time.time()

        # Set output directory and update system
        self.opt.out_dir = self.wrk_dir
        self.opt.system_path = self.system_path

        # Set and run command
        cmd = self.opt.set_command(system=self.sys)
        self.logger.info(f'Running: {" ".join(cmd)}')
        subprocess.run(cmd)

        self.logger.info(" pred time--- %.2f seconds ---" % (time.time() - start_time))