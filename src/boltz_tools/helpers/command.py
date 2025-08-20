import multiprocessing
import logging
from ..helpers import utils

logger = logging.getLogger('boltz-tools.helpers')

class Command:
    def __init__(self, options=None, options_path=None):
        """
        Initialize Command object from dictionary or YAML file.

        Args:
            options (dict, optional): Pre-loaded options dictionary.
            options_path (str, optional): Path to YAML file with options.
        """
        self.logger = logging.getLogger('boltz-tools.helpers.command.Command')

        if options and options_path:
            raise ValueError("Provide either 'options' or 'options_path', not both.")

        if options:
            self.options = options
        elif options_path:
            self.logger.debug(f"Loading options YAML from {options_path}")
            self.options = utils.read_yaml(path=options_path)
        else:
            raise ValueError("Either 'options' or 'options_path' must be provided.")
    
    def set_command(self, system):
        """
        Build the 'boltz predict' command based on system and options.

        Args:
            system (System): System object for querying settings.

        Returns:
            list: Complete command ready for subprocess execution.
        """
        cmd = ["boltz", "predict", self.system_path, "--out_dir", self.out_dir]

        # Include MSA server option if not defined in the system
        if system.find_value(key="msa") is None:
            cmd.append("--use_msa_server")

        for item in self.options.get("options", []):
            for key, value in item.items():
                if key == "use_msa_server" or value in (None, "None", "False"):
                    continue
                if value == "True":
                    cmd.append(f"--{key}")
                elif value == "multiprocessing.cpu_count()":
                    cmd.extend([f"--{key}", str(multiprocessing.cpu_count())])
                else:
                    cmd.extend([f"--{key}", str(value)])

        return cmd