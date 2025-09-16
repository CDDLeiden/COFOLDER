import multiprocessing
import logging
from ..helpers import utils

logger = logging.getLogger('boltz-tools.helpers')

import os
import subprocess

def download_cache(path: str):
    """
    Downloads the Boltz cache by running a minimal prediction.

    Parameters
    ----------
    path : str
        Path to the Boltz cache directory to use.
    """

    tmp_dir = "./tmp"
    os.makedirs(tmp_dir, exist_ok=True)

    fasta_path = os.path.join(tmp_dir, "tmp.fasta")
    # Write a single "A" residue
    with open(fasta_path, "w") as f:
        f.write(">A|protein|\nA\n")

    # Build the command
    cmd = [
        "boltz", "predict",
        fasta_path,
        "--cache", path,
        "--out_dir", tmp_dir,
        "--use_msa_server",
        "--recycling_steps", "1",
        "--sampling_steps", "10"
    ]

    # Run the command
    subprocess.run(cmd, check=True)

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
    
    def find_value(self, key=None, path=None):
        """
        Retrieve a value by path or recursively search by key.

        Args:
            key (str, optional): Key to search recursively.
            path (list, optional): Specific path to the value.
        Returns:
            The value found, or None if not found.
        """
        if path:
            cur = self.options
            for p in path:
                if isinstance(cur, dict):
                    cur = cur[p]
                elif isinstance(cur, list):
                    cur = cur[int(p)]
                else:
                    raise ValueError(f"Path {path} is invalid")
            return cur

        if key:
            def search(d):
                results = []
                if isinstance(d, dict):
                    for k, v in d.items():
                        if k == key:
                            results.append(v)
                        results.extend(search(v))
                elif isinstance(d, list):
                    for item in d:
                        results.extend(search(item))
                return results

            found = search(self.options)
            if not found:
                logger.info(f"Key not found in system: {key}")
                return None
            if len(found) > 1:
                raise ValueError(f"Key '{key}' appears multiple times; use path instead.")
            return found[0]
