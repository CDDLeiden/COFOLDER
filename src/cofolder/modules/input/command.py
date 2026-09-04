"""Command building and options management for Boltz predictions.

This module provides the Command class for managing Boltz command-line
options and building subprocess commands for predictions.
"""

import multiprocessing
import logging
import os
from cofolder.modules.utils import helpers, read, write

logger = logging.getLogger(__name__)

def download_cache(path: str):
    """Download the Boltz cache by running a minimal prediction.

    Triggers cache download by executing a minimal Boltz prediction with
    a single amino acid. This ensures all required model weights and
    data files are downloaded.

    Parameters
    ----------
    path : str
        Path to the Boltz cache directory to use.

    Notes
    -----
    Creates a temporary directory and FASTA file for the minimal prediction.
    Uses fast settings (minimal recycling and sampling steps) to reduce time.
    """

    tmp_dir = "./tmp"
    os.makedirs(tmp_dir, exist_ok=True)

    fasta_path = os.path.join(tmp_dir, "tmp.fasta")
    # Write a single "A" residue
    write.write_fasta(fasta_path, sequence="A", header="A|protein|")

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
    from cofolder.modules.runners.boltz_runner import run_boltz

    run_boltz(cmd)

class Command:
    """Manage Boltz command-line options and build prediction commands.

    The Command class handles loading, updating, and querying Boltz options,
    and builds complete command-line arguments for subprocess execution.

    Parameters
    ----------
    options : dict, optional
        Pre-loaded options dictionary containing Boltz parameters.
    options_path : str, optional
        Path to YAML file with options configuration.

    Raises
    ------
    ValueError
        If both options and options_path are provided, or if neither is provided.

    Examples
    --------
    >>> # Load from dictionary
    >>> opts = {"options": [{"cache": "~/.boltz"}]}
    >>> cmd = Command(options=opts)

    >>> # Load from YAML file
    >>> cmd = Command(options_path="options.yaml")
    """
    def __init__(self, options=None, options_path=None):
        self.logger = logging.getLogger(__name__)
        if options and options_path:
            raise ValueError("Provide either 'options' or 'options_path', not both.")

        if options:
            self.options = options
        elif options_path:
            self.logger.debug(f"Loading options YAML from {options_path}")
            self.options = read.read_yaml(path=options_path)
        else:
            raise ValueError("Either 'options' or 'options_path' must be provided.")
    
    def set_command(self, system):
        """Build the complete 'boltz predict' command for subprocess execution.

        Constructs the command-line arguments by combining system path,
        output directory, and all configured options. Automatically adds
        MSA server flag if no MSA is defined in the system.

        Parameters
        ----------
        system : System
            System object for querying molecular system settings.

        Returns
        -------
        list of str
            Complete command as a list ready for subprocess.run().

        Notes
        -----
        - Filters out None, "None", and "False" values
        - Handles boolean flags (value="True" becomes just --flag)
        - Supports multiprocessing.cpu_count() for auto-detection
        - Automatically adds --use_msa_server if system has no MSA defined
        """
        cmd = ["boltz", 
               "predict", 
               str(self.system_path), 
               "--out_dir", 
               str(self.out_dir),
               "--seed",
               str(getattr(self, "seed", 0))]

        # Include the MSA server option when at least one protein entity is
        # unresolved. Checking entities individually also supports mixed
        # precomputed/generated-MSA systems.
        sequences = system.find_value(key="sequences") or []
        missing_protein_msa = any(
            isinstance(entry, dict)
            and isinstance(entry.get("protein"), dict)
            and not str(entry["protein"].get("msa") or "").strip()
            for entry in sequences
        )
        if missing_protein_msa:
            cmd.append("--use_msa_server")
            logger.info(
                "At least one protein MSA is unresolved; "
                "adding --use_msa_server to command."
            )

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
        """Retrieve a value from options by path or recursive key search.

        Parameters
        ----------
        key : str, optional
            Key name to search for recursively throughout the options.
        path : list, optional
            Specific path to the value as a list of keys/indices.

        Returns
        -------
        any or None
            The value found at the specified location, or None if not found.

        Raises
        ------
        ValueError
            If path navigation fails or if the key appears multiple times.
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
        
    def update_options(self, value, path=None, parent_key=None, sub_key=None):
        """Update options configuration at a specific path or key.

        Parameters
        ----------
        value : any
            The value to set at the specified location.
        path : list, optional
            Path to the target location as a list of keys/indices.
        parent_key : str, optional
            Top-level or nested key to search for and update.
        sub_key : str, optional
            Sub-key within the parent_key dictionary to update.

        Raises
        ------
        ValueError
            If the path traversal fails due to type mismatch.
        KeyError
            If the parent_key is not found in the options.
        """
        if path:
            d = self.options
            for key in path[:-1]:
                if isinstance(d, dict):
                    d = d.setdefault(key, {})
                elif isinstance(d, list):
                    idx = int(key)
                    d = d[idx]
                else:
                    raise ValueError(f"Cannot traverse into object at {key} in path {path}")

            last_key = path[-1]
            if isinstance(d, dict):
                d[last_key] = value
            elif isinstance(d, list):
                d[int(last_key)] = value
            else:
                raise ValueError(f"Cannot set value at path {path}")

        elif parent_key:
            def set_key(d):
                if isinstance(d, dict):
                    if parent_key in d:
                        if sub_key:
                            d[parent_key] = d.get(parent_key, {})
                            d[parent_key][sub_key] = value
                        else:
                            d[parent_key] = value
                        return True
                    return any(set_key(v) for v in d.values())
                if isinstance(d, list):
                    return any(set_key(i) for i in d)
                return False

            if not set_key(self.options):
                raise KeyError(f"Parent key '{parent_key}' not found in the options.")
