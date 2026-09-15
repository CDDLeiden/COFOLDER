"""Command building and options management for Boltz predictions.

This module provides the Command class for managing Boltz command-line
options and building subprocess commands for predictions.
"""

import logging
import multiprocessing
import tempfile
from pathlib import Path

from cofolder.modules.utils import read, write

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

    with tempfile.TemporaryDirectory(prefix="cofolder-boltz-cache-") as tmp:
        tmp_dir = Path(tmp)
        fasta_path = tmp_dir / "setup.fasta"
        write.write_fasta(fasta_path, sequence="A", header="A|protein|")

        cmd = [
            "boltz", "predict",
            str(fasta_path),
            "--cache", str(path),
            "--out_dir", str(tmp_dir),
            "--use_msa_server",
            "--recycling_steps", "1",
            "--sampling_steps", "10"
        ]

        from cofolder.modules.runners.boltz_runner import run_boltz

        run_boltz(cmd)


_BOLTZ2_BUILTIN_CCD_MARKERS = ("ALA.pkl", "GLY.pkl")


def _missing_boltz2_cache_components(path: str | Path) -> tuple[str, ...]:
    """Return missing or empty components of the default Boltz2 cache."""
    cache_path = Path(path)
    missing = []
    for filename in ("boltz2_conf.ckpt", "mols.tar"):
        candidate = cache_path / filename
        if not candidate.is_file() or candidate.stat().st_size == 0:
            missing.append(filename)

    mols_dir = cache_path / "mols"
    if not mols_dir.is_dir():
        missing.append("mols/")
    else:
        for filename in _BOLTZ2_BUILTIN_CCD_MARKERS:
            candidate = mols_dir / filename
            if not candidate.is_file() or candidate.stat().st_size == 0:
                missing.append(f"mols/{filename}")
    return tuple(missing)


def _ensure_boltz2_cache(path: str | Path) -> None:
    """Initialize an incomplete default Boltz2 cache exactly once."""
    import fcntl

    cache_path = Path(path)
    cache_path.mkdir(parents=True, exist_ok=True)
    lock_path = cache_path / ".cofolder_setup.lock"
    with lock_path.open("a+b") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        missing = _missing_boltz2_cache_components(cache_path)
        if not missing:
            return

        # Upstream extraction checks directory existence. Removing only an empty
        # directory lets it repair the unreachable-branch artifact without touching
        # any valid or custom CCD files.
        mols_dir = cache_path / "mols"
        if mols_dir.is_dir() and not any(mols_dir.iterdir()):
            mols_dir.rmdir()

        logger.info(
            "Boltz2 cache at %s is incomplete (%s); running default setup.",
            cache_path,
            ", ".join(missing),
        )
        download_cache(cache_path)
        missing = _missing_boltz2_cache_components(cache_path)
        if missing:
            raise RuntimeError(
                "Boltz2 cache setup did not produce required components at "
                f"{cache_path}: {', '.join(missing)}"
            )

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
        if options is not None and options_path is not None:
            raise ValueError("Provide either 'options' or 'options_path', not both.")

        if options is not None:
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
                if (
                    key == "use_msa_server"
                    or value is None
                    or value is False
                    or value in ("None", "False")
                ):
                    continue
                if value is True or value == "True":
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
