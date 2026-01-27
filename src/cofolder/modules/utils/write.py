import os
import pandas as pd

from pathlib import Path
import matplotlib.pyplot as plt
from rdkit import Chem
import pickle
import yaml

import logging

logger = logging.getLogger(__name__)

def write_csv(
    df: pd.DataFrame,
    output_path: str,
    index: bool = False
) -> None:
    """Write a DataFrame to CSV, creating parent directories if needed.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to write.
    output_path : str
        Path to output CSV file.
    index : bool, default=False
        Whether to include the DataFrame index in the CSV.

    Returns
    -------
    None
    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    try:
        df.to_csv(output_path, index=index)
        logger.info("CSV written: %s", output_path)
    except Exception as e:
        logger.error("Failed to write CSV %s: %s", output_path, e)
        raise

def write_fasta(fasta_path: str, sequence: str, header: str = "A|protein|") -> None:
    """
    Write a single-sequence FASTA file.

    Parameters
    ----------
    fasta_path : str
        Path to the output FASTA file.
    sequence : str
        Amino acid or nucleotide sequence to write.
    header : str, default "A|protein|"
        FASTA header line (without the leading '>').

    Notes
    -----
    Creates parent directories if they do not exist.
    """
    fasta_path = Path(fasta_path)
    fasta_path.parent.mkdir(parents=True, exist_ok=True)

    with open(fasta_path, "w") as f:
        f.write(f">{header}\n{sequence}\n")
    
    logger.info(f"FASTA file written to: {fasta_path}")

def write_sdf(mols, path: str):
    """Write RDKit molecules to an SDF file."""
    writer = Chem.SDWriter(path)
    count = 0
    for mol in mols:
        if mol is not None:
            writer.write(mol)
            count += 1
    writer.close()
    logger.info("Wrote %d molecules to SDF: %s", count, path)

def write_pickle(obj, path: str):
    """Write a Python object to disk using pickle."""
    with open(path, "wb") as f:
        pickle.dump(obj, f)
    logger.debug("Pickle written: %s", path)

def write_yaml(self, path):
    """Save the system configuration to a YAML file.

    Writes the current system dictionary to a YAML file, preserving
    the order of keys.

    Parameters
    ----------
    path : str
        Output file path for the YAML file.

    Notes
    -----
    The YAML is written with sort_keys=False to preserve insertion order.
    """
    with open(path, "w") as file:
        yaml.dump(self.system, file, sort_keys=False, default_flow_style=False)

def save_plot(output_path: str | Path):
    """
    Save the current matplotlib figure to disk and close it.

    Parameters
    ----------
    output_path : str or Path
        Destination file path for the plot.

    Notes
    -----
    - Parent directories are created automatically
    - The figure is always closed after saving
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    plt.savefig(output_path)
    plt.close()

    logger.info("Plot saved to %s", output_path)

def delete_last_line(file_path):
    """Delete the last line from a file in-place.

    Reads the entire file, removes the last line, and writes the
    modified content back to the same file.

    Parameters
    ----------
    file_path : str
        Path to the file to modify.

    Examples
    --------
    >>> delete_last_line("output.txt")
    >>> delete_last_line("results.log")

    Notes
    -----
    This operation modifies the file in-place. If the file is empty,
    no action is taken. For large files, consider using alternative
    methods to avoid loading the entire file into memory.
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()
    if lines:
        with open(file_path, 'w') as f:
            f.writelines(lines[:-1])