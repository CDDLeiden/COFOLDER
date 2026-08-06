import pandas as pd

import yaml
import json
from rdkit import Chem

import logging

logger = logging.getLogger(__name__)

def read_csv(path, columns):
    """Read a CSV file and validate expected columns.

    Loads a CSV file into a pandas DataFrame and checks for the presence
    of specified columns, logging warnings for any missing columns.

    Parameters
    ----------
    path : str
        Path to the CSV file to read.
    columns : list of str
        List of column names to validate in the CSV file.

    Returns
    -------
    pd.DataFrame or None
        DataFrame containing the CSV data, or None if file not found.

    Examples
    --------
    >>> df = read_csv("compounds.csv", ["smiles", "id", "affinity"])
    >>> df = read_csv("results.csv", ["name", "value"])
    """
    try:
        df = pd.read_csv(path)
        logger.info(f"Read {path} containing {len(df)} entries")

        for col in columns:
            if col and col in df.columns:
                logger.info(f"Column {col} in {path}")
            else:
                logger.warning(f"Column {col} not in {path}")

        return df

    except FileNotFoundError:
        logging.error(f"File not found: {path}")

def read_json(path):
    """Read and parse a JSON file.

    Loads a JSON file and logs its presence at INFO level.
    Full contents are logged only at DEBUG level.

    Parameters
    ----------
    path : str
        Path to the JSON file to read.

    Returns
    -------
    dict or None
        Dictionary containing the parsed JSON data, or None if an error occurs.

    Raises
    ------
    FileNotFoundError
        If the specified file does not exist.
    json.JSONDecodeError
        If the JSON file is malformed or cannot be parsed.

    Examples
    --------
    >>> config = read_json("config.json")
    >>> results = read_json("confidence_model_0.json")
    """
    try:
        with open(path, "r") as file:
            data = json.load(file)
            logging.info("%s loaded successfully.", path)

            if isinstance(data, dict):
                logger.debug("%s contents:", path)
                for key, value in data.items():
                    logger.debug("\t%s: %s", key, value)
            else:
                logger.debug("\t<Non-dict JSON root of type %s>", type(data).__name__)

            return data

    except FileNotFoundError:
        logging.error("File not found: %s", path)
    except json.JSONDecodeError as e:
        logging.error("Error parsing JSON: %s", e)

def read_yaml(path):
    """Read and parse a YAML file.

    Loads a YAML file and logs its contents for debugging purposes.

    Parameters
    ----------
    path : str
        Path to the YAML file to read.

    Returns
    -------
    dict or None
        Dictionary containing the parsed YAML data, or None if an error occurs.

    Raises
    ------
    FileNotFoundError
        If the specified file does not exist.
    yaml.YAMLError
        If the YAML file is malformed or cannot be parsed.

    Examples
    --------
    >>> config = read_yaml("config.yaml")
    >>> system = read_yaml("system.yaml")
    """
    try:
        with open(path, 'r') as file:
            data = yaml.safe_load(file)
            logging.info(f"{path} loaded successfully. Contents:")
            for key, value in data.items():
                logging.info(f"\t{key}: {value}")

            return data

    except FileNotFoundError:
        logging.error(f"File not found: {path}")
    except yaml.YAMLError as e:
        logging.error(f"Error parsing YAML: {e}")

def read_sdf(path):
    """Read molecules from an SDF file.

    Loads molecules from an SDF (Structure Data File) using RDKit,
    filtering out invalid molecules.

    Parameters
    ----------
    path : str
        Path to the SDF file to read.

    Returns
    -------
    list of rdkit.Chem.Mol or None
        List of valid RDKit molecule objects, or None if an error occurs.

    Examples
    --------
    >>> mols = read_sdf("ligands.sdf")
    >>> mols = read_sdf("compounds.sdf")

    Notes
    -----
    Invalid molecules (those that RDKit cannot parse) are automatically
    filtered out and not included in the returned list.
    """
    try:
        suppl = Chem.SDMolSupplier(path)
        mols = [mol for mol in suppl if mol is not None]

        logging.info(f"{path} loaded successfully. Total valid molecules: {len(mols)}")

        return mols

    except FileNotFoundError:
        logging.error(f"File not found: {path}")
    except Exception as e:
        logging.error(f"Error reading SDF: {e}")
