# Script containing general functions
import os
import pandas as pd
import random
from typing import Optional, List

import logging

logger = logging.getLogger(__name__)

def create_dir(path: str):
    """Ensure that a directory exists, creating it if necessary.

    Parameters
    ----------
    path : str
        The directory path to check and/or create.

    Examples
    --------
    >>> create_dir("./output")
    >>> create_dir("/tmp/my_project/results")
    """
    logger.debug(f"Checking directory: '{path}'")

    if not os.path.isdir(path):
        logger.debug(f"Directory does not exist. Creating: '{path}'")
        os.makedirs(path, exist_ok=True)
        logger.info(f"Created working directory '{path}'")
    else:
        logger.debug(f"Directory already exists: '{path}'")

def set_dir(path):
    """Create a directory if it doesn't exist.

    Parameters
    ----------
    path : str
        The directory path to create.

    Notes
    -----
    This is a legacy function. Consider using `create_dir()` instead,
    which has more detailed logging and type hints.

    Examples
    --------
    >>> set_dir("./experiments")
    """
    if not (os.path.isdir(path)):
        os.makedirs(path)
        logger.info("Created working dir {0}".format(path))

def parse_list_as_str(
    list_as_str: str,
    separator: str = ",",
    item_type=str,
    expected_length: int = None
):
    """Convert a separated string into a typed Python list.

    Parses a command-line string value containing delimited items and
    converts each item to the specified type. Useful for parsing CLI
    arguments like device lists, seeds, or other numeric sequences.

    Parameters
    ----------
    list_as_str : str
        The raw string containing a list (e.g., "1,2,3").
    separator : str, default=","
        Delimiter used to separate items in the string.
    item_type : type, default=str
        The expected type of each element (str, int, float, etc.).
    expected_length : int, optional
        If provided, the parsed list must match this exact length.

    Returns
    -------
    list
        The processed and type-casted list.

    Raises
    ------
    ValueError
        If type conversion fails or if the parsed list length doesn't
        match expected_length.

    Examples
    --------
    >>> parse_list_as_str("1,2,3", item_type=int)
    [1, 2, 3]

    >>> parse_list_as_str("a;b;c", separator=";")
    ['a', 'b', 'c']

    >>> parse_list_as_str("1.5,2.5", item_type=float, expected_length=2)
    [1.5, 2.5]
    """
    logger.debug(f"Parsing list: raw='{list_as_str}', "
                 f"separator='{separator}', item_type={item_type.__name__}, "
                 f"expected_length={expected_length}")

    # Split into items
    items = [item.strip() for item in list_as_str.split(separator) if item.strip()]

    # Convert each item to specified type
    try:
        typed_items = [item_type(item) for item in items]
    except Exception:
        raise ValueError(
            f"Failed to convert items to type {item_type.__name__}: {items}"
        )

    # Check list length (if included)
    if expected_length is not None and len(typed_items) != expected_length:
        raise ValueError(
            f"List length must be {expected_length}, but got {len(typed_items)}"
        )
    
    logger.debug(f"Finished parsing. Result: {typed_items}")

    return typed_items

def get_seeds(
    repeats: int,
    seed: Optional[int] = None,
    logger: Optional[logging.Logger] = None
) -> (int, List[int]):
    """
    Generate a global seed and run seeds for repeated workflow runs.

    Parameters
    ----------
    repeats : int
        Number of repeats for which seeds should be generated.
    seed : int, optional
        Global seed to use. If None, a random global seed is generated.
    logger : logging.Logger, optional
        Logger for informational/debug messages. If None, logging is skipped.

    Returns
    -------
    global_seed : int
        The global seed used for generating run seeds.
    run_seeds : list of int
        List of seeds for each run (length == repeats).

    Notes
    -----
    - If repeats == 1, the run seed list contains only the global seed.
    - If repeats > 1, run seeds are generated deterministically from the global seed.
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    # Generate global seed if not provided
    if seed is None:
        global_seed = generate_seeds(num_seeds=1, seed=None)[0]
        logger.info("No global seed provided. Generated random global seed: %d", global_seed)
    else:
        global_seed = seed
        logger.info("Using provided global seed: %d", global_seed)

    # Generate run seeds based on repeats
    if repeats == 1:
        run_seeds = [global_seed]
        logger.debug("Single repeat: using global seed as run seed: %s", run_seeds)
    else:
        run_seeds = generate_seeds(num_seeds=repeats, seed=global_seed)
        logger.debug(
            "Multiple repeats: %d run seeds generated from global seed %d: %s",
            repeats,
            global_seed,
            run_seeds
        )

    logger.info("Run seeds to be used for this workflow: %s", run_seeds)
    return global_seed, run_seeds

def generate_seeds(num_seeds: int, seed: Optional[int] = None) -> List[int]:
    """
    Generate a list of random integer seeds for reproducibility.

    Parameters
    ----------
    num_seeds : int
        Number of seeds to generate (must be >= 1).
    seed : int or None
        Optional global seed for reproducibility. If None, a random global
        seed is generated.

    Returns
    -------
    List[int]
        List of generated seeds.
    """
    if num_seeds < 1:
        raise ValueError(f"num_seeds must be >= 1, got {num_seeds}")

    # Generate or use the global seed
    global_seed = seed if seed is not None else random.randint(0, 2**32 - 1)
    rng = random.Random(global_seed)
    
    return [rng.randint(0, 2**32 - 1) for _ in range(num_seeds)]

def drop_and_log_nans(df: pd.DataFrame, cols: list, context: str = "") -> pd.DataFrame:
    """Drop rows with NaN values and log the number of rows removed.

    Removes rows containing NaN values in any of the specified columns
    and logs the count of removed rows for debugging and quality control.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to process.
    cols : list of str
        Column names to check for NaN values.
    context : str, default=""
        Optional context string to include in the log message for clarity.

    Returns
    -------
    pd.DataFrame
        DataFrame with rows containing NaN in specified columns removed.

    Examples
    --------
    >>> import pandas as pd
    >>> import numpy as np
    >>> df = pd.DataFrame({'a': [1, 2, np.nan], 'b': [4, np.nan, 6]})
    >>> drop_and_log_nans(df, ['a', 'b'], context="correlation analysis")
         a    b
    0  1.0  4.0
    """
    nan_count = df[cols].isna().any(axis=1).sum()
    if nan_count > 0:
        logging.info(f"Removed {nan_count} rows with NaN in {cols} {f'for {context}' if context else ''}.")
    return df.dropna(subset=cols)

