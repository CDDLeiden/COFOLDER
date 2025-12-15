import os
import logging

logger = logging.getLogger("evaluate")

def create_dir(path: str):
    """
    Ensure that a directory exists. If it does not exist, create it.

    Args:
        path (str): The directory path to check and/or create.
    """
    logger.debug(f"Checking directory: '{path}'")

    if not os.path.isdir(path):
        logger.debug(f"Directory does not exist. Creating: '{path}'")
        os.makedirs(path, exist_ok=True)
        logger.info(f"Created working directory '{path}'")
    else:
        logger.debug(f"Directory already exists: '{path}'")

def parse_list_as_str(
    list_as_str: str,
    separator: str = ",",
    item_type=str,
    expected_length: int = None
):
    """
    Convert a comma-separated command-line value into a typed Python list.

    Args:
        list_as_str (str): The raw string containing a list (e.g. "1,2,3").
        separator (str): Separator used in the list (default: comma).
        item_type (type): The expected type of each element (str, int, float, etc.).
        expected_length (int, optional): If provided, list must match this length.
        

    Returns:
        list: The processed and type-casted list.

    Raises:
        ValueError: If type conversion fails or length is incorrect.
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