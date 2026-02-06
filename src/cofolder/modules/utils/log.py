"""Logging configuration utilities for cofolder.

This module provides utilities for setting up consistent logging across
the cofolder package, including formatters and handlers for both console
and file output.
"""

import logging
import sys
from pathlib import Path
from typing import Optional

DEFAULT_FORMAT = (
    "%(asctime)s | "
    "%(levelname)-8s | "
    "%(filename)s:%(lineno)d | "
    "%(funcName)s | "
    "%(message)s"
)
"""str: Default log format for all cofolder loggers.

Format includes timestamp, log level, logger name, file location,
function name, and the actual message, separated by pipes for readability.

Example output:
    2025-01-15 10:30:45,123 | INFO     | cofolder.predict | predict.py:42 | run | Starting prediction
"""

def setup_root_logger(
    level: int = logging.INFO,
    log_file: Optional[Path] = None
) -> None:
    """Configure the root logger with console and optional file handlers.

    Sets up the root logger with a standardized format for consistent logging
    across the entire application. If the root logger already has handlers,
    this function returns immediately to prevent duplicate handler registration.
    This makes it safe to call multiple times from different modules.

    Parameters
    ----------
    level : int, default=logging.INFO
        The logging level to set for the root logger. Common values:
        - logging.DEBUG (10): Detailed information for debugging
        - logging.INFO (20): General informational messages
        - logging.WARNING (30): Warning messages
        - logging.ERROR (40): Error messages
        - logging.CRITICAL (50): Critical errors
    log_file : str or Path, optional
        If provided, adds a FileHandler that writes logs to this path.
        The directory must exist or be creatable.

    Returns
    -------
    None


    See Also
    --------
    DEFAULT_FORMAT : The standard log format used by all handlers
    """
    root = logging.getLogger()

    if root.handlers:
        return

    formatter = logging.Formatter(DEFAULT_FORMAT)

    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(formatter)

    root.setLevel(level)
    root.addHandler(console)

    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        root.addHandler(file_handler)