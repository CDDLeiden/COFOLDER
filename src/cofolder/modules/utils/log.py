"""Logging configuration utilities for cofolder.

This module provides utilities for setting up consistent logging across
the cofolder package, including formatters and handlers for both console
and file output.
"""

import logging
import sys
from pathlib import Path
from typing import Optional

CONCISE_FORMAT = "%(levelname)s | %(message)s"
"""str: Concise format used for normal console progress."""

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
    2025-01-15 10:30:45,123 | INFO     | cofolder.validate | validate.py:42 | run | Starting validation
"""

def setup_root_logger(
    level: int = logging.INFO,
    log_file: Optional[Path] = None
) -> None:
    """Configure the root logger with console and optional file handlers.

    Sets up the root logger with a standardized format for consistent logging
    across the entire application. Repeated calls update the root log level and
    ensure the requested console/file handlers exist without duplicating them.
    Only handlers created by COFOLDER are updated; unrelated application handlers
    remain untouched.

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

    detailed_formatter = logging.Formatter(DEFAULT_FORMAT)
    console_formatter = logging.Formatter(
        DEFAULT_FORMAT if level <= logging.DEBUG else CONCISE_FORMAT
    )

    root.setLevel(level)

    console = next(
        (
            handler
            for handler in root.handlers
            if getattr(handler, "_cofolder_console_handler", False)
        ),
        None,
    )
    if console is None:
        console = logging.StreamHandler(sys.stdout)
        console._cofolder_console_handler = True
        root.addHandler(console)
    console.setFormatter(console_formatter)

    if log_file:
        log_path = str(Path(log_file).resolve())
        file_handler = next(
            (
                handler
                for handler in root.handlers
                if isinstance(handler, logging.FileHandler)
                and getattr(handler, "_cofolder_file_handler", False)
                and Path(handler.baseFilename).resolve() == Path(log_path)
            ),
            None,
        )
        if file_handler is None:
            file_handler = logging.FileHandler(log_path)
            file_handler._cofolder_file_handler = True
            root.addHandler(file_handler)
        file_handler.setFormatter(detailed_formatter)
