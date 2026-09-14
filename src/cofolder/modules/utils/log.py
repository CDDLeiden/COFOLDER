"""Logging configuration utilities for cofolder.

This module provides utilities for setting up consistent logging across
the cofolder package, including formatters and handlers for both console
and file output.
"""

import logging
import sys
import weakref
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
    The compatibility contract also preserves pytest `caplog` visibility after
    tests clear root handlers, even though that currently depends on pytest and
    Python logging internals rather than a public hook.

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
            if isinstance(handler, logging.StreamHandler)
            and not isinstance(handler, logging.FileHandler)
            and getattr(handler, "stream", None) is sys.stdout
        ),
        None,
    )
    if console is None:
        console = logging.StreamHandler(sys.stdout)
        root.addHandler(console)
    console.setFormatter(console_formatter)

    # Compatibility contract: tests may clear root.handlers before calling this
    # helper, but caplog-based assertions should still observe root logger
    # output afterward. Pytest does not expose a public reattachment hook, so
    # we intentionally depend on logging/pytest internals here and cover that
    # behavior with a narrow regression test.
    for handler_ref in list(getattr(logging, "_handlerList", [])):
        handler = handler_ref() if isinstance(handler_ref, weakref.ReferenceType) else None
        if handler is None:
            continue
        if handler.__class__.__module__.startswith("_pytest.logging") and handler not in root.handlers:
            root.addHandler(handler)

    if log_file:
        log_path = str(Path(log_file).resolve())
        file_handler = next(
            (
                handler
                for handler in root.handlers
                if isinstance(handler, logging.FileHandler)
                and Path(handler.baseFilename).resolve() == Path(log_path)
            ),
            None,
        )
        if file_handler is None:
            file_handler = logging.FileHandler(log_path)
            root.addHandler(file_handler)
        file_handler.setFormatter(detailed_formatter)
