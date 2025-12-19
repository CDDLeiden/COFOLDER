# src/boltz_lab/modules/utils/log.py

import logging
import sys
from pathlib import Path
from typing import Optional

DEFAULT_FORMAT = (
    "%(asctime)s | "
    "%(levelname)-8s | "
    "%(name)s | "
    "%(filename)s:%(lineno)d | "
    "%(funcName)s | "
    "%(message)s"
)

def setup_root_logger(
    level: int = logging.INFO,
    log_file: Optional[Path] = None
) -> None:
    """
    Configure the root logger exactly once.
    Safe for CLI, notebooks, and programmatic use.
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