# Script containing general functions

import logging
import os

def set_dir(path):
    if path:
        os.makedirs(path, exist_ok=True)
        logging.info(f"Created or verified existence of directory: {path}")
    else:
        logging.warning("No path ound in wrapper entry")

    return path