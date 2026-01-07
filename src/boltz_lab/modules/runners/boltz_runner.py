import subprocess
import time
import logging
from typing import Sequence

logger = logging.getLogger(__name__)

def run_boltz(
    cmd: Sequence[str],
    check: bool = True
) -> subprocess.CompletedProcess:
    """
    Execute a Boltz command via subprocess.

    Parameters
    ----------
    cmd : Sequence[str]
        Fully constructed Boltz command.
    check : bool, optional
        Raise exception on non-zero exit code.

    Returns
    -------
    subprocess.CompletedProcess
    """
    start_time = time.time()
    logger.info("Running: %s", " ".join(cmd))

    result = subprocess.run(cmd, check=check)

    logger.info(
        "Boltz execution completed in %.2f seconds",
        time.time() - start_time
    )

    return result