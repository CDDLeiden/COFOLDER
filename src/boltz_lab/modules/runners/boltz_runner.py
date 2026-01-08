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
    Execute a Boltz command via subprocess and log stdout/stderr in real-time.

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

    # Use Popen to capture stdout/stderr line by line
    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        universal_newlines=True
    )

    # Log output in real-time
    output_lines = []
    for line in process.stdout:
        line = line.rstrip()
        output_lines.append(line)
        logger.info(line)

    process.stdout.close()
    retcode = process.wait()

    logger.info(
        "Boltz execution completed in %.2f seconds",
        time.time() - start_time
    )

    if check and retcode != 0:
        raise subprocess.CalledProcessError(retcode, cmd, "\n".join(output_lines))
    
    return subprocess.CompletedProcess(args=cmd, returncode=retcode, stdout="\n".join(output_lines))