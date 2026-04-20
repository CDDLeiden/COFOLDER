from __future__ import annotations

import logging
import re
import subprocess
from time import perf_counter
from typing import Sequence

from cofolder.modules.utils.timing import DebugTimingCollector

logger = logging.getLogger(__name__)

_ANSI_ESCAPE_RE = re.compile(r"\x1b\[[0-9;?]*[ -/]*[@-~]")


def _normalize_boltz_output_line(line: str) -> str:
    cleaned = _ANSI_ESCAPE_RE.sub("", str(line or ""))
    return cleaned.replace("\r", "").strip()


def _parse_boltz_stage_timings(
    timed_lines: Sequence[tuple[float, str]],
    total_elapsed: float,
) -> dict[str, float]:
    """Parse stage timings from timestamped Boltz stdout lines."""

    def finalize_msa_window() -> None:
        nonlocal active_msa_start, active_msa_complete, msa_total
        if (
            active_msa_start is not None
            and active_msa_complete is not None
            and active_msa_complete >= active_msa_start
        ):
            msa_total += active_msa_complete - active_msa_start
        active_msa_start = None
        active_msa_complete = None

    msa_total = 0.0
    active_msa_start: float | None = None
    active_msa_complete: float | None = None
    affinity_start: float | None = None

    for offset, raw_line in timed_lines:
        line = _normalize_boltz_output_line(raw_line)
        if not line:
            continue

        if line.startswith("Calling MSA server for target"):
            finalize_msa_window()
            active_msa_start = offset
            active_msa_complete = None
            continue

        if active_msa_start is not None and "COMPLETE: 100%" in line:
            active_msa_complete = offset

        if line.startswith("Predicting property: affinity") or line.startswith(
            "Running affinity prediction for"
        ):
            finalize_msa_window()

        if line.startswith("Running affinity prediction for") and affinity_start is None:
            affinity_start = offset

    finalize_msa_window()

    timings: dict[str, float] = {}
    if msa_total > 0.0:
        timings["boltz.msa"] = msa_total
    if affinity_start is not None and total_elapsed >= affinity_start:
        timings["boltz.affinity_prediction"] = total_elapsed - affinity_start
    return timings


def _record_boltz_timings(
    timings: DebugTimingCollector | None,
    label_prefix: str | None,
    total_elapsed: float,
    timed_lines: Sequence[tuple[float, str]],
) -> None:
    if timings is None:
        return

    active_logger = timings.logger or logger
    prefix = f"{label_prefix}." if label_prefix else ""
    timings.record(f"{prefix}boltz.total", total_elapsed, logger=active_logger)
    for label, elapsed in _parse_boltz_stage_timings(timed_lines, total_elapsed).items():
        timings.record(f"{prefix}{label}", elapsed, logger=active_logger)


def run_boltz(
    cmd: Sequence[str],
    check: bool = True,
    timings: DebugTimingCollector | None = None,
    label_prefix: str | None = None,
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
    start_time = perf_counter()
    logger.info("Running: %s", " ".join(cmd))

    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        universal_newlines=True,
    )

    output_lines: list[str] = []
    timed_lines: list[tuple[float, str]] = []
    assert process.stdout is not None
    for line in process.stdout:
        line = line.rstrip()
        output_lines.append(line)
        timed_lines.append((perf_counter() - start_time, line))
        logger.info(line)

    process.stdout.close()
    retcode = process.wait()
    total_elapsed = perf_counter() - start_time

    logger.info("Boltz execution completed in %.2f seconds", total_elapsed)
    _record_boltz_timings(
        timings=timings,
        label_prefix=label_prefix,
        total_elapsed=total_elapsed,
        timed_lines=timed_lines,
    )

    if check and retcode != 0:
        raise subprocess.CalledProcessError(retcode, cmd, "\n".join(output_lines))

    return subprocess.CompletedProcess(
        args=cmd,
        returncode=retcode,
        stdout="\n".join(output_lines),
    )
