"""Debug-only timing helpers for benchmark-friendly pipeline logging."""

from __future__ import annotations

import logging
from contextlib import contextmanager
from time import perf_counter
from typing import Iterator


class DebugTimingCollector:
    """Collect and emit stable debug timing log lines."""

    def __init__(self, logger: logging.Logger | None = None):
        self.logger = logger
        self._entries: list[tuple[str, float]] = []

    def is_enabled(self, logger: logging.Logger | None = None) -> bool:
        active_logger = logger or self.logger
        return bool(active_logger and active_logger.isEnabledFor(logging.DEBUG))

    @contextmanager
    def measure(
        self,
        label: str,
        logger: logging.Logger | None = None,
    ) -> Iterator[None]:
        active_logger = logger or self.logger
        if not self.is_enabled(active_logger):
            yield
            return

        start = perf_counter()
        try:
            yield
        finally:
            self.record(label, perf_counter() - start, logger=active_logger)

    def record(
        self,
        label: str,
        elapsed: float,
        logger: logging.Logger | None = None,
    ) -> None:
        active_logger = logger or self.logger
        if not self.is_enabled(active_logger):
            return

        elapsed = float(elapsed)
        self._entries.append((label, elapsed))
        active_logger.debug("TIMER | %s | %.3fs", label, elapsed)

    def log_summary(self, logger: logging.Logger | None = None) -> None:
        active_logger = logger or self.logger
        if not self.is_enabled(active_logger):
            return

        if not self._entries:
            active_logger.debug("TIMER SUMMARY | no timing entries recorded")
            return

        total = sum(sec for _, sec in self._entries)
        active_logger.debug(
            "TIMER SUMMARY | entries=%d | total=%.3fs",
            len(self._entries),
            total,
        )
        for label, sec in self._entries:
            active_logger.debug("TIMER SUMMARY | %s | %.3fs", label, sec)
