"""Internal executable discovery with inspectable resolution diagnostics."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass
import logging
import os
from pathlib import Path
import shutil


LOGGER = logging.getLogger(__name__)
_AUTODETECT = object()


@dataclass(frozen=True)
class ExecutableCandidate:
    """One executable candidate considered during resolution."""

    source: str
    value: str
    resolved_path: str | None
    outcome: str

    def describe(self) -> str:
        location = self.resolved_path or self.value
        return f"{self.source}={location!r}: {self.outcome}"


@dataclass(frozen=True)
class ExecutableResolution:
    """Selected executable and the ordered candidates considered for it."""

    path: str | None
    source: str | None
    candidates: tuple[ExecutableCandidate, ...]

    def diagnostic_text(self) -> str:
        if not self.candidates:
            return "no executable candidates were available"
        return "; ".join(candidate.describe() for candidate in self.candidates)


def _source_checkout_root() -> Path | None:
    candidate = Path(__file__).resolve().parents[4]
    if (candidate / "pyproject.toml").is_file() and (
        candidate / "src/cofolder"
    ).is_dir():
        return candidate
    return None


def _is_path_value(value: str) -> bool:
    expanded = os.path.expanduser(value)
    return (
        Path(expanded).is_absolute()
        or value.startswith((".", "~"))
        or os.sep in value
        or (os.altsep is not None and os.altsep in value)
    )


def _inspect_path(source: str, value: str, path: Path) -> ExecutableCandidate:
    expanded = path.expanduser()
    if not expanded.exists():
        return ExecutableCandidate(source, value, None, "missing")
    if not expanded.is_file():
        return ExecutableCandidate(
            source, value, str(expanded.absolute()), "not a regular file"
        )
    resolved = str(expanded.resolve())
    if not os.access(expanded, os.X_OK):
        return ExecutableCandidate(source, value, resolved, "not executable")
    return ExecutableCandidate(source, value, resolved, "selected")


def _inspect_candidate(
    source: str,
    value: str,
    which: Callable[[str], str | None],
) -> ExecutableCandidate:
    expanded = Path(value).expanduser()
    if expanded.exists() or _is_path_value(value):
        return _inspect_path(source, value, expanded)

    found = which(value)
    if found is None:
        return ExecutableCandidate(source, value, None, "unresolved command")
    return _inspect_path(source, value, Path(found))


def resolve_mmseqs_executable(
    mmseqs_bin: str | None = None,
    *,
    environ: Mapping[str, str] | None = None,
    home: Path | None = None,
    source_root: Path | None | object = _AUTODETECT,
    which: Callable[[str], str | None] = shutil.which,
) -> ExecutableResolution:
    """Resolve MMseqs2 using COFOLDER's stable, ordered lookup policy.

    The optional environment, home, source-root and ``which`` arguments make the
    installed and source layouts testable without changing process-global state.
    They are internal implementation seams rather than public configuration.
    """
    environment = os.environ if environ is None else environ
    user_home = Path.home() if home is None else Path(home)
    checkout = _source_checkout_root() if source_root is _AUTODETECT else source_root

    ordered: list[tuple[str, str]] = []
    if mmseqs_bin:
        ordered.append(("explicit argument", str(mmseqs_bin)))
    env_value = environment.get("COFOLDER_MMSEQS_BIN")
    if env_value:
        ordered.append(("COFOLDER_MMSEQS_BIN", env_value))
    ordered.append(
        (
            "user vendor",
            str(user_home / ".cofolder/vendor/mmseqs/bin/mmseqs"),
        )
    )
    if checkout is not None:
        source_vendor = Path(checkout) / "vendor/mmseqs"
        ordered.extend(
            [
                ("source vendor bin", str(source_vendor / "bin/mmseqs")),
                ("source vendor root", str(source_vendor / "mmseqs")),
            ]
        )
    ordered.append(("PATH", "mmseqs"))

    attempts: list[ExecutableCandidate] = []
    seen: set[str] = set()
    for source, value in ordered:
        candidate = _inspect_candidate(source, value, which)
        identity = candidate.resolved_path or str(Path(value).expanduser().absolute())
        if identity in seen:
            continue
        seen.add(identity)
        attempts.append(candidate)
        if candidate.outcome == "selected":
            resolution = ExecutableResolution(
                path=candidate.resolved_path,
                source=source,
                candidates=tuple(attempts),
            )
            for attempt in attempts[:-1]:
                LOGGER.debug("Rejected MMseqs candidate: %s", attempt.describe())
            LOGGER.debug(
                "Selected MMseqs executable from %s: %s",
                source,
                candidate.resolved_path,
            )
            return resolution

    resolution = ExecutableResolution(None, None, tuple(attempts))
    for attempt in attempts:
        LOGGER.debug("Rejected MMseqs candidate: %s", attempt.describe())
    LOGGER.debug("No MMseqs executable selected")
    return resolution

