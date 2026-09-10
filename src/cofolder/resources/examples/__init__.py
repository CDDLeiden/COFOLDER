"""Access to the runnable example bundle shipped with COFOLDER."""

from __future__ import annotations

from importlib.resources import files
from importlib.resources.abc import Traversable
from pathlib import Path
import shutil


def examples_root() -> Traversable:
    """Return the packaged example collection without assuming a filesystem path."""

    return files(__package__)


def copy_examples(
    destination: str | Path,
    *,
    overwrite: bool = False,
) -> Path:
    """Copy all packaged examples into *destination*.

    Existing managed files are rejected before anything is written unless
    ``overwrite`` is true. Files in the destination that are not part of the
    example bundle are never removed.
    """

    destination_path = Path(destination).expanduser()
    resources = sorted(
        (
            item
            for item in examples_root().iterdir()
            if item.is_file() and item.name != "__init__.py"
        ),
        key=lambda item: item.name,
    )
    conflicts = [destination_path / item.name for item in resources]
    conflicts = [path for path in conflicts if path.exists()]
    if conflicts and not overwrite:
        names = ", ".join(path.name for path in conflicts)
        raise FileExistsError(
            f"Example files already exist in '{destination_path}': {names}. "
            "Pass overwrite=True to replace them."
        )
    non_files = [path for path in conflicts if not path.is_file()]
    if non_files:
        names = ", ".join(path.name for path in non_files)
        raise FileExistsError(
            f"Cannot overwrite non-file example targets in '{destination_path}': "
            f"{names}."
        )

    destination_path.mkdir(parents=True, exist_ok=True)
    for resource in resources:
        with resource.open("rb") as source, (destination_path / resource.name).open(
            "wb"
        ) as target:
            shutil.copyfileobj(source, target)
    return destination_path.resolve()


__all__ = ["copy_examples", "examples_root"]
