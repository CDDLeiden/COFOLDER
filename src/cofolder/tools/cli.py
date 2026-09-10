"""Dispatcher for installed COFOLDER setup and data-preparation tools."""

from __future__ import annotations

import argparse
from collections.abc import Callable, Sequence
import importlib
from pathlib import Path
import sys

from cofolder.resources.examples import copy_examples


COMMANDS = {
    "copy-examples": "Copy the packaged runnable examples to a workspace.",
    "setup-openfold3": "Prepare the OpenFold3 cache and checkpoints.",
    "install-mmseqs": "Install a verified vendored MMseqs2 executable.",
    "fetch-bias-training-data": "Fetch CCD and optional MMseqs2 reference data.",
    "build-bias-training-data": "Build bias-training reference tables.",
}

MODULES = {
    "setup-openfold3": "cofolder.tools.setup_openfold3",
    "install-mmseqs": "cofolder.tools.install_mmseqs",
    "fetch-bias-training-data": "cofolder.tools.fetch_bias_training_data",
    "build-bias-training-data": "cofolder.tools.build_bias_training_data",
}


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cofolder-tools",
        description="Installed COFOLDER example, setup, and data-preparation tools.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("command", nargs="?", choices=COMMANDS)
    parser.epilog = "Commands:\n  " + "\n  ".join(
        f"{name:<26} {description}" for name, description in COMMANDS.items()
    )
    return parser


def _copy_examples_main(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="cofolder-tools copy-examples",
        description=COMMANDS["copy-examples"],
    )
    parser.add_argument("destination", type=Path)
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace managed example files that already exist.",
    )
    args = parser.parse_args(argv)
    try:
        destination = copy_examples(args.destination, overwrite=args.overwrite)
    except FileExistsError as exc:
        parser.error(str(exc))
    print(f"[done] copied_examples={destination}")
    return 0


def _load_main(module_name: str) -> Callable[[Sequence[str] | None], int]:
    module = importlib.import_module(module_name)
    return module.main


def main(argv: Sequence[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    parser = _parser()
    if not args:
        parser.print_help()
        return 0
    if args[0] in {"-h", "--help"}:
        parser.print_help()
        return 0

    command = args.pop(0)
    if command not in COMMANDS:
        parser.error(f"argument command: invalid choice: '{command}'")
    if command == "copy-examples":
        return _copy_examples_main(args)
    return _load_main(MODULES[command])(args)


if __name__ == "__main__":
    raise SystemExit(main())
