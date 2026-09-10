"""Installed entry point for the existing bias-training builder."""

from __future__ import annotations

from collections.abc import Sequence


def main(argv: Sequence[str] | None = None) -> int:
    """Delegate to the retained analytics implementation during R09 migration."""

    from cofolder.modules.analytics.build_bias_training_data import main as builder_main

    return builder_main(argv, prog="cofolder-tools build-bias-training-data")


if __name__ == "__main__":
    raise SystemExit(main())
