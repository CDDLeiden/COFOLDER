#!/usr/bin/env python3
"""Compatibility wrapper for the packaged bias-data fetcher."""

from cofolder.tools.fetch_bias_training_data import main


if __name__ == "__main__":
    raise SystemExit(main())
