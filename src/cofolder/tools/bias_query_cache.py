"""Inspect or safely invalidate the shared bias-query cache."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from datetime import datetime, timezone
import fcntl
import json
from pathlib import Path

from cofolder.modules.analytics.bias_database import (
    BIAS_QUERY_CACHE_SCHEMA_VERSION,
    DEFAULT_BIAS_QUERY_CACHE,
    bias_query_cache_readiness,
)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="cofolder-tools bias-query-cache")
    parser.add_argument("action", choices=("inspect", "invalidate"))
    parser.add_argument("--cache_path", type=Path, default=DEFAULT_BIAS_QUERY_CACHE)
    parser.add_argument("--kind", choices=("protein", "ligand", "all"), default="all")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    root, ready, state = bias_query_cache_readiness(args.cache_path)
    if args.action == "inspect":
        counts = {
            kind: sum(1 for _ in (root / kind).glob("*/*/manifest.json"))
            if (root / kind).is_dir() else 0
            for kind in ("protein", "ligand")
        }
        print(json.dumps({"path": str(root), "ready": ready, "state": state, "entries": counts}))
        return 0 if ready else 1
    marker = root / "cache.json"
    if not marker.is_file():
        raise ValueError(f"Refusing to invalidate an unmarked cache directory: {root}")
    payload = json.loads(marker.read_text(encoding="utf-8"))
    if payload.get("schema_version") != BIAS_QUERY_CACHE_SCHEMA_VERSION:
        raise ValueError(f"Refusing to invalidate unsupported cache schema: {marker}")
    lock_path = root / ".invalidate.lock"
    with lock_path.open("a+", encoding="utf-8") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        recovery = root.parent / (
            f"{root.name}.invalidated-" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        )
        recovery.mkdir()
        kinds = ("protein", "ligand") if args.kind == "all" else (args.kind,)
        for kind in kinds:
            source = root / kind
            if source.exists():
                source.replace(recovery / kind)
        print(f"[done] invalidated={','.join(kinds)} recovery={recovery}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
