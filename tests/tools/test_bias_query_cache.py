from __future__ import annotations

from pathlib import Path

import pandas as pd

from cofolder.modules.analytics import bias_database
from cofolder.tools import bias_query_cache


def test_inspects_and_recoverably_invalidates_marked_cache(temp_dir, capsys):
    root = temp_dir / "cache"
    identity = {"fixture": True}
    key = bias_database._cache_key("protein", identity)
    bias_database._cached_frame(
        root=root, kind="protein", key=key, identity=identity,
        builder=lambda: pd.DataFrame([{"value": 1}]),
    )
    assert bias_query_cache.main(["inspect", "--cache_path", str(root)]) == 0
    assert '"protein": 1' in capsys.readouterr().out
    assert bias_query_cache.main([
        "invalidate", "--cache_path", str(root), "--kind", "protein"
    ]) == 0
    output = capsys.readouterr().out
    recovery = Path(output.split("recovery=", 1)[1].strip())
    assert (recovery / "protein").is_dir()
    assert not (root / "protein").exists()
