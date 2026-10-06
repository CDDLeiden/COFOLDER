"""Tests for reusable protein-MSA helpers."""

import json
from pathlib import Path

import pytest

from cofolder.modules.input.system import System
from cofolder.modules.runners.msa import (
    capture_generated_msas,
    inject_cached_msas,
    msa_cache_key,
    resolve_declared_msa_paths,
)


def test_capture_and_inject_generated_csv_msa(temp_dir):
    repeat_dir = temp_dir / "repeat"
    generated_dir = repeat_dir / "boltz_results_screen_system" / "msa"
    generated_dir.mkdir(parents=True)
    (generated_dir / "screen_system_0.csv").write_text(
        "key,sequence\n-1,MKRAAT\n-1,MKRAAA\n",
        encoding="utf-8",
    )
    cache_dir = temp_dir / "shared" / "msa"
    source = System(
        system={"sequences": [{"protein": {"id": "A", "sequence": "MKRAAT"}}]}
    )

    assert (
        capture_generated_msas(
            source,
            generated_dir=repeat_dir,
            cache_dir=cache_dir,
            runner_name="boltz2",
        )
        == 1
    )

    target = System(
        system={"sequences": [{"protein": {"id": "Z", "sequence": "MKRAAT"}}]}
    )
    assert inject_cached_msas(target, cache_dir) == 1
    msa_path = target.system["sequences"][0]["protein"]["msa"]
    assert msa_path.startswith(str(cache_dir.resolve()))
    manifest = json.loads((cache_dir / "manifest.json").read_text())
    assert manifest["runner"] == "boltz2"
    assert manifest["version"] == 2


def test_cache_identity_includes_msa_generation_settings(temp_dir):
    generated_dir = temp_dir / "repeat" / "msa"
    generated_dir.mkdir(parents=True)
    (generated_dir / "generated.csv").write_text(
        "key,sequence\n-1,MKRAAT\n",
        encoding="utf-8",
    )
    cache_dir = temp_dir / "shared" / "msa" / "boltz2"
    source = System(
        system={"sequences": [{"protein": {"id": "A", "sequence": "MKRAAT"}}]}
    )
    greedy = {"msa_pairing_strategy": "greedy", "max_msa_seqs": 8192}
    complete = {"msa_pairing_strategy": "complete", "max_msa_seqs": 8192}

    assert (
        capture_generated_msas(
            source,
            generated_dir=generated_dir,
            cache_dir=cache_dir,
            runner_name="boltz2",
            settings=greedy,
        )
        == 1
    )

    incompatible = System(
        system={"sequences": [{"protein": {"id": "A", "sequence": "MKRAAT"}}]}
    )
    assert inject_cached_msas(incompatible, cache_dir, settings=complete) == 0
    assert "msa" not in incompatible.system["sequences"][0]["protein"]

    compatible = System(
        system={"sequences": [{"protein": {"id": "A", "sequence": "MKRAAT"}}]}
    )
    assert inject_cached_msas(compatible, cache_dir, settings=greedy) == 1
    entry = next(
        iter(json.loads((cache_dir / "manifest.json").read_text())["proteins"].values())
    )
    assert entry["settings"] == greedy
    assert entry["settings_sha256"]


def test_mixed_precomputed_and_generated_msas_only_injects_missing_entity(temp_dir):
    supplied = temp_dir / "supplied.a3m"
    supplied.write_text(">query\nAAAA\n", encoding="utf-8")
    generated_dir = temp_dir / "repeat" / "msa"
    generated_dir.mkdir(parents=True)
    (generated_dir / "screen_system_1.csv").write_text(
        "key,sequence\n-1,BBBB\n",
        encoding="utf-8",
    )
    cache_dir = temp_dir / "shared" / "msa"
    target = System(
        system={
            "sequences": [
                {"protein": {"id": "A", "sequence": "AAAA", "msa": str(supplied)}},
                {"protein": {"id": "B", "sequence": "BBBB"}},
            ]
        }
    )

    assert (
        capture_generated_msas(
            target,
            generated_dir=generated_dir,
            cache_dir=cache_dir,
            runner_name="boltz2",
        )
        == 1
    )
    assert inject_cached_msas(target, cache_dir) == 1
    assert target.system["sequences"][0]["protein"]["msa"] == str(supplied)
    assert Path(target.system["sequences"][1]["protein"]["msa"]).is_file()


def test_corrupt_or_sequence_mismatched_cached_msa_is_not_injected(temp_dir):
    cache_dir = temp_dir / "msa"
    cache_dir.mkdir()
    cached = cache_dir / "protein_bad.csv"
    cached.write_text("key,sequence\n-1,DIFFERENT\n", encoding="utf-8")
    (cache_dir / "manifest.json").write_text(
        json.dumps(
            {
                "version": 1,
                "proteins": {msa_cache_key("MKRAAT"): {"path": cached.name}},
            }
        ),
        encoding="utf-8",
    )
    target = System(
        system={"sequences": [{"protein": {"id": "A", "sequence": "MKRAAT"}}]}
    )

    with pytest.raises(ValueError, match="does not match its protein sequence"):
        inject_cached_msas(target, cache_dir)
    assert "msa" not in target.system["sequences"][0]["protein"]


def test_relative_precomputed_msa_is_resolved_from_source_yaml_directory(temp_dir):
    msa_path = temp_dir / "inputs" / "protein.a3m"
    msa_path.parent.mkdir()
    msa_path.write_text(">query\nMKRAAT\n", encoding="utf-8")
    target = System(
        system={
            "sequences": [
                {"protein": {"id": "A", "sequence": "MKRAAT", "msa": "protein.a3m"}}
            ]
        }
    )

    assert resolve_declared_msa_paths(target, base_dir=msa_path.parent) == 1
    assert target.system["sequences"][0]["protein"]["msa"] == str(msa_path.resolve())


def test_explicit_empty_msa_sentinel_is_preserved(temp_dir):
    target = System(
        system={
            "sequences": [
                {"protein": {"id": "A", "sequence": "MKRAAT", "msa": "empty"}}
            ]
        }
    )

    assert resolve_declared_msa_paths(target, base_dir=temp_dir) == 0
    assert target.system["sequences"][0]["protein"]["msa"] == "empty"
