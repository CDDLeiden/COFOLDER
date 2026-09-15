"""Cross-revision semantic parity coverage for the S6 responsibility extracts."""

from __future__ import annotations

import ast
import csv
import hashlib
import importlib
import json
import subprocess
import sys
from collections import Counter
from contextlib import ExitStack
from dataclasses import asdict, replace
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest
import yaml

from cofolder.modules.contracts import (
    BackendVersionStatus,
    OutputIdentity,
    RepeatSeedProvenance,
    RunnerBackendIdentity,
    SeedOrigin,
    SeedPlan,
    WorkflowKind,
    bundle_from_frames,
    write_public_bundle,
)
from cofolder.modules.runners import PlannedExecution, RunnerExecutionPlan
from cofolder.recipes.bias import Bias
from cofolder.recipes.oracle import Oracle, OracleGate, OracleGatePolicy
from cofolder.recipes.screen import Screen
from cofolder.recipes.validate import Validate
from tests.recipes.test_validate import _patch_validate_pipeline


EXPECTATIONS_PATH = (
    Path(__file__).parents[1] / "fixtures" / "structural_parity_s6.json"
)
TEXT_SUFFIXES = {".cif", ".csv", ".json", ".jsonl", ".pdb", ".txt", ".yaml"}


def _write_inputs(root: Path) -> tuple[Path, Path]:
    system_path = root / "system.yaml"
    options_path = root / "options.yaml"
    system_path.write_text(
        yaml.safe_dump(
            {
                "sequences": [
                    {"protein": {"id": "A", "sequence": "MKRAAT"}},
                    {"ligand": {"id": "B", "smiles": "CCO"}},
                ]
            }
        ),
        encoding="utf-8",
    )
    options_path.write_text(
        yaml.safe_dump(
            {
                "version": 1,
                "runtime": {"diffusion_samples": 2},
                "runner": {"recycling_steps": 3},
            }
        ),
        encoding="utf-8",
    )
    return system_path, options_path


def _fixed_plan() -> RunnerExecutionPlan:
    repeats = (
        RepeatSeedProvenance(1, 17, 17, 17, 17, SeedOrigin.USER_SPECIFIED),
        RepeatSeedProvenance(2, 17, 17, 23, 23, SeedOrigin.USER_SPECIFIED),
    )
    return RunnerExecutionPlan(
        backend=RunnerBackendIdentity(
            runner_name="boltz2",
            backend_name="boltz",
            version="2.2.1",
            version_status=BackendVersionStatus.DETECTED,
        ),
        seed_plan=SeedPlan(17, 17, SeedOrigin.USER_SPECIFIED, repeats),
        executions=tuple(
            PlannedExecution(repeat_id, "model-a", sample_id, seed)
            for repeat_id, seed in ((1, 17), (2, 23))
            for sample_id in (0, 1)
        ),
    )


def _set_run_id(recipe, run_id: str) -> None:
    recipe.run_id = run_id
    if getattr(recipe, "output_identity", None) is not None:
        recipe.output_identity = replace(recipe.output_identity, run_id=run_id)


def _normalize(value, root: Path):
    if isinstance(value, dict):
        return {
            key: _normalize(item, root)
            for key, item in value.items()
            if key != "created_at"
        }
    if isinstance(value, list):
        return [_normalize(item, root) for item in value]
    if isinstance(value, str):
        return value.replace(str(root), "<WORKDIR>")
    return value


def _semantic_tree(root: Path, *, normalization_root: Path) -> dict[str, object]:
    files: dict[str, object] = {}
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        relative = path.relative_to(root).as_posix()
        if path.suffix not in TEXT_SUFFIXES:
            files[relative] = {"binary": True}
            continue
        if path.suffix == ".csv":
            with path.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.reader(handle))
            timestamp_index = (
                rows[0].index("created_at")
                if rows and "created_at" in rows[0]
                else None
            )
            files[relative] = [
                [
                    _normalize(cell, normalization_root)
                    for index, cell in enumerate(row)
                    if index != timestamp_index
                ]
                for row in rows
            ]
        elif path.suffix == ".json":
            files[relative] = _normalize(
                json.loads(path.read_text(encoding="utf-8")), normalization_root
            )
        elif path.suffix == ".jsonl":
            files[relative] = [
                _normalize(json.loads(line), normalization_root)
                for line in path.read_text(encoding="utf-8").splitlines()
            ]
        elif path.suffix == ".yaml":
            files[relative] = _normalize(
                yaml.safe_load(path.read_text(encoding="utf-8")), normalization_root
            )
        else:
            files[relative] = _normalize(
                path.read_text(encoding="utf-8"), normalization_root
            )
    return files


def _digest(tree: dict[str, object]) -> str:
    payload = json.dumps(tree, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _record_highlights(results_dir: Path) -> dict[str, object]:
    manifest = json.loads((results_dir / "manifest.json").read_text(encoding="utf-8"))
    records = [
        json.loads(line)
        for line in (results_dir / "records.jsonl").read_text(encoding="utf-8").splitlines()
    ]
    identities = sorted(
        {
            (
                item.get("compound_id"),
                item.get("effective_seed"),
                item.get("repeat_id"),
                item.get("model_id"),
                item.get("sample_id"),
                item.get("entity_id"),
                item.get("chain_id"),
            )
            for item in records
        },
        key=lambda value: tuple("" if item is None else str(item) for item in value),
    )
    metric_states = Counter(
        (item["metric_name"], item["status"])
        for item in records
        if item["record_kind"] == "metric"
    )
    computed_values = sorted(
        {
            (item["metric_name"], json.dumps(item.get("value"), sort_keys=True))
            for item in records
            if item["record_kind"] == "metric"
            and item["status"] == "computed"
            and item.get("value") is not None
        }
    )
    return {
        "status": manifest["status"],
        "record_counts": manifest["record_counts"],
        "requested_metrics": manifest["requested_metrics"],
        "artifacts": [item["relative_path"] for item in manifest["artifacts"]],
        "identity_axes": identities,
        "metric_state_counts": [
            [name, state, count]
            for (name, state), count in sorted(metric_states.items())
        ],
        "computed_values": [[name, json.loads(value)] for name, value in computed_values],
        "failures": [
            [item["stage"], item["error_code"]]
            for item in records
            if item["record_kind"] == "failure"
        ],
    }


def _validate_snapshot(root: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, object]:
    system_path, options_path = _write_inputs(root)
    _patch_validate_pipeline(monkeypatch, system_name="system")
    monkeypatch.setattr(
        "cofolder.recipes.validate.scaffold_reproduction_metrics",
        lambda **kwargs: (_ for _ in ()).throw(RuntimeError("fixture analytics failure")),
    )

    system_rows = []
    chain_rows = []
    for repeat_id, seed in ((1, 17), (2, 23)):
        for sample_id in (0, 1):
            system_rows.append(
                {
                    "model_name": "model-a",
                    "repeat": repeat_id,
                    "diffusion_sample": sample_id,
                    "confidence_score": 0.9 - sample_id * 0.1,
                    "ptm": 0.8,
                    "iptm": 0.7,
                }
            )
            for chain_id, entity_type, value in (
                ("A", "protein", 0.85),
                ("B", "ligand", 0.65),
            ):
                chain_rows.append(
                    {
                        "CHAIN_ID": chain_id,
                        "ENTITY_ID": f"entity:{0 if chain_id == 'A' else 1}",
                        "ENTITY_TYPE": entity_type,
                        "model_name": "model-a",
                        "repeat": repeat_id,
                        "diffusion_sample": sample_id,
                        "chains_ptm": value,
                        **(
                            {"affinity_pred_value": 6.1 + sample_id}
                            if chain_id == "B"
                            else {}
                        ),
                    }
                )
    monkeypatch.setattr(
        "cofolder.recipes.validate.gather.merge_runner_results",
        lambda **kwargs: (pd.DataFrame(system_rows), pd.DataFrame(chain_rows), []),
    )

    def gather_structures(**kwargs):
        structures = root / "results" / "structures"
        structures.mkdir(parents=True, exist_ok=True)
        (structures / "model-a.cif").write_text("data_model_a\n", encoding="utf-8")

    monkeypatch.setattr(
        "cofolder.recipes.validate.gather.gather_structures", gather_structures
    )
    validator = Validate(
        wrk_dir=str(root),
        system_path=str(system_path),
        options_path=str(options_path),
        repeats=2,
        seed=17,
        scoring_functions=["confidence_metrics", "affinity_metrics"],
        assess_robustness=False,
        execution_plan=_fixed_plan(),
    )
    _set_run_id(validator, "parity-validate")
    validator.run()
    tree = _semantic_tree(root / "results", normalization_root=root)
    highlights = _record_highlights(root / "results")
    return {
        "digest": _digest(tree),
        "files": list(tree),
        "status": highlights["status"],
        "record_counts": highlights["record_counts"],
        "failures": highlights["failures"],
        "identity_axes": highlights["identity_axes"],
        "metric_state_counts": highlights["metric_state_counts"],
        "computed_values": highlights["computed_values"],
        "seed_plan": asdict(validator.seed_plan),
    }


def _screen_snapshot(root: Path) -> dict[str, object]:
    system_path, options_path = _write_inputs(root)
    library = root / "library.csv"
    pd.DataFrame(
        {
            "compound_id": ["hit", "failed", "miss"],
            "smiles": ["CC", "CCC", "CCO"],
            "series": ["one", "one", "two"],
        }
    ).to_csv(library, index=False)
    vectors = {"compound_000001": "[1, 1, 0]", "compound_000003": "[0, 0, 1]"}

    def validate_run(validator):
        if validator.wrk_dir.name == "compound_000002":
            raise RuntimeError("fixture backend failure")
        results = validator.wrk_dir / "results"
        child_seed = RepeatSeedProvenance(
            1, 17, 17, 17, 17, SeedOrigin.USER_SPECIFIED
        )
        child_backend = RunnerBackendIdentity(
            "boltz2", "boltz", "2.2.1", BackendVersionStatus.DETECTED
        )
        identity = OutputIdentity(
            workflow=WorkflowKind.VALIDATE,
            run_id=f"child-{validator.wrk_dir.name}",
            system_id="screen_system",
            runner_id="boltz2",
            runner_version="2.2.1",
            backend_name="boltz",
            backend_version_status=BackendVersionStatus.DETECTED,
            effective_seed=17,
            repeat_id=1,
            model_id="boltz2",
            sample_id=0,
        )
        bundle = bundle_from_frames(
            pd.DataFrame(
                [{"model_name": "boltz2", "repeat": 1, "diffusion_sample": 0}]
            ),
            pd.DataFrame(
                [
                    {
                        "CHAIN_ID": "A",
                        "ENTITY_TYPE": "protein",
                        "model_name": "boltz2",
                        "repeat": 1,
                        "diffusion_sample": 0,
                    },
                    {
                        "CHAIN_ID": "B",
                        "ENTITY_TYPE": "ligand",
                        "model_name": "boltz2",
                        "repeat": 1,
                        "diffusion_sample": 0,
                        "ifp_distance": vectors[validator.wrk_dir.name],
                    },
                ]
            ),
            identity=identity,
            requested_metrics={"ifp_distance"},
            backend=child_backend,
            seed_plan=SeedPlan(
                17, 17, SeedOrigin.USER_SPECIFIED, (child_seed,)
            ),
        )
        write_public_bundle(bundle, results)

    with patch("cofolder.recipes.screen.Validate.run", autospec=True, side_effect=validate_run):
        screen = Screen(
            wrk_dir=str(root / "screen"),
            system_path=str(system_path),
            options_path=str(options_path),
            ligand_chain="B",
            library=str(library),
            smiles_column="smiles",
            col_id="compound_id",
            merge_data="series",
            seed=17,
            ifp_filter_threshold=0.5,
            pocket_coverage_reference="110",
            cluster_ifps=True,
        )
        _set_run_id(screen, "parity-screen")
        returned = screen.run()
    results_dir = root / "screen" / "results"
    tree = _semantic_tree(results_dir, normalization_root=root)
    highlights = _record_highlights(results_dir)
    return {
        "digest": _digest(tree),
        "files": list(tree),
        "decisions": returned[
            ["compound_id", "status", "ifp_filter_status", "ifp_cluster_id"]
        ].where(pd.notna(returned), None).values.tolist(),
        "status": highlights["status"],
        "record_counts": highlights["record_counts"],
        "artifacts": highlights["artifacts"],
        "failures": highlights["failures"],
        "identity_axes": highlights["identity_axes"],
    }


def _oracle_snapshot(root: Path) -> dict[str, object]:
    system_path, options_path = _write_inputs(root)

    def validate_run(validator):
        results = validator.wrk_dir / "results"
        results.mkdir(parents=True, exist_ok=True)
        pd.DataFrame([{"confidence_score": 0.4}]).to_csv(
            results / "system_metrics.csv", index=False
        )
        pd.DataFrame(
            [{"CHAIN_ID": "B", "ENTITY_TYPE": "ligand", "affinity_pred_value": 6.0}]
        ).to_csv(results / "chain_metrics.csv", index=False)

    with patch("cofolder.recipes.oracle.Validate.run", autospec=True, side_effect=validate_run):
        oracle = Oracle(
            wrk_dir=str(root / "oracle"),
            system_path=str(system_path),
            options_path=str(options_path),
            input_smiles="CCO",
            score_components={
                "ligand_B__affinity_pred_value": -1.0,
                "system__confidence_score": 2.0,
            },
            score_gates=[OracleGate("system__confidence_score", "ge", 0.5)],
            gate_policy=OracleGatePolicy("fixed_penalty", -3.0),
            scoring_functions=["affinity_metrics", "confidence_metrics"],
        )
        _set_run_id(oracle, "parity-oracle")
        value = oracle.run()
    results_dir = root / "oracle" / "results"
    tree = _semantic_tree(results_dir, normalization_root=root)
    highlights = _record_highlights(results_dir)
    return {
        "digest": _digest(tree),
        "files": list(tree),
        "value": value,
        "status": highlights["status"],
        "record_counts": highlights["record_counts"],
        "requested_metrics": highlights["requested_metrics"],
        "computed_values": highlights["computed_values"],
    }


def _bias_snapshot(root: Path) -> dict[str, object]:
    system_path = root / "system.yaml"
    system_path.write_text(
        yaml.safe_dump(
            {
                "sequences": [
                    {"protein": {"id": "A", "sequence": "MAAA"}},
                    {"protein": {"id": "C", "sequence": "MBBB"}},
                    {"ligand": {"id": "B", "smiles": "CCO"}},
                    {"ligand": {"id": "D", "smiles": "CCN"}},
                ]
            }
        ),
        encoding="utf-8",
    )
    structure = root / "custom.pdb"
    structure.write_text("HEADER PARITY\n", encoding="utf-8")
    protein = root / "protein.csv"
    ligand = root / "ligand.csv"
    custom_protein = root / "custom_protein.csv"
    custom_ligand = root / "custom_ligand.csv"
    protein.write_text(
        "pdb_id,release_date,sequence,sequence_similarity\n"
        "1AAA,2022-01-01,MAAA,100.0\n"
        "2CCC,2022-01-01,MBBB,100.0\n"
        "3ONE,2022-01-01,MAAA,90.0\n",
        encoding="utf-8",
    )
    ligand.write_text(
        "pdb_id,release_date,ligand_id,smiles,ecfp_similarity\n"
        "1AAA,2022-01-01,LIGA,CCO,1.0\n"
        "2CCC,2022-01-01,LIGD,CCN,1.0\n"
        "3ONE,2022-01-01,LOW,CCCC,0.1\n",
        encoding="utf-8",
    )
    custom_protein.write_text(
        "sequence,dataset_name,source_structure_path\n"
        f"MAAA,custom_set,{structure.name}\n",
        encoding="utf-8",
    )
    custom_ligand.write_text(
        "smiles,dataset_name,source_reference_path\n"
        f"CCO,custom_set,{structure.name}\n",
        encoding="utf-8",
    )
    from cofolder.modules.analytics import bias as bias_module

    enrichment_functions = (
        bias_module._enrich_mixed_bias_dataset_with_pdb_backfill,
        bias_module._enrich_same_type_ligand_pair_dataset_with_pdb_backfill,
    )
    owners = {
        importlib.import_module(function.__module__) for function in enrichment_functions
    }
    checkpoint_calls = []
    with ExitStack() as stack:
        for owner in owners:
            original = owner._write_progress_checkpoint

            def capture_checkpoint(*args, _original=original, **kwargs):
                checkpoint_calls.append(
                    [
                        kwargs["progress_label"],
                        kwargs["completed_lookups"],
                        kwargs["total_lookups"],
                        Path(kwargs["output_path"]).name,
                    ]
                )
                return _original(*args, **kwargs)

            stack.enter_context(
                patch.object(owner, "_write_progress_checkpoint", capture_checkpoint)
            )
        bias = Bias(
            wrk_dir=str(root / "bias"),
            system_path=str(system_path),
            protein_training_data_path=str(protein),
            ligand_training_data_path=str(ligand),
            custom_protein_reference_path=str(custom_protein),
            custom_ligand_reference_path=str(custom_ligand),
        )
        _set_run_id(bias, "parity-bias")
        system_df, chain_df = bias.run()
    results_dir = root / "bias" / "results"
    tree = _semantic_tree(results_dir, normalization_root=root)
    highlights = _record_highlights(results_dir)
    return {
        "digest": _digest(tree),
        "files": list(tree),
        "system_columns": system_df.columns.tolist(),
        "chain_ids": chain_df["CHAIN_ID"].tolist(),
        "status": highlights["status"],
        "record_counts": highlights["record_counts"],
        "artifacts": highlights["artifacts"],
        "computed_values": highlights["computed_values"],
        "checkpoint_calls": checkpoint_calls,
    }


def collect_structural_parity(root: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, object]:
    snapshots = {}
    for name, collector in (
        ("validate", lambda path: _validate_snapshot(path, monkeypatch)),
        ("screen", _screen_snapshot),
        ("oracle", _oracle_snapshot),
        ("bias", _bias_snapshot),
    ):
        path = root / name
        path.mkdir(parents=True)
        snapshots[name] = collector(path)
    return json.loads(json.dumps(snapshots, allow_nan=False))


def test_s6_workflow_semantics_match_pre_refactor_fixtures(temp_dir, monkeypatch):
    expected = json.loads(EXPECTATIONS_PATH.read_text(encoding="utf-8"))
    assert collect_structural_parity(temp_dir, monkeypatch) == expected


def test_extracted_modules_import_without_runtime_side_effects(temp_dir):
    modules = (
        "cofolder.recipes._execution",
        "cofolder.recipes._diagnostics",
        "cofolder.recipes._results",
        "cofolder.recipes._screen_postprocess",
        "cofolder.modules.analytics.aggregation",
        "cofolder.modules.analytics.bias",
        "cofolder.modules.analytics.bias._artifacts",
        "cofolder.modules.analytics.bias._datasets",
        "cofolder.modules.analytics.bias._enrichment",
        "cofolder.modules.analytics.bias._references",
        "cofolder.modules.analytics.bias._similarity",
        "cofolder.modules.analytics.bias._workflow",
        "cofolder.modules.runners._ligand_preparation",
    )
    code = f"""
import sys
def deny(event, args):
    if event in {{'subprocess.Popen', 'socket.connect'}}:
        raise RuntimeError(event)
sys.addaudithook(deny)
for name in {modules!r}:
    __import__(name)
"""
    before = tuple(temp_dir.iterdir())
    completed = subprocess.run(
        [sys.executable, "-c", code],
        cwd=temp_dir,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert tuple(temp_dir.iterdir()) == before


def test_extracted_foundations_do_not_import_workflows():
    source_root = Path(__file__).parents[2] / "src" / "cofolder"
    paths = [
        source_root / "modules" / "analytics" / "aggregation.py",
        source_root / "modules" / "runners" / "_ligand_preparation.py",
        *(source_root / "modules" / "analytics" / "bias").glob("*.py"),
    ]
    offenders = []
    for path in paths:
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            module = node.module if isinstance(node, ast.ImportFrom) else None
            names = [item.name for item in node.names] if isinstance(node, ast.Import) else []
            if (module and module.startswith("cofolder.recipes")) or any(
                name.startswith("cofolder.recipes") for name in names
            ):
                offenders.append(f"{path.relative_to(source_root)}:{node.lineno}")
    assert offenders == []
