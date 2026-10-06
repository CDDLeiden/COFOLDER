"""Run a reproducible public/custom bias cutoff comparison matrix."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
import json
from pathlib import Path

import yaml


def _resolve(base: Path, value: str | None) -> str | None:
    if not value:
        return None
    path = Path(value).expanduser()
    return str(path.resolve() if path.is_absolute() else (base / path).resolve())


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="cofolder-tools run-bias-matrix")
    parser.add_argument("config", type=Path)
    parser.add_argument("--preflight_only", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    config_path = args.config.resolve()
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or config.get("schema_version") != 1:
        raise ValueError("Bias matrix config must use schema_version: 1.")
    workflow = str(config.get("workflow", "bias")).lower()
    if workflow not in {"bias", "validate"}:
        raise ValueError("Bias matrix workflow must be 'bias' or 'validate'.")
    cutoff = str(config.get("release_cutoff", "2023-06-01"))
    if cutoff != "2023-06-01":
        raise ValueError(
            "The S5.1 comparison matrix release_cutoff must be 2023-06-01."
        )
    base = config_path.parent
    assessment_value = _resolve(base, config.get("assessment_root"))
    if not assessment_value:
        raise ValueError("Bias matrix requires assessment_root.")
    assessment_root = Path(assessment_value)
    custom = _resolve(base, config.get("custom_bias_reference_path"))
    if not custom:
        raise ValueError("Bias matrix requires custom_bias_reference_path for plus-custom scenarios.")
    custom_source: dict[str, str] = {"path": custom}
    try:
        from cofolder.modules.analytics.bias_database import (
            validate_custom_bias_reference_bundle,
        )

        custom_bundle = validate_custom_bias_reference_bundle(custom)
        custom_source["fingerprint"] = custom_bundle.fingerprint
        custom_source["dataset_name"] = str(
            custom_bundle.manifest.get("dataset_name", "")
        )
    except ValueError:
        # The public recipe/preflight reports the actionable validation failure.
        pass
    common = {
        "system_path": _resolve(base, config.get("system_path")),
        "bias_training_data_protein_path": _resolve(base, config.get("bias_training_data_protein_path")),
        "bias_training_data_ligand_path": _resolve(base, config.get("bias_training_data_ligand_path")),
        "bias_query_cache_path": _resolve(base, config.get("bias_query_cache_path")),
        "bias_chains": config.get("bias_chains"),
        "bias_protein_similarity_threshold": float(
            config.get("bias_protein_similarity_threshold", 0.25)
        ),
        "bias_ligand_similarity_threshold": float(
            config.get("bias_ligand_similarity_threshold", 0.35)
        ),
    }
    if not common["system_path"]:
        raise ValueError("Bias matrix requires system_path.")
    if not common["bias_training_data_protein_path"] or not common["bias_training_data_ligand_path"]:
        raise ValueError("Bias matrix requires both public bias source bundles.")
    scenarios = (
        {
            "name": "pre-2023-06-01",
            "release_cutoff": cutoff,
            "custom": None,
            "public": True,
        },
        {
            "name": "whole-snapshot",
            "release_cutoff": "whole",
            "custom": None,
            "public": True,
        },
        {
            "name": "pre-2023-06-01-plus-custom",
            "release_cutoff": cutoff,
            "custom": custom,
            "public": True,
        },
        {
            "name": "custom-only",
            "release_cutoff": "whole",
            "custom": custom,
            "public": False,
        },
    )
    records: list[dict[str, object]] = []
    failed = False
    for scenario in scenarios:
        name = str(scenario["name"])
        supplement = scenario["custom"]
        scenario_cutoff = str(scenario["release_cutoff"])
        destination = assessment_root / name
        kwargs = {
            **common,
            "wrk_dir": str(destination),
            "bias_release_cutoff": scenario_cutoff,
            "custom_bias_reference_path": supplement,
        }
        if not scenario["public"]:
            kwargs["bias_training_data_protein_path"] = None
            kwargs["bias_training_data_ligand_path"] = None
        try:
            if workflow == "bias":
                from cofolder.recipes.bias import Bias
                recipe = Bias(**kwargs)
            else:
                from cofolder.recipes.validate import Validate
                options = _resolve(base, config.get("options_path"))
                if not options:
                    raise ValueError("Validate matrix requires options_path.")
                recipe = Validate(
                    options_path=options,
                    runner=str(config.get("runner", "boltz2")),
                    assess_bias=True,
                    scoring_functions=config.get("scoring_functions", []),
                    **kwargs,
                )
            if args.preflight_only:
                report = recipe.preflight()
                status = "ready" if report.ready else "not-ready"
                if not report.ready:
                    failed = True
            else:
                recipe.run()
                status = "success"
            record = {
                "name": name,
                "release_cutoff": scenario_cutoff,
                "custom": bool(supplement),
                "public": bool(scenario["public"]),
                "custom_source": custom_source if supplement else None,
                "status": status,
                "output_dir": str(destination / "results"),
            }
            reference_manifest = destination / "results/bias_train/reference_manifest.json"
            if reference_manifest.is_file():
                provenance = json.loads(reference_manifest.read_text(encoding="utf-8"))
                record["reference_manifest"] = str(reference_manifest)
                record["sources"] = provenance.get("request", {}).get("sources")
                record["queries"] = provenance.get("request", {}).get("queries")
                record["release_policy"] = provenance.get("request", {}).get("release_policy")
            records.append(record)
        except Exception as exc:
            failed = True
            records.append(
                {
                    "name": name,
                    "release_cutoff": scenario_cutoff,
                    "custom": bool(supplement),
                    "public": bool(scenario["public"]),
                    "custom_source": custom_source if supplement else None,
                    "status": "failed",
                    "error": str(exc),
                    "output_dir": str(destination / "results"),
                }
            )
    if not args.preflight_only:
        assessment_root.mkdir(parents=True, exist_ok=True)
        temporary_manifest = assessment_root / ".matrix_manifest.json.tmp"
        temporary_manifest.write_text(
            json.dumps({"schema_version": 1, "workflow": workflow, "scenarios": records}, indent=2),
            encoding="utf-8",
        )
        temporary_manifest.replace(assessment_root / "matrix_manifest.json")
    else:
        print(json.dumps({"workflow": workflow, "scenarios": records}, indent=2))
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
