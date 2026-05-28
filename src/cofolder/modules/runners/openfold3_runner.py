from __future__ import annotations

import copy
import json
import logging
import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter
from typing import Any
from typing import Mapping
from typing import Sequence

import pandas as pd
import yaml

from cofolder.modules.runners.base import BaseRunner
from cofolder.modules.runners.contracts import (
    RunnerCompanionArtifact,
    RunnerExecutionRequest,
    RunnerExecutionResult,
    RunnerMetricOutcome,
    RunnerPreparationResult,
    RunnerRuntime,
)
from cofolder.modules.utils import read
from cofolder.modules.utils.timing import DebugTimingCollector

logger = logging.getLogger(__name__)

_SAMPLE_MODEL_RE = re.compile(
    r"(?P<query>.+)_seed_(?P<seed>\d+)_sample_(?P<sample>\d+)_model(?P<suffix>\.cif|\.mmcif|\.pdb)$"
)
_OPENFOLD3_CONFIDENCE_COLUMNS = (
    "system_metrics.ptm",
    "system_metrics.iptm",
    "system_metrics.avg_plddt",
    "system_metrics.gpde",
    "system_metrics.disorder",
    "system_metrics.has_clash",
    "system_metrics.sample_ranking_score",
    "chain_metrics.chain_ptm",
    "chain_metrics.prefix:chain_pair_iptm_",
    "chain_metrics.prefix:bespoke_iptm_",
)
_OPENFOLD3_CONFIDENCE_ARTIFACTS = (
    "companion_artifacts.plddt",
    "companion_artifacts.pae",
    "companion_artifacts.pde",
)
_OPENFOLD3_FULL_CONFIDENCE_LABELS = ("plddt", "pae", "pde")
_OPENFOLD3_DEFAULT_CACHE = Path.home() / ".openfold3"


def resolve_openfold3_cache_root(env: Mapping[str, str] | None = None) -> Path:
    environment = os.environ if env is None else env
    configured = environment.get("OPENFOLD_CACHE")
    if configured:
        return Path(configured).expanduser()
    return _OPENFOLD3_DEFAULT_CACHE


def check_openfold3_setup_ready(env: Mapping[str, str] | None = None) -> tuple[bool, str | None]:
    cache_root = resolve_openfold3_cache_root(env)
    ckpt_root_path = cache_root / "ckpt_root"
    if not ckpt_root_path.exists():
        return (
            False,
            (
                "The selected 'openfold3' runner is installed, but its setup is incomplete. "
                f"Expected setup marker '{ckpt_root_path}'. Run `scripts/setup_openfold3.sh` "
                "after installing the OpenFold3 extra."
            ),
        )

    ckpt_root_value = ckpt_root_path.read_text(encoding="utf-8").strip()
    if not ckpt_root_value:
        return (
            False,
            (
                "The selected 'openfold3' runner is installed, but its setup marker is empty. "
                f"Re-run `scripts/setup_openfold3.sh` for OPENFOLD_CACHE='{cache_root}'."
            ),
        )

    checkpoint_root = Path(ckpt_root_value).expanduser()
    if not checkpoint_root.is_absolute():
        checkpoint_root = (cache_root / checkpoint_root).resolve()
    if not checkpoint_root.exists():
        return (
            False,
            (
                "The selected 'openfold3' runner is installed, but its configured checkpoint root "
                f"'{checkpoint_root}' does not exist. Re-run `scripts/setup_openfold3.sh`."
            ),
        )
    return True, None


@dataclass(slots=True)
class OpenFold3Options:
    settings: dict[str, Any]

    def find_value(self, key: str | None = None, path: Sequence[str | int] | None = None) -> Any:
        if path:
            current: Any = self.settings
            for token in path:
                if isinstance(current, dict):
                    current = current[token]
                elif isinstance(current, list):
                    current = current[int(token)]
                else:
                    raise ValueError(f"Path {list(path)!r} is invalid for OpenFold3 options.")
            return current

        if key is None:
            return None

        def search(value: Any) -> list[Any]:
            found: list[Any] = []
            if isinstance(value, dict):
                for child_key, child_value in value.items():
                    if child_key == key:
                        found.append(child_value)
                    found.extend(search(child_value))
            elif isinstance(value, list):
                for child in value:
                    found.extend(search(child))
            return found

        matches = search(self.settings)
        if not matches:
            return None
        if len(matches) > 1:
            raise ValueError(f"Key {key!r} appears multiple times in OpenFold3 options.")
        return matches[0]

    @property
    def cache_path(self) -> str | None:
        value = self.find_value(key="cache_path")
        return None if value is None else str(value)

    @property
    def diffusion_samples(self) -> int:
        for key in ("samples_per_seed", "diffusion_samples", "num_samples"):
            value = self.find_value(key=key)
            if value is not None:
                return int(value)
        return 1

    @property
    def executable(self) -> str:
        value = self.find_value(key="executable")
        return str(value) if value is not None else "run_openfold"

    @property
    def subcommand(self) -> str:
        value = self.find_value(key="subcommand")
        return str(value) if value is not None else "predict"

    @property
    def extra_args(self) -> list[str]:
        value = self.find_value(key="extra_args")
        if value is None:
            return []
        if not isinstance(value, list):
            raise ValueError("OpenFold3 option 'extra_args' must be a list when provided.")
        return [str(item) for item in value]

    def to_config_payload(self) -> dict[str, Any]:
        return copy.deepcopy(self.settings)

    def to_runner_yaml_settings(self, *, seed: int) -> dict[str, Any]:
        settings = copy.deepcopy(self.settings)
        for local_key in (
            "cache_path",
            "diffusion_samples",
            "samples_per_seed",
            "num_samples",
            "executable",
            "subcommand",
            "extra_args",
        ):
            settings.pop(local_key, None)
        experiment_settings = settings.get("experiment_settings")
        if experiment_settings is None:
            experiment_settings = {}
            settings["experiment_settings"] = experiment_settings
        elif not isinstance(experiment_settings, dict):
            raise ValueError(
                "OpenFold3 option 'experiment_settings' must be a mapping when provided."
            )
        experiment_settings["seeds"] = [int(seed)]
        return settings


def run_openfold3(
    *,
    query_json_path: Path,
    runner_yaml_path: Path,
    output_dir: Path,
    diffusion_samples: int,
    options: OpenFold3Options,
    timings: DebugTimingCollector | None = None,
    label_prefix: str | None = None,
) -> subprocess.CompletedProcess[str]:
    start_time = perf_counter()
    cmd = [
        options.executable,
        options.subcommand,
        f"--query-json={query_json_path}",
        f"--output-dir={output_dir}",
        f"--num-diffusion-samples={diffusion_samples}",
        f"--runner-yaml={runner_yaml_path}",
    ]
    cmd.extend(options.extra_args)
    logger.info("Running: %s", " ".join(str(token) for token in cmd))
    completed = subprocess.run(
        cmd,
        check=True,
        capture_output=True,
        text=True,
    )
    if timings is not None:
        label = f"{label_prefix}.openfold3.total" if label_prefix else "openfold3.total"
        timings.record(label, perf_counter() - start_time, logger=logger)
    if completed.stdout:
        logger.info("%s", completed.stdout.strip())
    if completed.stderr:
        logger.info("%s", completed.stderr.strip())
    return completed


class OpenFold3Runner(BaseRunner):
    name = "openfold3"
    capabilities = {"confidence_metrics"}

    def check_availability(self) -> tuple[bool, str | None]:
        return self.check_distribution_available(
            distribution_name="openfold3",
            missing_message=(
                "The selected 'openfold3' runner is not installed. Install it with "
                "`python -m pip install -e \".[openfold3]\"` from a COFOLDER checkout."
            ),
        )

    def load_options(self, options_path: Path) -> OpenFold3Options:
        settings = read.read_yaml(path=options_path)
        if not isinstance(settings, dict):
            raise ValueError("OpenFold3 options YAML must contain a mapping at the document root.")
        return OpenFold3Options(settings=settings)

    def prepare_system(
        self,
        system_obj: Any,
        options_obj: OpenFold3Options,
        wrk_dir: Path,
        conformers: str | None,
        sdf_file: Path | None,
        logger: logging.Logger | None,
    ) -> RunnerPreparationResult:
        return RunnerPreparationResult(
            system_obj=system_obj,
            options_obj=options_obj,
            runtime=RunnerRuntime(
                cache_path=options_obj.cache_path,
                diffusion_samples=options_obj.diffusion_samples,
                model_name=self.name,
            ),
        )

    def run(self, request: RunnerExecutionRequest) -> RunnerExecutionResult:
        options_obj = request.options_obj
        if not isinstance(options_obj, OpenFold3Options):
            raise TypeError(
                "OpenFold3Runner requires OpenFold3Options from load_options()."
            )

        raw_output_dir = request.repeat_dir / "openfold3"
        raw_output_dir.mkdir(parents=True, exist_ok=True)
        query_json_path = raw_output_dir / "inference_query_set.json"
        runner_yaml_path = raw_output_dir / "runner.yaml"
        query_json_path.write_text(
            json.dumps(self._build_query_payload(request.system_name, request.system_obj), indent=2),
            encoding="utf-8",
        )
        runner_yaml_path.write_text(
            yaml.safe_dump(
                options_obj.to_runner_yaml_settings(seed=request.seed),
                sort_keys=False,
            ),
            encoding="utf-8",
        )

        run_openfold3(
            query_json_path=query_json_path,
            runner_yaml_path=runner_yaml_path,
            output_dir=raw_output_dir,
            diffusion_samples=options_obj.diffusion_samples,
            options=options_obj,
            timings=request.timings,
            label_prefix=request.label_prefix,
        )

        normalized_dir = request.repeat_dir / "normalized"
        structures_dir = normalized_dir / "structures"
        structures_dir.mkdir(parents=True, exist_ok=True)

        seed_dir = raw_output_dir / request.system_name / f"seed_{request.seed}"
        sample_numbers = self._discover_sample_numbers(seed_dir)
        sample_records = self._copy_structures(
            seed_dir=seed_dir,
            target_dir=structures_dir,
            repeat=request.repeat,
            system_name=request.system_name,
            sample_numbers=sample_numbers,
        )
        system_df, chain_df, confidence_payloads, confidence_issues = self._normalize_metrics(
            seed_dir=seed_dir,
            request=request,
            sample_numbers=sample_numbers,
        )
        companion_artifacts = self._copy_confidence_artifacts(
            normalized_dir=normalized_dir,
            confidence_payloads=confidence_payloads,
        )

        system_metrics_path = normalized_dir / "system_metrics.csv"
        chain_metrics_path = normalized_dir / "chain_metrics.csv"
        system_df.to_csv(system_metrics_path, index=False)
        chain_df.to_csv(chain_metrics_path, index=False)

        runtime = RunnerRuntime(
            cache_path=options_obj.cache_path,
            diffusion_samples=len(sample_numbers),
            model_name=self.name,
        )
        metric_outcomes = self._build_metric_outcomes(
            system_df=system_df,
            chain_df=chain_df,
            companion_artifacts=companion_artifacts,
            confidence_issues=confidence_issues,
        )

        manifest_path = normalized_dir / "manifest.json"
        manifest_path.write_text(
            json.dumps(
                {
                    "runner": self.name,
                    "capabilities": sorted(self.capabilities),
                    "repeat": request.repeat,
                    "raw_output_dir": str(raw_output_dir),
                    "normalized_dir": str(normalized_dir),
                    "system_metrics_path": str(system_metrics_path),
                    "chain_metrics_path": str(chain_metrics_path),
                    "structures_dir": str(structures_dir),
                    "diffusion_samples": len(sample_numbers),
                    "runtime_context": {
                        "cache_path": runtime.cache_path,
                        "diffusion_samples": runtime.diffusion_samples,
                        "openfold3_model": runtime.model_name,
                        "timing_path": str(seed_dir / "timing.json"),
                    },
                    "warnings": [],
                    "sample_records": sample_records,
                    "companion_artifacts": [
                        {
                            "label": artifact.label,
                            "relative_path": artifact.relative_path,
                            "kind": artifact.kind,
                            "description": artifact.description,
                        }
                        for artifact in companion_artifacts
                    ],
                },
                indent=2,
            ),
            encoding="utf-8",
        )

        return RunnerExecutionResult(
            runner_name=self.name,
            raw_output_dir=raw_output_dir,
            normalized_dir=normalized_dir,
            structures_dir=structures_dir,
            system_metrics_path=system_metrics_path,
            chain_metrics_path=chain_metrics_path,
            manifest_path=manifest_path,
            diffusion_samples=len(sample_numbers),
            capabilities=set(self.capabilities),
            warnings=[],
            runtime=runtime,
            sample_records=sample_records,
            metric_outcomes=metric_outcomes,
            companion_artifacts=companion_artifacts,
        )

    def _build_query_payload(self, system_name: str, system_obj: Any) -> dict[str, Any]:
        sequences = system_obj.find_value(key="sequences") or []
        chains: list[dict[str, Any]] = []
        for index, entry in enumerate(sequences):
            if not isinstance(entry, dict):
                raise ValueError(f"Sequence entry {index} must be a mapping, got {type(entry)!r}.")
            if "protein" in entry:
                chains.append(self._build_protein_chain(entry["protein"], index=index))
                continue
            if "ligand" in entry:
                chains.append(self._build_ligand_chain(entry["ligand"], index=index))
                continue
            raise ValueError(
                f"Sequence entry {index} must contain either 'protein' or 'ligand' for OpenFold3."
            )
        return {"queries": {system_name: {"chains": chains}}}

    @staticmethod
    def _build_protein_chain(protein: Any, *, index: int) -> dict[str, Any]:
        if not isinstance(protein, dict):
            raise ValueError(f"Protein entry {index} must be a mapping.")
        chain_ids = OpenFold3Runner._normalize_chain_ids(protein.get("id"), entry_label=f"protein[{index}]")
        sequence = protein.get("sequence") or protein.get("fasta")
        if not sequence:
            raise ValueError(
                f"Protein entry {index} must define either 'sequence' or 'fasta' for OpenFold3."
            )
        return {
            "molecule_type": "protein",
            "chain_ids": chain_ids,
            "sequence": str(sequence),
        }

    @staticmethod
    def _build_ligand_chain(ligand: Any, *, index: int) -> dict[str, Any]:
        if not isinstance(ligand, dict):
            raise ValueError(f"Ligand entry {index} must be a mapping.")
        chain_ids = OpenFold3Runner._normalize_chain_ids(ligand.get("id"), entry_label=f"ligand[{index}]")
        smiles = ligand.get("smiles")
        ccd_codes = ligand.get("ccd_codes")
        if ccd_codes is None and ligand.get("ccd") is not None:
            ccd_codes = [ligand.get("ccd")]
        elif isinstance(ccd_codes, str):
            ccd_codes = [ccd_codes]

        if ccd_codes:
            return {
                "molecule_type": "ligand",
                "chain_ids": chain_ids,
                "ccd_codes": [str(code) for code in ccd_codes],
            }
        if smiles:
            return {
                "molecule_type": "ligand",
                "chain_ids": chain_ids,
                "smiles": str(smiles),
            }
        raise ValueError(
            f"Ligand entry {index} must define either 'smiles' or 'ccd'/'ccd_codes' for OpenFold3."
        )

    @staticmethod
    def _normalize_chain_ids(value: Any, *, entry_label: str) -> list[str]:
        if value is None:
            raise ValueError(f"{entry_label} is missing required chain id(s).")
        if isinstance(value, list):
            chain_ids = [str(item) for item in value if str(item).strip()]
        else:
            chain_ids = [str(value)]
        if not chain_ids:
            raise ValueError(f"{entry_label} must contain at least one non-empty chain id.")
        return chain_ids

    @staticmethod
    def _discover_sample_numbers(seed_dir: Path) -> list[int]:
        sample_numbers: set[int] = set()
        for path in seed_dir.glob("*_model.*"):
            match = _SAMPLE_MODEL_RE.match(path.name)
            if match is None:
                continue
            sample_numbers.add(int(match.group("sample")))
        if not sample_numbers:
            raise ValueError(f"OpenFold3 produced no sample structures under {seed_dir!s}.")
        return sorted(sample_numbers)

    @staticmethod
    def _copy_structures(
        *,
        seed_dir: Path,
        target_dir: Path,
        repeat: int,
        system_name: str,
        sample_numbers: Sequence[int],
    ) -> list[dict[str, Any]]:
        records: list[dict[str, Any]] = []
        for sample_number in sample_numbers:
            source = OpenFold3Runner._find_sample_file(
                seed_dir,
                sample_number=sample_number,
                suffix="_model",
                allowed_suffixes=(".cif", ".mmcif", ".pdb"),
            )
            if source is None:
                raise ValueError(
                    f"OpenFold3 sample {sample_number} is missing a structure file under {seed_dir!s}."
                )
            target_name = f"{repeat}_{system_name}_model_{sample_number - 1}{source.suffix}"
            shutil.copy2(source, target_dir / target_name)
            records.append(
                {
                    "repeat": repeat,
                    "diffusion_sample": sample_number - 1,
                    "cif_file": target_name,
                }
            )
        return records

    def _normalize_metrics(
        self,
        *,
        seed_dir: Path,
        request: RunnerExecutionRequest,
        sample_numbers: Sequence[int],
    ) -> tuple[
        pd.DataFrame,
        pd.DataFrame,
        dict[str, list[tuple[str, Any]]],
        list[str],
    ]:
        system_rows: list[dict[str, Any]] = []
        chain_rows: list[dict[str, Any]] = []
        confidence_payloads: dict[str, list[tuple[str, Any]]] = {
            label: [] for label in _OPENFOLD3_FULL_CONFIDENCE_LABELS
        }
        confidence_issues: list[str] = []
        chain_order = self._chain_order(request.system_obj)

        for sample_number in sample_numbers:
            structure_path = self._find_sample_file(
                seed_dir,
                sample_number=sample_number,
                suffix="_model",
                allowed_suffixes=(".cif", ".mmcif", ".pdb"),
            )
            if structure_path is None:
                continue

            aggregated_path = seed_dir / (
                f"{request.system_name}_seed_{request.seed}_sample_{sample_number}_confidences_aggregated.json"
            )
            aggregated = read.read_json(aggregated_path) if aggregated_path.exists() else {}
            if not aggregated_path.exists():
                confidence_issues.append(
                    f"sample {sample_number} is missing confidences_aggregated.json"
                )
            system_row = {
                "cif_file": f"{request.repeat}_{request.system_name}_model_{sample_number - 1}{structure_path.suffix}",
                "model_name": request.system_name,
                "repeat": request.repeat,
                "diffusion_sample": sample_number - 1,
            }
            for key, value in aggregated.items():
                if key in {"chain_ptm", "chain_pair_iptm", "bespoke_iptm"}:
                    continue
                if not isinstance(value, (dict, list)):
                    system_row[key] = value
            system_rows.append(system_row)

            missing_scalar_fields = [
                key
                for key in (
                    "ptm",
                    "iptm",
                    "avg_plddt",
                    "gpde",
                    "disorder",
                    "has_clash",
                    "sample_ranking_score",
                )
                if key not in aggregated
            ]
            if missing_scalar_fields:
                confidence_issues.append(
                    f"sample {sample_number} is missing aggregated confidence fields "
                    f"{missing_scalar_fields!r}"
                )

            chain_ptm = aggregated.get("chain_ptm") or {}
            chain_pair_iptm = aggregated.get("chain_pair_iptm") or {}
            bespoke_iptm = aggregated.get("bespoke_iptm") or {}
            if not chain_ptm:
                confidence_issues.append(
                    f"sample {sample_number} is missing chain_ptm values"
                )
            if not any((chain_pair_iptm.get(chain_id) or {}) for chain_id in chain_order):
                confidence_issues.append(
                    f"sample {sample_number} is missing chain_pair_iptm values"
                )
            if not any((bespoke_iptm.get(chain_id) or {}) for chain_id in chain_order):
                confidence_issues.append(
                    f"sample {sample_number} is missing bespoke_iptm values"
                )
            for conf_chain_id, chain_id in enumerate(chain_order):
                row = {
                    "conf_chain_id": conf_chain_id,
                    "cif_file": system_row["cif_file"],
                    "model_name": request.system_name,
                    "repeat": request.repeat,
                    "diffusion_sample": sample_number - 1,
                }
                if chain_id in chain_ptm:
                    row["chain_ptm"] = chain_ptm[chain_id]
                for other_chain_id, value in (chain_pair_iptm.get(chain_id) or {}).items():
                    row[f"chain_pair_iptm_{chain_id}_{other_chain_id}"] = value
                for other_chain_id, value in (bespoke_iptm.get(chain_id) or {}).items():
                    row[f"bespoke_iptm_{chain_id}_{other_chain_id}"] = value
                chain_rows.append(row)

            full_confidence_path = seed_dir / (
                f"{request.system_name}_seed_{request.seed}_sample_{sample_number}_confidences.json"
            )
            if not full_confidence_path.exists():
                confidence_issues.append(
                    f"sample {sample_number} is missing confidences.json"
                )
                continue

            full_confidence = read.read_json(full_confidence_path) or {}
            artifact_base_name = (
                f"{request.system_name}_seed_{request.seed}_sample_{sample_number}"
            )
            for label in _OPENFOLD3_FULL_CONFIDENCE_LABELS:
                if label not in full_confidence:
                    confidence_issues.append(
                        f"sample {sample_number} is missing full-confidence field {label!r}"
                    )
                    continue
                confidence_payloads[label].append((artifact_base_name, full_confidence[label]))

        return (
            pd.DataFrame(system_rows),
            pd.DataFrame(chain_rows),
            confidence_payloads,
            confidence_issues,
        )

    @staticmethod
    def _copy_confidence_artifacts(
        *,
        normalized_dir: Path,
        confidence_payloads: dict[str, list[tuple[str, Any]]],
    ) -> list[RunnerCompanionArtifact]:
        artifacts: list[RunnerCompanionArtifact] = []
        artifacts_dir = normalized_dir / "artifacts"
        for label in _OPENFOLD3_FULL_CONFIDENCE_LABELS:
            payloads = confidence_payloads.get(label) or []
            if not payloads:
                continue
            label_dir = artifacts_dir / label
            label_dir.mkdir(parents=True, exist_ok=True)
            for artifact_base_name, payload in payloads:
                artifact_path = label_dir / f"{artifact_base_name}_{label}.json"
                artifact_path.write_text(
                    json.dumps(payload),
                    encoding="utf-8",
                )
            artifacts.append(
                RunnerCompanionArtifact(
                    label=label,
                    relative_path=f"artifacts/{label}",
                    kind="directory",
                    description=f"OpenFold3 normalized {label} outputs for all diffusion samples.",
                )
            )
        return artifacts

    @staticmethod
    def _build_metric_outcomes(
        *,
        system_df: pd.DataFrame,
        chain_df: pd.DataFrame,
        companion_artifacts: Sequence[RunnerCompanionArtifact],
        confidence_issues: Sequence[str],
    ) -> dict[str, RunnerMetricOutcome]:
        if confidence_issues:
            confidence_outcome = RunnerMetricOutcome(
                state="failed",
                message=(
                    "OpenFold3 confidence outputs were incomplete: "
                    + "; ".join(dict.fromkeys(str(issue) for issue in confidence_issues))
                    + "."
                ),
            )
        elif not system_df.empty and not chain_df.empty and companion_artifacts:
            confidence_outcome = RunnerMetricOutcome(
                state="computed",
                required_columns=_OPENFOLD3_CONFIDENCE_COLUMNS,
                required_artifacts=_OPENFOLD3_CONFIDENCE_ARTIFACTS,
            )
        else:
            confidence_outcome = RunnerMetricOutcome(
                state="missing",
                required_columns=_OPENFOLD3_CONFIDENCE_COLUMNS,
                required_artifacts=_OPENFOLD3_CONFIDENCE_ARTIFACTS,
                message="OpenFold3 confidence outputs were not found in the seed/sample tree.",
            )

        return {
            "confidence_metrics": confidence_outcome,
            "affinity_metrics": RunnerMetricOutcome(state="unsupported"),
            "affinity_metrics_ext": RunnerMetricOutcome(state="unsupported"),
        }

    @staticmethod
    def _find_sample_file(
        seed_dir: Path,
        *,
        sample_number: int,
        suffix: str,
        allowed_suffixes: Sequence[str],
    ) -> Path | None:
        for extension in allowed_suffixes:
            matches = sorted(seed_dir.glob(f"*_sample_{sample_number}{suffix}{extension}"))
            if matches:
                return matches[0]
        return None

    @staticmethod
    def _chain_order(system_obj: Any) -> list[str]:
        chain_ids: list[str] = []
        for entry in system_obj.find_value(key="sequences") or []:
            if not isinstance(entry, dict):
                continue
            payload = entry.get("protein") or entry.get("ligand")
            if not isinstance(payload, dict):
                continue
            ids = payload.get("id")
            if isinstance(ids, list):
                chain_ids.extend(str(value) for value in ids)
            elif ids is not None:
                chain_ids.append(str(ids))
        return chain_ids


RUNNER = OpenFold3Runner()
