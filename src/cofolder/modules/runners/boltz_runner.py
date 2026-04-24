from __future__ import annotations

import copy
import json
import logging
import re
import shutil
import subprocess
from pathlib import Path
from time import perf_counter
from typing import Any, Sequence

import pandas as pd

from cofolder.modules.analytics import stats
from cofolder.modules.input.command import Command
from cofolder.modules.runners.base import BaseRunner, RunnerPreparation, RunnerRequest, RunnerResult
from cofolder.modules.utils import read
from cofolder.modules.utils.timing import DebugTimingCollector

logger = logging.getLogger(__name__)

_ANSI_ESCAPE_RE = re.compile(r"\x1b\[[0-9;?]*[ -/]*[@-~]")


def _normalize_boltz_output_line(line: str) -> str:
    cleaned = _ANSI_ESCAPE_RE.sub("", str(line or ""))
    return cleaned.replace("\r", "").strip()


def _parse_boltz_stage_timings(
    timed_lines: Sequence[tuple[float, str]],
    total_elapsed: float,
) -> dict[str, float]:
    """Parse stage timings from timestamped Boltz stdout lines."""

    def finalize_msa_window() -> None:
        nonlocal active_msa_start, active_msa_complete, msa_total
        if (
            active_msa_start is not None
            and active_msa_complete is not None
            and active_msa_complete >= active_msa_start
        ):
            msa_total += active_msa_complete - active_msa_start
        active_msa_start = None
        active_msa_complete = None

    msa_total = 0.0
    active_msa_start: float | None = None
    active_msa_complete: float | None = None
    affinity_start: float | None = None

    for offset, raw_line in timed_lines:
        line = _normalize_boltz_output_line(raw_line)
        if not line:
            continue

        if line.startswith("Calling MSA server for target"):
            finalize_msa_window()
            active_msa_start = offset
            active_msa_complete = None
            continue

        if active_msa_start is not None and "COMPLETE: 100%" in line:
            active_msa_complete = offset

        if line.startswith("Predicting property: affinity") or line.startswith(
            "Running affinity prediction for"
        ):
            finalize_msa_window()

        if line.startswith("Running affinity prediction for") and affinity_start is None:
            affinity_start = offset

    finalize_msa_window()

    timings: dict[str, float] = {}
    if msa_total > 0.0:
        timings["boltz.msa"] = msa_total
    if affinity_start is not None and total_elapsed >= affinity_start:
        timings["boltz.affinity_prediction"] = total_elapsed - affinity_start
    return timings


def _record_boltz_timings(
    timings: DebugTimingCollector | None,
    label_prefix: str | None,
    total_elapsed: float,
    timed_lines: Sequence[tuple[float, str]],
) -> None:
    if timings is None:
        return

    active_logger = timings.logger or logger
    prefix = f"{label_prefix}." if label_prefix else ""
    timings.record(f"{prefix}boltz.total", total_elapsed, logger=active_logger)
    for label, elapsed in _parse_boltz_stage_timings(timed_lines, total_elapsed).items():
        timings.record(f"{prefix}{label}", elapsed, logger=active_logger)


def run_boltz(
    cmd: Sequence[str],
    check: bool = True,
    timings: DebugTimingCollector | None = None,
    label_prefix: str | None = None,
) -> subprocess.CompletedProcess:
    """Execute a Boltz command via subprocess and log stdout/stderr in real-time."""
    start_time = perf_counter()
    logger.info("Running: %s", " ".join(cmd))

    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        universal_newlines=True,
    )

    output_lines: list[str] = []
    timed_lines: list[tuple[float, str]] = []
    assert process.stdout is not None
    for line in process.stdout:
        line = line.rstrip()
        output_lines.append(line)
        timed_lines.append((perf_counter() - start_time, line))
        logger.info(line)

    process.stdout.close()
    retcode = process.wait()
    total_elapsed = perf_counter() - start_time

    logger.info("Boltz execution completed in %.2f seconds", total_elapsed)
    _record_boltz_timings(
        timings=timings,
        label_prefix=label_prefix,
        total_elapsed=total_elapsed,
        timed_lines=timed_lines,
    )

    if check and retcode != 0:
        raise subprocess.CalledProcessError(retcode, cmd, "\n".join(output_lines))

    return subprocess.CompletedProcess(
        args=cmd,
        returncode=retcode,
        stdout="\n".join(output_lines),
    )


class BoltzRunner(BaseRunner):
    """Shared parent for Boltz-family runners."""

    name = "boltz"
    capabilities = {
        "confidence_metrics",
        "affinity_metrics",
        "affinity_metrics_ext",
    }
    model_name: str | None = "boltz2"

    def load_options(self, options_path: Path) -> Command:
        command = Command(options_path=str(options_path))
        self._set_command_option(command, "model", self.model_name)
        return command

    def prepare_system(
        self,
        system_obj: Any,
        options_obj: Command,
        wrk_dir: Path,
        conformers: str | None,
        sdf_file: Path | None,
        logger: logging.Logger,
    ) -> RunnerPreparation:
        warnings: list[str] = []
        runtime_context = {
            "cache_path": options_obj.find_value(key="cache") or "~/.boltz",
            "diffusion_samples": options_obj.find_value(key="diffusion_samples") or 1,
            "boltz_model": self.model_name,
        }

        if conformers:
            from cofolder.modules.entities import ligand

            ligand.handle_conformers(
                sys_obj=system_obj,
                opt_obj=options_obj,
                wrk_dir=wrk_dir,
                conformers=conformers,
                sdf_file=sdf_file,
                logger=logger,
            )

        return RunnerPreparation(
            system_obj=system_obj,
            options_obj=options_obj,
            warnings=warnings,
            runtime_context=runtime_context,
        )

    def run(self, request: RunnerRequest) -> RunnerResult:
        options_obj = Command(options=copy.deepcopy(request.options_obj.options))
        options_obj.seed = request.seed
        options_obj.out_dir = request.repeat_dir
        options_obj.system_path = request.system_path

        cmd = options_obj.set_command(system=request.system_obj)
        run_boltz(
            cmd,
            timings=request.timings,
            label_prefix=request.label_prefix,
        )

        raw_output_dir = (
            request.repeat_dir
            / f"boltz_results_{request.system_name}"
            / "predictions"
            / request.system_name
        )
        normalized_dir = request.repeat_dir / "normalized"
        structures_dir = normalized_dir / "structures"
        structures_dir.mkdir(parents=True, exist_ok=True)

        diffusion_samples = int(options_obj.find_value(key="diffusion_samples") or 1)
        sample_records = self._copy_structures(
            raw_output_dir=raw_output_dir,
            target_dir=structures_dir,
            repeat=request.repeat,
            system_name=request.system_name,
            diffusion_samples=diffusion_samples,
        )
        system_df, chain_df = self._normalize_metrics(
            raw_output_dir=raw_output_dir,
            request=request,
            diffusion_samples=diffusion_samples,
        )

        system_metrics_path = normalized_dir / "system_metrics.csv"
        chain_metrics_path = normalized_dir / "chain_metrics.csv"
        system_df.to_csv(system_metrics_path, index=False)
        chain_df.to_csv(chain_metrics_path, index=False)

        runtime_context = {
            "cache_path": options_obj.find_value(key="cache") or "~/.boltz",
            "diffusion_samples": diffusion_samples,
            "boltz_model": self.model_name,
        }
        manifest_path = normalized_dir / "manifest.json"
        manifest = {
            "runner": self.name,
            "capabilities": sorted(self.capabilities),
            "repeat": request.repeat,
            "raw_output_dir": str(raw_output_dir),
            "normalized_dir": str(normalized_dir),
            "system_metrics_path": str(system_metrics_path),
            "chain_metrics_path": str(chain_metrics_path),
            "structures_dir": str(structures_dir),
            "diffusion_samples": diffusion_samples,
            "runtime_context": runtime_context,
            "warnings": [],
            "sample_records": sample_records,
        }
        manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

        return RunnerResult(
            runner_name=self.name,
            raw_output_dir=raw_output_dir,
            normalized_dir=normalized_dir,
            structures_dir=structures_dir,
            system_metrics_path=system_metrics_path,
            chain_metrics_path=chain_metrics_path,
            manifest_path=manifest_path,
            diffusion_samples=diffusion_samples,
            capabilities=set(self.capabilities),
            warnings=[],
            runtime_context=runtime_context,
            sample_records=sample_records,
        )

    def _copy_structures(
        self,
        raw_output_dir: Path,
        target_dir: Path,
        repeat: int,
        system_name: str,
        diffusion_samples: int,
    ) -> list[dict[str, Any]]:
        records: list[dict[str, Any]] = []
        for sample_idx in range(diffusion_samples):
            copied_name = None
            for suffix in (".cif", ".mmcif", ".pdb"):
                file_path = raw_output_dir / f"{system_name}_model_{sample_idx}{suffix}"
                if not file_path.exists():
                    continue

                copied_name = f"{repeat}_{system_name}_model_{sample_idx}{suffix}"
                shutil.copy2(file_path, target_dir / copied_name)
                break

            if copied_name is None:
                copied_name = f"{repeat}_{system_name}_model_{sample_idx}.cif"

            records.append(
                {
                    "repeat": repeat,
                    "diffusion_sample": sample_idx,
                    "cif_file": copied_name,
                }
            )
        return records

    def _normalize_metrics(
        self,
        raw_output_dir: Path,
        request: RunnerRequest,
        diffusion_samples: int,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        system_rows: list[dict[str, Any]] = []
        chain_rows: list[dict[str, Any]] = []
        affinity_payload = self._read_affinity_payload(raw_output_dir, request.system_name)
        binder_chain_id = self._resolve_binder_chain_id(request.system_obj)

        for sample_idx in range(diffusion_samples):
            conf_path = raw_output_dir / f"confidence_{request.system_name}_model_{sample_idx}.json"
            conf_data = read.read_json(conf_path) if conf_path.exists() else {}
            system_row = {
                "cif_file": f"{request.repeat}_{request.system_name}_model_{sample_idx}.cif",
                "model_name": request.system_name,
                "repeat": request.repeat,
                "diffusion_sample": sample_idx,
            }
            chain_metrics = {}

            if conf_data:
                for key, value in conf_data.items():
                    if key in {"chains_ptm", "pair_chains_iptm"}:
                        continue
                    if not isinstance(value, (dict, list)):
                        system_row[key] = value

                chain_metrics = conf_data.get("chains_ptm", {}) or {}
                pair_chain_metrics = conf_data.get("pair_chains_iptm", {}) or {}
            else:
                pair_chain_metrics = {}

            system_rows.append(system_row)

            chain_ids = sorted(
                set(chain_metrics.keys()) | set(pair_chain_metrics.keys()),
                key=int,
            )
            for conf_chain_id in chain_ids:
                row = {
                    "conf_chain_id": int(conf_chain_id),
                    "cif_file": f"{request.repeat}_{request.system_name}_model_{sample_idx}.cif",
                    "model_name": request.system_name,
                    "repeat": request.repeat,
                    "diffusion_sample": sample_idx,
                }
                if conf_chain_id in chain_metrics:
                    row["chains_ptm"] = chain_metrics[conf_chain_id]

                for other_chain_id, value in (pair_chain_metrics.get(conf_chain_id) or {}).items():
                    row[f"pair_chains_iptm_{other_chain_id}"] = value

                if (
                    binder_chain_id is not None
                    and affinity_payload is not None
                    and int(conf_chain_id) == binder_chain_id
                ):
                    row.update(affinity_payload)

                chain_rows.append(row)

        system_df = pd.DataFrame(system_rows)
        chain_df = pd.DataFrame(chain_rows)
        return system_df, chain_df

    @staticmethod
    def _read_affinity_payload(raw_output_dir: Path, system_name: str) -> dict[str, Any] | None:
        affinity_path = raw_output_dir / f"affinity_{system_name}.json"
        if not affinity_path.exists():
            return None

        affinity_data = read.read_json(affinity_path) or {}
        affinity_pred_value = affinity_data.get("affinity_pred_value")
        affinity_probability_binary = affinity_data.get("affinity_probability_binary")
        if affinity_pred_value is None or affinity_probability_binary is None:
            return None

        pIC50, IC50_M = stats.affinity_to_pic50_and_ic50(affinity_pred_value)
        pIC50_kcal_per_mol = stats.affinity_to_pic50_kcal_per_mol(affinity_pred_value)

        return {
            "affinity_pred_value": affinity_pred_value,
            "affinity_probability_binary": affinity_probability_binary,
            "pIC50": pIC50,
            "IC50_M": IC50_M,
            "pIC50_kcal_per_mol": pIC50_kcal_per_mol,
        }

    @staticmethod
    def _set_command_option(command: Command, key: str, value: Any | None) -> None:
        options = command.options.setdefault("options", [])
        found = False
        cleaned_options = []

        for item in options:
            if not isinstance(item, dict):
                cleaned_options.append(item)
                continue

            if key in item:
                found = True
                if value is not None:
                    item[key] = value
                else:
                    item = {k: v for k, v in item.items() if k != key}

            if item:
                cleaned_options.append(item)

        if not found and value is not None:
            cleaned_options.append({key: value})

        command.options["options"] = cleaned_options

    @staticmethod
    def _resolve_binder_chain_id(system_obj: Any) -> int | None:
        properties = system_obj.find_value(key="properties") or []
        for prop in properties:
            if not isinstance(prop, dict) or "affinity" not in prop:
                continue
            binder = prop["affinity"].get("binder")
            if binder is None:
                return None
            chain_ids = BoltzRunner._ordered_chain_ids(system_obj)
            binder = str(binder)
            if binder in chain_ids:
                return chain_ids.index(binder)
            try:
                return int(binder)
            except ValueError:
                return None
        return None

    @staticmethod
    def _ordered_chain_ids(system_obj: Any) -> list[str]:
        ordered_chain_ids: list[str] = []
        sequences = system_obj.find_value(key="sequences") or []
        for seq_entry in sequences:
            if not isinstance(seq_entry, dict):
                continue
            entity_type = next(iter(seq_entry))
            entity_data = seq_entry[entity_type]
            chain_ids = entity_data.get("id")
            if chain_ids is None:
                continue
            if not isinstance(chain_ids, list):
                chain_ids = [chain_ids]
            ordered_chain_ids.extend([str(cid) for cid in chain_ids])
        return ordered_chain_ids
