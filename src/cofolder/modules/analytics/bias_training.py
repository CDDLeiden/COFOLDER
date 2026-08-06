from __future__ import annotations

import csv
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from time import perf_counter

from cofolder.modules.utils.timing import DebugTimingCollector


_PROTEIN_STAGE_MARKERS = (
    "[info] using protein query sequence",
    "[info] reusing existing protein CSV:",
    "[warn] failed to read existing protein CSV, recomputing with MMseqs:",
    "[info] skipping MMseqs protein expansion",
)

def _resolve_mmseqs_bin() -> str | None:
    explicit_env = os.environ.get("COFOLDER_MMSEQS_BIN")
    candidates = []
    if explicit_env:
        candidates.append(explicit_env)
    candidates.extend(
        [
            str(Path.home() / ".cofolder/vendor/mmseqs/bin/mmseqs"),
            str(Path(__file__).resolve().parents[4] / "vendor/mmseqs/bin/mmseqs"),
            str(Path(__file__).resolve().parents[4] / "vendor/mmseqs/mmseqs"),
        ]
    )
    from_path = shutil.which("mmseqs")
    if from_path:
        candidates.append(from_path)

    for cand in candidates:
        p = Path(cand).expanduser()
        if p.exists() and p.is_file():
            return str(p)
        found = shutil.which(cand)
        if found:
            return found
    return None


def run_build_bias_training_data(
    system_path: Path,
    components_cif: Path,
    output_protein_csv: Path,
    output_ligand_csv: Path,
    release_cutoff: str,
    ligand_similarity_threshold: float = 0.35,
    overwrite: bool = True,
    skip_bias_csv: bool = True,
    skip_protein_mmseqs: bool | None = None,
    skip_ligand_ecfp: bool = False,
    ligand_chains: set[str] | None = None,
    timings: DebugTimingCollector | None = None,
) -> None:
    """Run the packaged bias-training builder in a subprocess.

    This keeps validate() in control of when and how bias training references
    are built, while reusing the existing builder implementation.
    """
    resolved_mmseqs = _resolve_mmseqs_bin()
    if skip_protein_mmseqs is None:
        skip_protein_mmseqs = resolved_mmseqs is None

    threshold = float(ligand_similarity_threshold)
    no_hits_msg = "No CCD hits above threshold"
    no_pre_cutoff_msg = "No pre-cutoff CCD-linked PDB entries"
    context = f"system={system_path} ligand_csv={output_ligand_csv}"
    ligand_elapsed_total = 0.0
    protein_elapsed_total = 0.0

    while threshold >= 0.0:
        cmd: list[str] = [
            sys.executable,
            "-u",
            "-m",
            "cofolder.modules.analytics.build_bias_training_data",
            "--system_path",
            str(system_path),
            "--components_cif",
            str(components_cif),
            "--output_protein_csv",
            str(output_protein_csv),
            "--output_ligand_csv",
            str(output_ligand_csv),
            "--release_cutoff",
            str(release_cutoff),
            "--ligand_similarity_threshold",
            str(threshold),
        ]
        if overwrite:
            cmd.append("--overwrite")
        if skip_bias_csv:
            cmd.append("--skip_bias_csv")
        if skip_protein_mmseqs:
            cmd.append("--skip_protein_mmseqs")
        elif resolved_mmseqs is not None:
            cmd.extend(["--mmseqs_bin", resolved_mmseqs])
        if skip_ligand_ecfp:
            cmd.append("--skip_ligand_ecfp")
        if ligand_chains:
            cmd.append("--ligand_chains")
            cmd.extend(sorted(ligand_chains))

        returncode, combined, timed_lines = _run_bias_training_subprocess(cmd)
        ligand_elapsed, protein_elapsed = _split_bias_build_phase_timings(
            timed_lines=timed_lines,
        )
        ligand_elapsed_total += ligand_elapsed
        protein_elapsed_total += protein_elapsed
        pre_cutoff_count = _extract_pre_cutoff_count(combined)
        if pre_cutoff_count == 0:
            next_threshold = round(threshold - 0.05, 2)
            if next_threshold < 0.0:
                _write_empty_training_csvs(
                    output_protein_csv=output_protein_csv,
                    output_ligand_csv=output_ligand_csv,
                )
                _record_bias_build_phase_timings(
                    timings=timings,
                    ligand_elapsed=ligand_elapsed_total,
                    protein_elapsed=protein_elapsed_total,
                )
                print(
                    f"[warn] {no_pre_cutoff_msg} ({context}) after threshold backoff; "
                    "wrote empty training CSVs and continuing without ligand bias hits.",
                    file=sys.stderr,
                )
                return
            print(
                f"[warn] {no_pre_cutoff_msg} ({context}) at threshold={threshold:.2f}; "
                f"retrying with {next_threshold:.2f}",
                file=sys.stderr,
            )
            threshold = next_threshold
            continue

        if returncode == 0:
            _record_bias_build_phase_timings(
                timings=timings,
                ligand_elapsed=ligand_elapsed_total,
                protein_elapsed=protein_elapsed_total,
            )
            return

        if _looks_like_invalid_smiles_error(combined):
            print(
                f"[warn] Invalid SMILES detected during bias-training build ({context}); "
                "continuing with ligand ECFP skipped for this run.",
                file=sys.stderr,
            )
            fallback_cmd = list(cmd)
            if "--skip_ligand_ecfp" not in fallback_cmd:
                fallback_cmd.append("--skip_ligand_ecfp")
            fallback_returncode, fallback_output, fallback_timed_lines = _run_bias_training_subprocess(
                fallback_cmd
            )
            fallback_ligand_elapsed, fallback_protein_elapsed = _split_bias_build_phase_timings(
                timed_lines=fallback_timed_lines,
            )
            ligand_elapsed_total += fallback_ligand_elapsed
            protein_elapsed_total += fallback_protein_elapsed
            if fallback_returncode == 0:
                _record_bias_build_phase_timings(
                    timings=timings,
                    ligand_elapsed=ligand_elapsed_total,
                    protein_elapsed=protein_elapsed_total,
                )
                return
            raise subprocess.CalledProcessError(
                fallback_returncode,
                fallback_cmd,
                output=fallback_output,
                stderr="",
            )

        if no_hits_msg in combined:
            next_threshold = round(threshold - 0.05, 2)
            if next_threshold < 0.0:
                _write_empty_training_csvs(
                    output_protein_csv=output_protein_csv,
                    output_ligand_csv=output_ligand_csv,
                )
                _record_bias_build_phase_timings(
                    timings=timings,
                    ligand_elapsed=ligand_elapsed_total,
                    protein_elapsed=protein_elapsed_total,
                )
                print(
                    f"[warn] No CCD hits found after threshold backoff ({context}); "
                    "wrote empty training CSVs and continuing without ligand bias hits.",
                    file=sys.stderr,
                )
                return
            print(
                f"[warn] {no_hits_msg} ({context}) at threshold={threshold:.2f}; "
                f"retrying with {next_threshold:.2f}",
                file=sys.stderr,
            )
            threshold = next_threshold
            continue

        raise subprocess.CalledProcessError(
            returncode,
            cmd,
            output=combined,
            stderr="",
        )


def _run_bias_training_subprocess(
    cmd: list[str],
) -> tuple[int, str, list[tuple[float, str]]]:
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    assert proc.stdout is not None

    start = perf_counter()
    timed_lines: list[tuple[float, str]] = []
    output_parts: list[str] = []

    try:
        for line in proc.stdout:
            output_parts.append(line)
            timed_lines.append((perf_counter() - start, line.rstrip("\n")))
            print(line, end="")
    finally:
        proc.stdout.close()

    returncode = proc.wait()
    combined = "".join(output_parts)
    return returncode, combined, timed_lines


def _split_bias_build_phase_timings(
    timed_lines: list[tuple[float, str]],
) -> tuple[float, float]:
    if not timed_lines:
        return 0.0, 0.0

    total_elapsed = max(timed_lines[-1][0], 0.0)
    protein_start = None

    for elapsed, line in timed_lines:
        if line.startswith(_PROTEIN_STAGE_MARKERS):
            protein_start = elapsed
            break

    if protein_start is None:
        return round(total_elapsed, 10), 0.0

    ligand_elapsed = round(max(protein_start, 0.0), 10)
    protein_elapsed = round(max(total_elapsed - protein_start, 0.0), 10)
    return ligand_elapsed, protein_elapsed


def _record_bias_build_phase_timings(
    timings: DebugTimingCollector | None,
    ligand_elapsed: float,
    protein_elapsed: float,
) -> None:
    if timings is None:
        return
    timings.record("bias.training_data.build.ligand", ligand_elapsed)
    timings.record("bias.training_data.build.protein", protein_elapsed)


def _extract_pre_cutoff_count(text: str) -> int | None:
    matches = re.findall(r"selected_pdb_entries_pre_cutoff=(\d+)", text or "")
    if not matches:
        return None
    return int(matches[-1])


def _looks_like_invalid_smiles_error(text: str) -> bool:
    lowered = (text or "").lower()
    patterns = [
        "smiles parse error",
        "failed parsing smiles",
        "failed to parse smiles",
        "invalid smiles",
    ]
    return any(p in lowered for p in patterns)


def _write_empty_training_csvs(output_protein_csv: Path, output_ligand_csv: Path) -> None:
    output_protein_csv.parent.mkdir(parents=True, exist_ok=True)
    output_ligand_csv.parent.mkdir(parents=True, exist_ok=True)

    with output_protein_csv.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "query_chain_id",
                "pdb_id",
                "release_date",
                "sequence_similarity",
                "sequence",
            ]
        )

    with output_ligand_csv.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "query_chain_id",
                "pdb_id",
                "release_date",
                "ligand_id",
                "ecfp_similarity",
                "smiles",
            ]
        )
