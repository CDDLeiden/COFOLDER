from __future__ import annotations

import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path


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
) -> None:
    """Run scripts/build_bias_training_data.py from package code.

    This keeps validate() in control of when and how bias training references
    are built, while reusing the existing builder implementation.
    """
    repo_root = Path(__file__).resolve().parents[4]
    script_path = repo_root / "scripts" / "build_bias_training_data.py"
    if not script_path.exists():
        raise FileNotFoundError(f"build script not found: {script_path}")

    resolved_mmseqs = _resolve_mmseqs_bin()
    if skip_protein_mmseqs is None:
        skip_protein_mmseqs = resolved_mmseqs is None

    threshold = float(ligand_similarity_threshold)
    no_hits_msg = "No CCD hits above threshold"
    context = f"system={system_path} ligand_csv={output_ligand_csv}"

    while threshold >= 0.0:
        cmd: list[str] = [
            sys.executable,
            str(script_path),
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

        proc = subprocess.run(cmd, check=False, capture_output=True, text=True)
        if proc.stdout:
            print(proc.stdout, end="")
        if proc.stderr:
            print(proc.stderr, end="", file=sys.stderr)

        if proc.returncode == 0:
            return

        combined = f"{proc.stdout}\n{proc.stderr}"
        if _looks_like_invalid_smiles_error(combined):
            print(
                f"[warn] Invalid SMILES detected during bias-training build ({context}); "
                "continuing with ligand ECFP skipped for this run.",
                file=sys.stderr,
            )
            fallback_cmd = list(cmd)
            if "--skip_ligand_ecfp" not in fallback_cmd:
                fallback_cmd.append("--skip_ligand_ecfp")
            fallback = subprocess.run(
                fallback_cmd,
                check=False,
                capture_output=True,
                text=True,
            )
            if fallback.stdout:
                print(fallback.stdout, end="")
            if fallback.stderr:
                print(fallback.stderr, end="", file=sys.stderr)
            if fallback.returncode == 0:
                return
            raise subprocess.CalledProcessError(
                fallback.returncode,
                fallback_cmd,
                output=fallback.stdout,
                stderr=fallback.stderr,
            )

        if no_hits_msg in combined:
            next_threshold = round(threshold - 0.05, 2)
            if next_threshold < 0.0:
                _write_empty_training_csvs(
                    output_protein_csv=output_protein_csv,
                    output_ligand_csv=output_ligand_csv,
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
            proc.returncode,
            cmd,
            output=proc.stdout,
            stderr=proc.stderr,
        )


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
