from __future__ import annotations

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

    subprocess.run(cmd, check=True)
