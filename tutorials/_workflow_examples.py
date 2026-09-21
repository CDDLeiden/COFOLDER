"""Shared, testable command builders for the public workflow tutorials."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml


@dataclass(frozen=True, slots=True)
class BiasExampleInputs:
    system_path: Path
    protein_reference_path: Path
    ligand_reference_path: Path


def prepare_bias_example(workspace: Path) -> BiasExampleInputs:
    """Create the small, offline bias example used by the notebook and tests."""

    workspace.mkdir(parents=True, exist_ok=True)
    system_path = workspace / "bias_system.yaml"
    system_path.write_text(
        yaml.safe_dump(
            {
                "version": 1,
                "sequences": [
                    {"protein": {"id": "A", "sequence": "MKRAAT"}},
                    {"ligand": {"id": "B", "smiles": "CCO"}},
                ],
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    reference_structure = workspace / "custom_reference.pdb"
    reference_structure.write_text("HEADER CUSTOM\n", encoding="utf-8")
    protein_reference_path = workspace / "custom_protein.csv"
    protein_reference_path.write_text(
        "sequence,dataset_name,source_structure_path\n"
        f"MKRAAT,private_proteins,{reference_structure.name}\n",
        encoding="utf-8",
    )
    ligand_reference_path = workspace / "custom_ligand.csv"
    ligand_reference_path.write_text(
        "smiles,dataset_name,source_reference_path\n"
        f"CCO,private_ligands,{reference_structure.name}\n",
        encoding="utf-8",
    )
    return BiasExampleInputs(
        system_path,
        protein_reference_path,
        ligand_reference_path,
    )


def build_bias_command(inputs: BiasExampleInputs, output_dir: Path) -> list[str]:
    return [
        "cofolder",
        "bias",
        "--system_path",
        str(inputs.system_path),
        "--wrk_dir",
        str(output_dir),
        "--custom_protein_reference_path",
        str(inputs.protein_reference_path),
        "--custom_ligand_reference_path",
        str(inputs.ligand_reference_path),
    ]


def build_validate_command(
    examples_dir: Path,
    output_dir: Path,
    *,
    system_name: str = "system.yaml",
    conformer_sdf: bool = False,
) -> list[str]:
    command = [
        "cofolder",
        "validate",
        "-s",
        str(examples_dir / system_name),
        "-o",
        str(examples_dir / "options.yaml"),
        "-w",
        str(output_dir),
        "--runner",
        "boltz2",
    ]
    if conformer_sdf:
        command.extend(
            ["--conformers", "sdf", "--sdf_file", str(examples_dir / "ethanol.sdf")]
        )
    return command


def build_screen_command(
    examples_dir: Path,
    output_dir: Path,
    *,
    library_name: str = "ligand_screen.csv",
) -> list[str]:
    command = [
        "cofolder",
        "screen",
        "-s",
        str(examples_dir / "system_screen.yaml"),
        "-o",
        str(examples_dir / "options.yaml"),
        "-c",
        str(examples_dir / library_name),
        "--ligand_chain",
        "B",
        "-w",
        str(output_dir),
        "--runner",
        "boltz2",
    ]
    if library_name.endswith(".csv"):
        command.extend(
            [
                "--col_id",
                "Name",
                "--smiles_column",
                "SMILES",
                "--merge_data",
                "pIC50",
            ]
        )
    else:
        command.extend(["--id_property", "ID"])
    return command


def build_oracle_command(examples_dir: Path, output_dir: Path) -> list[str]:
    return [
        "cofolder",
        "oracle",
        "-s",
        str(examples_dir / "system.yaml"),
        "-o",
        str(examples_dir / "options.yaml"),
        "--input_smiles",
        "CCO",
        "--output_metric",
        "system__confidence_score",
        "--scoring_functions",
        "confidence_metrics",
        "--aggregate",
        "first",
        "-w",
        str(output_dir),
        "--runner",
        "boltz2",
    ]


def with_preflight(command: list[str]) -> list[str]:
    """Return a copy of a workflow command configured for no-service validation."""

    return [*command, "--preflight_only"]
