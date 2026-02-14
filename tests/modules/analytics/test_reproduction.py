"""Tests for reproduction metric scaffolding and phase-2 RMSD metrics."""

import pandas as pd

from cofolder.modules.analytics.reproduction import scaffold_reproduction_metrics


def _write_reference_pdb(path):
    path.write_text(
        """\
ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00 20.00           N
ATOM      2  CA  ALA A   1      11.000  10.000  10.000  1.00 20.00           C
ATOM      3  C   ALA A   1      12.000  10.000  10.000  1.00 20.00           C
ATOM      4  O   ALA A   1      13.000  10.000  10.000  1.00 20.00           O
HETATM    5  C1  LIG L   1       5.000   5.000   5.000  1.00 20.00           C
TER
END
"""
    )


def _write_reference_with_duplicate_ligands_pdb(path):
    path.write_text(
        """\
ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00 20.00           N
ATOM      2  CA  ALA A   1      11.000  10.000  10.000  1.00 20.00           C
ATOM      3  C   ALA A   1      12.000  10.000  10.000  1.00 20.00           C
ATOM      4  O   ALA A   1      13.000  10.000  10.000  1.00 20.00           O
HETATM    5  C1  LIG L   1       5.000   5.000   5.000  1.00 20.00           C
HETATM    6  C1  LIG M   1      15.000   5.000   5.000  1.00 20.00           C
TER
END
"""
    )


def _write_predicted_pdb(path):
    path.write_text(
        """\
ATOM      1  N   ALA A   1      10.700  10.000  10.000  1.00 20.00           N
ATOM      2  CA  ALA A   1      11.000  10.000  10.000  1.00 20.00           C
ATOM      3  C   ALA A   1      12.300  10.000  10.000  1.00 20.00           C
ATOM      4  O   ALA A   1      13.400  10.000  10.000  1.00 20.00           O
HETATM    5  C1  LIG Z   1       6.500   5.000   5.000  1.00 20.00           C
TER
END
"""
    )


def test_scaffold_reproduction_metrics_adds_schema_and_rmsd(temp_dir):
    wrk_dir = temp_dir
    structures_dir = wrk_dir / "results" / "structures"
    structures_dir.mkdir(parents=True, exist_ok=True)

    reference_path = wrk_dir / "reference.pdb"
    predicted_path = structures_dir / "1_demo_model_0.pdb"

    _write_reference_pdb(reference_path)
    _write_predicted_pdb(predicted_path)

    system_df = pd.DataFrame(
        [
            {
                "idx": 0,
                "cif_file": "1_demo_model_0.pdb",
                "model_name": "demo",
                "repeat": 1,
                "diffusion_sample": 0,
            }
        ]
    )

    chain_df = pd.DataFrame(
        [
            {
                "idx": 0,
                "CHAIN_ID": "A",
                "ENTITY_TYPE": "protein",
                "conf_chain_id": 0,
                "cif_file": "1_demo_model_0.pdb",
                "model_name": "demo",
                "repeat": 1,
                "diffusion_sample": 0,
            },
            {
                "idx": 1,
                "CHAIN_ID": "Z",
                "ENTITY_TYPE": "ligand",
                "conf_chain_id": 1,
                "cif_file": "1_demo_model_0.pdb",
                "model_name": "demo",
                "repeat": 1,
                "diffusion_sample": 0,
            },
        ]
    )

    out_system_df, out_chain_df = scaffold_reproduction_metrics(
        system_df=system_df,
        chain_df=chain_df,
        reference_path=reference_path,
        wrk_dir=wrk_dir,
    )

    assert "protein_rmsd_ref" in out_chain_df.columns
    assert "ligand_rmsd_ref" in out_chain_df.columns
    assert "protein_rmsd_ref_mean" in out_system_df.columns
    assert "ligand_rmsd_ref_mean" in out_system_df.columns

    protein_rows = out_chain_df[out_chain_df["ENTITY_TYPE"] != "ligand"]
    ligand_rows = out_chain_df[out_chain_df["ENTITY_TYPE"] == "ligand"]

    assert protein_rows["protein_rmsd_ref"].notna().any()
    assert ligand_rows["ligand_rmsd_ref"].notna().any()

    assert out_system_df["protein_rmsd_ref_mean"].notna().any()
    assert out_system_df["ligand_rmsd_ref_mean"].notna().any()


def test_ligand_uses_closest_reference_instance_when_duplicates_exist(temp_dir):
    wrk_dir = temp_dir
    structures_dir = wrk_dir / "results" / "structures"
    structures_dir.mkdir(parents=True, exist_ok=True)

    reference_path = wrk_dir / "reference_dup.pdb"
    predicted_path = structures_dir / "1_demo_model_0.pdb"

    _write_reference_with_duplicate_ligands_pdb(reference_path)
    _write_predicted_pdb(predicted_path)

    system_df = pd.DataFrame(
        [
            {
                "idx": 0,
                "cif_file": "1_demo_model_0.pdb",
                "model_name": "demo",
                "repeat": 1,
                "diffusion_sample": 0,
            }
        ]
    )

    chain_df = pd.DataFrame(
        [
            {
                "idx": 0,
                "CHAIN_ID": "A",
                "ENTITY_TYPE": "protein",
                "conf_chain_id": 0,
                "cif_file": "1_demo_model_0.pdb",
                "model_name": "demo",
                "repeat": 1,
                "diffusion_sample": 0,
            },
            {
                "idx": 1,
                "CHAIN_ID": "Z",
                "ENTITY_TYPE": "ligand",
                "conf_chain_id": 1,
                "cif_file": "1_demo_model_0.pdb",
                "model_name": "demo",
                "repeat": 1,
                "diffusion_sample": 0,
            },
        ]
    )

    _, out_chain_df = scaffold_reproduction_metrics(
        system_df=system_df,
        chain_df=chain_df,
        reference_path=reference_path,
        wrk_dir=wrk_dir,
    )

    ligand_rmsd = float(
        out_chain_df.loc[out_chain_df["ENTITY_TYPE"] == "ligand", "ligand_rmsd_ref"].iloc[0]
    )

    # Predicted ligand x=6.5 is closer to reference ligand at x=5.0 than x=15.0.
    assert abs(ligand_rmsd - 1.5) < 1e-6
