"""Tests for reference-based reproduction metrics."""

import json
import logging

import pandas as pd

from cofolder.modules.analytics.reproduction import scaffold_reproduction_metrics


def _write_reference_pdb(path):
    path.write_text(
        """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00 20.00           C
ATOM      3  C   ALA A   1       2.000   0.000   0.000  1.00 20.00           C
ATOM      4  O   ALA A   1       3.000   0.000   0.000  1.00 20.00           O
ATOM      5  N   GLY A   2      10.000   0.000   0.000  1.00 20.00           N
ATOM      6  CA  GLY A   2      11.000   0.000   0.000  1.00 20.00           C
ATOM      7  C   GLY A   2      12.000   0.000   0.000  1.00 20.00           C
ATOM      8  O   GLY A   2      13.000   0.000   0.000  1.00 20.00           O
HETATM    9  C1  LIG L   1       1.200   1.200   0.000  1.00 20.00           C
HETATM   10  O1  LIG L   1       2.600   1.200   0.000  1.00 20.00           O
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


def _write_predicted_pdb(path, ligand_coords):
    c1, o1 = ligand_coords
    path.write_text(
        f"""\
ATOM      1  N   ALA A   1       0.000   0.000   0.100  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.000   0.000   0.100  1.00 20.00           C
ATOM      3  C   ALA A   1       2.000   0.000   0.100  1.00 20.00           C
ATOM      4  O   ALA A   1       3.000   0.000   0.100  1.00 20.00           O
ATOM      5  N   GLY A   2      10.000   0.000   0.100  1.00 20.00           N
ATOM      6  CA  GLY A   2      11.000   0.000   0.100  1.00 20.00           C
ATOM      7  C   GLY A   2      12.000   0.000   0.100  1.00 20.00           C
ATOM      8  O   GLY A   2      13.000   0.000   0.100  1.00 20.00           O
HETATM    9  C1  LIG Z   1      {c1[0]:7.3f}  {c1[1]:7.3f}  {c1[2]:7.3f}  1.00 20.00           C
HETATM   10  O1  LIG Z   1      {o1[0]:7.3f}  {o1[1]:7.3f}  {o1[2]:7.3f}  1.00 20.00           O
TER
END
"""
    )


def _write_predicted_single_atom_pdb(path):
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


def _make_system_df(cif_name):
    return pd.DataFrame(
        [
            {
                "idx": 0,
                "cif_file": cif_name,
                "model_name": "demo",
                "repeat": 1,
                "diffusion_sample": 0,
            }
        ]
    )


def _make_chain_df(cif_name, ifp_vector):
    return pd.DataFrame(
        [
            {
                "idx": 0,
                "CHAIN_ID": "A",
                "ENTITY_TYPE": "protein",
                "conf_chain_id": 0,
                "cif_file": cif_name,
                "model_name": "demo",
                "repeat": 1,
                "diffusion_sample": 0,
            },
            {
                "idx": 1,
                "CHAIN_ID": "Z",
                "ENTITY_TYPE": "ligand",
                "conf_chain_id": 1,
                "cif_file": cif_name,
                "model_name": "demo",
                "repeat": 1,
                "diffusion_sample": 0,
                "ifp_distance": json.dumps(ifp_vector),
            },
        ]
    )


def test_scaffold_reproduction_metrics_adds_schema_and_metrics(temp_dir):
    wrk_dir = temp_dir
    structures_dir = wrk_dir / "results" / "structures"
    structures_dir.mkdir(parents=True, exist_ok=True)

    reference_path = wrk_dir / "reference.pdb"
    predicted_path = structures_dir / "1_demo_model_0.pdb"

    _write_reference_pdb(reference_path)
    _write_predicted_pdb(predicted_path, ligand_coords=((1.3, 1.2, 0.0), (2.7, 1.2, 0.0)))

    out_system_df, out_chain_df = scaffold_reproduction_metrics(
        system_df=_make_system_df("1_demo_model_0.pdb"),
        chain_df=_make_chain_df("1_demo_model_0.pdb", ifp_vector=[1, 0]),
        reference_path=reference_path,
        wrk_dir=wrk_dir,
    )

    assert "pocket_coverage_ref" in out_chain_df.columns
    assert "sucos_ref" in out_chain_df.columns
    assert "sucos_shape_ref" in out_chain_df.columns
    assert "sucos_feature_ref" in out_chain_df.columns
    assert "pocket_coverage_ref_mean" in out_system_df.columns
    assert "sucos_ref_mean" in out_system_df.columns

    lig_row = out_chain_df[out_chain_df["ENTITY_TYPE"] == "ligand"].iloc[0]
    assert lig_row["ligand_pose_overlap_ref"] == lig_row["sucos_ref"]
    assert 0.0 <= float(lig_row["pocket_coverage_ref"]) <= 1.0
    assert 0.0 <= float(lig_row["sucos_ref"]) <= 1.0


def test_reference_dependent_metrics_remain_nan_without_reference(temp_dir):
    wrk_dir = temp_dir
    out_system_df, out_chain_df = scaffold_reproduction_metrics(
        system_df=_make_system_df("1_demo_model_0.pdb"),
        chain_df=_make_chain_df("1_demo_model_0.pdb", ifp_vector=[1, 0]),
        reference_path=None,
        wrk_dir=wrk_dir,
    )

    lig_row = out_chain_df[out_chain_df["ENTITY_TYPE"] == "ligand"].iloc[0]
    assert pd.isna(lig_row["pocket_coverage_ref"])
    assert pd.isna(lig_row["sucos_ref"])
    assert pd.isna(out_system_df["pocket_coverage_ref_mean"].iloc[0])
    assert pd.isna(out_system_df["sucos_ref_mean"].iloc[0])


def test_pocket_coverage_exact_overlap_is_one(temp_dir):
    wrk_dir = temp_dir
    structures_dir = wrk_dir / "results" / "structures"
    structures_dir.mkdir(parents=True, exist_ok=True)

    reference_path = wrk_dir / "reference.pdb"
    predicted_path = structures_dir / "1_demo_model_0.pdb"
    _write_reference_pdb(reference_path)
    _write_predicted_pdb(predicted_path, ligand_coords=((1.3, 1.2, 0.0), (2.7, 1.2, 0.0)))

    _, out_chain_df = scaffold_reproduction_metrics(
        system_df=_make_system_df("1_demo_model_0.pdb"),
        chain_df=_make_chain_df("1_demo_model_0.pdb", ifp_vector=[1, 0]),
        reference_path=reference_path,
        wrk_dir=wrk_dir,
        reproduction_metrics=["pocket_coverage"],
    )

    value = float(out_chain_df[out_chain_df["ENTITY_TYPE"] == "ligand"]["pocket_coverage_ref"].iloc[0])
    assert value == 1.0


def test_pocket_coverage_partial_overlap(temp_dir):
    wrk_dir = temp_dir
    structures_dir = wrk_dir / "results" / "structures"
    structures_dir.mkdir(parents=True, exist_ok=True)

    reference_path = wrk_dir / "reference.pdb"
    predicted_path = structures_dir / "1_demo_model_0.pdb"
    _write_reference_pdb(reference_path)
    _write_predicted_pdb(predicted_path, ligand_coords=((1.3, 1.2, 0.0), (2.7, 1.2, 0.0)))

    _, out_chain_df = scaffold_reproduction_metrics(
        system_df=_make_system_df("1_demo_model_0.pdb"),
        chain_df=_make_chain_df("1_demo_model_0.pdb", ifp_vector=[0, 1]),
        reference_path=reference_path,
        wrk_dir=wrk_dir,
        reproduction_metrics=["pocket_coverage"],
    )

    value = float(out_chain_df[out_chain_df["ENTITY_TYPE"] == "ligand"]["pocket_coverage_ref"].iloc[0])
    assert value == 0.0


def test_pocket_coverage_empty_reference_ifp_warns_and_nan(temp_dir, caplog):
    wrk_dir = temp_dir
    structures_dir = wrk_dir / "results" / "structures"
    structures_dir.mkdir(parents=True, exist_ok=True)

    reference_path = wrk_dir / "reference_far.pdb"
    predicted_path = structures_dir / "1_demo_model_0.pdb"
    _write_reference_pdb(reference_path)
    _write_predicted_pdb(predicted_path, ligand_coords=((50.0, 50.0, 50.0), (51.4, 50.0, 50.0)))

    caplog.set_level(logging.WARNING)
    _, out_chain_df = scaffold_reproduction_metrics(
        system_df=_make_system_df("1_demo_model_0.pdb"),
        chain_df=_make_chain_df("1_demo_model_0.pdb", ifp_vector=[1, 0]),
        reference_path=reference_path,
        wrk_dir=wrk_dir,
        reproduction_metrics=["pocket_coverage"],
        logger=logging.getLogger("test_reproduction"),
    )

    value = out_chain_df[out_chain_df["ENTITY_TYPE"] == "ligand"]["pocket_coverage_ref"].iloc[0]
    assert pd.isna(value)
    assert "Reference IFP has no active bits for pocket coverage" in caplog.text


def test_sucos_is_higher_for_near_identical_than_displaced(temp_dir):
    wrk_dir = temp_dir
    structures_dir = wrk_dir / "results" / "structures"
    structures_dir.mkdir(parents=True, exist_ok=True)

    reference_path = wrk_dir / "reference.pdb"
    _write_reference_pdb(reference_path)

    near_name = "1_demo_model_0.pdb"
    far_name = "1_demo_model_1.pdb"
    _write_predicted_pdb(
        structures_dir / near_name,
        ligand_coords=((1.25, 1.25, 0.0), (2.65, 1.25, 0.0)),
    )
    _write_predicted_pdb(
        structures_dir / far_name,
        ligand_coords=((8.0, 8.0, 8.0), (9.4, 8.0, 8.0)),
    )

    system_df = pd.DataFrame(
        [
            {"idx": 0, "cif_file": near_name, "model_name": "demo", "repeat": 1, "diffusion_sample": 0},
            {"idx": 1, "cif_file": far_name, "model_name": "demo", "repeat": 1, "diffusion_sample": 1},
        ]
    )

    chain_df = pd.DataFrame(
        [
            {"idx": 0, "CHAIN_ID": "A", "ENTITY_TYPE": "protein", "conf_chain_id": 0, "cif_file": near_name, "model_name": "demo", "repeat": 1, "diffusion_sample": 0},
            {"idx": 1, "CHAIN_ID": "Z", "ENTITY_TYPE": "ligand", "conf_chain_id": 1, "cif_file": near_name, "model_name": "demo", "repeat": 1, "diffusion_sample": 0, "ifp_distance": json.dumps([1, 0])},
            {"idx": 2, "CHAIN_ID": "A", "ENTITY_TYPE": "protein", "conf_chain_id": 0, "cif_file": far_name, "model_name": "demo", "repeat": 1, "diffusion_sample": 1},
            {"idx": 3, "CHAIN_ID": "Z", "ENTITY_TYPE": "ligand", "conf_chain_id": 1, "cif_file": far_name, "model_name": "demo", "repeat": 1, "diffusion_sample": 1, "ifp_distance": json.dumps([1, 0])},
        ]
    )

    _, out_chain_df = scaffold_reproduction_metrics(
        system_df=system_df,
        chain_df=chain_df,
        reference_path=reference_path,
        wrk_dir=wrk_dir,
        reproduction_metrics=["sucos"],
    )

    near_sucos = float(out_chain_df[out_chain_df["cif_file"] == near_name]["sucos_ref"].dropna().iloc[0])
    far_sucos = float(out_chain_df[out_chain_df["cif_file"] == far_name]["sucos_ref"].dropna().iloc[0])
    assert near_sucos > far_sucos
    assert 0.0 <= near_sucos <= 1.0
    assert 0.0 <= far_sucos <= 1.0


def test_ligand_uses_closest_reference_instance_when_duplicates_exist(temp_dir):
    wrk_dir = temp_dir
    structures_dir = wrk_dir / "results" / "structures"
    structures_dir.mkdir(parents=True, exist_ok=True)

    reference_path = wrk_dir / "reference_dup.pdb"
    predicted_path = structures_dir / "1_demo_model_0.pdb"

    _write_reference_with_duplicate_ligands_pdb(reference_path)
    _write_predicted_single_atom_pdb(predicted_path)

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
                "ifp_distance": json.dumps([1]),
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
    assert abs(ligand_rmsd - 1.5) < 1e-6
