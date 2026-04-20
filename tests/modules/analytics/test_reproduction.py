"""Tests for reference-based reproduction metrics."""

import json
import logging

import pandas as pd

from cofolder.modules.analytics.reproduction import scaffold_reproduction_metrics


def _write_reference_pdb(path, ligand_coords=((1.2, 1.2, 0.0), (2.6, 1.2, 0.0))):
    c1, o1 = ligand_coords
    path.write_text(
        f"""\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00 20.00           C
ATOM      3  C   ALA A   1       2.000   0.000   0.000  1.00 20.00           C
ATOM      4  O   ALA A   1       3.000   0.000   0.000  1.00 20.00           O
ATOM      5  N   GLY A   2      10.000   0.000   0.000  1.00 20.00           N
ATOM      6  CA  GLY A   2      11.000   0.000   0.000  1.00 20.00           C
ATOM      7  C   GLY A   2      12.000   0.000   0.000  1.00 20.00           C
ATOM      8  O   GLY A   2      13.000   0.000   0.000  1.00 20.00           O
HETATM    9  C1  LIG L   1    {c1[0]:8.3f}{c1[1]:8.3f}{c1[2]:8.3f}  1.00 20.00           C
HETATM   10  O1  LIG L   1    {o1[0]:8.3f}{o1[1]:8.3f}{o1[2]:8.3f}  1.00 20.00           O
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


def _write_predicted_pdb(path, ligand_coords, protein_offset=(0.0, 0.0, 0.1)):
    c1, o1 = ligand_coords
    px, py, pz = protein_offset
    path.write_text(
        f"""\
ATOM      1  N   ALA A   1    {0.0 + px:8.3f}{0.0 + py:8.3f}{0.0 + pz:8.3f}  1.00 20.00           N
ATOM      2  CA  ALA A   1    {1.0 + px:8.3f}{0.0 + py:8.3f}{0.0 + pz:8.3f}  1.00 20.00           C
ATOM      3  C   ALA A   1    {2.0 + px:8.3f}{0.0 + py:8.3f}{0.0 + pz:8.3f}  1.00 20.00           C
ATOM      4  O   ALA A   1    {3.0 + px:8.3f}{0.0 + py:8.3f}{0.0 + pz:8.3f}  1.00 20.00           O
ATOM      5  N   GLY A   2    {10.0 + px:8.3f}{0.0 + py:8.3f}{0.0 + pz:8.3f}  1.00 20.00           N
ATOM      6  CA  GLY A   2    {11.0 + px:8.3f}{0.0 + py:8.3f}{0.0 + pz:8.3f}  1.00 20.00           C
ATOM      7  C   GLY A   2    {12.0 + px:8.3f}{0.0 + py:8.3f}{0.0 + pz:8.3f}  1.00 20.00           C
ATOM      8  O   GLY A   2    {13.0 + px:8.3f}{0.0 + py:8.3f}{0.0 + pz:8.3f}  1.00 20.00           O
HETATM    9  C1  LIG Z   1    {c1[0]:8.3f}{c1[1]:8.3f}{c1[2]:8.3f}  1.00 20.00           C
HETATM   10  O1  LIG Z   1    {o1[0]:8.3f}{o1[1]:8.3f}{o1[2]:8.3f}  1.00 20.00           O
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


def _write_reference_sequence_fallback_pdb(path):
    path.write_text(
        """\
ATOM      1  N   ALA A 670       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A 670       1.000   0.000   0.000  1.00 20.00           C
ATOM      3  C   ALA A 670       2.000   1.000   0.000  1.00 20.00           C
ATOM      4  O   ALA A 670       2.500   2.000   0.000  1.00 20.00           O
ATOM      5  N   GLY A 671      10.000   0.000   1.000  1.00 20.00           N
ATOM      6  CA  GLY A 671      11.000   0.000   1.000  1.00 20.00           C
ATOM      7  C   GLY A 671      12.000   1.000   1.000  1.00 20.00           C
ATOM      8  O   GLY A 671      12.500   2.000   1.000  1.00 20.00           O
ATOM      9  N   SER A 672      10.000  10.000   2.000  1.00 20.00           N
ATOM     10  CA  SER A 672      11.000  10.000   2.000  1.00 20.00           C
ATOM     11  C   SER A 672      12.000  11.000   2.000  1.00 20.00           C
ATOM     12  O   SER A 672      12.500  12.000   2.000  1.00 20.00           O
HETATM   13  C1  LIG L   1      11.200  10.800   2.500  1.00 20.00           C
HETATM   14  O1  LIG L   1      12.600  10.800   2.500  1.00 20.00           O
TER
END
"""
    )


def _write_predicted_sequence_fallback_pdb(path):
    path.write_text(
        """\
ATOM      1  N   MET X   1      70.000  40.000  30.000  1.00 20.00           N
ATOM      2  CA  MET X   1      71.000  40.000  30.000  1.00 20.00           C
ATOM      3  C   MET X   1      72.000  41.000  30.000  1.00 20.00           C
ATOM      4  O   MET X   1      72.500  42.000  30.000  1.00 20.00           O
ATOM      5  N   THR X   2      80.000  45.000  35.000  1.00 20.00           N
ATOM      6  CA  THR X   2      81.000  45.000  35.000  1.00 20.00           C
ATOM      7  C   THR X   2      82.000  46.000  35.000  1.00 20.00           C
ATOM      8  O   THR X   2      82.500  47.000  35.000  1.00 20.00           O
ATOM      9  N   ALA X   3     100.000 100.000 100.000  1.00 20.00           N
ATOM     10  CA  ALA X   3     101.000 100.000 100.000  1.00 20.00           C
ATOM     11  C   ALA X   3     102.000 101.000 100.000  1.00 20.00           C
ATOM     12  O   ALA X   3     102.500 102.000 100.000  1.00 20.00           O
ATOM     13  N   GLY X   4     110.000 100.000 101.000  1.00 20.00           N
ATOM     14  CA  GLY X   4     111.000 100.000 101.000  1.00 20.00           C
ATOM     15  C   GLY X   4     112.000 101.000 101.000  1.00 20.00           C
ATOM     16  O   GLY X   4     112.500 102.000 101.000  1.00 20.00           O
ATOM     17  N   SER X   5     110.000 110.000 102.000  1.00 20.00           N
ATOM     18  CA  SER X   5     111.000 110.000 102.000  1.00 20.00           C
ATOM     19  C   SER X   5     112.000 111.000 102.000  1.00 20.00           C
ATOM     20  O   SER X   5     112.500 112.000 102.000  1.00 20.00           O
HETATM   21  C1  LIG Z   1     111.200 110.800 102.500  1.00 20.00           C
HETATM   22  O1  LIG Z   1     112.600 110.800 102.500  1.00 20.00           O
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
    assert "pocket_coverage_custom" in out_chain_df.columns
    assert "sucos_ref" in out_chain_df.columns
    assert "sucos_shape_ref" in out_chain_df.columns
    assert "sucos_feature_ref" in out_chain_df.columns
    assert "pocket_coverage_ref_mean" in out_system_df.columns
    assert "pocket_coverage_custom_mean" in out_system_df.columns
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


def test_pocket_coverage_can_be_computed_from_custom_reference_without_reference(temp_dir):
    wrk_dir = temp_dir

    out_system_df, out_chain_df = scaffold_reproduction_metrics(
        system_df=_make_system_df("1_demo_model_0.pdb"),
        chain_df=_make_chain_df("1_demo_model_0.pdb", ifp_vector=[1, 1]),
        reference_path=None,
        wrk_dir=wrk_dir,
        pocket_coverage_reference="0100000000",
        reproduction_metrics=["pocket_coverage"],
    )

    lig_row = out_chain_df[out_chain_df["ENTITY_TYPE"] == "ligand"].iloc[0]
    assert float(lig_row["pocket_coverage_custom"]) == 1.0
    assert pd.isna(lig_row["pocket_coverage_ref"])
    assert pd.isna(lig_row["sucos_ref"])
    assert float(out_system_df["pocket_coverage_custom_mean"].iloc[0]) == 1.0
    assert pd.isna(out_system_df["pocket_coverage_ref_mean"].iloc[0])


def test_custom_reference_residue_formats_with_reference(temp_dir):
    wrk_dir = temp_dir
    structures_dir = wrk_dir / "results" / "structures"
    structures_dir.mkdir(parents=True, exist_ok=True)

    reference_path = wrk_dir / "reference.pdb"
    predicted_path = structures_dir / "1_demo_model_0.pdb"
    _write_reference_pdb(reference_path)
    _write_predicted_pdb(predicted_path, ligand_coords=((1.3, 1.2, 0.0), (2.7, 1.2, 0.0)))

    # residue order for protein chain A is [1, 2]; choose residue 1 as reference pocket
    for custom_ref in ["1", "A1", "A1 S1"]:
        _, out_chain_df = scaffold_reproduction_metrics(
            system_df=_make_system_df("1_demo_model_0.pdb"),
            chain_df=_make_chain_df("1_demo_model_0.pdb", ifp_vector=[1, 0]),
            reference_path=reference_path,
            wrk_dir=wrk_dir,
            pocket_coverage_reference=custom_ref,
            reproduction_metrics=["pocket_coverage"],
        )
        value = float(
            out_chain_df[out_chain_df["ENTITY_TYPE"] == "ligand"]["pocket_coverage_custom"].iloc[0]
        )
        assert value == 1.0


def test_custom_reference_residue_label_mismatch_warns(temp_dir, caplog):
    wrk_dir = temp_dir
    structures_dir = wrk_dir / "results" / "structures"
    structures_dir.mkdir(parents=True, exist_ok=True)

    reference_path = wrk_dir / "reference.pdb"
    predicted_path = structures_dir / "1_demo_model_0.pdb"
    _write_reference_pdb(reference_path)
    _write_predicted_pdb(predicted_path, ligand_coords=((1.3, 1.2, 0.0), (2.7, 1.2, 0.0)))

    caplog.set_level(logging.WARNING)
    _, out_chain_df = scaffold_reproduction_metrics(
        system_df=_make_system_df("1_demo_model_0.pdb"),
        chain_df=_make_chain_df("1_demo_model_0.pdb", ifp_vector=[0, 1]),
        reference_path=reference_path,
        wrk_dir=wrk_dir,
        pocket_coverage_reference="A2",  # residue 2 is GLY in this fixture
        reproduction_metrics=["pocket_coverage"],
        logger=logging.getLogger("test_reproduction"),
    )

    value = float(
        out_chain_df[out_chain_df["ENTITY_TYPE"] == "ligand"]["pocket_coverage_custom"].iloc[0]
    )
    assert value == 1.0
    assert "Custom pocket reference residue label mismatch" in caplog.text


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
    _write_reference_pdb(
        reference_path,
        ligand_coords=((50.0, 50.0, 50.0), (51.4, 50.0, 50.0)),
    )
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


def test_ligand_rmsd_is_computed_after_alignment_to_reference(temp_dir):
    wrk_dir = temp_dir
    structures_dir = wrk_dir / "results" / "structures"
    structures_dir.mkdir(parents=True, exist_ok=True)

    reference_path = wrk_dir / "reference.pdb"
    predicted_path = structures_dir / "1_demo_model_0.pdb"

    _write_reference_pdb(reference_path)
    _write_predicted_pdb(
        predicted_path,
        ligand_coords=((101.2, 101.2, 100.0), (102.6, 101.2, 100.0)),
        protein_offset=(100.0, 100.0, 100.0),
    )

    _, out_chain_df = scaffold_reproduction_metrics(
        system_df=_make_system_df("1_demo_model_0.pdb"),
        chain_df=_make_chain_df("1_demo_model_0.pdb", ifp_vector=[1, 0]),
        reference_path=reference_path,
        wrk_dir=wrk_dir,
        reproduction_metrics=["protein_rmsd", "ligand_rmsd"],
    )

    ligand_rmsd = float(
        out_chain_df.loc[out_chain_df["ENTITY_TYPE"] == "ligand", "ligand_rmsd_ref"].iloc[0]
    )
    protein_rmsd = float(
        out_chain_df.loc[out_chain_df["ENTITY_TYPE"] == "protein", "protein_rmsd_ref"].iloc[0]
    )

    assert ligand_rmsd < 1e-4
    assert protein_rmsd < 1e-6


def test_sequence_fallback_alignment_handles_chain_and_residue_id_mismatch(temp_dir):
    wrk_dir = temp_dir
    structures_dir = wrk_dir / "results" / "structures"
    structures_dir.mkdir(parents=True, exist_ok=True)

    reference_path = wrk_dir / "reference_seq_fallback.pdb"
    predicted_name = "1_demo_model_0.pdb"
    predicted_path = structures_dir / predicted_name

    _write_reference_sequence_fallback_pdb(reference_path)
    _write_predicted_sequence_fallback_pdb(predicted_path)

    system_df = _make_system_df(predicted_name)
    chain_df = pd.DataFrame(
        [
            {
                "idx": 0,
                "CHAIN_ID": "X",
                "ENTITY_TYPE": "protein",
                "conf_chain_id": 0,
                "cif_file": predicted_name,
                "model_name": "demo",
                "repeat": 1,
                "diffusion_sample": 0,
            },
            {
                "idx": 1,
                "CHAIN_ID": "Z",
                "ENTITY_TYPE": "ligand",
                "conf_chain_id": 1,
                "cif_file": predicted_name,
                "model_name": "demo",
                "repeat": 1,
                "diffusion_sample": 0,
                "ifp_distance": json.dumps([1, 0, 0]),
            },
        ]
    )

    _, out_chain_df = scaffold_reproduction_metrics(
        system_df=system_df,
        chain_df=chain_df,
        reference_path=reference_path,
        wrk_dir=wrk_dir,
        reproduction_metrics=["protein_rmsd", "ligand_rmsd"],
    )

    ligand_rmsd = float(
        out_chain_df.loc[out_chain_df["ENTITY_TYPE"] == "ligand", "ligand_rmsd_ref"].iloc[0]
    )
    protein_rmsd = float(
        out_chain_df.loc[out_chain_df["ENTITY_TYPE"] == "protein", "protein_rmsd_ref"].iloc[0]
    )

    assert ligand_rmsd < 1e-4
    assert protein_rmsd < 1e-4


def test_sequence_fallback_alignment_logs_strategy(temp_dir, caplog):
    wrk_dir = temp_dir
    structures_dir = wrk_dir / "results" / "structures"
    structures_dir.mkdir(parents=True, exist_ok=True)

    reference_path = wrk_dir / "reference_seq_fallback.pdb"
    predicted_name = "1_demo_model_0.pdb"
    predicted_path = structures_dir / predicted_name

    _write_reference_sequence_fallback_pdb(reference_path)
    _write_predicted_sequence_fallback_pdb(predicted_path)

    logger = logging.getLogger("test_reproduction")
    caplog.set_level(logging.DEBUG)

    scaffold_reproduction_metrics(
        system_df=_make_system_df(predicted_name),
        chain_df=pd.DataFrame(
            [
                {
                    "idx": 0,
                    "CHAIN_ID": "X",
                    "ENTITY_TYPE": "protein",
                    "conf_chain_id": 0,
                    "cif_file": predicted_name,
                    "model_name": "demo",
                    "repeat": 1,
                    "diffusion_sample": 0,
                },
                {
                    "idx": 1,
                    "CHAIN_ID": "Z",
                    "ENTITY_TYPE": "ligand",
                    "conf_chain_id": 1,
                    "cif_file": predicted_name,
                    "model_name": "demo",
                    "repeat": 1,
                    "diffusion_sample": 0,
                    "ifp_distance": json.dumps([1, 0, 0]),
                },
            ]
        ),
        reference_path=reference_path,
        wrk_dir=wrk_dir,
        reproduction_metrics=["protein_rmsd", "ligand_rmsd"],
        logger=logger,
    )

    assert "using sequence_fallback" in caplog.text
    assert "using naive_first_n" not in caplog.text
