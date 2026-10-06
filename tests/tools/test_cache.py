from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import pickle
from types import SimpleNamespace

import pytest
from rdkit import Chem

from cofolder.tools import cache


def _write_sdf(path: Path, records):
    writer = Chem.SDWriter(str(path))
    for identifier, smiles in records:
        molecule = Chem.MolFromSmiles(smiles)
        if identifier is not None:
            molecule.SetProp("ID", identifier)
        writer.write(molecule)
    writer.close()


def _complete(path: Path):
    path.mkdir(parents=True, exist_ok=True)
    (path / "boltz2_conf.ckpt").write_bytes(b"checkpoint")
    (path / "mols.tar").write_bytes(b"archive")
    (path / "mols").mkdir(exist_ok=True)
    for name in ("ALA.pkl", "GLY.pkl"):
        (path / "mols" / name).write_bytes(b"ccd")


def test_population_validates_before_writes_and_never_runs_setup(temp_dir):
    root = temp_dir / "cache"
    with pytest.raises(FileNotFoundError):
        cache.populate_ccd_cache_from_sdf(temp_dir / "missing.sdf", "ID", cache_path=root)
    assert not root.exists()

    source = temp_dir / "ligands.sdf"
    _write_sdf(source, [("ETH", "CCO")])
    with pytest.raises(RuntimeError, match="setup-boltz2-cache"):
        cache.populate_ccd_cache_from_sdf(source, "ID", cache_path=root)
    assert not root.exists()


def test_population_rejects_invalid_ids_duplicates_and_conflict_mode(temp_dir):
    root = temp_dir / "cache"
    source = temp_dir / "ligands.sdf"
    _write_sdf(source, [("ETH", "CCO")])
    with pytest.raises(ValueError, match="Invalid on_conflict mode"):
        cache.populate_ccd_cache_from_sdf(
            source, "ID", cache_path=root, on_conflict="keep"
        )
    missing = temp_dir / "missing-id.sdf"
    _write_sdf(missing, [(None, "CCO")])
    with pytest.raises(ValueError, match="missing required property"):
        cache.populate_ccd_cache_from_sdf(missing, "ID", cache_path=root)
    duplicate = temp_dir / "duplicate.sdf"
    _write_sdf(duplicate, [("ETH", "CCO"), ("ETH", "CCN")])
    with pytest.raises(ValueError, match="Duplicate molecule IDs.*ETH"):
        cache.populate_ccd_cache_from_sdf(duplicate, "ID", cache_path=root)
    invalid = temp_dir / "invalid-id.sdf"
    _write_sdf(invalid, [("TOO_LONG", "CCO")])
    with pytest.raises(ValueError, match="Invalid CCD identifiers"):
        cache.populate_ccd_cache_from_sdf(invalid, "ID", cache_path=root)
    assert not root.exists()


def test_population_preserves_conflict_modes(monkeypatch, temp_dir):
    root = temp_dir / "cache"
    _complete(root)
    (root / "mols" / "ETH.pkl").write_bytes(b"existing")
    source = temp_dir / "ligands.sdf"
    _write_sdf(source, [("ETH", "CCO"), ("ETN", "CCN")])
    converted = []
    monkeypatch.setattr(
        "cofolder.modules.entities.ligand.mol_to_ccd",
        lambda identifier, molecule, boltz_path: converted.append(identifier),
    )
    cache.populate_ccd_cache_from_sdf(
        source, "ID", cache_path=root, on_conflict="use_cache"
    )
    assert converted == ["ETN"]
    converted.clear()
    cache.populate_ccd_cache_from_sdf(source, "ID", cache_path=root)
    assert converted == ["ETH", "ETN"]


def test_population_creates_readable_ccd_artifacts_for_multiple_records(
    monkeypatch, temp_dir
):
    root = temp_dir / "cache"
    _complete(root)
    source = temp_dir / "ligands.sdf"
    _write_sdf(source, [("ETH", "CCO"), ("ETN", "CCN")])
    parsed = SimpleNamespace(
        rdkit_bounds_constraints=[],
        chiral_atom_constraints=[],
        stereo_bond_constraints=[],
        planar_ring_5_constraints=[],
        planar_ring_6_constraints=[],
        planar_bond_constraints=[],
    )
    monkeypatch.setattr(
        "cofolder.modules.entities.ligand._load_parse_ccd_residue",
        lambda: lambda *args: parsed,
    )

    cache.populate_ccd_cache_from_sdf(source, "ID", cache_path=root)

    for identifier, expected_smiles in (("ETH", "CCO"), ("ETN", "CCN")):
        artifact = root / "mols" / f"{identifier}.pkl"
        assert artifact.is_file()
        with artifact.open("rb") as stream:
            molecule = pickle.load(stream)
        assert Chem.MolToSmiles(molecule, isomericSmiles=True) == expected_smiles
        assert pickle.loads(bytes.fromhex(molecule.GetProp("MOL_NAME"))) == identifier


def test_population_surfaces_conversion_failures_with_identifier(
    monkeypatch, temp_dir
):
    root = temp_dir / "cache"
    _complete(root)
    source = temp_dir / "ligands.sdf"
    _write_sdf(source, [("ETH", "CCO")])
    monkeypatch.setattr(
        "cofolder.modules.entities.ligand.mol_to_ccd",
        lambda *args, **kwargs: (_ for _ in ()).throw(ValueError("bad conversion")),
    )

    with pytest.raises(RuntimeError, match="ETH: bad conversion"):
        cache.populate_ccd_cache_from_sdf(source, "ID", cache_path=root)

    assert not (root / "mols" / "ETH.pkl").exists()


def test_explicit_setup_is_locked_and_verified(monkeypatch, temp_dir):
    root = temp_dir / "cache"
    calls = []

    def complete(path):
        calls.append(path)
        _complete(path)

    monkeypatch.setattr(cache, "_download_boltz2_cache", complete)
    with ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(lambda _: cache.setup_boltz2_cache(root), range(2)))
    assert results == [root, root]
    assert calls == [root]


def test_setup_failure_is_explicit(monkeypatch, temp_dir):
    monkeypatch.setattr(cache, "_download_boltz2_cache", lambda path: None)
    with pytest.raises(RuntimeError, match="did not produce required components"):
        cache.setup_boltz2_cache(temp_dir / "cache")
