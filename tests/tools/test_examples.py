from __future__ import annotations

import pytest

from cofolder.resources.examples import copy_examples, examples_root


EXPECTED_NAMES = {
    "4HJO.cif",
    "4HJO.pdb",
    "README.md",
    "bias_matrix.yaml",
    "custom_bias_complexes.yaml",
    "ifp_clustering_demo.py",
    "ligand_screen.csv",
    "options.yaml",
    "system.yaml",
    "system_covalent.yaml",
    "system_nucleic_acid.yaml",
    "system_screen.yaml",
}


def _resource_bytes() -> dict[str, bytes]:
    return {
        item.name: item.read_bytes()
        for item in examples_root().iterdir()
        if item.is_file() and item.name in EXPECTED_NAMES
    }


def test_example_locator_and_copy_are_complete_and_byte_identical(temp_dir):
    expected = _resource_bytes()
    destination = temp_dir / "workspace with spaces" / "examples"

    copied_to = copy_examples(destination)

    assert copied_to == destination.resolve()
    assert expected.keys() == EXPECTED_NAMES
    assert {path.name for path in destination.iterdir()} == EXPECTED_NAMES
    assert {
        path.name: path.read_bytes() for path in destination.iterdir()
    } == expected


def test_copy_refuses_all_conflicts_before_writing(temp_dir):
    destination = temp_dir / "examples"
    destination.mkdir()
    (destination / "system.yaml").write_text("user content", encoding="utf-8")

    with pytest.raises(FileExistsError, match="system.yaml"):
        copy_examples(destination)

    assert {path.name for path in destination.iterdir()} == {"system.yaml"}
    assert (destination / "system.yaml").read_text(encoding="utf-8") == "user content"


def test_copy_overwrites_managed_files_and_preserves_unrelated_files(temp_dir):
    destination = temp_dir / "examples"
    destination.mkdir()
    (destination / "system.yaml").write_text("stale", encoding="utf-8")
    unrelated = destination / "notes.txt"
    unrelated.write_text("keep me", encoding="utf-8")

    copy_examples(destination, overwrite=True)

    assert (destination / "system.yaml").read_bytes() == _resource_bytes()["system.yaml"]
    assert unrelated.read_text(encoding="utf-8") == "keep me"


def test_copy_never_replaces_a_directory_with_a_managed_file(temp_dir):
    destination = temp_dir / "examples"
    (destination / "system.yaml").mkdir(parents=True)

    with pytest.raises(FileExistsError, match="non-file"):
        copy_examples(destination, overwrite=True)

    assert (destination / "system.yaml").is_dir()
    assert {path.name for path in destination.iterdir()} == {"system.yaml"}
