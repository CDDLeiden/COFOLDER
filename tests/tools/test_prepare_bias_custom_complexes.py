from __future__ import annotations

from pathlib import Path
import json

import pandas as pd
import pytest
import yaml

from cofolder.modules.analytics.bias_database import validate_custom_bias_reference_bundle
from cofolder.tools.prepare_bias_custom_complexes import prepare


REPOSITORY = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize("structure_name", ["4HJO.cif", "4HJO.pdb"])
def test_prepares_selected_structure_with_stable_pair_identity(
    temp_dir, structure_name
):
    structure = REPOSITORY / "src/cofolder/resources/examples" / structure_name
    manifest = temp_dir / "custom.yaml"
    manifest.write_text(
        yaml.safe_dump(
            {
                "schema_version": 1,
                "dataset_name": "two-complex-fixture",
                "complexes": [
                    {
                        "complex_id": "custom-4hjo",
                        "structure_path": str(structure),
                        "proteins": [{"chain_id": "A"}],
                        "ligands": [
                            {
                                "chain_id": "A",
                                "residue_name": "AQ4",
                                "residue_number": 1001,
                                "smiles": "C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1",
                            }
                        ],
                    },
                    {
                        "complex_id": "custom-4hjo-copy",
                        "structure_path": str(structure),
                        "proteins": [{"chain_id": "A"}],
                        "ligands": [{
                            "chain_id": "A", "residue_name": "AQ4",
                            "residue_number": 1001,
                            "smiles": "C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1",
                        }],
                    },
                ],
            }
        ),
        encoding="utf-8",
    )
    output = prepare(manifest, temp_dir / "bundle")
    bundle = validate_custom_bias_reference_bundle(output)
    proteins = pd.read_csv(bundle.protein_path)
    ligands = pd.read_csv(bundle.ligand_path)
    assert proteins["complex_id"].tolist() == ["custom-4hjo", "custom-4hjo-copy"]
    assert ligands["complex_id"].tolist() == ["custom-4hjo", "custom-4hjo-copy"]
    assert proteins["sequence"].str.len().iloc[0] > 100
    assert json.loads((output / "manifest.json").read_text())["complex_ids"] == ["custom-4hjo", "custom-4hjo-copy"]


def test_rejects_duplicate_complex_identity(temp_dir):
    manifest = temp_dir / "invalid.yaml"
    manifest.write_text(
        yaml.safe_dump(
            {
                "schema_version": 1,
                "dataset_name": "invalid",
                "complexes": [
                    {"complex_id": "same", "structure_path": "missing.pdb"},
                    {"complex_id": "same", "structure_path": "missing.pdb"},
                ],
            }
        ),
        encoding="utf-8",
    )
    with pytest.raises(ValueError):
        prepare(manifest, temp_dir / "bundle")
