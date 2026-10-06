# Ligand Handling Tutorial

COFOLDER 1.0 accepts native SMILES, SDF/MOL structural records, and preprocessed
Chemical Component Dictionary (CCD) entries. Copy the distributable fixtures before
working through the examples:

```bash
cofolder-tools copy-examples ./cofolder-examples
cd cofolder-examples
```

## Route 1: Native SMILES

Put SMILES directly in a system YAML when the selected runner supports native SMILES:

```yaml
version: 1
sequences:
  - protein:
      id: A
      sequence: ACDEFGHIK
  - ligand:
      id: B
      smiles: CCO
```

Validate and canonicalize SMILES through the shared public boundary when preparing
inputs programmatically:

```python
from cofolder.modules.input import LigandSourceIdentity, validate_smiles

ligand = validate_smiles(
    "C[C@H](O)F",
    source=LigandSourceIdentity("entity:1", ("B",)),
)
print(ligand.canonical_smiles)
```

The validator preserves defined stereochemistry in canonical isomeric SMILES. It does
not invent unspecified stereocentres or protonation/tautomer states.

## Route 2: SDF Conformers for Validate

The system still declares canonical chemistry as SMILES. `--conformers sdf` supplies
coordinates for that same molecule:

```bash
cofolder validate -s system_screen.yaml -o options.yaml \
  --conformers sdf --sdf_file ethanol.sdf \
  --runner boltz2 -w ./validate-sdf
```

The first valid SDF molecule must normalize to the same isomeric SMILES as the selected
system ligand. A chemistry mismatch, unreadable record, or empty file fails before
backend execution. Direct `sdf`, `mol`, and `pdb` fields are not system-YAML ligand
fields.

Use generated conformers when no external coordinates are required:

```bash
cofolder validate -s system_screen.yaml -o options.yaml --conformers 2D
cofolder validate -s system_screen.yaml -o options.yaml --conformers 3D
```

## Route 3: SDF and MOL Libraries for Screen

Screen consumes structural libraries directly:

```bash
cofolder screen -s system_screen.yaml -o options.yaml \
  -c ethanol.sdf --ligand_chain B --id_property ID \
  --runner boltz2 -w ./screen-sdf

cofolder screen -s system_screen.yaml -o options.yaml \
  -c ethanol.mol --ligand_chain B --id_property ID \
  --runner boltz2 -w ./screen-mol
```

- Each SDF `$$$$` block is one source record; MOL has exactly one record.
- `--id_property` selects the SDF property used as the compound ID. Missing values use
  deterministic `record_000001` identifiers.
- Duplicate IDs fail by default. Use `--duplicate_id_policy suffix` or `source_index`
  only when deterministic rewriting is intended.
- Malformed records become source-mapped failure records and do not hide valid members.
- Coordinates and defined stereochemistry from a structural record are retained for
  ligand preparation; chemistry is also normalized to canonical isomeric SMILES.

## Route 4: Standard and Custom CCD Inputs

Use a CCD code when exact atom names are required, especially for covalent constraints,
or when a backend-compatible component has already been prepared:

```yaml
version: 1
sequences:
  - protein:
      id: A
      sequence: ACDEFGHIK
  - ligand:
      id: B
      ccd: ATP
```

CCD identifiers must be nonempty alphanumeric codes no longer than five characters.
Atom names in constraints must exist in the cached component. A preprocessed CCD is
accepted as-is; COFOLDER does not rederive its named atoms and bonds from SMILES.

The packaged `system_custom_ccd.yaml` uses the custom code `ET5`. Initialize a Boltz2
cache, populate it from the packaged SDF, then run the CCD system with the same cache:

```bash
cofolder-tools setup-boltz2-cache --cache-path ./cache/.boltz
cofolder-tools populate-ccd-cache \
  --sdf ethanol.sdf --property-id ID --cache-path ./cache/.boltz
cofolder validate -s system_custom_ccd.yaml -o options.yaml \
  --runner boltz2 -w ./validate-custom-ccd
```

`setup-boltz2-cache` may download backend assets. `populate-ccd-cache` is explicit and
supports conflict policies; it never silently replaces an existing component.

### Convert SMILES through the public Python API

```python
from pathlib import Path

from rdkit import Chem
from cofolder.modules.entities.ligand import mol_to_ccd

cache = Path("cache/.boltz")
molecule = Chem.MolFromSmiles("CCO")
mol_to_ccd("ET5", molecule, boltz_path=cache)
```

This produces `cache/.boltz/mols/ET5.pkl`, including Boltz-compatible atom names,
connectivity, stereochemical properties, and geometric constraints. Initialize the
rest of the cache first when it will be used for prediction.

## Reusable Conversion Utilities

```python
from cofolder.modules.entities.ligand import (
    csv_to_sdf,
    generate_2d_conformers,
    generate_3d_conformers,
    iterate_sdf_records,
)

csv_to_sdf(
    csv_path="compounds.csv",
    smiles_col="smiles",
    output_sdf_path="compounds.sdf",
    property_cols=["compound_id"],
)
generate_2d_conformers("compounds.sdf", "compounds-2d.sdf")
generate_3d_conformers("compounds.sdf", "compounds-3d.sdf")

for source_index, compound_id, molblock in iterate_sdf_records(
    "compounds-3d.sdf", "compound_id"
):
    print(source_index, compound_id, len(molblock))
```

2D coordinates are intended for layout. 3D generation uses ETKDG followed by UFF and
can fail for molecules that cannot be embedded or optimized. Inspect generated
stereochemistry rather than assuming that 3D generation resolves unspecified centres.

## Failure Guide

| Failure | Resolution |
| --- | --- |
| Invalid SMILES | Parse and sanitize with `validate_smiles`; correct the source chemistry. |
| SDF chemistry mismatch | Use an SDF record for the exact selected system ligand, including stereochemistry. |
| Missing SDF identifier | Set the requested property or accept the deterministic source-record fallback. |
| Duplicate identifiers | Correct the library or select an explicit duplicate-ID policy. |
| Invalid CCD identifier | Use a unique alphanumeric code of at most five characters. |
| Missing custom CCD | Initialize and populate the same cache path configured in `options.yaml`. |
| CCD conflict | Select `overwrite` or `use_cache` deliberately; the installed command defaults to `overwrite`. |

## Not a COFOLDER 1.0 Workflow

Ligand soaking and standalone molecule generation are post-v1 ambitions. External
molecule generators can call Oracle with proposed SMILES, but COFOLDER does not claim a
standalone generator or soaking workflow in this release.

## Related

- [Virtual Screening](screening.md)
- [Configuration Guide](../getting-started/configuration.md)
- [Entities API](../api/modules/entities.md)
