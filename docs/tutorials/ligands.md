# Ligand Handling Tutorial

Comprehensive guide to working with ligands in Boltz-Lab.

## Overview

Boltz-Lab supports multiple ligand input formats and provides utilities for ligand preparation and manipulation.

## Ligand Input Formats

### SMILES Strings

The simplest format for small molecules:

```yaml
ligand:
  smiles: "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
  ccd: "IBP"
```

### SDF Files

For molecules with 3D coordinates:

```yaml
ligand:
  sdf: "path/to/ligand.sdf"
  ccd: "LIG"
```

### PDB/CIF Files

For ligands from crystal structures:

```yaml
ligand:
  pdb: "path/to/ligand.pdb"
  ccd: "ATP"
```

### CCD Identifiers

For standard residues in the Chemical Component Dictionary:

```yaml
ligand:
  ccd: "ATP"  # Adenosine triphosphate
```

## Conformer Generation

### 2D Conformers

Fast generation of 2D coordinates:

```bash
boltz-lab predict \
  -s system.yaml \
  -b options.yaml \
  --generate_conformers 2D
```

### 3D Conformers

Generate 3D conformers using ETKDG + UFF:

```bash
boltz-lab predict \
  -s system.yaml \
  -b options.yaml \
  --generate_conformers 3D
```

## Programmatic Ligand Handling

### Convert CSV to SDF

```python
from boltz_lab.modules.entities.ligand import csv_to_sdf

csv_to_sdf(
    csv_path="compounds.csv",
    smiles_col="smiles",
    output_sdf_path="compounds.sdf",
    property_cols=["compound_id", "mw", "logp"]
)
```

### Generate 2D Conformers

```python
from boltz_lab.modules.entities.ligand import generate_2d_conformers

generate_2d_conformers(
    sdf_in="input.sdf",
    sdf_out="output_2d.sdf"
)
```

### Generate 3D Conformers

```python
from boltz_lab.modules.entities.ligand import generate_3d_conformers

generate_3d_conformers(
    sdf_in="input.sdf",
    sdf_out="output_3d.sdf"
)
```

### Create CCD Entries

Convert molecules to Boltz CCD format:

```python
from rdkit import Chem
from boltz_lab.modules.entities.ligand import mol_to_ccd

mol = Chem.MolFromSmiles("CCO")
mol_to_ccd("ETH", mol, boltz_path="~/.boltz")
```

### Cache Molecules from SDF

```python
from boltz_lab.modules.entities.ligand import cache_mols_from_sdf

cache_mols_from_sdf(
    file_path="ligands.sdf",
    property_id="ID",
    on_conflict="overwrite",
    cache="~/.boltz/"
)
```

## Working with SDF Files

### Iterate Through SDF

```python
from boltz_lab.modules.entities.ligand import iterate_sdf_records

for idx, mol_id, molblock in iterate_sdf_records("library.sdf", "ID"):
    print(f"Processing {mol_id}...")
    # Do something with molblock
```

## Best Practices

### CCD Naming

- Keep identifiers ≤ 5 characters
- Use alphanumeric characters only
- Avoid conflicts with standard CCD codes
- Use descriptive names when possible

```python
from boltz_lab.modules.entities.ligand import sanitize_mol_id

safe_id = sanitize_mol_id("VERY_LONG_COMPOUND_ID")  # Truncates to 5 chars
```

### SMILES Validation

Always validate SMILES before processing:

```python
from rdkit import Chem

smiles = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
mol = Chem.MolFromSmiles(smiles)

if mol is None:
    print("Invalid SMILES!")
else:
    print("Valid SMILES")
```

### Stereochemistry

Ensure stereochemistry is properly defined:

```python
from rdkit import Chem
from rdkit.Chem import AllChem

mol = Chem.MolFromSmiles("C[C@H](O)c1ccccc1")
AllChem.AssignStereochemistryFrom3D(mol)
```

## Advanced Features

### Custom CCD Properties

Add custom properties to molecules:

```python
from boltz_lab.modules.entities.ligand import add_pickled_prop

mol = Chem.MolFromSmiles("CCO")
add_pickled_prop(mol, "custom_score", 0.95)
```

### Extract Constraints

Extract geometric constraints from molecules:

```python
from boltz_lab.modules.entities.ligand import extract_constraints

constraints = [...] # From parsed residue
bond_constraints = extract_constraints(constraints, 'is_bond')
```

## Common Issues

### Invalid SMILES

**Problem:** SMILES string cannot be parsed

**Solution:**
- Validate with RDKit before use
- Check for typos or invalid characters
- Use canonical SMILES

### CCD Conflicts

**Problem:** CCD identifier already exists

**Solution:**
```python
cache_mols_from_sdf(
    file_path="ligands.sdf",
    property_id="ID",
    on_conflict="overwrite"  # or "use_cache"
)
```

### 3D Generation Failures

**Problem:** 3D conformer generation fails

**Solution:**
- Try 2D conformers instead
- Check molecule for unusual features
- Simplify structure if possible

## Examples

### Complete Workflow

```python
from rdkit import Chem
from boltz_lab.modules.entities.ligand import (
    csv_to_sdf,
    generate_3d_conformers,
    cache_mols_from_sdf
)

# 1. Convert CSV to SDF
csv_to_sdf(
    csv_path="compounds.csv",
    smiles_col="smiles",
    property_cols=["id", "name"]
)

# 2. Generate 3D conformers
generate_3d_conformers(
    sdf_in="compounds.sdf",
    sdf_out="compounds_3d.sdf"
)

# 3. Cache to Boltz
cache_mols_from_sdf(
    file_path="compounds_3d.sdf",
    property_id="id"
)
```

## Related

- [API Reference: Entities](../api/modules/entities.md)
- [Virtual Screening Tutorial](screening.md)
- [Configuration Guide](../getting-started/configuration.md)
