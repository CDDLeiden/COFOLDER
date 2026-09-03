# Configuration

## System Configuration

The system YAML file defines proteins, ligands, DNA, RNA, and optional constraints.

### Basic Structure

```yaml
version: 1
sequences:
  - protein:
      id: "protein_1"
      sequence: "SEQUENCE..."
  - ligand:
      id: "ligand_1"
      smiles: "SMILES_STRING"
```

### Protein Specification

```yaml
protein:
  id: "my_protein"
  sequence: "MKTAYIAKQRQISFV..."
```

### Ligand Specification

Multiple formats are supported:

#### SMILES Format
```yaml
ligand:
  id: L
  smiles: "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
```

#### CCD Identifier Only
```yaml
ligand:
  id: L
  ccd: "ATP"  # Use pre-existing CCD entry
```

Specify either `smiles` or `ccd`, not both. Conformer preparation replaces a
SMILES representation with a generated CCD before validation and execution.

### DNA and RNA

```yaml
- dna:
    id: D
    sequence: ATGC
- rna:
    id: R
    sequence: AUGC
```

Chain IDs must be non-empty and unique across every entity. A list-valued `id`
creates identical copies in that exact order.

### Runner input support

| Runner | Protein | Ligand | DNA | RNA | Bond | Pocket | Contact |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `boltz1` | yes | yes | yes | yes | yes | one, distance 6 | no |
| `boltz2` | yes | yes | yes | yes | yes | yes | yes |
| `boltz-community` | yes | yes | yes | yes | yes | yes | yes |
| `openfold3` | yes | yes | yes | yes | no | one | no |

Unsupported combinations fail during COFOLDER preflight, before the backend is
started. OpenFold3 bond inputs are rejected because OpenFold3 0.5.0 accepts a
`covalent_bonds` field but does not consume it during query construction.

### Constraints

Constraints are top-level entries. Residue positions are 1-based.

```yaml
constraints:
  - bond:
      atom1: [A, 145, SG]
      atom2: [L, 1, C12]
  - pocket:
      binder: L
      contacts: [[A, 140], [A, 145]]
      max_distance: 6.0
  - contact:
      token1: [A, 140]
      token2: [D, 2]
      max_distance: 6.0
```

Bond endpoints are `[chain_id, residue_id, atom_name]`. Pocket/contact tokens
use a residue number for polymer chains and an atom name for ligand chains.
Referenced chains, residues, and atoms are checked against the prepared system.
OpenFold3 translates one pocket constraint to its native query schema and does
not accept `force` on that constraint.

`screen` preserves constraints not involving a replaced field. If replacement
or SMILES/SDF-to-CCD conversion removes a referenced ligand atom, that row fails
with the invalid chain/residue/atom reference instead of submitting stale input.

### Multiple Chains

```yaml
version: 1
sequences:
  - protein:
      id: "chain_A"
      sequence: "SEQUENCE_A..."
  - protein:
      id: "chain_B"
      sequence: "SEQUENCE_B..."
  - ligand:
      id: L
      smiles: "SMILES..."
```

## Runner Options

The options YAML file controls runner prediction parameters.

### Basic Options

```yaml
out_dir: output          # Output directory
devices: [0]            # GPU devices to use
num_models: 1           # Number of models to run
recycling_steps: 3      # Number of recycling iterations
diffusion_samples: 1    # Number of diffusion samples
```

### Advanced Options

```yaml
out_dir: output
devices: [0, 1]         # Multi-GPU support
num_models: 5           # Ensemble prediction
recycling_steps: 5      # More recycling for accuracy
diffusion_samples: 10   # Multiple samples for diversity
sampling_steps: 200     # Diffusion sampling steps
diffusion_temperature: 1.0  # Temperature for sampling
```

## Command-Line Options

### Common Arguments

- `-w, --wrk_dir`: Working directory (default: current directory)
- `-s, --system_path`: Path to system YAML file (required)
- `-o, --options_path`: Path to runner options YAML file (required)
- `-d, --debug`: Enable debug logging

### Conformer Generation

```bash
--conformers {2D,3D,sdf}
```

Generate conformers for SMILES input:
- `2D`: Generate 2D coordinates
- `3D`: Generate 3D conformers using ETKDG + UFF
- `sdf`: Reuse conformers from an input SDF file (requires `--sdf_file`)

### Screen-Specific Options

```bash
-v, --variable "sequences,1,ligand,smiles"
-c, --variable_csv compounds.csv
--col_variable smiles
--col_id compound_id
```

### Validate-Specific Options

```bash
--repeats 5                    # Number of validation repeats
--seed 42                      # Base seed for reproducibility
--conformers {2D,3D,sdf}       # Optional conformer handling
--sdf_file ligands.sdf         # SDF input when using --conformers sdf
```

## Environment Variables

### Boltz Cache

By default, Boltz stores CCD files in `~/.boltz/`. You can override this:

```bash
export BOLTZ_CACHE=/path/to/cache
```

### CUDA Configuration

```bash
export CUDA_VISIBLE_DEVICES=0,1  # Specify GPU devices
```

## Best Practices

1. **Use absolute paths**: Avoid relative paths in YAML files
2. **Validate SMILES**: Ensure SMILES strings are valid before screening
3. **CCD naming**: Keep CCD identifiers to 5 characters or less
4. **GPU memory**: Monitor GPU memory usage for large systems
5. **Reproducibility**: Use fixed seeds for reproducible results

## Example Configurations

See the `examples/` directory in the repository for complete configuration examples:

- `examples/system.yaml` - Basic protein-ligand system
- `examples/system_screen.yaml` - Screening configuration
- `examples/system_covalent.yaml` - Covalent binding system
- `examples/system_nucleic_acid.yaml` - DNA/RNA system with a pocket constraint
- `examples/options.yaml` - Standard Boltz options
