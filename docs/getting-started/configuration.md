# Configuration

## System Configuration

The system YAML file defines the molecular system to be co-folded. It specifies proteins, ligands, and other molecular entities.

### Basic Structure

```yaml
version: 1
sequences:
  - protein:
      id: "protein_1"
      fasta: "SEQUENCE..."
  - ligand:
      smiles: "SMILES_STRING"
      ccd: "RESIDUE_NAME"
```

### Protein Specification

```yaml
protein:
  id: "my_protein"
  fasta: "MKTAYIAKQRQISFV..."
  # or
  pdb: "path/to/protein.pdb"
```

### Ligand Specification

Multiple formats are supported:

#### SMILES Format
```yaml
ligand:
  smiles: "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
  ccd: "IBP"  # Chemical Component Dictionary identifier
```

#### PDB/CIF Format
```yaml
ligand:
  pdb: "path/to/ligand.pdb"
  ccd: "LIG"
```

#### CCD Identifier Only
```yaml
ligand:
  ccd: "ATP"  # Use pre-existing CCD entry
```

### Multiple Chains

```yaml
version: 1
sequences:
  - protein:
      id: "chain_A"
      fasta: "SEQUENCE_A..."
  - protein:
      id: "chain_B"
      fasta: "SEQUENCE_B..."
  - ligand:
      smiles: "SMILES..."
      ccd: "LIG"
```

## Boltz Options

The options YAML file controls Boltz prediction parameters.

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
- `-b, --boltz_options_path`: Path to Boltz options YAML file (required)
- `-d, --debug`: Enable debug logging

### Conformer Generation

```bash
--generate_conformers {2D,3D}
```

Generate conformers for SMILES input:
- `2D`: Generate 2D coordinates
- `3D`: Generate 3D conformers using ETKDG + UFF

### Screen-Specific Options

```bash
-v, --variable "sequences,0,ligand,smiles"
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
- `examples/options.yaml` - Standard Boltz options
