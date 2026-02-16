# Quick Start Guide

This guide will walk you through your first COFOLDER workflow.

## Basic Workflow

### 1. Prepare Your System

Create a system YAML file (`system.yaml`) defining your protein-ligand system:

```yaml
version: 1
sequences:
  - protein:
      id: "protein_1"
      fasta: "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL"
  - ligand:
      smiles: "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
      ccd: "IBP"
```

### 2. Configure Boltz Options

Create a Boltz options file (`options.yaml`):

```yaml
out_dir: output
devices: [0]
num_models: 1
recycling_steps: 3
diffusion_samples: 1
```

### 3. Run Validation

Execute the validation workflow:

```bash
cofolder validate -s system.yaml -b options.yaml -w ./output
```

## Output Files

After successful execution, you'll find:

- `predictions.cif` - Predicted structure in CIF format
- `confidence_model_0.json` - Confidence metrics for the prediction
- Logs and additional metadata

## Common Commands

### Validate a Single System

```bash
cofolder validate -s system.yaml -b options.yaml
```

### Screen a Library

```bash
cofolder screen -s system.yaml -b options.yaml -c compounds.csv --col_id compound_id --variable sequences,0,ligand,smiles --col_variable smiles
```

### Use as Oracle

```bash
cofolder oracle -s system.yaml -b options.yaml
```

## Getting Help

For detailed information about any command:

```bash
cofolder validate --help
cofolder screen --help
cofolder oracle --help
```

## Next Steps

- Explore the [User Guide](../user-guide/overview.md) for detailed command documentation
- Check out [Tutorials](../tutorials/basic.md) for step-by-step examples
- Review [Configuration](configuration.md) for advanced options
