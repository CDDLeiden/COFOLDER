# Quick Start Guide

Copy the packaged inputs and inspect a complete Screen plan without running a
backend or service:

```bash
cofolder-tools copy-examples ./cofolder-example
cofolder screen --preflight_only \
  --system_path ./cofolder-example/system_screen.yaml \
  --options_path ./cofolder-example/options.yaml \
  --library ./cofolder-example/ligand_screen.csv \
  --col_id Name --smiles_column SMILES
```

The example contains one ligand entity, so `--ligand_chain` is inferred.

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
      id: B
      smiles: "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
```

### 2. Configure Runner Options

Create a runner options file (`options.yaml`):

```yaml
version: 1
runtime:
  cache_path: ./cache/.boltz
  diffusion_samples: 1
runner:
  devices: 1
  recycling_steps: 3
```

### 3. Run Validation

Execute the validation workflow:

```bash
cofolder validate -s system.yaml -o options.yaml -w ./output
```

## Output Files

After successful execution, you'll find:

- `raw/` - runner-owned raw execution artifacts
- `results/records.jsonl` - authoritative versioned records
- `results/metrics.csv` - long-form metric view
- `results/successes.csv` and `results/failures.csv` - outcome views
- `results/manifest.json` - invocation manifest
- `results/structures/` - gathered output structures
- logs and additional metadata in the working directory

## Common Commands

### Validate a Single System

```bash
cofolder validate -s system.yaml -o options.yaml
```

### Screen a Library

```bash
cofolder screen -s system.yaml -o options.yaml -c compounds.csv --col_id compound_id --ligand_chain B --smiles_column smiles
```

### Use as Oracle

```bash
cofolder oracle -s system.yaml -o options.yaml \
  --input_smiles "CCO" --output_metric system__confidence_score \
  --scoring_functions confidence_metrics
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
