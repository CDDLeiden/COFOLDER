# Validate Command

The `validate` command performs co-folding of a single protein-ligand system.

## Basic Usage

```bash
cofolder validate -s system.yaml -o options.yaml
```

## Arguments

### Required Arguments

- `-s, --system_path`: Path to system YAML file
- `-o, --options_path`: Path to runner options YAML file

### Optional Arguments

- `-w, --wrk_dir`: Working directory (default: current directory)
- `--conformers {2D,3D,sdf}`: Generate or load conformers for SMILES input
- `--sdf_file`: SDF path (required with `--conformers sdf`)
- `--assess_bias`: Add shared bias diagnostics under `results/bias_train/`
- `-d, --debug`: Enable debug logging

## Examples

### Basic Validation

```bash
cofolder validate \
  -s examples/system.yaml \
  -o examples/options.yaml \
  -w ./output
```

### With 3D Conformer Generation

```bash
cofolder validate \
  -s system.yaml \
  -o options.yaml \
  --conformers 3D
```

### With Debug Logging

```bash
cofolder validate \
  -s system.yaml \
  -o options.yaml \
  -d
```

When `--debug` is enabled, COFOLDER emits benchmark-friendly timing lines in the log:

- `TIMER | <label> | <seconds>s` during the run
- `TIMER SUMMARY | <label> | <seconds>s` at the end

This is intended for component benchmarking. For reproducible benchmark numbers, use
warmed Boltz cache/model files and fresh prediction outputs so one-time downloads or
cache-hit skips do not distort timings.

### With Shared Bias Diagnostics

```bash
cofolder validate \
  -s system.yaml \
  -o options.yaml \
  --assess_bias \
  --protein_training_data_path protein_training_data.csv \
  --ligand_training_data_path ligand_training_data.csv
```

This uses the same bias-reporting seam as the standalone [`bias`](bias.md) workflow and
writes those diagnostics to `<wrk_dir>/results/bias_train/` alongside the runner-backed
prediction outputs.

## System File Format

```yaml
version: 1
sequences:
  - protein:
      id: "protein_1"
      fasta: "MKTAYIAKQRQISFV..."
  - ligand:
      smiles: "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
      ccd: "IBP"
```

## Output

The command creates:

- `raw/` - runner-owned raw execution artifacts
- `results/records.jsonl` - authoritative success, metric, and failure records
- `results/successes.csv`, `results/metrics.csv`, and `results/failures.csv` - tabular views
- `results/manifest.json` - schema, invocation, evidence, and artifact metadata
- `results/structures/` - gathered output structures
- `results/bias_train/` - optional shared bias diagnostics when `--assess_bias` is enabled
- log files in the working directory

Bias outputs remain distinct from validation metrics, model-confidence outputs, and any
affinity predictions. They are reference-overlap diagnostics, not a proxy for binding
affinity or expected cofolding success.

## Tips

- Use absolute paths in YAML files for reproducibility
- Validate SMILES strings before prediction
- Check GPU memory availability for large systems
- Use debug mode to troubleshoot issues

## Related

- [Bias Command](bias.md)
- [Configuration Guide](../getting-started/configuration.md)
- [Validate API Reference](../api/recipes/validate.md)
- [Basic Tutorial](../tutorials/basic.md)
