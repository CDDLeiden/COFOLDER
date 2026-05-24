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
- `results/system_metrics.csv` - merged system-level metrics
- `results/chain_metrics.csv` - merged chain-level metrics
- `results/structures/` - gathered output structures
- log files in the working directory

## Tips

- Use absolute paths in YAML files for reproducibility
- Validate SMILES strings before prediction
- Check GPU memory availability for large systems
- Use debug mode to troubleshoot issues

## Related

- [Configuration Guide](../getting-started/configuration.md)
- [Validate API Reference](../api/recipes/validate.md)
- [Basic Tutorial](../tutorials/basic.md)
