# Predict Command

The `predict` command performs co-folding of a single protein-ligand system.

## Basic Usage

```bash
boltz-lab predict -s system.yaml -b options.yaml
```

## Arguments

### Required Arguments

- `-s, --system_path`: Path to system YAML file
- `-b, --boltz_options_path`: Path to Boltz options YAML file

### Optional Arguments

- `-w, --wrk_dir`: Working directory (default: current directory)
- `--generate_conformers {2D,3D}`: Generate conformers for SMILES input
- `-d, --debug`: Enable debug logging

## Examples

### Basic Prediction

```bash
boltz-lab predict \
  -s examples/system.yaml \
  -b examples/options.yaml \
  -w ./output
```

### With 3D Conformer Generation

```bash
boltz-lab predict \
  -s system.yaml \
  -b options.yaml \
  --generate_conformers 3D
```

### With Debug Logging

```bash
boltz-lab predict \
  -s system.yaml \
  -b options.yaml \
  -d
```

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

- `predictions.cif` - Predicted structure
- `confidence_model_0.json` - Confidence scores
- Log files in the working directory

## Tips

- Use absolute paths in YAML files for reproducibility
- Validate SMILES strings before prediction
- Check GPU memory availability for large systems
- Use debug mode to troubleshoot issues

## Related

- [Configuration Guide](../getting-started/configuration.md)
- [Predict API Reference](../api/recipes/predict.md)
- [Basic Tutorial](../tutorials/basic.md)
