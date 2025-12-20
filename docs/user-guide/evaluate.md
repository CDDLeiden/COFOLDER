# Evaluate Command

The `evaluate` command performs comprehensive validation of predictions against reference structures.

## Basic Usage

```bash
boltz-lab evaluate \
  -s system.yaml \
  -b options.yaml \
  -i reference.pdb \
  --repeats 5
```

## Arguments

### Required Arguments

- `-s, --system_path`: Path to system YAML file
- `-b, --boltz_options_path`: Path to Boltz options YAML file

### Optional Arguments

- `-w, --wrk_dir`: Working directory
- `--repeats N`: Number of prediction repeats (default: 1)
- `--seeds`: Comma-separated list of random seeds
- `-i, --input_pdb`: Reference structure for RMSD calculation
- `--ifp`: Interaction fingerprint specification
- `--generate_conformers {2D,3D}`: Generate conformers
- `--log_name`: Base name for log file (default: "log")
- `-d, --debug`: Enable debug logging

## Examples

### Basic Evaluation

```bash
boltz-lab evaluate \
  -s system.yaml \
  -b options.yaml \
  -i reference.pdb \
  --repeats 5
```

### With Fixed Seeds

```bash
boltz-lab evaluate \
  -s system.yaml \
  -b options.yaml \
  -i reference.pdb \
  --repeats 3 \
  --seeds "42,123,456"
```

### With Interaction Fingerprints

```bash
boltz-lab evaluate \
  -s system.yaml \
  -b options.yaml \
  -i reference.pdb \
  --ifp "true"
```

### Manual IFP Specification

```bash
boltz-lab evaluate \
  -s system.yaml \
  -b options.yaml \
  --ifp '{"A:123":"ARG", "B:45":"TYR"}'
```

## Metrics

The evaluate command calculates:

### RMSD Metrics
- Protein RMSD
- Ligand RMSD
- Overall system RMSD

### Interaction Fingerprints
- Interaction overlap with reference
- Residue-level contact analysis
- Binding mode similarity

### Confidence Scores
- Model confidence
- Per-residue confidence
- Ligand confidence

## IFP Specification

Three modes for interaction fingerprints:

1. **Auto-extract from crystal structure:**
   ```bash
   --ifp "true"  # Requires --input_pdb
   ```

2. **Disable IFP calculation:**
   ```bash
   --ifp "false"
   ```

3. **Manual specification:**
   ```bash
   --ifp '{"A:123":"ARG", "B:45":"TYR"}'
   ```

## Output

Creates comprehensive results including:

- Individual prediction structures
- RMSD statistics (mean, std, min, max)
- IFP overlap scores
- Confidence metrics
- Detailed logs

## Reproducibility

For reproducible results:

```bash
boltz-lab evaluate \
  -s system.yaml \
  -b options.yaml \
  --repeats 10 \
  --seeds "1,2,3,4,5,6,7,8,9,10"
```

## Performance Considerations

- More repeats = better statistics but longer runtime
- Use parallel jobs for independent evaluations
- Consider GPU memory for large systems

## Tips

- Use at least 5 repeats for reliable statistics
- Compare against experimental structures when available
- Document random seeds for reproducibility
- Check RMSD distributions for multimodal binding

## Related

- [Evaluate API Reference](../api/recipes/evaluate.md)
- [Configuration Guide](../getting-started/configuration.md)
