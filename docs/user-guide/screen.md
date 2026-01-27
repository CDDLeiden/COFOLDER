# Screen Command

The `screen` command performs high-throughput virtual screening of compound libraries.

## Basic Usage

```bash
cofolder screen \
  -s system.yaml \
  -b options.yaml \
  -v "sequences,0,ligand,smiles" \
  -c compounds.csv \
  --col_variable smiles \
  --col_id compound_id
```

## Arguments

### Required Arguments

- `-s, --system_path`: Path to system YAML file
- `-b, --boltz_options_path`: Path to Boltz options YAML file
- `-v, --variable`: Path in system YAML to update (comma-separated)

### Input Sources (choose one)

**CSV Input:**
- `-c, --variable_csv`: Path to CSV file
- `--col_variable`: Column containing variables (e.g., SMILES)
- `--col_id`: Column containing compound IDs

**SDF Input:**
- `--variable_sdf`: Path to SDF file
- `--property_id`: Property name for compound ID

### Optional Arguments

- `-w, --wrk_dir`: Working directory
- `--generate_conformers {2D,3D}`: Generate conformers
- `--merge_data`: Columns/properties to merge into output
- `-d, --debug`: Enable debug logging

## Examples

### Screen from CSV

```bash
cofolder screen \
  -s system.yaml \
  -b options.yaml \
  -v "sequences,0,ligand,smiles" \
  -c library.csv \
  --col_variable smiles \
  --col_id compound_id \
  --merge_data "mw,logp,tpsa"
```

### Screen from SDF

```bash
cofolder screen \
  -s system.yaml \
  -b options.yaml \
  -v "sequences,0,ligand,smiles" \
  --variable_sdf library.sdf \
  --property_id ID \
  --generate_conformers 3D
```

## Variable Path

The `-v, --variable` argument specifies the nested path in the system YAML:

```yaml
# For -v "sequences,0,ligand,smiles"
sequences:
  - ligand:        # sequences[0]
      smiles: ...  # will be replaced
      ccd: "LIG"
```

## CSV Format

```csv
compound_id,smiles,mw,logp
CMPD001,CC(C)Cc1ccc(cc1)C(C)C(=O)O,206.28,3.5
CMPD002,CC(=O)OC1=CC=CC=C1C(=O)O,180.16,1.19
```

## Output

Creates a results CSV with:
- Compound IDs
- Prediction metrics
- Merged data from input
- File paths to structures

## Performance Tips

1. Split large libraries into batches
2. Use multiple GPUs with parallel jobs
3. Pre-validate SMILES to avoid failures
4. Use 2D conformers for faster screening
5. Monitor disk space for large libraries

## Related

- [Configuration Guide](../getting-started/configuration.md)
- [Screen API Reference](../api/recipes/screen.md)
- [Virtual Screening Tutorial](../tutorials/screening.md)
