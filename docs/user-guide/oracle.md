# Oracle Command

The `oracle` command runs `validate` once for a **single input ligand**
and returns **one scalar value** selected from validate outputs.

## Basic Usage

```bash
cofolder oracle \
  -s system.yaml \
  -b options.yaml \
  --input_smiles "CCO" \
  --output_metric affinity_pred_value \
  --aggregate first
```

## Required Arguments

- `-s, --system_path`: Path to system YAML file
- `-b, --boltz_options_path`: Path to Boltz options YAML file
- exactly one input:
  - `--input_smiles <smiles>`
  - `--input_mol_file <path>`
- `--output_metric <column_name>`: metric to extract (e.g. `affinity_pred_value`)

## Optional Arguments

- `--aggregate {first,mean,max,min,median}`: how to reduce multiple metric rows
- all common validate arguments are supported and forwarded
- `-w, --wrk_dir`: Working directory
- `-d, --debug`: Enable debug logging

## Scoring Function Checks

Oracle validates that your requested `--output_metric` is compatible
with enabled `--scoring_functions` before running.

Examples:

- `ifp_distance` requires `ifp_distance`
- `sasa_norm_heavy` requires `sasa_normalized`
- `affinity_pred_value` requires affinity metrics
- `bias_*` requires `--assess_bias`

## Output

- Scalar return value from `Oracle.run()`
- CSV file: `<wrk_dir>/oracle_result.csv`

## Related

- [Oracle API Reference](../api/recipes/oracle.md)
- [Validate Command](validate.md)
