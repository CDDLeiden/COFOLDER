# Screen Command

The `screen` command runs **`validate` once per row** in a CSV file.
For each row, it adapts one or more fields in `system.yaml`, runs validation,
and writes a screening summary file.

## Basic Usage

```bash
cofolder screen \
  -s system.yaml \
  -o options.yaml \
  --variable_csv compounds.csv \
  --col_id compound_id \
  --variable sequences,1,ligand,smiles --col_variable smiles
```

## Multi-Variable Mapping

Use repeated `--variable` / `--col_variable` pairs in the same order:

```bash
cofolder screen \
  -s system.yaml \
  -o options.yaml \
  --variable_csv compounds.csv \
  --col_id compound_id \
  --variable sequences,1,ligand,smiles --col_variable smiles \
  --variable sequences,1,ligand,ccd --col_variable ccd
```

## Required Arguments

- `-s, --system_path`: Path to system YAML file
- `-o, --options_path`: Path to runner options YAML file
- `-c, --variable_csv`: Path to CSV file
- `--col_id`: Column containing row IDs
- `--variable`: Repeatable YAML path to update (comma-separated)
- `--col_variable`: Repeatable CSV column mapped to each `--variable`

## Validation Rules

1. At least one `--variable/--col_variable` pair is required.
2. Number of `--variable` entries must equal number of `--col_variable` entries.
3. `--col_id` and each `--col_variable` must exist in the CSV.
4. Empty values in mapped columns fail that row (screening continues).

## Optional Arguments

- `-w, --wrk_dir`: Working directory
- `--merge_data`: Comma-separated metadata columns copied into summary
- All common validate options are supported and forwarded, including:
  - scoring functions
  - bias options (`--assess_bias`, `--bias_*`)
  - reproduction options (`--reference_path`, `--reproduction_metrics`)
  - robustness options

If you only need pre-cofolding bias diagnostics for one system, prefer the dedicated
[`bias`](bias.md) workflow instead of running `screen`.

## Output

`screen` writes:

- Per-row validate outputs under `<wrk_dir>/<index>_<id>/...`
- Summary CSV: `<wrk_dir>/screen_results.csv`
- Merged input+scores CSV: `<wrk_dir>/screen_results_with_scores.csv`

Summary columns include:

- `index`
- `<col_id>`
- `status` (`success` / `failed`)
- `error_message`
- `run_dir`
- mapped `col_variable` values
- optional `merge_data` values

`screen_results_with_scores.csv` includes:

- all original input CSV columns
- run metadata (`index`, `status`, `error_message`, `run_dir`)
- extracted system-level score columns prefixed as `system__...`
- extracted chain-level score columns prefixed as `<entity>_<chain_id>__...`

## Failure Behavior

If one row fails, screening continues for remaining rows.
Failures are recorded in `screen_results.csv`.

## Related

- [Bias Command](bias.md)
- [Configuration Guide](../getting-started/configuration.md)
- [Virtual Screening Tutorial](../tutorials/screening.md)
- [Validate Command](validate.md)
