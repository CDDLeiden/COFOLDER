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
- `--ifp_filter_threshold FLOAT`: annotate compounds using inclusive reference
  overlap (`0` to `1`). Omit this option to disable filtering.
- `--ifp_ligand_chain CHAIN`: ligand chain to evaluate. A single ligand is
  selected automatically; this option is required when multiple ligand chains
  are present.

Reference-overlap filtering requires `ifp_distance` scoring and
`--pocket_coverage_reference`. The reference can be a binary bitstring, residue
numbers or labels, or a file containing one of those forms:

```bash
cofolder screen \
  -s system.yaml -o options.yaml -c compounds.csv \
  --col_id compound_id \
  --variable sequences,1,ligand,smiles --col_variable smiles \
  --scoring_functions ifp_distance \
  --pocket_coverage_reference "A25 G48 Y51" \
  --ifp_filter_threshold 0.6
```

Overlap is the fraction of active reference-pocket bits also present in the
predicted distance IFP. A row passes when overlap is greater than or equal to
the threshold. Filtering only annotates results; it never deletes rejected
rows or their prediction artifacts.

If you only need pre-cofolding bias diagnostics for one system, prefer the dedicated
[`bias`](bias.md) workflow instead of running `screen`.

## Protein MSA Reuse

For Boltz-family runners, `screen` resolves each missing protein MSA once and stores
the reusable artifact under `<wrk_dir>/shared/msa/<runner>/`. The resolved MSA is injected into
later repeats and every subsequent ligand-specific system YAML, so the MSA server is
not called again for the fixed protein system.

If the original system YAML already supplies `msa` for a protein, that file is used
directly and the server is not called for that protein. Relative MSA paths are resolved
relative to the original system YAML before row-specific YAML files are written. A
system may mix supplied and missing MSAs; only missing protein MSAs are generated.

The shared cache is matched by protein sequence and includes a manifest. A missing,
corrupt, or sequence-mismatched cached artifact is not silently injected.

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
- filter audit columns: `ifp_filter_pass`, `ifp_filter_status`,
  `ifp_filter_reason`, `ifp_filter_overlap`, `ifp_filter_threshold`, and
  `ifp_filter_reference`

`screen_results_with_scores.csv` includes:

- all original input CSV columns
- run metadata (`index`, `status`, `error_message`, `run_dir`)
- extracted system-level score columns prefixed as `system__...`
- extracted chain-level score columns prefixed as `<entity>_<chain_id>__...`
- the same six filter audit columns as the summary CSV

Filter status is `accepted`, `rejected`, `not_evaluable`, or `not_applied`.
Missing or malformed IFPs, all-zero references, and unequal vector lengths are
`not_evaluable`; vectors are never truncated for filtering.

The consolidated schema always includes the manuscript-facing confidence, affinity,
SASA, distance-IFP, training-set-proximity, and pocket-coverage columns. Unsupported or
unconfigured metrics are represented by empty values rather than omitted columns.
Pairwise confidence retains the runner's numeric-chain column and also exposes a
chain-ID alias such as `ligand_B__pair_chains_iptm_A`.

When used from Python, `Screen.run()` returns the same merged results as a
`pandas.DataFrame` after writing the CSV files.

## Failure Behavior

If one row fails, screening continues for remaining rows.
Failures are recorded in `screen_results.csv`.

## Related

- [Bias Command](bias.md)
- [Configuration Guide](../getting-started/configuration.md)
- [Virtual Screening Tutorial](../tutorials/screening.md)
- [Validate Command](validate.md)
