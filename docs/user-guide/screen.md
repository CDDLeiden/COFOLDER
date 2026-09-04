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
- `--cluster_ifps`: cluster evaluable binary distance IFPs after all rows finish.
- `--ifp_cluster_similarity_threshold FLOAT`: inclusive Jaccard-similarity cut
  for clustering (`0` to `1`, default `0.5`).

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

## IFP Clustering Without a Reference

Clustering does not require a reference structure or pocket definition. It is useful
for discovering recurring predicted contact patterns, including a reproducible
all-zero “no contacts” pattern:

```bash
cofolder screen \
  -s system.yaml -o options.yaml -c compounds.csv \
  --col_id compound_id \
  --variable sequences,1,ligand,smiles --col_variable smiles \
  --scoring_functions ifp_distance \
  --cluster_ifps \
  --ifp_cluster_similarity_threshold 0.5
```

COFOLDER uses average-linkage agglomerative clustering over Jaccard distance and
cuts the tree at `distance <= 1 - similarity_threshold`. Stable IDs (`IFP001`,
`IFP002`, …) are assigned by each cluster's earliest input row. Missing, malformed,
or vector-length-incompatible IFPs receive `ifp_cluster_status=not_evaluable` and
no cluster ID. The first valid fingerprint defines the expected vector length.

## Full Diagnostic Screen

This command requests the structural, affinity, bias, distance-IFP, and pocket
diagnostics used for publication-facing review. Replace the reference and training
data paths with datasets appropriate for the target:

```bash
cofolder screen \
  -s system.yaml -o options.yaml -c compounds.csv \
  --col_id compound_id \
  --variable sequences,1,ligand,smiles --col_variable smiles \
  --repeats 3 \
  --scoring_functions confidence_metrics affinity_metrics affinity_metrics_ext \
    sasa sasa_normalized ifp_distance \
  --assess_bias \
  --protein_training_data_path protein_training_data.csv \
  --ligand_training_data_path ligand_training_data.csv \
  --reference_path reference_complex.cif \
  --reproduction_metrics protein_rmsd ligand_rmsd sucos pocket_coverage \
  --pocket_coverage_reference "A25 G48 Y51" \
  --ifp_filter_threshold 0.6 \
  --cluster_ifps
```

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
- Cluster summary (when enabled): `<wrk_dir>/ifp_cluster_summary.csv`

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
- clustering columns: `ifp_cluster_id` and `ifp_cluster_status`

`screen_results_with_scores.csv` includes:

- all original input CSV columns
- run metadata (`index`, `status`, `error_message`, `run_dir`)
- extracted system-level score columns prefixed as `system__...`
- extracted chain-level score columns prefixed as `<entity>_<chain_id>__...`
- the same six filter audit columns as the summary CSV
- the same two clustering columns as the summary CSV

Filter status is `accepted`, `rejected`, `not_evaluable`, or `not_applied`.
Missing or malformed IFPs, all-zero references, and unequal vector lengths are
`not_evaluable`; vectors are never truncated for filtering.

The consolidated schema always includes the manuscript-facing confidence, affinity,
SASA, distance-IFP, training-set-proximity, and pocket-coverage columns. Unsupported or
unconfigured metrics are represented by empty values rather than omitted columns.
Pairwise confidence retains the runner's numeric-chain column and also exposes a
chain-ID alias such as `ligand_B__pair_chains_iptm_A`.

### Stable manuscript-facing score columns

The merged output always declares the following columns for configured system chain
IDs, using empty values when a backend or run does not provide a metric:

| Scope | Stable columns |
| --- | --- |
| System | `system__confidence_score`, `system__ptm`, `system__iptm` |
| Protein proximity | `system__bias_prot_sim_train_max`, `system__bias_prot_sim_train_pairwise_max`, `protein_<ID>__bias_prot_sim_train`, `protein_<ID>__bias_prot_sim_train_pairwise` |
| Ligand proximity | `system__bias_lig_sim_train_max`, `ligand_<ID>__bias_lig_sim_train` |
| Affinity/binding | `ligand_<ID>__affinity_pred_value`, `ligand_<ID>__affinity_probability_binary`, `ligand_<ID>__pIC50`, `ligand_<ID>__IC50_M`, `ligand_<ID>__pIC50_kcal_per_mol` |
| Chain confidence | `protein_<ID>__chains_ptm`, `ligand_<ID>__chains_ptm`, numeric `pair_chains_iptm_<index>` fields when emitted, and semantic `ligand_<ID>__pair_chains_iptm_<protein_ID>` aliases |
| Exposure and contacts | `ligand_<ID>__sasa`, `ligand_<ID>__sasa_norm_heavy`, `ligand_<ID>__ifp_distance`, `ligand_<ID>__ifp_prolif` |
| Pocket/reference | `system__pocket_coverage_ref`, `system__pocket_coverage_ref_mean`, `system__pocket_coverage_custom`, `system__pocket_coverage_custom_mean`, `ligand_<ID>__pocket_coverage_ref`, `ligand_<ID>__pocket_coverage_custom`, `ligand_<ID>__ligand_pose_overlap_ref` |

`ifp_cluster_status` is `clustered`, `not_evaluable`, or `not_applied`.
`ifp_cluster_summary.csv` contains `ifp_cluster_id`, `size`, JSON `member_ids`,
`medoid_compound_id`, JSON `consensus_ifp`, and
`mean_within_cluster_jaccard_similarity`. Consensus bits are present in at least
half of cluster members; medoid ties resolve to the earliest input row. A singleton
cluster has mean within-cluster similarity `1.0`.

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
