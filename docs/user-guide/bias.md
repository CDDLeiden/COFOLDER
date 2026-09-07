# Bias Command

The `bias` command assesses protein and ligand reference overlap directly from
`system.yaml`. It is a standalone diagnostic workflow for pre-cofolding decision
support and does not require a backend runner or `options.yaml`.

## Basic Usage

### Public Reference CSVs

```bash
cofolder bias \
  --system_path system.yaml \
  --wrk_dir ./bias_output \
  --protein_training_data_path protein_training_data.csv \
  --ligand_training_data_path ligand_training_data.csv
```

### Custom-Only References

```bash
cofolder bias \
  --system_path system.yaml \
  --wrk_dir ./bias_output \
  --custom_protein_reference_path custom_protein.csv \
  --custom_ligand_reference_path custom_ligand.csv
```

### Build Public Training Data First

```bash
cofolder bias \
  --system_path system.yaml \
  --wrk_dir ./bias_output \
  --build_bias_training_data \
  --protein_training_data_path ./training/protein_training_data.csv
```

Build mode is the path that may require `mmseqs2`, depending on how you prepare your
public protein reference data. It also requires a readable `components.cif`, either at
`<protein_training_data_path parent>/ccd/components.cif` or via
`--bias_training_components_cif`.

## What It Requires

- `--system_path`: required
- One of these reference modes:
  - public protein and ligand training data paths
  - custom protein and custom ligand reference paths
  - a public-plus-custom combination
  - build mode with `--build_bias_training_data` and `--protein_training_data_path`

The command validates these inputs, including build-mode `components.cif` resolution,
before recipe execution. It does not look up a runner and does not require
`--options_path`.

## Output Layout

Successful runs publish the public contract under `<wrk_dir>/results/` and supporting
bias artifacts under `<wrk_dir>/results/bias_train/`:

```text
results/
├── records.jsonl
├── successes.csv
├── metrics.csv
├── failures.csv
├── manifest.json
└── bias_train/
├── protein_training_data.csv
├── ligand_training_data_<CHAIN>.csv
├── bias_training_data.csv
├── reference_landscape_summary.csv
├── bias_reference_overlap_scatter.png
└── bias_reference_overlap_scatter.pdf
```

These files are shared with `validate --assess_bias`, which writes the same bias artifact
layout under the same `results/bias_train/` directory.

Artifact notes:

- Ligand reference views are usually written per chain as `ligand_training_data_<CHAIN>.csv`.
- Some ligand-only compatibility paths still write `ligand_training_data.csv`.
- The scatter plot is written only when the merged plotting dataset has both protein and
  ligand similarity axes. Otherwise the workflow writes
  `bias_reference_overlap_scatter.skipped.txt` instead of `.png`/`.pdf`.

## Key Arguments

- `--bias_release_cutoff YYYY-MM-DD`: filters public reference rows by release date
- `--bias_ligand_similarity_threshold 0.0-1.0`: threshold used when building public ligand training data
- `--bias_chains A B`: restricts analysis to selected chains
- `--bias_training_components_cif PATH`: explicit `components.cif` for build mode
- `--custom_protein_reference_path PATH`: custom protein CSV
- `--custom_ligand_reference_path PATH`: custom ligand CSV or SDF

## Interpretation

Bias outputs are reference-overlap diagnostics:

- `metrics.csv` reports long-form bias-related overlap summaries
- `protein_training_data.csv` and `ligand_training_data_<CHAIN>.csv` show the reference
  rows used for each query chain
- `bias_training_data.csv` is the merged plotting dataset behind the bias scatter plot
- `reference_landscape_summary.csv` compares nearest overall, public, and custom references

Protein similarity fields are deliberately method-specific:

- `sequence_similarity` and `plot_sequence_similarity` contain only MMseqs `pident`
  (percent and 0-1 normalized, respectively).
- `sequence_similarity_pairwise` and `plot_sequence_similarity_pairwise` contain only
  Biopython PairwiseAligner identity (percent and 0-1 normalized, respectively).
- `sequence_similarity_method` is `mmseqs_pident`, `pairwise_aligner`, or `unavailable`.
  A protein row cannot contain both score types.
- `bias_prot_sim_train_max` summarizes only MMseqs values, while
  `bias_prot_sim_train_pairwise_max` summarizes only PairwiseAligner values.

The nearest-reference columns in `reference_landscape_summary.csv` include matching
`*_similarity_method` columns so consumers do not have to infer which scale produced a
reported similarity.

Keep these outputs separate from other COFOLDER surfaces:

- They are not binding-affinity predictions.
- They are not model-confidence scores.
- They are not structure-derived reproduction metrics.
- They do not predict expected cofolding success.

Custom references can change the nearest-reference landscape, but that does not establish
model quality or biological performance.

## Related

- [Installation](../getting-started/installation.md)
- [Validate Command](validate.md)
- [Bias API Reference](../api/recipes/bias.md)
