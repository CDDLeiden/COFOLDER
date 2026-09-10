# Installation

## Prerequisites

COFOLDER requires:

- Python 3.11 or 3.12
- CUDA-compatible GPU (for Boltz)
- Git (for installation from source)

## Installation from Source

The recommended method is to install directly from GitHub:

```bash
# Clone the repository
git clone https://github.com/CDDLeiden/COFOLDER.git
cd COFOLDER

# Install the base package
pip install -e .
```

This installs the base COFOLDER package and its core dependencies, including:

- `rdkit` - For molecular structure handling
- `matplotlib` - For plotting and visualization
- `seaborn` - For enhanced visualizations

Backend runners are installed separately through optional extras:

```bash
# Boltz 1 runner
pip install -e ".[boltz1]"

# Boltz 2 runner
pip install -e ".[boltz2]"

# Boltz Community runner
pip install -e ".[boltz-community]"

# OpenFold3 runner
pip install -e ".[openfold3]"
```

The integrated `openfold3` runner still has a second setup step after installation because the upstream model cache, checkpoints, and CCD need to be prepared:

```bash
# Install the COFOLDER package with the OpenFold3 backend extra
pip install -e ".[openfold3]"

# Keep the downloaded cache, parameters, and CCD in the repository cache folder
export OPENFOLD_CACHE="$PWD/cache/.openfold3-cache"
cofolder-tools setup-openfold3
```

Notes for OpenFold3:

- The COFOLDER `openfold3` extra installs the upstream `openfold3` package.
- The example above sets `OPENFOLD_CACHE` to the gitignored `cache/.openfold3-cache` directory. If the variable is omitted, `cofolder-tools setup-openfold3` retains the upstream `~/.openfold3` default.
- `cofolder-tools setup-openfold3` also answers the standard upstream setup prompts explicitly: it uses `OPENFOLD_CACHE` for both cache/checkpoint-root questions, selects parameter download choice `1`, and declines forced redownloads by default. Set `OPENFOLD3_FORCE_DOWNLOAD_PARAMETERS=yes` only when you intentionally want to replace an existing checkpoint.
- The tool answers the upstream integration-test prompt with `no` by default so bootstrap does not hang in an unexpected interactive test path. If you intentionally want those tests, run `OPENFOLD3_RUN_INTEGRATION_TESTS=yes cofolder-tools setup-openfold3`.
- If you want all published checkpoints instead of the default checkpoint only, run `OPENFOLD3_PARAMETER_CHOICE=2 cofolder-tools setup-openfold3`.
- If you want a different cache root, set `OPENFOLD_CACHE` before running the tool, for example `OPENFOLD_CACHE=/scratch/$USER/.openfold3 cofolder-tools setup-openfold3`.
- Upstream docs note that first inference can also download default model parameters into `$HOME/.openfold3`, but this project prefers the explicit setup script so readiness is established before prediction runs.
- The current OpenFold3 integration is confidence-only: it exposes OpenFold3-native confidence outputs and keeps affinity groups unsupported unless a future code change intentionally widens that capability.
- For clean-install manual acceptance in a fresh environment, see the [Backend Acceptance Tutorial](../tutorials/backend-acceptance.md).

For clean-install backend acceptance in fresh environments, see the [Backend Acceptance Tutorial](../tutorials/backend-acceptance.md).

## Verify Installation

After installation, verify that COFOLDER is correctly installed:

```bash
cofolder --version
cofolder-tools --help
```

To create a local, editable copy of the packaged example bundle from any working
directory:

```bash
cofolder-tools copy-examples ./cofolder-examples
```

You should see output showing the version number.

## Optional Dependencies

### Standalone Bias Workflow

The dedicated `cofolder bias` command is part of the base package. You do not need a
cofolding backend runner just to inspect reference-overlap diagnostics.

Base-package-only bias usage works for:

- custom-only reference inputs
- prebuilt public protein/ligand training CSVs
- mixed public-plus-custom reference analysis

Install `mmseqs2` only if you want to build public protein bias-training data from a
system definition in the same environment where you run `cofolder`:

```bash
conda install -c conda-forge -c bioconda mmseqs2
```

Alternative without conda-forge:

```bash
cofolder-tools install-mmseqs
export COFOLDER_MMSEQS_BIN="$HOME/.cofolder/vendor/mmseqs/bin/mmseqs"
```

Verify:

```bash
which mmseqs
mmseqs --version
```

Bias-only setup after installing the base package:

```bash
# 1) Optional: fetch CCD + mmseqs DB if you want to build public training CSVs
cofolder-tools fetch-bias-training-data \
  --output_root /path/to/training_data \
  --overwrite

# 2) Optional: build public protein/ligand training CSVs
cofolder-tools build-bias-training-data \
  --system_path /path/to/system.yaml \
  --components_cif /path/to/training_data/ccd/components.cif \
  --output_protein_csv /path/to/training_data/protein_training_data.csv \
  --output_ligand_csv /path/to/training_data/ligand_training_data.csv \
  --release_cutoff 2023-06-01 \
  --overwrite

# 3a) Run standalone bias from prebuilt public references
cofolder bias \
  --system_path /path/to/system.yaml \
  --wrk_dir /path/to/bias_run \
  --protein_training_data_path /path/to/training_data/protein_training_data.csv \
  --ligand_training_data_path /path/to/training_data/ligand_training_data.csv

# 3b) Or run standalone bias from custom-only references
cofolder bias \
  --system_path /path/to/system.yaml \
  --wrk_dir /path/to/bias_run \
  --custom_protein_reference_path /path/to/custom_protein.csv \
  --custom_ligand_reference_path /path/to/custom_ligand.csv
```

When the builder writes `bias_training_data.csv`, its `sequence_similarity` values are
MMseqs `pident` percentages. The accompanying `sequence_similarity_method` column is
`mmseqs_pident` when a hit exists and `unavailable` when MMseqs returned no result; an
unavailable protein similarity is left empty so it cannot pass a protein-similarity
threshold. The configured threshold filters the standalone `protein_training_data.csv`
but does not discard MMseqs hits needed to annotate PDBs in the combined table.

Downstream bias analysis never puts a PairwiseAligner result in `sequence_similarity`.
When PairwiseAligner is used for sequence-only custom references, its percentage is
written to `sequence_similarity_pairwise`, with method `pairwise_aligner`; the MMseqs
column remains empty.

`cofolder bias` writes diagnostic artifacts under `<wrk_dir>/results/bias_train/`.
These are reference-overlap outputs for pre-cofolding decision support. They are not
binding-affinity predictions, model-confidence scores, or guarantees of cofolding
success.

If you want the same bias diagnostics attached to a runner-backed validation run, use
`validate --assess_bias` with the same public reference CSVs:

```bash
cofolder validate ... \
  --assess_bias \
  --protein_training_data_path /path/to/training_data/protein_training_data.csv \
  --ligand_training_data_path /path/to/training_data/ligand_training_data.csv
```

### Documentation Tools

To build the documentation locally:

```bash
pip install -e ".[docs]"
```

This installs:

- MkDocs with Material theme
- mkdocstrings for API reference generation

### Tutorial Notebooks

To run the interactive marimo tutorials:

```bash
pip install -e ".[tutorials]"
```

Then launch a tutorial notebook, for example:

```bash
marimo edit tutorials/bias.py
```

## Development Installation

For development work, install with all optional dependencies:

```bash
pip install -e ".[docs,test,development]"
```

Add the backend extra you need in the same environment, for example:

```bash
pip install -e ".[boltz2]"
```

If you are developing against OpenFold3 instead, keep the same development environment and install the matching COFOLDER backend extra plus the setup script:

```bash
pip install -e ".[openfold3]"
cofolder-tools setup-openfold3
```

## Troubleshooting

### CUDA Issues

If you encounter CUDA-related errors, ensure:

1. You have a CUDA-compatible GPU
2. CUDA drivers are properly installed
3. PyTorch is installed with CUDA support

### RDKit Issues

If RDKit installation fails, try installing via conda:

```bash
conda install -c conda-forge rdkit
```

## Next Steps

Once installed, proceed to the [Quick Start Guide](quickstart.md) to begin using COFOLDER.
