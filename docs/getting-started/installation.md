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
```

For clean-install backend acceptance in fresh environments, see the [Backend Acceptance Tutorial](../tutorials/backend-acceptance.md).

## Verify Installation

After installation, verify that COFOLDER is correctly installed:

```bash
cofolder --version
```

You should see output showing the version number.

## Optional Dependencies

### Bias Assessment (MMseqs2)

If you want protein sequence-similarity bias metrics (`--assess_bias` with protein training data),
install `mmseqs2` in the same environment where you run `cofolder`:

```bash
conda install -c conda-forge -c bioconda mmseqs2
```

Alternative without conda-forge (vendor MMseqs2 binary):

```bash
chmod +x scripts/install_mmseqs_vendor.sh
scripts/install_mmseqs_vendor.sh
export COFOLDER_MMSEQS_BIN="$HOME/.cofolder/vendor/mmseqs/bin/mmseqs"
```

Verify:

```bash
which mmseqs
mmseqs --version
```

Step-by-step bias setup after installing the base package and the backend you plan to use:

```bash
# 1) fetch CCD + mmseqs DB
python scripts/fetch_bias_training_data.py \
  --output_root /path/to/training_data \
  --overwrite

# 2) build protein/ligand training CSVs
python scripts/build_bias_training_data.py \
  --system_path /path/to/system.yaml \
  --components_cif /path/to/training_data/ccd/components.cif \
  --output_protein_csv /path/to/training_data/protein_training_data.csv \
  --output_ligand_csv /path/to/training_data/ligand_training_data.csv \
  --release_cutoff 2023-06-01 \
  --overwrite

# 3) run validate with bias
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

## Development Installation

For development work, install with all optional dependencies:

```bash
pip install -e ".[docs,test]"
```

Add the backend extra you need in the same environment, for example:

```bash
pip install -e ".[boltz2]"
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
