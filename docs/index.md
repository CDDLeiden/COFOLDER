# COFOLDER Documentation

Welcome to the official documentation for **COFOLDER**, a collection of command-line
workflows for runner-backed protein-ligand cofolding and pre-cofolding bias diagnostics.

## Overview

COFOLDER provides four main workflows:

- **Bias Diagnostics**: Inspect training-reference overlap directly from `system.yaml`
- **Single System Predictions**: Co-fold individual protein-ligand systems
- **Virtual Screening**: Run validate-style prediction across a ligand library
- **Oracle Functions**: Return one scalar from a runner-backed validation

## Key Features

- **Easy-to-use CLI**: Intuitive command-line interface for all workflows
- **Flexible Input Formats**: Support for SMILES, SDF, PDB, and CIF formats
- **Comprehensive Analytics**: Built-in tools for bias, confidence, reproduction, and structure-derived diagnostics
- **Conformer Generation**: Automated 2D and 3D conformer generation for ligands
- **Batch Processing**: Efficient screening of large compound libraries

## Quick Start

```bash
# Install COFOLDER
git clone https://github.com/CDDLeiden/COFOLDER.git
cd COFOLDER
pip install -e .

# Optional for standalone public-reference bias-data builds
conda install -c conda-forge -c bioconda mmseqs2

# Optional alternative: vendor mmseqs2 binary
scripts/install_mmseqs_vendor.sh
export COFOLDER_MMSEQS_BIN="$HOME/.cofolder/vendor/mmseqs/bin/mmseqs"

# Run standalone bias without any backend runner
cofolder bias \
  --system_path system.yaml \
  --custom_protein_reference_path custom_protein.csv \
  --custom_ligand_reference_path custom_ligand.csv

# Install a backend only when you want runner-backed prediction workflows
pip install -e ".[boltz2]"

# Run a simple prediction
cofolder validate -s system.yaml -o options.yaml
```

## Navigation

- **[Getting Started](getting-started/installation.md)**: Installation and basic setup
- **[User Guide](user-guide/overview.md)**: Detailed guides for `bias`, `validate`, `screen`, and `oracle`
- **[Tutorials](tutorials/basic.md)**: Step-by-step tutorials and examples
- **[Adding New Runners](tutorials/runners.md)**: Authoritative guide for building a new backend runner and normalizing its outputs for COFOLDER analytics
- **[Backend Acceptance Tutorial](tutorials/backend-acceptance.md)**: Clean-install, real-backend validation lane for runner and CLI changes
- **[API Reference](api/recipes/bias.md)**: Complete API documentation

## Support

For issues, questions, or contributions, please visit our [GitHub repository](https://github.com/CDDLeiden/COFOLDER).

## Citation

If you use COFOLDER in your research, please cite the Boltz paper and this repository.
