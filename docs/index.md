# COFOLDER Documentation

Welcome to the official documentation for **COFOLDER**, a comprehensive collection of command-line tools designed for performing various co-folding tasks using Boltz.

## Overview

COFOLDER provides a suite of utilities that streamline protein-ligand co-folding workflows, including:

- **Single System Predictions**: Co-fold individual protein-ligand systems
- **Virtual Screening**: High-throughput screening of ligand libraries
- **Oracle Functions**: Use Boltz as a scoring function for molecular design

## Key Features

- **Easy-to-use CLI**: Intuitive command-line interface for all workflows
- **Flexible Input Formats**: Support for SMILES, SDF, PDB, and CIF formats
- **Comprehensive Analytics**: Built-in tools for RMSD calculation, interaction fingerprints, and visualization
- **Conformer Generation**: Automated 2D and 3D conformer generation for ligands
- **Batch Processing**: Efficient screening of large compound libraries

## Quick Start

```bash
# Install COFOLDER
git clone https://github.com/CDDLeiden/COFOLDER.git
cd COFOLDER
pip install -e .

# Optional for bias assessment: install mmseqs2
conda install -c conda-forge -c bioconda mmseqs2

# Optional alternative: vendor mmseqs2 binary (no conda-forge)
scripts/install_mmseqs_vendor.sh
export COFOLDER_MMSEQS_BIN="$HOME/.cofolder/vendor/mmseqs/bin/mmseqs"

# Verify mmseqs2
which mmseqs
mmseqs --version

# Run a simple prediction
cofolder validate -s system.yaml -b options.yaml
```

## Navigation

- **[Getting Started](getting-started/installation.md)**: Installation and basic setup
- **[User Guide](user-guide/overview.md)**: Detailed guides for each command
- **[Tutorials](tutorials/basic.md)**: Step-by-step tutorials and examples
- **[Adding New Runners](tutorials/runners.md)**: Authoritative guide for building a new backend runner and normalizing its outputs for COFOLDER analytics
- **[API Reference](api/recipes/validate.md)**: Complete API documentation

## Support

For issues, questions, or contributions, please visit our [GitHub repository](https://github.com/CDDLeiden/COFOLDER).

## Citation

If you use COFOLDER in your research, please cite the Boltz paper and this repository.
