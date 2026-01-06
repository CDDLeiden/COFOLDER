# Boltz-Lab Documentation

Welcome to the official documentation for **Boltz-Lab**, a comprehensive collection of command-line tools designed for performing various co-folding tasks using Boltz.

## Overview

Boltz-Lab provides a suite of utilities that streamline protein-ligand co-folding workflows, including:

- **Single System Predictions**: Co-fold individual protein-ligand systems
- **Virtual Screening**: High-throughput screening of ligand libraries
- **Oracle Functions**: Use Boltz as a scoring function for molecular design
- **Evaluation Tools**: Assess and validate co-folding results with detailed metrics

## Key Features

- **Easy-to-use CLI**: Intuitive command-line interface for all workflows
- **Flexible Input Formats**: Support for SMILES, SDF, PDB, and CIF formats
- **Comprehensive Analytics**: Built-in tools for RMSD calculation, interaction fingerprints, and visualization
- **Conformer Generation**: Automated 2D and 3D conformer generation for ligands
- **Batch Processing**: Efficient screening of large compound libraries

## Quick Start

```bash
# Install boltz-lab
git clone https://github.com/CDDLeiden/boltz-lab.git
cd boltz-lab
pip install -e .

# Run a simple prediction
boltz-lab predict -s system.yaml -b options.yaml
```

## Navigation

- **[Getting Started](getting-started/installation.md)**: Installation and basic setup
- **[User Guide](user-guide/overview.md)**: Detailed guides for each command
- **[Tutorials](tutorials/basic.md)**: Step-by-step tutorials and examples
- **[API Reference](api/recipes/predict.md)**: Complete API documentation

## Support

For issues, questions, or contributions, please visit our [GitHub repository](https://github.com/CDDLeiden/boltz-lab).

## Citation

If you use Boltz-Lab in your research, please cite the Boltz paper and this repository.
