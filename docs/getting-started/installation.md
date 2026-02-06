# Installation

## Prerequisites

COFOLDER requires:

- Python 3.10, 3.11, or 3.12
- CUDA-compatible GPU (for Boltz)
- Git (for installation from source)

## Installation from Source

The recommended method is to install directly from GitHub:

```bash
# Clone the repository
git clone https://github.com/CDDLeiden/cofolder.git
cd cofolder

# Install in editable mode
pip install -e .
```

This will install COFOLDER along with all required dependencies, including:

- `boltz[cuda]` - The core Boltz package with CUDA support
- `rdkit` - For molecular structure handling
- `matplotlib` - For plotting and visualization
- `seaborn` - For enhanced visualizations

## Verify Installation

After installation, verify that COFOLDER is correctly installed:

```bash
cofolder --version
```

You should see output showing the version number.

## Optional Dependencies

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
