# COFOLDER

COFOLDER is a collection of command-line workflows for protein-ligand co-folding and pre-cofolding diagnostics.

The supported public entry points are:

- `bias`: inspect protein and ligand reference overlap directly from `system.yaml`
- `validate`: run a single runner-backed co-folding workflow
- `screen`: run validate-style workflows across a ligand library
- `oracle`: return one scalar metric from a runner-backed run

COFOLDER supports a backend-free base install for the standalone `bias` workflow, plus optional backend extras for prediction workflows. Structural-confidence outputs, validation metrics, and affinity-related outputs are treated as distinct concepts; structural confidence should not be described as a proxy for binding affinity.

## Start Here

If you are new to the repository, use this path:

1. Install COFOLDER with the backend you need.
   See [docs/getting-started/installation.md](docs/getting-started/installation.md).
2. Run your first command.
   See [docs/getting-started/quickstart.md](docs/getting-started/quickstart.md).
3. Pick the workflow that matches your task.
   See [docs/user-guide/overview.md](docs/user-guide/overview.md).
4. Move to worked examples and hands-on tutorials.
   Use [examples/README.md](examples/README.md), [docs/tutorials/index.md](docs/tutorials/index.md), and [tutorials/README.md](tutorials/README.md).

## Installation

COFOLDER requires Python 3.11 or 3.12.

Base install:

```bash
git clone https://github.com/CDDLeiden/COFOLDER.git
cd COFOLDER
python -m pip install -e .
```

Recommended prediction install with the default backend:

```bash
python -m pip install -e ".[boltz2]"
```

OpenFold3 install:

```bash
python -m pip install -e ".[openfold3]"
scripts/setup_openfold3.sh
```

The `bias` workflow is part of the base package and does not require a co-folding backend.

## Common Commands

```bash
cofolder [-h] [-v] {bias,validate,screen,oracle}
```

Examples:

```bash
cofolder bias \
  --system_path system.yaml \
  --wrk_dir ./bias_out \
  --custom_protein_reference_path custom_protein.csv \
  --custom_ligand_reference_path custom_ligand.csv

cofolder validate \
  -s examples/system.yaml \
  -o examples/options.yaml \
  -w ./validate_out

cofolder screen \
  -s examples/system_screen.yaml \
  -o examples/options.yaml \
  -c examples/ligand_screen.csv \
  --col_id compound_id \
  --ligand_chain B --smiles_column smiles \
  -w ./screen_out

cofolder oracle \
  -s examples/system.yaml \
  -o examples/options.yaml \
  --input_smiles "CCO" \
  --output_metric affinity_pred_value \
  --aggregate first
```

For detailed setup, workflow selection, and command options, use:

- [docs/getting-started/installation.md](docs/getting-started/installation.md)
- [docs/getting-started/quickstart.md](docs/getting-started/quickstart.md)
- [docs/user-guide/overview.md](docs/user-guide/overview.md)

## Repository Guide

These top-level paths now have distinct roles:

- [`docs/`](docs/): authoritative public documentation for installation, workflow selection, tutorials, and API reference
- [`examples/`](examples/): copyable input files and small fixtures used by quickstarts and docs
- [`tutorials/`](tutorials/): interactive notebook tutorials for hands-on exploration
- [`scripts/`](scripts/): helper utilities for optional setup and bias-data preparation
- [`legacy/`](legacy/): archival material preserved for traceability, not the recommended public path
- [`run_ui.sh`](run_ui.sh): optional launcher for the Streamlit UI in an environment where the UI extra is already installed
- [`LICENSE`](LICENSE) and [`THIRD_PARTY_SOFTWARE.md`](THIRD_PARTY_SOFTWARE.md): repository licensing and third-party attribution

For a short tour of how these pieces fit together, see [docs/getting-started/repository-tour.md](docs/getting-started/repository-tour.md).

## Optional UI

COFOLDER also includes an optional Streamlit interface:

```bash
python -m pip install -e ".[ui]"
./run_ui.sh
```

The UI launcher no longer installs packages implicitly. This keeps setup reproducible and makes the required environment explicit.

## Folder-Specific Notes

- [examples/README.md](examples/README.md) explains which example files to use first.
- [tutorials/README.md](tutorials/README.md) gives a recommended notebook order.
- [scripts/README.md](scripts/README.md) separates supported helpers from maintainer-oriented utilities.
- [legacy/README.md](legacy/README.md) explains what is archived there and why it should not be the first stop for new users.

## Documentation

The MkDocs site is the best public-facing documentation surface. To build it locally:

```bash
python -m pip install -e ".[docs]"
mkdocs serve
```

The docs site complements, rather than replaces, the repository tutorials:

- `docs/` is the curated reading path
- `examples/` is the copy-and-run path
- `tutorials/` is the interactive notebook path

## Licenses

COFOLDER is released under the MIT License. Third-party attributions for adapted or bundled components are listed in [THIRD_PARTY_SOFTWARE.md](THIRD_PARTY_SOFTWARE.md).
