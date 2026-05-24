# COFOLDER

## Introduction
COFOLDER is a collection of command-line utilities for performing co-folding workflows through pluggable backends. The default built-in backend is the `boltz2` runner, which powers the current `validate`, `screen`, and `oracle` workflows in this repository.

COFOLDER can be installed without a backend so the package remains flexible. Additional backends are possible, but they require backend-specific runners and integration code.

Third-party license attributions for vendored code are listed in `THIRD_PARTY_LICENSES.md`.


## Installation
Install directly from GitHub for newest updates. On Python 3.11 or 3.12, the recommended setup installs the default `boltz2` backend:
```
git clone https://github.com/CDDLeiden/COFOLDER.git
cd COFOLDER
pip install -e ".[boltz2]"
```

This is the recommended installation for current prediction workflows on supported backend Python versions.

COFOLDER itself requires Python 3.11+. The current built-in backend extras split by package line:

- `boltz1`: packaged for Python `<3.13`
- `boltz2`: packaged for Python `<3.13`
- `boltz-community`: installed from Git and not currently version-gated in `pyproject.toml`

### Backend-Agnostic Base Install
If you only want the base COFOLDER package without a backend:
```bash
pip install -e .
```

This installs COFOLDER without a co-folding backend. On Python 3.11 or 3.12, you can add the default `boltz2` backend later with:
```bash
pip install -e ".[boltz2]"
```

Available backend extras are `boltz1`, `boltz2`, and `boltz-community`.

### Backend Notes
- **Default backend**: `boltz2`
- **Backend-free install**: supported
- **Backend Python support**:
  - `boltz1`: Python `<3.13`
  - `boltz2`: Python `<3.13`
  - `boltz-community`: no Python version marker currently declared in `pyproject.toml`
- **Other backends**: possible, but require backend-specific runners and integration work before COFOLDER commands can use them

### Optional: Bias-Assessment Setup (MMseqs2)
For protein sequence-similarity bias metrics, install `mmseqs2` in the same environment where you run `cofolder`:
```bash
conda install -c conda-forge -c bioconda mmseqs2
```

## Usage
The main command is `cofolder`, which supports several subcommands:
```
cofolder [-h] [-v] {validate,screen,oracle}
```

### Subcommands
- **validate**: Co-fold and validate a single system using the default `boltz2` backend.
- **screen**: Co-fold a library using the default `boltz2` backend for virtual screening.
- **oracle**: Run single-input oracle scoring (`--input_smiles` or `--input_mol_file`) with the default `boltz2` backend and return one metric value.

Use the -h flag with any command to see detailed usage:
```
cofolder -h
cofolder validate -h
cofolder screen -h
cofolder oracle -h
```

Validate example:
```bash
cofolder validate \
  -s system.yaml \
  -o options.yaml \
  -w ./validate_out
```

Screening example (CSV -> per-row validate wrapper):
```bash
cofolder screen \
  -s system.yaml \
  -o options.yaml \
  -c compounds.csv \
  --col_id compound_id \
  --variable sequences,1,ligand,smiles --col_variable smiles
```
Outputs include `screen_results.csv` and `screen_results_with_scores.csv`.

Oracle example:
```bash
cofolder oracle \
  -s system.yaml \
  -o options.yaml \
  --input_smiles "CCO" \
  --output_metric affinity_pred_value \
  --aggregate first
```

## Package Structure

This is a high-level snapshot, not the authoritative runner authoring reference. For runner contract details and examples, use [docs/tutorials/runners.md](docs/tutorials/runners.md).

```
src/cofolder/
├── __main__.py
├── cli.py
├── recipes/
│   ├── validate.py
│   ├── screen.py
│   └── oracle.py
├── modules/
│   ├── analytics/
│   │   ├── align.py
│   │   ├── bias.py
│   │   ├── bias_training.py
│   │   ├── dataset.py
│   │   ├── plots.py
│   │   ├── reproduction.py
│   │   ├── stats.py
│   │   ├── structure.py
│   │   └── sucos.py
│   ├── entities/
│   │   └── ligand.py
│   ├── input/
│   │   ├── command.py
│   │   └── system.py
│   ├── runners/
│   │   ├── __init__.py
│   │   ├── base.py
│   │   ├── boltz1_runner.py
│   │   ├── boltz2_runner.py
│   │   ├── boltz_community_runner.py
│   │   ├── boltz_runner.py
│   │   ├── contracts.py
│   │   └── validators.py
│   └── utils/
│       ├── gather.py
│       ├── helpers.py
│       ├── read.py
│       └── write.py
└── ui/
    └── app.py

scripts/
├── fetch_bias_training_data.py
├── build_bias_training_data.py
├── install_mmseqs_vendor.sh
└── fetch_bias_training_data_tmp_mmseqs_env.sh

docs/
tests/
tutorials/
examples/
LICENSE
pyproject.toml
README.md
THIRD_PARTY_LICENSES.md
```
