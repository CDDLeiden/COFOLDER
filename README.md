# COFOLDER

## Introduction
COFOLDER is a collection of command-line utilities for performing co-folding workflows through pluggable backends. The default backend is Boltz, which powers the current `validate`, `screen`, and `oracle` workflows in this repository.

COFOLDER can be installed without a backend so the package remains flexible. Additional backends are possible, but they require backend-specific runners and integration code.

Third-party license attributions for vendored code are listed in `THIRD_PARTY_LICENSES.md`.


## Installation
Install directly from GitHub for newest updates. The default setup installs the Boltz backend:
```
git clone https://github.com/CDDLeiden/COFOLDER.git
cd COFOLDER
pip install -e ".[boltz]"
```

This is the recommended installation for current prediction workflows.

COFOLDER currently supports Python 3.11 and 3.12.

### Backend-Agnostic Base Install
If you only want the base COFOLDER package without a backend:
```bash
pip install -e .
```

This installs COFOLDER without Boltz. To run the default prediction workflows later, add the Boltz backend with:
```bash
pip install -e ".[boltz]"
```

The `boltz` extra currently installs `boltz[cuda]`.

### Backend Notes
- **Default backend**: Boltz
- **Backend-free install**: supported
- **Other backends**: possible, but require backend-specific runners and integration work before COFOLDER commands can use them

### Optional: Bias-Assessment Setup (MMseqs2)
For protein sequence-similarity bias metrics, install `mmseqs2` in the same environment where you run `cofolder`:
```bash
conda install -c conda-forge -c bioconda mmseqs2
```

Alternative without conda-forge (vendor a release binary):
```bash
chmod +x scripts/install_mmseqs_vendor.sh
scripts/install_mmseqs_vendor.sh
export PATH="$HOME/.cofolder/vendor/mmseqs/bin:$PATH"
```

Verify:
```bash
which mmseqs
mmseqs --help
```

Step-by-step:
```bash
# 1) install COFOLDER with the default Boltz backend
git clone https://github.com/CDDLeiden/COFOLDER.git
cd COFOLDER
pip install -e ".[boltz]"

# Alternative: install backend-free base package only
# pip install -e .

# 2) optional for bias assessment: install mmseqs2
conda install -c conda-forge -c bioconda mmseqs2

# 2b) optional alternative: vendor mmseqs2 binary (no conda-forge)
chmod +x scripts/install_mmseqs_vendor.sh
scripts/install_mmseqs_vendor.sh
export PATH="$HOME/.cofolder/vendor/mmseqs/bin:$PATH"

# 3) fetch bias training assets (CCD + mmseqs DB)
python scripts/fetch_bias_training_data.py \
  --output_root /path/to/training_data \
  --overwrite

# 4) run validate with bias
cofolder validate ... \
  --assess_bias \
  --protein_training_data_path /path/to/training_data/protein_training_data.csv \
  --ligand_training_data_path /path/to/training_data/ligand_training_data.csv
```

## Usage
The main command is `cofolder`, which supports several subcommands:
```
cofolder [-h] [-v] {validate,screen,oracle}
```

### Subcommands
- **validate**: Co-fold and validate a single system using the default Boltz backend.
- **screen**: Co-fold a library using the default Boltz backend for virtual screening.
- **oracle**: Run single-input oracle scoring (`--input_smiles` or `--input_mol_file`) with the default Boltz backend and return one metric value.

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
  -b options.yaml \
  -w ./validate_out
```

Screening example (CSV -> per-row validate wrapper):
```bash
cofolder screen \
  -s system.yaml \
  -b options.yaml \
  -c compounds.csv \
  --col_id compound_id \
  --variable sequences,1,ligand,smiles --col_variable smiles
```
Outputs include `screen_results.csv` and `screen_results_with_scores.csv`.

Oracle example:
```bash
cofolder oracle \
  -s system.yaml \
  -b options.yaml \
  --input_smiles "CCO" \
  --output_metric affinity_pred_value \
  --aggregate first
```

## Package Structure
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
│   │   └── boltz_runner.py
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
