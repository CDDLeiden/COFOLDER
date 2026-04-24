# COFOLDER

## Introduction
COFOLDER is a collection of command-line utilities for performing co-folding workflows through pluggable runners. The current built-in runner is `boltz`, and the `validate`, `screen`, and `oracle` workflows dispatch through a generic runner interface rather than calling one backend directly.

COFOLDER can be installed without a backend so the package remains flexible. Additional backends are possible, but they require backend-specific runner implementations and normalization into COFOLDER's canonical output bundle.

Third-party license attributions for vendored code are listed in `THIRD_PARTY_LICENSES.md`.


## Installation
Install directly from GitHub for newest updates. The default setup installs the `boltz` runner:
```
git clone https://github.com/CDDLeiden/COFOLDER.git
cd COFOLDER
pip install -e ".[boltz]"
```

This is the recommended installation for current prediction workflows, because `boltz` is the only runner shipped in the current phase of the model-agnostic architecture.

`pyproject.toml` currently requires Python 3.11+.

### Runner-Agnostic Base Install
If you only want the base COFOLDER package without any runner dependency:
```bash
pip install -e .
```

This installs COFOLDER without `boltz`. To run the current prediction workflows later, add the `boltz` runner with:
```bash
pip install -e ".[boltz]"
```

The `boltz` extra currently installs `boltz[cuda]`.

### Runner Notes
- **Shipped runner**: `boltz`
- **CLI runner selection**: `--runner <name>`
- **Backend-free install**: supported
- **Other backends**: possible, but require backend-specific runners before COFOLDER commands can use them

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
# 1) install COFOLDER with the default shipped runner
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
  --runner boltz \
  --assess_bias \
  --protein_training_data_path /path/to/training_data/protein_training_data.csv \
  --ligand_training_data_path /path/to/training_data/ligand_training_data.csv
```

## Usage
The main command is `cofolder`, which supports several subcommands:
```bash
cofolder [-h] [-v] {validate,screen,oracle} ...
```

### Subcommands
- **validate**: Co-fold and validate a single system using the selected runner.
- **screen**: Co-fold a library using the selected runner for virtual screening.
- **oracle**: Run single-input oracle scoring (`--input_smiles` or `--input_mol_file`) with the selected runner and return one metric value.

### Common Runner Arguments
- `--runner`: select the cofolding runner. In the current release this is `boltz`.
- `-o` / `--options_path`: path to the runner options YAML.
- `--scoring_functions`: generic metric groups and analytics, such as `confidence_metrics`, `affinity_metrics`, `affinity_metrics_ext`, `ifp_distance`, `ifp_prolif`, `sasa`, and `sasa_normalized`.

If a requested metric group is not supported by the selected runner, COFOLDER warns and continues, leaving the unsupported output columns empty.

Use the -h flag with any command to see detailed usage:
```bash
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
  --runner boltz \
  --scoring_functions confidence_metrics affinity_metrics sasa \
  -w ./validate_out
```

Screening example (CSV -> per-row validate wrapper):
```bash
cofolder screen \
  -s system.yaml \
  -o options.yaml \
  --runner boltz \
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
  --runner boltz \
  --input_smiles "CCO" \
  --output_metric affinity_pred_value \
  --aggregate first
```

### Canonical Outputs
Runner-specific raw prediction artifacts remain in the per-repeat `raw/` tree, but COFOLDER normalizes them into a shared result bundle consumed by analytics:

- `results/system_metrics.csv`
- `results/chain_metrics.csv`
- `results/structures/`

Each repeat also writes normalized runner metadata under `raw/repeat_<n>/normalized/`, including a manifest describing the runner, capabilities, normalized sample records, and raw artifact locations.

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
│   │   ├── base.py
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
