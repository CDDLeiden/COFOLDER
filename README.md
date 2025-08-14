# boltz_tools

## Introduction 
Collection of helpful tools for performing various co-folding tasks using Boltz.

## Installation

Install directly from GitHub for newest updates:
```
git clone https://github.com/rlvandenbroek/boltz-tools.git
cd boltz-tools; pip install -e .
```

## Package Structure
```
boltz_tools/
│
├── __init__.py
├── __main__.py
├── cli.py                     # Root CLI
│
├── tools/                     # Tools / subcommands
│   ├── __init__.py
│   ├── predict.py
│   ├── validate.py
│   ├── screen.py
│   └── oracle.py
│
├── helpers/                   # General utilities
│   ├── __init__.py
│   ├── command.py
│   ├── helpers.py
│   └── system.py
│
├── analytics/                 # Data analysis utilities
│   ├── __init__.py
│   ├── stats.py
│   └── plotting.py
│
tests/
pyproject.toml
```

## How to run
from boltz_wrapper/templates copy options.yaml into your cwd (TODO: make this adjustable when calling the wrapper).
In options.yaml, specify wrapper options.
  - run_dir: directory where results will be stored
  - system: path to YAML system file
  - ligands:
    - lig_csv: path to CSV file containing SMILES
    - smiles_col: SMILES column name
    - id_col: Compound identifier column name
Additionally boltz options can be adjusted in this file.

for the YAML system file, as for now, keep the following format (i.e. first define protein and ligand sequences, only then all other sytem features):

sequences:
  - protein:
      id: [A]
      sequence: {your FASTA sequence}
  - ligand:
      id: [B]
      smiles: {your SMILES}

Run the screening.py with the following command: python {your_path}/boltz_wrapper/screening.py

## screening.py - WORK IN PROGRESS
Script for performing virtual screening using Boltz.

Input: 
  - [OPTIONAL] SMILES or FASTA sequence to inject into YAML system
  - [OPTIONAL] CSV file containing SMILES or FASTA sequences to inject into YAML system
  - [OPTIONAL] SMILES column when iterating through CSV file
  - [OPTIONAL] ID column when iterating through CSV file (if not specified, index will be used)
  - YAML system file containing: 
      - FASTA sequence of protein(s)
      - [OPTIONAL] Co-factors
      - [OPTIONAL] Boltz Constraints and/or Templates
Output:
  - output confidence metric (confidence_score, ptm, ...)    
  - output affinity (affinity_pred_value, affinity_probability_binary, ...)
  - RMSD of diffusion samples

### Usage: Merging Additional Columns from Input CSV to Output CSV
You can specify extra columns from your input CSV to be included in the output CSV using the `--merge_columns` argument:

```bash
python screening.py options.yaml --merge_columns "column1,column2,extra_info"
```
This will ensure that `column1`, `column2`, and `extra_info` from your input CSV are present in the output CSV, in addition to the default columns.

#### Output CSV Columns
The output CSV will always include the following columns at the front:
- `index`: Row number in the input CSV (starting from 1)
- `id`: Value from your specified ID column
- `basename`: Unique basename for each row
- `smiles`: Value from your specified SMILES column
- Any columns specified in `--merge_columns`

## oracle.py - NON-OPERATIONAL
Script for running Boltz as an oracle function.

Input: 
  - SMILES that varies per oracle call.
  - YAML system file containing:
      - FASTA sequence of protein
      - [OPTIONAL] Co-factors
      - [OPTIONAL] Contraints and/or Templates
  - choice of output metric:
      - output confidence metric (confidence_score, ptm, ...)    
      - output affinity (affinity_pred_value, affinity_probability_binary, ...)
      - RMSD of diffusion samples
Output:
  - Metric of choice   

## validate.py - NON-OPERATIONAL
Script for validating Boltz system based on ability to recreate poses / find 
correct binding pocket in reference PDB structure

Input: 
  - PDB with (ligand-)protein complex
  - [OPTIONAL] RESID of ligand in PDB
  - YAML system file containing: 
      - FASTA sequence of protein(s)
      - [OPTIONAL] SMILES of ligand
      - [OPTIONAL] Co-factors
      - [OPTIONAL] Boltz Constraints and/or Templates
  - YAML option file containing Boltz parameters
Output: 
  - RMSD co-folded protein to reference protein
  - RMSD co-folded ligand to reference ligand
