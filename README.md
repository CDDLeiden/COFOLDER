# boltz_wrapper

## Installation
Make sure to have an operational version of boltz working in your conda environment. 
pip install boltz=2.0.3

You can check if the istall is functional using the following command:
boltz predict --help

Note, if you have not downloaded the cache yet, the first run will download this, which may take a few hours.

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