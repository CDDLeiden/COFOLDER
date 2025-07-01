# boltz_wrapper

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