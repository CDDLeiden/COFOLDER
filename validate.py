# Script for validating Boltz system based on ability to recreate poses / find 
# correct binding pocket in reference PDB structure

# Input: 
#   - PDB with (ligand-)protein complex
#   - [OPTIONAL] RESID of ligand in PDB
#   - YAML system file containing: 
#       - FASTA sequence of protein(s)
#       - [OPTIONAL] SMILES of ligand
#       - [OPTIONAL] Co-factors
#       - [OPTIONAL] Boltz Constraints and/or Templates
#   - YAML option file containing Boltz parameters
# Output: 
#   - RMSD co-folded protein to reference protein
#   - RMSD co-folded ligand to reference ligand

# This script should co-fold the complex using YAML system file. Afterwards
# the system is aligned to the reference protein. Co-factors (e.g. ATP) can 
# be added to the system to guide co-folding, but results are only calculated 
# for the protein (and ligand).

# Output would be the RMSD of the diffusion samples to the reference protein(s)
# (and ligand).

# EXTRA IDEAS:
#   - Extract SMILES directly from PDB using RESID.
#   - Allow for Grid search of Boltz parameters, co-factors and Boltz Constraints/Templates
#   - Add IFP profiling for the co-folded and reference ligands

# TODO: everything

if __name__ == "__main__":
    pass