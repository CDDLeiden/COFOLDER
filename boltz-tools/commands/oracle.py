# Script for running Boltz as an oracle function.

# Input: 
#   - SMILES that varies per oracle call.
#   - YAML system file containing:
#       - FASTA sequence of protein
#       - [OPTIONAL] Co-factors
#       - [OPTIONAL] Contraints and/or Templates
#   - choice of output metric:
#       - output confidence metric (confidence_score, ptm, ...)    
#       - output affinity (affinity_pred_value, affinity_probability_binary, ...)
#       - RMSD of diffusion samples
# Output:
#   - Metric of choice   

# This script should take in a singular SMILES and return a singular float
# value. 

# EXTRA IDEAS:
#   - Add FASTA sequence as variable input.

# TODO: everything


if __name__ == "__main__":
    pass