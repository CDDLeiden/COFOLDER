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

import os

def add_arguments(parser):
    """Add validate-specific CLI arguments."""
    parser.add_argument('-w', '--wrk_dir',
                        type=str,
                        dest='wrk_dir',
                        help='Set Working directory if different from CWD.',
                        default=os.getcwd())

    parser.add_argument('-s', '--system_path',
                        type=str,
                        dest='system_path',
                        help='Path to system YAML file.',
                        required=True)

    parser.add_argument('-b', '--boltz_options_path',
                        type=str,
                        dest='options_path',
                        help='Path to Boltz options YAML file.',
                        required=True)
    
    parser.add_argument('-i', '--input_pdb',
                        type=str,
                        dest='input_pdb',
                        help='Path to reference PDB file containing the (ligand-)protein complex.',
                        required=True)

    parser.add_argument('--generate_conformers',
                        choices=['2D', '3D'],
                        default=None,
                        dest='generate_conformers',
                        help='Generate 2D or 3D conformers for CCD input. If not specified, \
                            original SMILES (csv) or MolBlock (sdf) are used as system input. \
                            Note: This only works for SMILES, not any other variable type.')
    
    parser.add_argument('-d', '--debug',
                        action='store_true',
                        help='Enable debug logging')
def main(args):
    print(f"[PREDICT] Running in {args.wrk_dir} with system={args.yaml_system}")