# Script for validating Boltz system based on ability to recreate poses / find 
# correct binding pocket in reference PDB structure

import os
import logging
from Bio import PDB

from ..helpers import command, conformers, system, utils

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

    parser.add_argument('--repeats',
                        type=int,
                        default=1,
                        dest='repeats',
                        help=('Number of repeats to run the calculation with different seeds. '
                              'Each repeat will use a different random seed if unspecified. Must be >=1.'))

    parser.add_argument('--seeds',
                        type=str,
                        default=None,
                        dest='seeds',
                        help=('Optional, List of integer seeds to use for repeats. '
                              'Must be the same length as --repeats, e.g. "[42, 123, 999]".'))
    
    parser.add_argument('-i', '--input_pdb',
                        type=str,
                        dest='input_pdb',
                        help='Path to reference CIF/PDB file containing a (ligand-)protein system. \
                            By adding a system, protein and ligand RMSD calculations will be \
                            performed between the co-folded and reference structures.')

    parser.add_argument('--ifp',
                        type=str,
                        dest='ifp',
                        help=('Interaction fingerprint (IFP) specification. Accepts: \n'
                              '  - "true": extract from crystal structure (requires input_pdb)\n'
                              '  - "false": skip IFP overlap calculation\n'
                              '  - A dictionary of interacting residues for manual specification, e.g.\n'
                              '      \'{"A:123":"ARG", "B:45":"TYR"}\''))

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
    """Run the validate system tool."""
    utils.set_dir(args.wrk_dir) 

    logger = logging.getLogger('boltz-tools')
    initiate_logger(logger, debug=args.debug, wrk_dir=args.wrk_dir)
    logger.info("Boltz-tools validation started.")

    validate = Validate(
        wrk_dir=args.wrk_dir,
        system_path=args.system_path,
        options_path=args.options_path,
        repeats=args.repeats,
        seeds=args.seeds,
        input_pdb=args.input_pdb,
        ifp=args.ifp,
        generate_conformers=args.generate_conformers
    )

def initiate_logger(logger, debug, wrk_dir):
    log_file = os.path.join(wrk_dir, 'boltz-tools.log')
    open(log_file, 'w+').close()

    fh = logging.FileHandler(log_file)
    ch = logging.StreamHandler()

    fh.setLevel(logging.DEBUG if debug else logging.INFO)
    ch.setLevel(logging.WARNING)

    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d,%H:%M:%S'
    )
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    logger.setLevel(logging.DEBUG if debug else logging.INFO)
    logger.addHandler(fh)
    logger.addHandler(ch)

class Validate(object):
    """
    Class to validate Boltz system based on ability to recreate poses / find 
    correct binding pocket in reference PDB structure.
    """
    def __init__(self, **kwargs):
        self.wrk_dir = kwargs.get("wrk_dir")
        self.system_path = kwargs.get("system_path")
        self.options_path = kwargs.get("options_path")

        self.repeats = kwargs.get("repeats")
        self._seeds = kwargs.get("seeds", [])
        self.seeds = [int(v.strip()) if v.strip().isdigit() else v.strip() for v in self._seeds.split(",")] if self._seeds else []

        # Validate seeds length if provided
        if self.seeds is not None and len(self.seeds) != self.repeats:
            raise ValueError(f"Number of seeds ({len(self.seeds)}) must match --repeats ({self.repeats}).")

        self.input_pdb = kwargs.get("input_pdb")
        self.ifp = kwargs.get("ifp")
        self.generate_conformers = kwargs.get("generate_conformers")

        self.logger = logging.getLogger('boltz-tools.validate.Validate')
        self.logger.debug(f"Validate args: {kwargs}")

        # Set command and system objects
        self._options = utils.read_yaml(path=self.options_path)
        self.opt = command.Command(options=self._options)
        self._system = utils.read_yaml(path=self.system_path)
        self.sys = system.System(system=self._system)

        self.run()

    def run(self):
        reference_structure = self._load_reference_structure()

##### Generated code from here, to be implemented #####
        self._run_predict()
#### ^ WIP ^ ####
        
        self._extract_confidence_affinity()
        self._sequence_similarity()
        self._ligand_similarity()

        if reference_structure:
            self._calculate_rmsd_and_pocket(reference_structure)
            self._extract_ifp(reference_structure)
            self._calculate_ifp_overlap()

        self._aggregate_results()
        self.logger.info("Validation complete.")

    def _run_predict(self):
        self.logger.info(f"Running Boltz predictions for {self.repeats} repeats")
        self.results = []

        for i in range(self.repeats):
            seed = self.seeds[i] if self.seeds else None
            self.logger.info(f"Repeat {i+1}/{self.repeats}, using seed={seed}")
            self.opt.update_options(value=seed, path="options,seed")
            
            # Run Boltz predict (assuming self.opt.predict returns dict with results)
            result = self.opt.predict(system=self.sys, seed=seed)
            self.results.append(result)

    def _extract_confidence_affinity(self):
        for i, result in enumerate(self.results):
            confidence = result.get('confidence', None)
            affinity = result.get('affinity', None)
            self.logger.info(f"Repeat {i+1}: confidence={confidence}, affinity={affinity}")
            result['confidence'] = confidence
            result['affinity'] = affinity
    
    def _sequence_similarity(self):
        sequence = self.sys.get_sequence()  # assuming System class has this
        self.logger.info(f"Protein sequence: {sequence}")
        # placeholder for similarity calculation
        self.sequence_similarity_score = utils.compare_sequence_to_database(sequence)
  
    def _ligand_similarity(self):
        ligands = self.sys.get_ligands()
        self.ligand_similarity_scores = []
        for ligand in ligands:
            score = utils.compare_ligand_to_database(ligand)
            self.ligand_similarity_scores.append(score)
            self.logger.info(f"Ligand similarity: {score}")

    def _load_reference_structure(self):
        if not self.input_pdb:
            self.logger.warning("No input structure provided. Skipping template-based RMSD/IFP calculations.")
            return None

        file_ext = os.path.splitext(self.input_pdb)[1].lower()


        try:
            if file_ext in ['.cif', '.mmcif']:
                parser = PDB.MMCIFParser(QUIET=True)
            elif file_ext == '.pdb':
                parser = PDB.PDBParser(QUIET=True)
            else:
                raise RuntimeError(f"Unsupported file type: {file_ext}")

            structure = parser.get_structure('reference', self.input_pdb)
            self.logger.info(f"Loaded reference structure: {self.input_pdb}")
            return structure

        except Exception as e:
            raise RuntimeError(f"Failed to read structure file: {e}")

    def _calculate_rmsd_and_pocket(self, reference_structure):
        from Bio.SVDSuperimposer import SVDSuperimposer  # or custom utils
        # Placeholder function calls:
        self.protein_rmsd = utils.calculate_protein_rmsd(self.sys, reference_structure)
        self.ligand_rmsd = utils.calculate_ligand_rmsd(self.sys, reference_structure)
        self.binding_pocket_overlap = utils.calculate_binding_pocket_overlap(self.sys, reference_structure)

        self.logger.info(f"Protein RMSD: {self.protein_rmsd}")
        self.logger.info(f"Ligand RMSD: {self.ligand_rmsd}")
        self.logger.info(f"Binding pocket overlap: {self.binding_pocket_overlap}")

    def _extract_ifp(self, reference_structure):
        if not self.ifp or self.ifp.lower() == 'false':
            self.ifp_data = None
            return

        if self.ifp.lower() == 'true':
            self.ifp_data = utils.extract_ifp_from_pdb(reference_structure)
        else:
            # assume dictionary provided
            self.ifp_data = self.ifp
        self.logger.info(f"IFP data extracted: {self.ifp_data}") 

    def _calculate_ifp_overlap(self):
        if not self.ifp_data:
            self.ifp_overlap_score = None
            return
        self.ifp_overlap_score = utils.calculate_ifp_overlap(self.sys, self.ifp_data)
        self.logger.info(f"IFP overlap score: {self.ifp_overlap_score}")
  
    def _aggregate_results(self):
        # Assuming self.results contains Boltz prediction dicts
        self.logger.info("Aggregating results over repeats")
        self.aggregated_results = utils.aggregate_boltz_results(self.results)

if __name__ == "__main__":
    pass