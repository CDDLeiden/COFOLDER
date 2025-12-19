import os

from boltz_lab.modules.utils.helpers import parse_list_as_str

import logging
logger = logging.getLogger(__name__)

# legacy imports - to be striped and placed with modules
from Bio import PDB
import subprocess
from boltz_lab.modules.input import command, system
from boltz_lab.modules.utils import helpers  

class Evaluate(object):
    """High-level orchestrator for evaluation workflow."""
    def __init__(
        self,
        wrk_dir: str,
        system_path: str,
        options_path: str,
        repeats: int = 1,
        seeds: str | None = None,
        input_pdb: str | None = None,
        ifp: str | None = None,
        generate_conformers: str | None = None
    ):
        self.wrk_dir = wrk_dir
        self.system_path = system_path
        self.options_path = options_path
        self.repeats = repeats
        self._seeds = seeds
        self.seeds = parse_list_as_str(self._seeds, separator=',', item_type=int, expected_length=self.repeats) if self._seeds else []
        self.input_pdb = input_pdb
        self.ifp = ifp
        self.generate_conformers = generate_conformers

        logger.debug("Initializing Evaluate with parameters: %s", {
            "wrk_dir": wrk_dir,
            "system_path": system_path,
            "options_path": options_path,
            "repeats": repeats,
            "seeds": seeds,
            "input_pdb": input_pdb,
            "ifp": ifp,
            "generate_conformers": generate_conformers
        })

        self.opt = None
        self.sys = None

        logger.debug("Evaluate initialization complete.")

    def run(self):
        pass

    def _run(self):
        self.load_config()      # helper function that loads YAML config files
        # Load system and options
        #self._options = utils.read_yaml(path=self.options_path)
        #self.opt = command.Command(options=self._options)

        #self._system = utils.read_yaml(path=self.system_path)
        #self.sys = system.System(system=self._system)

        reference_structure = self._load_reference_structure()

        self.run_boltz()
        ##### Generated code from here, to be implemented #####
        
        self._extract_confidence_affinity()
        self._sequence_similarity()
        self._ligand_similarity()

        if reference_structure:
            self._calculate_rmsd_and_pocket(reference_structure)
            self._extract_ifp(reference_structure)
            self._calculate_ifp_overlap()

        self._aggregate_results()
        self.logger.info("Validation complete.")

    def run_boltz(self):
        self.logger.info(f"Running Boltz predictions for {self.repeats} repeats")
        self.results = []

        self.opt.system_path = self.system_path
        self.name = os.path.splitext(os.path.basename(self.system_path))[0]

        for i in range(self.repeats):
            self.opt.out_dir = os.path.join(self.wrk_dir, f'{self.name}_{str(i+1)}')
            seed = self.seeds[i] if self.seeds else None
            self.logger.info(f"Repeat {i+1}/{self.repeats}, using seed={seed}")
            self.opt.update_options(value=seed, path=["options", 0, "seed"]) 

            print(self.opt.options)
        
            cmd = self.opt.set_command(system=self.sys)
            self.logger.info(f'Running: {" ".join(cmd)}')
            subprocess.run(cmd)

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
        self.sequence_similarity_score = helpers.compare_sequence_to_database(sequence)
  
    def _ligand_similarity(self):
        ligands = self.sys.get_ligands()
        self.ligand_similarity_scores = []
        for ligand in ligands:
            score = helpers.compare_ligand_to_database(ligand)
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
        self.protein_rmsd = helpers.calculate_protein_rmsd(self.sys, reference_structure)
        self.ligand_rmsd = helpers.calculate_ligand_rmsd(self.sys, reference_structure)
        self.binding_pocket_overlap = helpers.calculate_binding_pocket_overlap(self.sys, reference_structure)
        self.logger.info(f"Protein RMSD: {self.protein_rmsd}")
        self.logger.info(f"Ligand RMSD: {self.ligand_rmsd}")
        self.logger.info(f"Binding pocket overlap: {self.binding_pocket_overlap}")

    def _extract_ifp(self, reference_structure):
        if not self.ifp or self.ifp.lower() == 'false':
            self.ifp_data = None
            return

        if self.ifp.lower() == 'true':
            self.ifp_data = helpers.extract_ifp_from_pdb(reference_structure)
        else:
            # assume dictionary provided
            self.ifp_data = self.ifp
        self.logger.info(f"IFP data extracted: {self.ifp_data}") 

    def _calculate_ifp_overlap(self):
        if not self.ifp_data:
            self.ifp_overlap_score = None
            return
        self.ifp_overlap_score = helpers.calculate_ifp_overlap(self.sys, self.ifp_data)
        self.logger.info(f"IFP overlap score: {self.ifp_overlap_score}")
  
    def _aggregate_results(self):
        # Assuming self.results contains Boltz prediction dicts
        self.logger.info("Aggregating results over repeats")
        self.aggregated_results = helpers.aggregate_boltz_results(self.results)

    def _parse_seeds(self, seeds_str: str | None) -> list[int]:
        """Parse seeds string into a list of integers (or strings if non-digit)."""
        if not seeds_str:
            return []
        seeds_list = [
            int(s.strip()) if s.strip().isdigit() else s.strip()
            for s in seeds_str.split(",")
        ]
        return seeds_list
