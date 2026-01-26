# Standard library
from contextlib import contextmanager
import json
import logging
from pathlib import Path
import pickle
import re
import tempfile

# Third-party libraries
import numpy as np
import pandas as pd
import gemmi
from rdkit import Chem
from rdkit.Chem import AllChem

# Project-specific / external tools
from Bio.PDB import MMCIFParser
from Bio.PDB.SASA import ShrakeRupley
import MDAnalysis as mda
import prolif as plf
from pdb2pqr.main import run_pdb2pqr

# Logger
logger = logging.getLogger(__name__)


class Structure:
    """
    Structure-level utilities for CIF-based interaction analysis.
    """

    def __init__(
        self,
        wrk_dir: Path,
        chain_df: pd.DataFrame,
        cif_folder: Path,
        chain_id: str = "CHAIN_ID",
        entity_type: str = "ENTITY_TYPE",
        cif_file: str = "cif_file",
    ):
        """
        Parameters
        ----------
        chain_df : pd.DataFrame
            Chain-level DataFrame.
        cif_folder : Path
            Base directory containing CIF files.
        chain_id : str
            Column name for chain identifiers.
        entity_type : str
            Column name for entity type (e.g. protein, ligand).
        cif_file : str
            Column name for CIF file names.
        """
        self.wrk_dir = Path(wrk_dir)
        self.chain_df = chain_df
        self.cif_folder = Path(cif_folder)

        self.chain_id_col = chain_id
        self.entity_type_col = entity_type
        self.cif_file_col = cif_file

        self.receptor_chain_id = self._derive_receptor_chain()


    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def add_ifp_distance(
        self,
        cutoff: float = 5.0,
        column_name: str = "ifp_distance",
    ) -> pd.DataFrame:
        """
        Add distance-based interaction fingerprints to chain_df.

        - Only ligand chains receive fingerprints
        - Protein rows are left as None
        - Fingerprints are stored as JSON strings (CSV-safe)
        """
        if self.receptor_chain_id is None:
            logger.warning(
                "Skipping interaction fingerprints: invalid receptor definition."
            )
            self.chain_df[column_name] = None
            return self.chain_df

        if column_name not in self.chain_df.columns:
            self.chain_df[column_name] = None

        for idx, row in self.chain_df.iterrows():
            if row[self.entity_type_col] != "ligand":
                continue

            cif_path = self.cif_folder / row[self.cif_file_col]
            if not cif_path.exists():
                logger.warning("Missing CIF file: %s", cif_path)
                continue

            fingerprint = self._distance_interaction_fingerprint(
                cif_path=cif_path,
                receptor_chain=self.receptor_chain_id,
                ligand_chain=str(row[self.chain_id_col]),
                cutoff=cutoff,
            )

            self.chain_df.at[idx, column_name] = fingerprint

        return self.chain_df
    
    def add_sasa(
        self,
        *,
        absolute: bool = False,
        normalized: bool = False,
        heavy_atoms_only: bool = True,
    ) -> pd.DataFrame:
        """
        Add solvent-accessible surface area (SASA) metrics per chain.

        Parameters
        ----------
        absolute : bool
            If True, compute absolute SASA.
        normalized : bool
            If True, compute normalized SASA.
        heavy_atoms_only : bool
            If True, normalize by heavy atoms only.
            Ignored if normalized=False.
        """

        if not absolute and not normalized:
            logger.debug("No SASA metrics requested; skipping.")
            return self.chain_df

        logger.info(
            "Starting SASA calculation | absolute=%s | normalized=%s",
            absolute,
            normalized,
        )

        if absolute and "sasa" not in self.chain_df.columns:
            self.chain_df["sasa"] = None

        if normalized and "sasa_norm_heavy" not in self.chain_df.columns:
            self.chain_df["sasa_norm_heavy"] = None

        parser = MMCIFParser(QUIET=True)
        sr = ShrakeRupley()

        for idx, row in self.chain_df.iterrows():
            cif_path = self.cif_folder / row[self.cif_file_col]
            chain_id = str(row[self.chain_id_col])

            if not cif_path.exists():
                logger.warning("Missing CIF file: %s", cif_path)
                continue

            try:
                structure = parser.get_structure("struct", str(cif_path))
                sr.compute(structure, level="A")

                model = structure[0]
                if chain_id not in model:
                    logger.warning(
                        "Chain %s not found in CIF %s",
                        chain_id,
                        cif_path.name,
                    )
                    continue

                chain = model[chain_id]

                total_sasa = 0.0
                atom_count = 0

                for residue in chain:
                    for atom in residue:
                        if not hasattr(atom, "sasa"):
                            continue

                        total_sasa += atom.sasa

                        if normalized:
                            if heavy_atoms_only:
                                if atom.element != "H":
                                    atom_count += 1
                            else:
                                atom_count += 1

                if absolute:
                    self.chain_df.at[idx, "sasa"] = float(round(total_sasa, 3))

                if normalized:
                    if atom_count == 0:
                        logger.warning(
                            "Normalized SASA requested but atom count is zero | "
                            "CIF=%s | chain=%s",
                            cif_path.name,
                            chain_id,
                        )
                        self.chain_df.at[idx, "sasa_norm_heavy"] = None
                    else:
                        value = total_sasa / atom_count
                        self.chain_df.at[idx, "sasa_norm_heavy"] = float(round(value, 5))

            except Exception:
                logger.exception(
                    "SASA calculation failed | CIF=%s | chain=%s",
                    cif_path.name,
                    chain_id,
                )

                if absolute:
                    self.chain_df.at[idx, "sasa"] = None
                if normalized:
                    self.chain_df.at[idx, "sasa_norm_heavy"] = None

        logger.info(
            "Completed SASA calculation | absolute=%s | normalized=%s",
            absolute,
            normalized,
        )

        return self.chain_df


    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    def _derive_receptor_chain(self) -> str | None:
        """
        Determine the receptor chain from chain_df.

        Rules:
        - Exactly one protein chain must exist
        """
        protein_chains = (
            self.chain_df
            .loc[self.chain_df[self.entity_type_col] == "protein", self.chain_id_col]
            .astype(str)
            .unique()
        )

        if len(protein_chains) != 1:
            logger.warning(
                "Expected exactly one protein chain, found %d: %s",
                len(protein_chains),
                protein_chains.tolist(),
            )
            return None

        return protein_chains[0]

    def _distance_interaction_fingerprint(
        self,
        cif_path: Path,
        receptor_chain: str,
        ligand_chain: str,
        cutoff: float,
    ) -> str:
        """
        Compute a residue-level distance-based interaction fingerprint.

        Returns
        -------
        str
            JSON-serialized list of 0/1 values (CSV-safe).
        """
        structure = gemmi.read_structure(str(cif_path))
        model = structure[0]

        try:
            chain_rec = model[receptor_chain]
            chain_lig = model[ligand_chain]
        except KeyError:
            logger.warning(
                "Missing receptor (%s) or ligand (%s) chain in CIF: %s",
                receptor_chain,
                ligand_chain,
                cif_path,
            )
            return json.dumps([])

        rec_residues = list(chain_rec)
        lig_atoms = [
            atom
            for res in chain_lig
            for atom in res
            if atom.element.name != "H"
        ]

        bitvector = np.zeros(len(rec_residues), dtype=int)

        for i, res in enumerate(rec_residues):
            rec_atoms = [atom for atom in res if atom.element.name != "H"]

            for ra in rec_atoms:
                for la in lig_atoms:
                    if ra.pos.dist(la.pos) <= cutoff:
                        bitvector[i] = 1
                        break
                if bitvector[i]:
                    break

        return json.dumps(bitvector.tolist())
    
    def add_ifp_prolif(
        self,
        column_name: str = "ifp_prolif",
        save_folder: Path = Path("results/ifp/prolif"),
        overwrite: bool = False,
    ) -> pd.DataFrame:
        """
        Generate ProLIF fingerprints for all ligand chains in chain_df.
        Each fingerprint is saved as a pickle file, and the filename
        is stored in the specified column.

        Args:
            overwrite: If False (default), skip rows whose pickle already exists.
                    If True, recompute and overwrite existing pickles.
        """
        save_folder = self.wrk_dir / save_folder
        save_folder.mkdir(parents=True, exist_ok=True)

        pdb_output_folder = self.wrk_dir / "results/structures/pdb"
        pdb_output_folder.mkdir(parents=True, exist_ok=True)

        for idx, row in self.chain_df.iterrows():
            try:
                if row[self.entity_type_col] != "ligand":
                    continue

                cif_path = self.cif_folder / row[self.cif_file_col]
                if not cif_path.exists():
                    logger.warning("Missing CIF file: %s", cif_path)
                    continue

                chain_id = row[self.chain_id_col]

                # Deterministic pickle name
                pickle_name = f"{cif_path.stem}_ligand_chain_{chain_id}_ifp.pkl"
                pickle_path = save_folder / pickle_name

                # ---- EARLY EXIT UNLESS OVERWRITE ----
                if pickle_path.exists() and not overwrite:
                    logger.info("IFP exists, skipping (overwrite=False): %s", pickle_path)
                    self.chain_df.at[idx, column_name] = pickle_name
                    continue

                pdb_file = self.cif_to_pdb(cif_path, pdb_output_folder)
                if pdb_file is None:
                    continue

                pickle_name = self._generate_ifp(
                    pdb_file,
                    chain_id,
                    save_folder,
                )
                self.chain_df.at[idx, column_name] = pickle_name

            except Exception as e:
                logger.warning("Failed processing row %d: %s", idx, e)
                logger.debug("Row data: %s", row)

        return self.chain_df


    def _generate_ifp(self, pdb_file: Path, chain_id: str, save_folder: Path) -> str:
        """
        Build ProLIF-compatible protein and ligand molecules from PDB,
        generate fingerprint, save as pickle, and return filename.
        Entire MDA and ProLIF lifetimes are contained in this function.
        """
        pdb_output_folder = pdb_file.parent

        # ---------- Load universe ----------
        universe = mda.Universe(str(pdb_file))

        # ---------- Protein ----------
        protein_pdb = pdb_output_folder / f"{pdb_file.stem}_protein_chain_{chain_id}.pdb"
        prot_atoms = universe.select_atoms("protein and not resname HOH")
        prot_atoms.write(str(protein_pdb))

        # Add hydrogens
        protein_pdb = self.add_hydrogens_with_pdb2pqr(protein_pdb, protein_pdb)

        # Sanitize protein for ProLIF
        sanitized_protein_pdb = self._sanitize_protein(protein_pdb)

        # ---------- Ligand ----------
        ligand_pdb = self._sanitize_ligand(universe, chain_id, pdb_output_folder, pdb_file)

        # ---------- Generate fingerprint ----------
        fp = plf.Fingerprint()
        ifp = fp.generate(
            plf.Molecule.from_mda(
                mda.Universe(str(ligand_pdb))
        ), 
            plf.Molecule.from_mda(
                mda.Universe(str(sanitized_protein_pdb)),
        ), 
            metadata=True
        )

        # Save fingerprint
        pickle_name = f"{pdb_file.stem}_ligand_chain_{chain_id}_ifp.pkl"
        pickle_path = save_folder / pickle_name
        with open(pickle_path, "wb") as f:
            pickle.dump(ifp, f)

        # Cleanup
        del universe, fp, ifp

        return pickle_name


    def cif_to_pdb(self, cif_file: Path, output_folder: Path) -> Path | None:
        """Convert CIF to PDB with deterministic naming."""
        if not cif_file.exists():
            logger.error("CIF not found: %s", cif_file)
            return None

        output_folder.mkdir(parents=True, exist_ok=True)
        pdb_file = output_folder / f"{cif_file.stem}.pdb"

        if pdb_file.exists():
            logger.info("PDB exists, skipping conversion: %s", pdb_file)
        else:
            structure = gemmi.read_structure(str(cif_file))
            structure.write_pdb(str(pdb_file))
            logger.info("Converted CIF -> PDB: %s", pdb_file)

        return pdb_file

    def _sanitize_protein(self, protein_pdb: Path, max_attempts: int = 100) -> Path:
        """
        Load a protein PDB, iteratively remove problematic hydrogens until
        RDKit/ProLIF can build a molecule without valence errors, and
        write a sanitized PDB.

        Input:
            protein_pdb : Path
                Protein PDB with hydrogens.

        Output:
            Path
                Path to sanitized protein PDB.
        """
        protein_pdb = Path(protein_pdb)

        # Load original Universe
        u_protein = mda.Universe(str(protein_pdb))

        # Deterministic output path (matches original behavior)
        sanitized_pdb = protein_pdb.parent / f"{protein_pdb.stem}_sanitized.pdb"

        # Run Universe-based sanitizer (new logic)
        self._sanitize_protein_universe(
            u_protein,
            sanitized_pdb=sanitized_pdb,
            max_attempts=max_attempts,
        )

        return sanitized_pdb

    def _sanitize_protein_universe(self, u_protein, sanitized_pdb: Path, max_attempts: int = 100):
        """
        Universe-based protein sanitization.
        See previous implementation for full documentation.
        """
        keep_mask = np.ones(len(u_protein.atoms), dtype=bool)

        for attempt in range(max_attempts):
            try:
                _ = plf.Molecule.from_mda(
                    u_protein.atoms[keep_mask]
                )

                logger.info(
                    "_sanitize_protein | Protein sanitized successfully after %d attempts.",
                    attempt,
                )

                u_protein.atoms[keep_mask].write(sanitized_pdb)

                return 

            except Chem.rdchem.AtomValenceException as e:
                msg = str(e)
                match = re.search(r"atom # (\d+) (\w+),", msg)
                if not match:
                    raise RuntimeError(f"Cannot parse RDKit error: {msg}") from e

                rdkit_idx = int(match.group(1))
                active_indices = np.flatnonzero(keep_mask)
                mda_idx = active_indices[rdkit_idx]
                atom = u_protein.atoms[mda_idx]

                if atom.name.startswith("H"):
                    keep_mask[mda_idx] = False
                else:
                    attached_hs = u_protein.select_atoms(
                        f"resid {atom.resid} and name H* and around 1.2 index {atom.index}"
                    )
                    attached_hs = attached_hs[keep_mask[attached_hs.indices]]

                    if not attached_hs:
                        raise RuntimeError(
                            f"No hydrogens removable from overbonded atom {atom.name}"
                        )

                    keep_mask[attached_hs[0].index] = False

        raise RuntimeError(
            f"_sanitize_protein | Failed to sanitize protein after {max_attempts} attempts"
        )

    def _sanitize_protein_obs(self, protein_pdb: Path, max_attempts: int = 100) -> Path:
        """
        Load protein PDB and remove problematic hydrogens to make it ProLIF-compatible.
        Returns sanitized PDB path. Memory-safe.
        """ 
        u_protein = mda.Universe(str(protein_pdb))
        keep_mask = np.ones(len(u_protein.atoms), dtype=bool)

        for attempt in range(max_attempts):
            try:
                # Write temporary PDB
                temp_pdb = protein_pdb.parent / f"{protein_pdb.stem}_sanitized.pdb"
                u_protein.atoms[keep_mask]
                u_protein.write(temp_pdb)

                plf.Molecule.from_mda(
                    mda.Universe(str(temp_pdb)),
                )

                return temp_pdb 

            except Chem.rdchem.AtomValenceException as e:
                # Identify offending hydrogen
                msg = str(e)
                match = re.search(r"atom # (\d+) (\w+),", msg)
                if not match:
                    raise RuntimeError(f"Cannot parse RDKit error: {msg}")
                rdkit_idx = int(match.group(1))
                active_indices = np.flatnonzero(keep_mask)
                mda_idx = active_indices[rdkit_idx]
                atom = u_protein.atoms[mda_idx]

                # Remove hydrogen(s)
                if atom.name.startswith("H"):
                    keep_mask[atom.index] = False
                else:
                    attached_hs = u_protein.select_atoms(
                        f"resid {atom.resid} and name H* and around 1.2 index {atom.index}"
                    )
                    attached_hs = attached_hs[keep_mask[attached_hs.indices]]
                    if not attached_hs:
                        raise RuntimeError("No hydrogens to remove from overbonded atom")
                    keep_mask[attached_hs[0].index] = False

        raise RuntimeError(f"Failed to sanitize protein after {max_attempts} attempts")


    def _sanitize_ligand(self, universe, chain_id: str, output_folder: Path, pdb_file: Path) -> Path:
        """
        Generate a ProLIF-compatible ligand PDB.
        Memory-safe: RDKit objects are deleted immediately after use.
        """
        canonical_chain = chain_id
        canonical_resid = 1
        canonical_resname = "LIG"

        lig_atoms = universe.select_atoms(f"segid {chain_id} or chainID {chain_id}")
        for res in lig_atoms.residues:
            res.resid = canonical_resid
            res.resname = canonical_resname
            res.segment.segid = canonical_chain
        for atom in lig_atoms:
            atom.chainID = canonical_chain

        output_folder.mkdir(parents=True, exist_ok=True)
        ligand_pdb_path = output_folder / f"{pdb_file.stem}_ligand_chain_{chain_id}.pdb"
        lig_atoms.write(str(ligand_pdb_path))

        # Load ligand with RDKit
        lig_mol = Chem.MolFromPDBFile(str(ligand_pdb_path), removeHs=False)
        n_h = sum(1 for atom in lig_mol.GetAtoms() if atom.GetAtomicNum() == 1)

        if n_h == 0:
            # Add hydrogens and optimize
            lig_mol_h = Chem.AddHs(lig_mol, addCoords=True)
            ff = AllChem.UFFGetMoleculeForceField(lig_mol_h)
            heavy_idx = [atom.GetIdx() for atom in lig_mol_h.GetAtoms() if atom.GetAtomicNum() > 1]
            for idx in heavy_idx:
                ff.AddFixedPoint(idx)
            ff.Minimize()

            # Write PDB
            pdb_block = Chem.MolToPDBBlock(lig_mol_h)
            with open(ligand_pdb_path, "w") as f:
                for line in pdb_block.splitlines():
                    if line.startswith(("HETATM", "ATOM")):
                        line = (
                            line[:17] + f"{canonical_resname:<3}" +
                            line[20:21] + canonical_chain +
                            f"{canonical_resid:>4}" +
                            line[26:]
                        )
                    f.write(line + "\n")

            del lig_mol_h, ff, heavy_idx, pdb_block

        del lig_mol
        return ligand_pdb_path

    def add_hydrogens_with_pdb2pqr(self, input_pdb: Path, output_pdb: Path, ph: float = 7) -> Path:
        """
        Add hydrogens to a protein PDB using PDB2PQR in a minimal, deterministic way.

        Args:
            input_pdb: Path to input protein PDB.
            output_pdb: Path to output PDB with hydrogens added.
            ph: pH for protonation (default 7).

        Returns:
            Path to the output PDB with hydrogens added.
        """
        output_pdb = Path(output_pdb)
        input_pdb = Path(input_pdb)
        
        # Ensure output folder exists
        output_pdb.parent.mkdir(parents=True, exist_ok=True)

        # Skip if already exists
        if output_pdb.exists() and output_pdb.resolve() != input_pdb.resolve():
            logger.info("Hydrogen-added PDB already exists, skipping PDB2PQR: %s", output_pdb)
            return output_pdb

        # Construct PDB2PQR command
        args = [
            str(input_pdb),
            str(output_pdb),
            "--ff", "PARSE",            # force field
            "--with-ph", str(ph),       # set pH
            "--pdb-output", str(output_pdb)
        ]

        try:
            # Execute PDB2PQR (blocking call)
            run_pdb2pqr(args)
            logger.info("Added hydrogens to protein via PDB2PQR: %s", output_pdb)
        except Exception as e:
            logger.error("PDB2PQR failed for %s: %s", input_pdb, str(e))
            raise RuntimeError(f"PDB2PQR failed for {input_pdb}") from e

        return output_pdb