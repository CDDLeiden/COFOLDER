# Standard library
import json
import logging
from pathlib import Path
import pickle
import re

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

    def add_ifp_prolif(
        self,
        column_name: str = "ifp_prolif",
        save_folder: Path = Path("results/ifp/prolif"),
    ) -> pd.DataFrame:
        """
        Add ProLIF interaction fingerprints to chain_df.

        - Must have exactly one protein defined in system
        - Converts CIF -> PDB on‑the‑fly for ProLIF compatibility
        - Saves fingerprint DataFrame per chain as CSV
        """
        save_folder = self.wrk_dir / save_folder
        save_folder.mkdir(parents=True, exist_ok=True)

        for idx, row in self.chain_df.iterrows():
            if row[self.entity_type_col] != "ligand":
                continue

            cif_path = self.cif_folder / row[self.cif_file_col]
            if not cif_path.exists():
                logger.warning("Missing CIF file: %s", cif_path)
                continue

            # Deterministic pickle name
            pickle_name = f"{cif_path.stem}_ligand_chain_{row[self.chain_id_col]}_ifp.pkl"
            pickle_path = save_folder / pickle_name

            # ---- EARLY EXIT UNLESS OVERWRITE ----
            if pickle_path.exists():
                logger.info("IFP exists, skipping (overwrite=False): %s", pickle_path)
                self.chain_df.at[idx, column_name] = pickle_name
                continue

            # Step 1: CIF -> PDB
            pdb_file = self.cif_to_pdb(
                cif_file=cif_path,
                output_folder=self.wrk_dir / "results/structures/pdb"
            )

            # Step 2: Split PDB
            split_files = self.split_pdb(
                pdb_file=pdb_file,
                chain_id=row[self.chain_id_col],
                output_folder=self.wrk_dir / "results/structures/pdb"
            )

            protein_file = Path(split_files.get("protein"))
            ligand_file = Path(split_files.get("ligand"))

            # ---- Load protein ----
            u_protein = mda.Universe(str(protein_file))
            u_protein = self._sanitize_protein(u_protein)
            protein_mol = plf.Molecule.from_mda(u_protein)

            logger.info(
                "Protein loaded: %s | residues=%d | atoms=%d",
                protein_file,
                protein_mol.n_residues,
                protein_mol.GetNumAtoms(),
            )

            # ---- Load ligand ----
            u_ligand = mda.Universe(str(ligand_file))
            ligand_mol = plf.Molecule.from_mda(u_ligand)

            logger.info(
                "Ligand loaded: %s | residues=%d | atoms=%d",
                ligand_file,
                ligand_mol.n_residues,
                ligand_mol.GetNumAtoms(),
            )

            # ---- Generate ProLIF fingerprints ----
            fp = plf.Fingerprint()
            ifp = fp.generate(ligand_mol, protein_mol, metadata=True)
            df = plf.to_dataframe({0: ifp}, fp.interactions)
            logger.info("ProLIF fingerprint matrix:\n%s", df.T)

            # ---- Save IFP pickle ----
            pickle_name = ligand_file.stem + "_ifp.pkl"
            pickle_path = save_folder / pickle_name

            with open(pickle_path, "wb") as f:
                pickle.dump(ifp, f)

            logger.info("IFP saved to %s", pickle_path)

            self.chain_df.at[idx, column_name] = pickle_name

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

    def cif_to_pdb(self, cif_file: Path, output_folder: Path) -> Path | None:
        """
        Convert CIF file to a PDB file with a deterministic name.

        Output name:
            <cif_filename>.pdb

        Args:
            cif_file: Path to CIF file.
            output_folder: Folder where the PDB will be written.

        Returns:
            Path to generated PDB file.
        """
        if not cif_file.exists():
            logger.error("CIF file not found: %s", cif_file)
            return None

        output_folder.mkdir(parents=True, exist_ok=True)

        pdb_file = output_folder / f"{cif_file.stem}.pdb"

        structure = gemmi.read_structure(str(cif_file))
        structure.write_pdb(str(pdb_file))

        logger.info("Converted CIF -> PDB: %s", pdb_file)

        return pdb_file


    def split_pdb(self, pdb_file: Path, chain_id: str, output_folder: Path) -> dict:
        """
        Split a PDB file into protein and ligand PDBs, ensuring hydrogens are present.
        Hydrogens are added using RDKit if none are found.

        Args:
            pdb_file: Path to input PDB.
            chain_id: Ligand chain ID.
            output_folder: Folder to save split PDBs.

        Returns:
            dict with keys 'protein' and 'ligand' and their PDB paths.
        """
        output_folder.mkdir(parents=True, exist_ok=True)
        universe = mda.Universe(str(pdb_file))

        # ---- Protein ----
        protein_pdb_path = output_folder / f"{pdb_file.stem}_protein.pdb"

        if protein_pdb_path.exists():
            logger.warning("Protein PDB already exists at %s. Skipping and not overwriting.", protein_pdb_path)
        else:
            # ---- Protein ----
            prot_atoms = universe.select_atoms("protein and not resname HOH")
            prot_atoms.write(str(protein_pdb_path))

            # Add hydrogens using PDB2PQR
            n_atoms_before = len(prot_atoms)
            protein_pdb_path = self.add_hydrogens_with_pdb2pqr(protein_pdb_path, protein_pdb_path)
            u_prot_h = mda.Universe(str(protein_pdb_path))
            n_atoms_after = len(u_prot_h.atoms)
            logger.info("Protein atoms: before H=%d, after H=%d", n_atoms_before, n_atoms_after)

            # Check proper protonation
            mol = Chem.MolFromPDBFile(str(protein_pdb_path), removeHs=False)
            if mol is None:
                logger.error("RDKit failed to parse protein PDB! Check hydrogens/valence.")
            else:
                logger.debug("Protein PDB parsed by RDKit: atoms=%d", mol.GetNumAtoms())

        # ---- Ligand ----
        ligand_pdb_path = self._sanitize_ligand(
            universe=universe,
            chain_id=chain_id,
            output_folder=output_folder,
            pdb_file=pdb_file
        )

        return {"protein": protein_pdb_path, "ligand": ligand_pdb_path}

    def add_hydrogens_with_pdb2pqr(self, input_pdb: Path, output_pdb: Path, ph: float = 7.4):
        """
        Add hydrogens to a protein PDB using PDB2PQR (minimal options).

        Args:
            input_pdb: Input protein PDB file
            output_pdb: Output PDB file with hydrogens added
            ph: pH for protonation (default 7.4)

        Returns:
            Path to output PDB with hydrogens
        """
        args = [
            str(input_pdb),
            str(output_pdb),
            "--ff", "PARSE",
            "--with-ph", str(ph),
            "--pdb-output", str(output_pdb),
        ]
        run_pdb2pqr(args)
        return output_pdb        

    def _remove_overbonding_hydrogens(self, u_protein):
        """Remove only the HZ3 hydrogen from Lys residues."""

        # Select all HZ3 hydrogens in Lys residues
        hs_to_remove = u_protein.select_atoms("resname LYS and name HZ3")

        if len(hs_to_remove) == 0:
            logger.info("_remove_overbonding_hydrogens | No Lys HZ3 atoms to remove.")
            return u_protein

        # Log which hydrogens are being removed
        removed_list = [f"{atom.resname}{atom.resid}-{atom.name}" for atom in hs_to_remove]
        logger.debug("_remove_overbonding_hydrogens | Removing Lys HZ3 atoms: %s", ", ".join(removed_list))

        # Remove the selected atoms
        u_clean = u_protein.atoms[np.setdiff1d(np.arange(len(u_protein.atoms)), hs_to_remove.indices)]

        logger.info("_remove_overbonding_hydrogens | Total HZ3 atoms removed: %d", len(hs_to_remove))
        return u_clean

    def _sanitize_protein(self, u_protein, max_attempts=100):
        """
        Iteratively remove hydrogens from overbonded atoms until
        RDKit can create a molecule without valence errors.

        No additional RDKit/plf calls are made inside the loop.
        Atom indices remain stable by using a persistent mask.
        """
        n_atoms = len(u_protein.atoms)
        keep_mask = np.ones(n_atoms, dtype=bool)

        attempts = 0

        while attempts < max_attempts:
            try:
                # Pass a VIEW to RDKit; original Universe untouched
                _ = plf.Molecule.from_mda(u_protein.atoms[keep_mask])

                if attempts == 0:
                    logger.info(
                        "_sanitize_protein | Protein sanitized successfully on first attempt."
                    )
                else:
                    logger.info(
                        "_sanitize_protein | Protein sanitized successfully after %d attempts.",
                        attempts,
                    )

                # Return a pruned Universe view
                return mda.Merge(u_protein.atoms[keep_mask])

            except Chem.rdchem.AtomValenceException as e:
                msg = str(e)
                match = re.search(r"atom # (\d+) (\w+),", msg)
                if not match:
                    logger.error(
                        "_sanitize_protein | Could not parse RDKit error: %s", msg
                    )
                    raise

                rdkit_index = int(match.group(1))
                atom_symbol = match.group(2)

                # RDKit index == index into u_protein.atoms[keep_mask]
                active_indices = np.flatnonzero(keep_mask)

                if rdkit_index >= len(active_indices):
                    logger.error(
                        "_sanitize_protein | RDKit atom index %d out of range",
                        rdkit_index,
                    )
                    raise

                mda_index = active_indices[rdkit_index]
                target_atom = u_protein.atoms[mda_index]

                logger.warning(
                    "_sanitize_protein | Atom %d (%s) overbonded, attempting to remove hydrogen.",
                    mda_index,
                    target_atom.name,
                )

                # ---- Case 1: overbonded atom is hydrogen ----
                if target_atom.name.startswith("H"):
                    logger.debug(
                        "_sanitize_protein | Removing overbonded hydrogen %s (atom id %d) from residue %s%d",
                        target_atom.name,
                        target_atom.index,
                        target_atom.resname,
                        target_atom.resid,
                    )
                    keep_mask[target_atom.index] = False

                # ---- Case 2: remove an attached hydrogen ----
                else:
                    target_ag = u_protein.atoms[target_atom.index : target_atom.index + 1]

                    attached_hs = u_protein.select_atoms(
                        "resid {} and name H* and around 1.2 group target_ag".format(
                            target_atom.resid
                        ),
                        target_ag=target_ag,
                    )

                    # Only consider still-active hydrogens
                    attached_hs = attached_hs[
                        keep_mask[attached_hs.indices]
                    ]

                    if len(attached_hs) == 0:
                        logger.error(
                            "_sanitize_protein | No hydrogens found to remove for overbonded atom %s%d",
                            target_atom.resname,
                            target_atom.resid,
                        )
                        raise

                    h = attached_hs[0]
                    logger.debug(
                        "_sanitize_protein | Removing attached hydrogen %s (atom id %d) from residue %s%d",
                        h.name,
                        h.index,
                        target_atom.resname,
                        target_atom.resid,
                    )
                    keep_mask[h.index] = False

                attempts += 1

        raise RuntimeError(
            f"_sanitize_protein | Failed to sanitize protein after {max_attempts} attempts"
        )
    
    def _sanitize_ligand(self, universe, chain_id: str, output_folder: Path, pdb_file: Path):
        canonical_chain = chain_id
        canonical_resid = 1
        canonical_resname = "LIG"

        # Select ligand atoms
        lig_atoms = universe.select_atoms(f"segid {chain_id} or chainID {chain_id}")

        # Update residues
        for res in lig_atoms.residues:
            res.resid = canonical_resid
            res.resname = canonical_resname
            res.segment.segid = canonical_chain

        # Update atom-level chainID
        for atom in lig_atoms:
            atom.chainID = canonical_chain

        # Write initial PDB
        output_folder.mkdir(parents=True, exist_ok=True)
        ligand_pdb_path = output_folder / f"{pdb_file.stem}_ligand_chain_{chain_id}.pdb"
        lig_atoms.write(str(ligand_pdb_path))
        n_atoms_before = len(lig_atoms)

        # Load PDB with RDKit
        lig_mol = Chem.MolFromPDBFile(str(ligand_pdb_path), removeHs=False)
        n_h = len([atom for atom in lig_mol.GetAtoms() if atom.GetAtomicNum() == 1])

        if n_h == 0:
            logger.warning("Ligand has no hydrogens; adding using RDKit.")
            lig_mol_h = Chem.AddHs(lig_mol, addCoords=True)
            ff = AllChem.UFFGetMoleculeForceField(lig_mol_h)
            heavy_idx = [atom.GetIdx() for atom in lig_mol_h.GetAtoms() if atom.GetAtomicNum() > 1]
            for idx in heavy_idx:
                ff.AddFixedPoint(idx)
            ff.Minimize()

            # Generate PDB block
            pdb_block = Chem.MolToPDBBlock(lig_mol_h)

            # Patch all HETATM/ATOM lines for canonical chain/resid/resname
            pdb_lines = []
            for line in pdb_block.splitlines():
                if line.startswith("HETATM") or line.startswith("ATOM"):
                    line = (
                        line[:17] + f"{canonical_resname:<3}" +  # resname
                        line[20:21] + canonical_chain +          # chainID
                        f"{canonical_resid:>4}" +                # resid
                        line[26:]
                    )
                pdb_lines.append(line)

            with open(ligand_pdb_path, "w") as f:
                f.write("\n".join(pdb_lines))

            n_atoms_after = lig_mol_h.GetNumAtoms()
        else:
            n_atoms_after = lig_mol.GetNumAtoms()

        logger.info("Ligand atoms: before H=%d, after H=%d", n_atoms_before, n_atoms_after)
        return ligand_pdb_path