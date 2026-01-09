import gemmi
import numpy as np
from pathlib import Path
import json
import logging
import pandas as pd

logger = logging.getLogger(__name__)


class Structure:
    """
    Structure-level utilities for CIF-based interaction analysis.
    """

    def __init__(
        self,
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