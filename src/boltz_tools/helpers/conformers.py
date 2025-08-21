import os
import pandas as pd
from typing import Optional, Union, List
import pickle

from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor

from ..helpers import command

import logging
ccd_logger = logging.getLogger('boltz-tools.helpers.ccd')

def iterate_sdf_records(sdf_path: str, id_property: str):
    """
    Yield (index, id, molblock) for each valid molecule in the SDF.
    """
    from rdkit import Chem
    suppl = Chem.SDMolSupplier(sdf_path)
    for i, mol in enumerate(suppl, 1):
        if mol is None:
            continue
        mol_id = mol.GetProp(id_property) if mol.HasProp(id_property) else f"mol_{i}"
        yield i, mol_id, Chem.MolToMolBlock(mol)

def _prepare_mol(mol):
    """Assign per-atom names and create a Mol copy via MolBlock."""
    for i, atom in enumerate(mol.GetAtoms()):
        atom.SetProp("name", f"{atom.GetSymbol()}{i+1}")
    return mol

def _save_mol(mol, mol_id, mols_dir):
    """Pickle the molecule to the specified directory."""
    out_path = os.path.join(mols_dir, f"{mol_id}.pkl")
    with open(out_path, "wb") as f:
        pickle.dump(mol, f)
    ccd_logger.info(f"Pickled molecule saved at {out_path} for ID: {mol_id}")


def cache_mols_from_sdf(file_path: str, property_id: str, cache: str = "~/.boltz/"):
    """
    Reads molecules from SDF, assigns atom names, and pickles each molecule.
    """
    file_path = os.path.expanduser(file_path)
    cache = os.path.expanduser(cache)

    # Ensure cache exists
    if not os.path.exists(cache):
        ccd_logger.info(f"Cache path {cache} does not exist. Downloading default cache...")
        command.download_cache(cache)

    mols_dir = os.path.join(cache, "mols")
    os.makedirs(mols_dir, exist_ok=True)

    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"Input file not found: {file_path}")

    suppl = Chem.SDMolSupplier(file_path)
    for mol in suppl:
        if mol is None:
            continue
        
        mol_id = mol.GetProp(property_id) if mol.HasProp(property_id) else None
        if not mol_id:
            ccd_logger.error(f"SDF molecule missing ID property '{property_id}'")
            continue
        try:
            mol = _prepare_mol(mol)
            _save_mol(mol, mol_id, mols_dir)
        except Exception as e:
            ccd_logger.error(f"Failed to process ID {mol_id}: {e}")

def csv_to_sdf(
    csv_path: str,
    smiles_col: str,
    output_sdf_path: Optional[str] = None,
    property_cols: Optional[Union[str, List[str]]] = None
) -> None:
    """
    Convert a CSV file with SMILES to an SDF file, writing specified columns as SDF properties.

    Parameters
    ----------
    csv_path : str
        Path to the input CSV file.
    smiles_col : str
        Name of the column containing SMILES strings.
    output_sdf_path : str, optional
        Path to the output SDF file. If None, replaces .csv with .sdf.
    property_cols : str or list of str, optional
        Column(s) to write as SDF properties for each molecule.
    """

    if output_sdf_path is None:
        output_sdf_path = os.path.splitext(csv_path)[0] + ".sdf"

    df = pd.read_csv(csv_path)
    if smiles_col not in df.columns:
        ccd_logger.error(f"SMILES column '{smiles_col}' not found in {csv_path}")
        return

    # Normalize property_cols to a list
    if property_cols is None:
        property_cols = []
    elif isinstance(property_cols, str):
        property_cols = [property_cols]
    else:
        property_cols = list(property_cols)

    # Warn for missing columns
    for col in property_cols:
        if col not in df.columns:
            ccd_logger.warning(f"Property column '{col}' not found in {csv_path}. It will be skipped.")
    property_cols = [col for col in property_cols if col in df.columns]

    writer = Chem.SDWriter(output_sdf_path)
    n_written = 0
    for idx, row in df.iterrows():
        smi = row[smiles_col]
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            ccd_logger.warning(f"Invalid SMILES at row {idx}: {smi}")
            continue
        # Set properties
        for col in property_cols:
            val = row[col]
            if pd.notnull(val):
                mol.SetProp(str(col), str(val))
        writer.write(mol)
        n_written += 1
    writer.close()
    ccd_logger.info(f"Wrote {n_written} molecules to {output_sdf_path}")

# -------------------------------
# Conformer Generation
# -------------------------------
def generate_2d_conformers(sdf_in: str, sdf_out: Optional[str] = None):
    """
    Generate 2D conformers for molecules in an SDF file.
    """
    sdf_out = sdf_out or sdf_in
    suppl = Chem.SDMolSupplier(sdf_in)
    mols = [mol for mol in suppl if mol is not None]

    for mol in mols:
        rdDepictor.Compute2DCoords(mol)

    writer = Chem.SDWriter(sdf_out)
    for mol in mols:
        writer.write(mol)
    writer.close()
    ccd_logger.info(f"Generated 2D conformers in {sdf_out}")


def generate_3d_conformers(sdf_in: str, sdf_out: Optional[str] = None):
    """
    Generate 3D conformers for molecules in an SDF file using ETKDG + UFF.
    """
    sdf_out = sdf_out or sdf_in
    suppl = Chem.SDMolSupplier(sdf_in)
    mols = [mol for mol in suppl if mol is not None]

    for mol in mols:
        try:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            AllChem.UFFOptimizeMolecule(mol)
        except Exception as e:
            ccd_logger.warning(f"3D conformer generation failed for molecule: {e}")

    writer = Chem.SDWriter(sdf_out)
    for mol in mols:
        writer.write(mol)
    writer.close()
    ccd_logger.info(f"Generated 3D conformers in {sdf_out}")
