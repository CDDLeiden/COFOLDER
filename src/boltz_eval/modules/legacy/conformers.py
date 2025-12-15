import os
import pickle
import logging
from collections import Counter
from pathlib import Path
from typing import Optional, Union, List

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor, rdmolops
from boltz.data.parse.mmcif_with_constraints import parse_ccd_residue

from ..modules import command

ccd_logger = logging.getLogger('boltz-tools.helpers.conformers')

def sanitize_mol_id(mol_id: str) -> str:
    """
    Ensures the molecule ID is at most 5 characters.
    If longer, truncates and prints a warning.
    """
    if len(mol_id) > 5:
        truncated = mol_id[:5]
        ccd_logger.warning(f"[WARNING] Molecule ID '{mol_id}' is longer than 5 characters. "
              f"Truncating to '{truncated}' to comply with CCD naming rules.")
        return truncated
    return mol_id

def add_pickled_prop(mol: Chem.Mol, name: str, value):
    """
    Pickle and set a property on the molecule.

    Parameters
    ----------
    mol : Chem.Mol
        RDKit molecule.
    name : str
        Property name.
    value : Any
        Value to pickle and set.
    """
    mol.SetProp(name, pickle.dumps(value).hex())

def extract_constraints(constraints, key, transpose=False):
    """
    Extract a list of values from constraint objects.

    Parameters
    ----------
    constraints : list
        List of constraint objects.
    key : str
        Attribute name to extract.
    transpose : bool, optional
        Whether to transpose the result (default False).

    Returns
    -------
    list or numpy.ndarray
        Extracted values, optionally transposed.
    """
    result = [getattr(obj, key) for obj in constraints]
    if transpose:
        import numpy as np
        return np.array(result).transpose()
    return result

def _prepare_mol(mol: Chem.Mol) -> Chem.Mol:
    """
    Assign per-atom names and return the molecule.

    Parameters
    ----------
    mol : Chem.Mol
        RDKit molecule.

    Returns
    -------
    Chem.Mol
        Molecule with atom names set.
    """
    for i, atom in enumerate(mol.GetAtoms()):
        atom.SetProp("name", f"{atom.GetSymbol()}{i+1}".upper())
    return mol

def generate_2d_conformers(sdf_in: str, sdf_out: Optional[str] = None):
    """
    Generate 2D conformers for molecules in an SDF file.

    Parameters
    ----------
    sdf_in : str
        Input SDF file path.
    sdf_out : str, optional
        Output SDF file path. If None, overwrites input file.
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

    Parameters
    ----------
    sdf_in : str
        Input SDF file path.
    sdf_out : str, optional
        Output SDF file path. If None, overwrites input file.
    """
    sdf_out = sdf_out or sdf_in
    suppl = Chem.SDMolSupplier(sdf_in)
    mols = [mol for mol in suppl if mol is not None]
    mols_3d = []

    for mol in mols:
        try:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            AllChem.UFFOptimizeMolecule(mol)
            mols_3d.append(mol)
        except Exception as e:
            ccd_logger.warning(f"3D conformer generation failed for molecule: {e}")

    writer = Chem.SDWriter(sdf_out)
    for mol in mols_3d:
        writer.write(mol)
    writer.close()
    ccd_logger.info(f"Generated 3D conformers in {sdf_out}")

def iterate_sdf_records(sdf_path: str, id_property: str):
    """
    Yield (index, id, molblock) for each valid molecule in the SDF file.

    Parameters
    ----------
    sdf_path : str
        Path to the SDF file.
    id_property : str
        Property name to use as molecule ID.

    Yields
    ------
    tuple
        (index, id, molblock) for each valid molecule.
    """
    suppl = Chem.SDMolSupplier(sdf_path)
    for i, mol in enumerate(suppl, 1):
        if mol is None:
            continue
        mol_id = mol.GetProp(id_property) if mol.HasProp(id_property) else f"mol_{i}"
        yield i, mol_id, Chem.MolToMolBlock(mol)

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

    if property_cols is None:
        property_cols = []
    elif isinstance(property_cols, str):
        property_cols = [property_cols]
    else:
        property_cols = list(property_cols)

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
        for col in property_cols:
            val = row[col]
            if pd.notnull(val):
                mol.SetProp(str(col), str(val))
        writer.write(mol)
        n_written += 1
    writer.close()
    ccd_logger.info(f"Wrote {n_written} molecules to {output_sdf_path}")

def _save_mol(mol: Chem.Mol, mol_id: str, mols_dir: str):
    """
    Pickle the molecule to the specified directory.

    Parameters
    ----------
    mol : Chem.Mol
        RDKit molecule.
    mol_id : str
        Molecule identifier.
    mols_dir : str
        Directory to save the pickle file.
    """
    out_path = os.path.join(mols_dir, f"{mol_id}.pkl")
    with open(out_path, "wb") as f:
        pickle.dump(mol, f)
    ccd_logger.info(f"Pickled molecule saved at {out_path} for ID: {mol_id}")

def mol_to_ccd(resname: str, mol: Chem.Mol, boltz_path: Union[str, os.PathLike] = os.path.expanduser("~/.boltz")):
    """
    Converts an RDKit Mol object to a CCD object and pickles it to the CCD directory.

    Parameters
    ----------
    resname : str
        Residue name (identifier for the molecule).cache_mols_from_sdf
    mol : Chem.Mol
        RDKit molecule to convert.
    boltz_path : str or Path, optional
        Path to boltz CCD directory (default ~/.boltz).

    Example
    -------
    >>> from rdkit import Chem
    >>> from boltz_tools.helpers.conformers import mol_to_ccd
    >>> smiles = "CCO"
    >>> mol = Chem.MolFromSmiles(smiles)
    >>> mol_to_ccd("CCO", mol)
    # This will save ~/.boltz/mols/CCO.pkl
    """
    if mol is None:
        raise Exception("Mol can not be null")

    mol = Chem.RemoveHs(mol)
    mol.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(mol)
    rdmolops.AssignStereochemistryFrom3D(mol)

    mol = _prepare_mol(mol)

    parsedResidue = parse_ccd_residue(resname, mol, 0)

    add_pickled_prop(mol, 'MOL_NAME', resname)
    add_pickled_prop(mol, 'symmetries', mol.GetSubstructMatches(mol, uniquify=False))
    add_pickled_prop(mol, 'pb_edge_index', extract_constraints(parsedResidue.rdkit_bounds_constraints, 'atom_idxs', transpose=True))
    add_pickled_prop(mol, 'pb_lower_bounds', extract_constraints(parsedResidue.rdkit_bounds_constraints, 'lower_bound'))
    add_pickled_prop(mol, 'pb_upper_bounds', extract_constraints(parsedResidue.rdkit_bounds_constraints, 'upper_bound'))
    add_pickled_prop(mol, 'pb_bond_mask', extract_constraints(parsedResidue.rdkit_bounds_constraints, 'is_bond'))
    add_pickled_prop(mol, 'pb_angle_mask', extract_constraints(parsedResidue.rdkit_bounds_constraints, 'is_angle'))
    add_pickled_prop(mol, 'chiral_atom_index', extract_constraints(parsedResidue.chiral_atom_constraints, 'atom_idxs', transpose=True))
    add_pickled_prop(mol, 'chiral_check_mask', extract_constraints(parsedResidue.chiral_atom_constraints, 'is_reference'))
    add_pickled_prop(mol, 'chiral_atom_orientations', extract_constraints(parsedResidue.chiral_atom_constraints, 'is_r'))
    add_pickled_prop(mol, 'stereo_bond_index', extract_constraints(parsedResidue.stereo_bond_constraints, 'atom_idxs', transpose=True))
    add_pickled_prop(mol, 'stereo_check_mask', extract_constraints(parsedResidue.stereo_bond_constraints, 'is_check'))
    add_pickled_prop(mol, 'stereo_bond_orientations', extract_constraints(parsedResidue.stereo_bond_constraints, 'is_e'))
    add_pickled_prop(mol, 'aromatic_5_ring_index', extract_constraints(parsedResidue.planar_ring_5_constraints, 'atom_idxs', transpose=True))
    add_pickled_prop(mol, 'aromatic_6_ring_index', extract_constraints(parsedResidue.planar_ring_6_constraints, 'atom_idxs', transpose=True))
    add_pickled_prop(mol, 'planar_double_bond_index', extract_constraints(parsedResidue.planar_bond_constraints, 'atom_idxs', transpose=True))

    Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)
    mols_dir = Path(boltz_path) / 'mols'
    mols_dir.mkdir(parents=True, exist_ok=True)
    
    ccd_logger.debug("Checking atom properties just before pickle:")
    atoms = list(mol.GetAtoms())
    ccd_logger.debug("Total atoms:", len(atoms))
    missing = []
    for atom in atoms:
        if not atom.HasProp("name"):
            missing.append(atom.GetIdx())
        # show first few atoms for inspection
        ccd_logger.debug(f"  idx {atom.GetIdx()} props: {atom.GetPropsAsDict()}")

    if missing:
        ccd_logger.warning(f"{len(missing)} atoms are missing 'name' property. Example indices: {missing[:10]}")
    else:
        ccd_logger.debug("All atoms have 'name' property.")
        
    with open(mols_dir / f'{resname}.pkl', 'wb') as f:
        pickle.dump(mol, f)

def cache_mols_from_sdf(file_path: str, property_id: str, on_conflict: str = 'overwrite', cache: str = "~/.boltz/"):
    """
    Reads molecules from SDF, assigns atom names, and pickles each molecule.

    Parameters
    ----------
    file_path : str
        Path to the SDF file.
    property_id : str
        Property name to use as molecule ID.
    on_conflict : {'use_cache', 'overwrite'}
        Behavior when CCD already exists in cache.
    cache : str, optional
        Path to cache directory (default ~/.boltz/).
    """
    file_path = os.path.expanduser(file_path)
    cache = os.path.expanduser(cache)
    mols_dir = os.path.join(cache, "mols")
    os.makedirs(mols_dir, exist_ok=True)

    if not os.path.exists(cache):
        ccd_logger.info(f"Cache path {cache} does not exist. Downloading default cache...")
        command.download_cache(cache)

    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"Input file not found: {file_path}")

    # Pre-check for duplicate IDs
    ids = []
    suppl = Chem.SDMolSupplier(file_path)
    for mol in suppl:
        if mol is not None and mol.HasProp(property_id):
            ids.append(mol.GetProp(property_id))
    duplicates = [i for i, c in Counter(ids).items() if c > 1]
    if duplicates:
        raise ValueError(f"Duplicate molecule IDs found in SDF: {duplicates}")

    # Existing cache
    cached_ids = {f.stem for f in Path(mols_dir).glob("*.pkl")}
    overlap_ids = set(ids) & cached_ids
    if overlap_ids:
        ccd_logger.info(f"Found {len(overlap_ids)} overlapping IDs with cache: {sorted(list(overlap_ids))[:10]}...")

    # Re-iterate to process molecules
    suppl = Chem.SDMolSupplier(file_path)
    n_success, n_fail, n_skipped = 0, 0, 0
    for mol in suppl:
        if mol is None:
            continue
        mol_id = mol.GetProp(property_id)
        mol_file = Path(mols_dir) / f"{mol_id}.pkl"

        if mol_file.exists():
            if on_conflict == "use_cache":
                ccd_logger.info(f"Using cached CCD for {mol_id}")
                n_skipped += 1
                continue
            elif on_conflict == "overwrite":
                ccd_logger.info(f"Overwriting cached CCD for {mol_id}")
            else:
                raise ValueError(f"Invalid on_conflict mode: {on_conflict}")

        try:
            mol_to_ccd(mol_id, mol)
            n_success += 1
        except Exception as e:
            ccd_logger.error(f"Failed to process ID {mol_id}: {e}")
            n_fail += 1

    ccd_logger.info(f"Caching complete — Success: {n_success}, Failed: {n_fail}, Skipped: {n_skipped}")

if __name__ == "__main__":
    print("Testing mol_to_ccd with ethanol (CCO)...")
    try:
        smiles = "CCO"
        mol = Chem.MolFromSmiles(smiles)
        print(Chem.MolToMolBlock(mol))
        print("Converting molecule to CCD...")
        mol = _prepare_mol(mol)
        mol_to_ccd("CCO", mol)
        print("Checking if CCD file was created...")
        out_path = Path(os.path.expanduser("~/.boltz/mols/CCO.pkl"))
        if out_path.exists():
            print(f"Success: CCD file created at {out_path}")
        else:
            print(f"Failure: CCD file not found at {out_path}")
    except Exception as e:
        print(f"Error during test: {e}")
