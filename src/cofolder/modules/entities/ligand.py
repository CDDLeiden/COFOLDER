"""Ligand processing and CCD cache management for Boltz predictions.

This module provides utilities for processing molecular structures, generating
conformers, and managing the Chemical Component Dictionary (CCD) cache used
by Boltz for ligand predictions.
"""
import importlib
import logging
import os
import pickle
from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, List, Optional, Union

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor, rdmolops

from cofolder.modules.utils import read, write

if TYPE_CHECKING:
    from cofolder.modules.input.ligand import LigandSourceIdentity, NormalizedLigand

logger = logging.getLogger(__name__)

def _load_parse_ccd_residue():
    """Load the Boltz-only CCD parser at the operation boundary."""
    try:
        module = importlib.import_module("boltz.data.parse.mmcif_with_constraints")
    except (ImportError, ModuleNotFoundError) as exc:
        raise ImportError(
            "Boltz CCD conversion requires a compatible Boltz backend. Install the "
            "matching COFOLDER backend extra, for example "
            "`pip install \"cofolder[boltz2]\"`."
        ) from exc
    return module.parse_ccd_residue


def sanitize_mol_id(mol_id: str) -> str:
    """Ensure molecule ID complies with CCD naming rules (max 5 characters).

    Parameters
    ----------
    mol_id : str
        Original molecule identifier.

    Returns
    -------
    str
        Sanitized molecule ID, truncated to 5 characters if necessary.

    Warnings
    --------
    Logs a warning if the ID is truncated.
    """
    if len(mol_id) > 5:
        truncated = mol_id[:5]
        logger.warning(f"[WARNING] Molecule ID '{mol_id}' is longer than 5 characters. "
              f"Truncating to '{truncated}' to comply with CCD naming rules.")
        return truncated
    return mol_id

def add_pickled_prop(mol: Chem.Mol, name: str, value):
    """Pickle and attach a property to an RDKit molecule.

    Serializes a value using pickle and stores it as a hex-encoded
    string property on the molecule. Used for storing complex data
    structures in CCD format.

    Parameters
    ----------
    mol : Chem.Mol
        RDKit molecule object to modify.
    name : str
        Property name to use for storage.
    value : any
        Value to pickle and attach (can be any picklable Python object).
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
    """Assign systematic atom names to all atoms in the molecule.

    Generates atom names in the format ElementSymbol + Index (e.g., C1, O2).

    Parameters
    ----------
    mol : Chem.Mol
        RDKit molecule object.

    Returns
    -------
    Chem.Mol
        The same molecule with 'name' property set on all atoms.
    """
    for i, atom in enumerate(mol.GetAtoms()):
        atom.SetProp("name", f"{atom.GetSymbol()}{i+1}".upper())
    return mol

def generate_2d_conformers(sdf_in: str, sdf_out: Optional[str] = None):
    """Generate 2D coordinate layouts for molecules in an SDF file.

    Computes 2D coordinates suitable for visualization and quick analysis.
    Useful for preparing ligands that only have SMILES or 1D representations.

    Parameters
    ----------
    sdf_in : str
        Path to input SDF file.
    sdf_out : str, optional
        Path to output SDF file. If None, overwrites the input file.

    Notes
    -----
    Uses RDKit's Compute2DCoords for layout generation.
    Invalid molecules (None) are skipped.
    """
    sdf_out = sdf_out or sdf_in
    mols = read.read_sdf(sdf_in)

    for mol in mols:
        rdDepictor.Compute2DCoords(mol)

    write.write_sdf(mols, sdf_out)
    logger.info(f"Generated 2D conformers in {sdf_out}")


def generate_3d_conformers(sdf_in: str, sdf_out: Optional[str] = None):
    """Generate 3D conformers for molecules using ETKDG and UFF optimization.

    Creates physically realistic 3D conformations by:
    1. Adding hydrogens
    2. Embedding with ETKDG (Experimental-Torsion Knowledge Distance Geometry)
    3. Optimizing geometry with UFF (Universal Force Field)

    Parameters
    ----------
    sdf_in : str
        Path to input SDF file.
    sdf_out : str, optional
        Path to output SDF file. If None, overwrites the input file.

    Warnings
    --------
    Logs a warning if 3D conformer generation fails for any molecule.
    Failed molecules are skipped and not written to output.

    Notes
    -----
    This is computationally more expensive than 2D conformer generation
    but produces geometries suitable for docking and structure prediction.
    """
    sdf_out = sdf_out or sdf_in
    mols = read.read_sdf(sdf_in)
    mols_3d = []

    for mol in mols:
        try:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            AllChem.UFFOptimizeMolecule(mol)
            mols_3d.append(mol)
        except Exception as e:
            logger.warning(f"3D conformer generation failed for molecule: {e}")

    write.write_sdf(mols_3d, sdf_out)
    logger.info(f"Generated 3D conformers in {sdf_out}")

@dataclass(frozen=True, slots=True)
class RawMolblockRecord:
    index: int
    source_record_id: str
    molblock: str
    title: str | None
    properties: Mapping[str, str]
    first_line: int


def _raw_molblock_record(text: str, *, index: int, first_line: int) -> RawMolblockRecord:
    import re

    lines = text.splitlines()
    title = lines[0].strip() if lines and lines[0].strip() else None
    properties: dict[str, str] = {}
    property_header = re.compile(r"^>\s*<([^>]+)>")
    position = 0
    while position < len(lines):
        match = property_header.match(lines[position].strip())
        if match:
            position += 1
            values: list[str] = []
            while position < len(lines) and lines[position].strip():
                values.append(lines[position])
                position += 1
            properties[match.group(1)] = "\n".join(values).strip()
        position += 1
    return RawMolblockRecord(
        index=index,
        source_record_id=f"record_{index:06d}",
        molblock=text.rstrip("\r\n"),
        title=title,
        properties=properties,
        first_line=first_line,
    )


def iter_sdf_records(path: str | Path) -> Iterator[RawMolblockRecord]:
    """Yield every delimiter-defined SDF record without asking RDKit to filter it."""
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines(
        keepends=True
    )
    buffered: list[str] = []
    first_line = 1
    index = 0
    for line_number, line in enumerate(lines, 1):
        if line.strip() == "$$$$":
            index += 1
            yield _raw_molblock_record(
                "".join(buffered), index=index, first_line=first_line
            )
            buffered = []
            first_line = line_number + 1
        else:
            buffered.append(line)
    if "".join(buffered).strip():
        index += 1
        yield _raw_molblock_record(
            "".join(buffered), index=index, first_line=first_line
        )


def read_molblock_record(path: str | Path) -> RawMolblockRecord:
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    return _raw_molblock_record(text, index=1, first_line=1)


def parse_molblock(
    record: RawMolblockRecord,
    *,
    source: "LigandSourceIdentity",
    source_path: str | Path,
) -> "NormalizedLigand":
    """Parse one raw molblock, preserving parse and sanitization as separate failures."""
    from cofolder.modules.input.compound_library import MolblockRecordParseError
    from cofolder.modules.input.ligand import validate_molecule

    structure_lines: list[str] = []
    for line in record.molblock.splitlines():
        structure_lines.append(line)
        if line.startswith("M  END"):
            break
    structure = "\n".join(structure_lines)
    try:
        mol = Chem.MolFromMolBlock(
            structure, sanitize=False, removeHs=False, strictParsing=True
        )
    except Exception as exc:
        mol = None
        parse_exc = exc
    else:
        parse_exc = None
    if mol is None:
        raise MolblockRecordParseError(
            f"RDKit could not parse molblock record {record.source_record_id}.",
            source_path=source_path,
            field_path=(record.index,),
            line=record.first_line,
            entity_id=source.entity_id,
            chain_id=source.chain_ids[0] if source.chain_ids else None,
            source_record_id=source.source_record_id,
        ) from parse_exc
    return validate_molecule(
        mol,
        source=source,
        source_path=Path(source_path),
        field_path=(record.index,),
        line=record.first_line,
    )


def iterate_sdf_records(sdf_path: str, id_property: str):
    """Compatibility iterator that now retains record indexes while yielding valid molblocks."""
    for record in iter_sdf_records(sdf_path):
        source_id = record.title if id_property == "_Name" else record.properties.get(id_property)
        mol_id = source_id or f"mol_{record.index}"
        try:
            mol = Chem.MolFromMolBlock(record.molblock, removeHs=False)
        except Exception:
            mol = None
        if mol is not None:
            yield record.index, mol_id, Chem.MolToMolBlock(mol)

def smiles_to_sdf(
    data: Union[str, List[str]],
    smiles_col: Optional[str] = None,
    output_sdf_path: Optional[str] = None,
    property_cols: Optional[Union[dict, List[dict], str, List[str]]] = None,
) -> None:
    """
    Convert either a single SMILES / list of SMILES or a CSV file to SDF format.

    Parameters
    ----------
    data : str or list of str
        Single SMILES string, list of SMILES, or path to a CSV file.
    smiles_col : str, optional
        Column name containing SMILES. Required if `data` is a CSV file.
    output_sdf_path : str, optional
        Path to output SDF file. Defaults to "output.sdf" for SMILES or
        CSV file base name for CSV input.
    property_cols : dict, list of dict, str, or list of str, optional
        For SMILES: dict mapping property names to values (or list for multiple SMILES)
        For CSV: column names to preserve as SDF properties
    """

    # --- Case 2: CSV input ---
    if isinstance(data, str) and os.path.exists(data):
        if smiles_col is None:
            raise ValueError("smiles_col must be provided for CSV input")

        df = pd.read_csv(data)
        if smiles_col not in df.columns:
            logger.error(f"SMILES column '{smiles_col}' not found in {data}")
            return

        if output_sdf_path is None:
            output_sdf_path = os.path.splitext(data)[0] + ".sdf"

        if property_cols is None:
            property_cols_list = []
        elif isinstance(property_cols, str):
            property_cols_list = [property_cols]
        else:
            property_cols_list = list(property_cols)

        for col in property_cols_list:
            if col not in df.columns:
                logger.warning(f"Property column '{col}' not found in {data}. Skipping.")
        property_cols_list = [c for c in property_cols_list if c in df.columns]

        mols = mols_from_dataframe(df, smiles_col, property_cols_list)
        write.write_sdf(mols, output_sdf_path)
        logger.info(f"Wrote {len(mols)} molecules to {output_sdf_path}")
        return

    # --- Case 1: SMILES input ---
    if isinstance(data, str) or (isinstance(data, list) and all(isinstance(s, str) for s in data)):
        smiles_list = [data] if isinstance(data, str) else list(data)
        if output_sdf_path is None:
            output_sdf_path = "output.sdf"

        # Normalize property_cols for SMILES
        if property_cols is None:
            props_list = [{}] * len(smiles_list)
        elif isinstance(property_cols, dict):
            props_list = [property_cols] * len(smiles_list)
        elif isinstance(property_cols, list) and all(isinstance(p, dict) for p in property_cols):
            props_list = property_cols
            if len(props_list) < len(smiles_list):
                props_list.extend([{}] * (len(smiles_list) - len(props_list)))
        else:
            # Single list of property names (values will be None)
            props_list = [{str(k): None for k in property_cols}] * len(smiles_list)

        mols = mols_from_smiles(smiles_list, props_list)
        write.write_sdf(mols, output_sdf_path)
        logger.info(f"Wrote {len(mols)} molecules to {output_sdf_path}")
        return

    raise ValueError("Invalid input for smiles_or_csv_to_sdf: must be SMILES string, list of SMILES, or CSV file path.")


def csv_to_sdf(
    csv_path: str,
    smiles_col: str,
    output_sdf_path: str | None = None,
    property_cols: Optional[Union[dict, List[dict], str, List[str]]] = None,
) -> None:
    """Compatibility wrapper that preserves the older CSV-specific API."""
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"CSV file not found: {csv_path}")

    return smiles_to_sdf(
        data=csv_path,
        smiles_col=smiles_col,
        output_sdf_path=output_sdf_path,
        property_cols=property_cols,
    )

# Adapted from jacktday/boltztools (MIT License).
# See THIRD_PARTY_SOFTWARE.md for attribution and full license text.
def mol_to_ccd(resname: str, mol: Chem.Mol, boltz_path: Union[str, os.PathLike] = os.path.expanduser("~/.boltz")):
    """Convert an RDKit molecule to Boltz CCD format and cache it.

    Processes an RDKit molecule by:
    1. Sanitizing and preparing the structure
    2. Parsing geometric and stereochemical constraints
    3. Pickling all properties for Boltz compatibility
    4. Saving to the CCD cache directory

    Parameters
    ----------
    resname : str
        Residue name (identifier) for the molecule. Must be 5 characters or less.
    mol : Chem.Mol
        RDKit molecule object to convert.
    boltz_path : str or Path, default="~/.boltz"
        Path to Boltz cache directory where CCD files are stored.

    Raises
    ------
    Exception
        If the molecule is None.

    Examples
    --------
    >>> from rdkit import Chem
    >>> from cofolder.modules.entities.ligand import mol_to_ccd
    >>> smiles = "CCO"
    >>> mol = Chem.MolFromSmiles(smiles)
    >>> mol_to_ccd("ETH", mol)
    # Saves to ~/.boltz/mols/ETH.pkl

    Notes
    -----
    This function contains adapted logic from jacktday/boltztools (MIT).

    The CCD format includes:
    - Atomic symmetries
    - Bond length and angle constraints
    - Chiral center definitions
    - Stereochemical bond configurations
    - Aromatic ring planarity constraints
    """
    if mol is None:
        raise Exception("Mol can not be null")

    mol = Chem.RemoveHs(mol)
    mol.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(mol)
    rdmolops.AssignStereochemistryFrom3D(mol)

    mol = _prepare_mol(mol)

    parse_ccd_residue = _load_parse_ccd_residue()
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
    
    logger.debug("Checking atom properties just before pickle:")
    atoms = list(mol.GetAtoms())
    logger.debug("Total atoms: %d", len(atoms))
    missing = []
    for atom in atoms:
        if not atom.HasProp("name"):
            missing.append(atom.GetIdx())
        # show first few atoms for inspection
        logger.debug(f"  idx {atom.GetIdx()} props: {atom.GetPropsAsDict()}")

    if missing:
        logger.warning(f"{len(missing)} atoms are missing 'name' property. Example indices: {missing[:10]}")
    else:
        logger.debug("All atoms have 'name' property.")
        
    write.write_pickle(mol, mols_dir / f"{resname}.pkl")

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

def mols_from_smiles(
    smiles_list,
    props_list=None,
):
    """Create RDKit molecules from SMILES with optional properties.

    Parameters
    ----------
    smiles_list : list[str]
        SMILES strings.
    props_list : list[dict], optional
        Per-molecule property dictionaries.

    Returns
    -------
    list[Chem.Mol]
        Valid RDKit molecules.
    """
    mols = []

    if props_list is None:
        props_list = [{}] * len(smiles_list)

    for i, smi in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            logger.warning("Invalid SMILES at index %d: %s", i, smi)
            continue

        for k, v in props_list[i].items():
            if v is not None:
                mol.SetProp(str(k), str(v))

        mols.append(mol)

    return mols

def mols_from_dataframe(
    df: pd.DataFrame,
    smiles_col: str,
    property_cols: list[str] = None
) -> list:
    """Convert a DataFrame of SMILES and optional properties to RDKit molecules.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe containing SMILES strings.
    smiles_col : str
        Column name containing SMILES.
    property_cols : list of str, optional
        Columns to store as molecule properties.

    Returns
    -------
    list of Chem.Mol
        RDKit molecules.
    """
    mols = []
    property_cols = property_cols or []

    for idx, row in df.iterrows():
        smi = row[smiles_col]
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            logger.warning("Invalid SMILES at row %s: %s", idx, smi)
            continue
        for col in property_cols:
            val = row[col]
            if pd.notnull(val):
                mol.SetProp(str(col), str(val))
        mols.append(mol)

    return mols
