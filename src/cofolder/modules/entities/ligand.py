"""Ligand processing and CCD cache management for Boltz predictions.

This module provides utilities for processing molecular structures, generating
conformers, and managing the Chemical Component Dictionary (CCD) cache used
by Boltz for ligand predictions.
"""
import fcntl
import hashlib
import importlib
import json
import logging
import os
import pickle
from collections import Counter
from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, List, Optional, Union

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor, rdmolops

from cofolder.modules.entities import ligand
from cofolder.modules.input import command
from cofolder.modules.utils import read, write

if TYPE_CHECKING:
    from cofolder.modules.input.ligand import LigandSourceIdentity, NormalizedLigand

logger = logging.getLogger(__name__)

_CCD_ALPHABET = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
_CCD_SPACE = len(_CCD_ALPHABET) ** 5


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


def _base36_fixed5(value: int) -> str:
    """Encode integer to 5-char base36 token."""
    value = value % _CCD_SPACE
    chars = []
    for _ in range(5):
        value, rem = divmod(value, len(_CCD_ALPHABET))
        chars.append(_CCD_ALPHABET[rem])
    return "".join(reversed(chars))


def _canonical_smiles(smiles: str) -> str:
    """Return canonical SMILES when possible; fallback to raw string."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    except Exception:
        pass
    return str(smiles).strip()


def _allocate_ccd_resname(smiles: str, boltz_path: Path) -> str:
    """Allocate stable 5-char CCD ID using shared cache map with file lock.

    This avoids collisions when multiple workers run in parallel and would
    otherwise all write to the same default ID (e.g. ``0_B``).
    """
    mols_dir = Path(boltz_path) / "mols"
    mols_dir.mkdir(parents=True, exist_ok=True)

    map_path = mols_dir / ".cofolder_ccd_map.json"
    lock_path = mols_dir / ".cofolder_ccd_map.lock"
    canonical = _canonical_smiles(smiles)

    with open(lock_path, "a+", encoding="utf-8") as lock_file:
        fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX)

        mapping = {"smiles_to_resname": {}, "resname_to_smiles": {}}
        if map_path.exists():
            try:
                with open(map_path, "r", encoding="utf-8") as fh:
                    loaded = json.load(fh)
                if isinstance(loaded, dict):
                    mapping["smiles_to_resname"] = dict(loaded.get("smiles_to_resname", {}))
                    mapping["resname_to_smiles"] = dict(loaded.get("resname_to_smiles", {}))
            except Exception:
                logger.warning("Failed reading CCD map at %s. Rebuilding mapping in memory.", map_path)

        smiles_to_resname = mapping["smiles_to_resname"]
        resname_to_smiles = mapping["resname_to_smiles"]

        # Fast path for reruns.
        existing = smiles_to_resname.get(canonical)
        if existing:
            return existing

        # Prevent accidental overwrite of any already-cached IDs.
        used_ids = set(resname_to_smiles.keys())
        used_ids.update(p.stem for p in mols_dir.glob("*.pkl"))

        seed = int(hashlib.sha1(canonical.encode("utf-8")).hexdigest(), 16) % _CCD_SPACE
        chosen = None
        for step in range(_CCD_SPACE):
            candidate = _base36_fixed5(seed + step)
            owner = resname_to_smiles.get(candidate)
            if owner == canonical:
                chosen = candidate
                break
            if owner is None and candidate not in used_ids:
                chosen = candidate
                break

        if chosen is None:
            raise RuntimeError("Unable to allocate unique 5-char CCD ID.")

        smiles_to_resname[canonical] = chosen
        resname_to_smiles[chosen] = canonical

        tmp_path = map_path.with_suffix(".json.tmp")
        with open(tmp_path, "w", encoding="utf-8") as fh:
            json.dump(mapping, fh, indent=2, sort_keys=True)
        os.replace(tmp_path, map_path)

        return chosen

def handle_conformers(
    sys_obj,
    opt_obj,
    wrk_dir: str,
    conformers: Optional[str] = None,
    sdf_file: Optional[Union[str, Path]] = None,
    logger: Optional[logging.Logger] = None,
) -> Optional[str]:
    """
    Prepare ligand conformers in SDF format, optionally generate 2D/3D conformers,
    convert the first molecule to CCD, and update the system object.

    Parameters
    ----------
    sys_obj : System
        The system object containing ligand information (must support find_value and update_system).
    opt_obj : Command
        Boltz options object (must support find_value).
    wrk_dir : str
        Working directory where temporary SDF may be written.
    conformers : {"2D", "3D", "sdf"}, optional
        Type of conformer generation requested.
    sdf_file : str or Path, optional
        Existing SDF file to use (required if conformers="sdf").
    global_seed : int, optional
        Global seed used as fallback for residue name if none found.
    logger : logging.Logger, optional
        Logger to use. Defaults to root logger if None.

    Returns
    -------
    resname : str
        mapping of ligand id -> CCD residue name..
    """
    logger.info("Creating conformer: %s", conformers)
    ccd_map = {}

    sequences = sys_obj.find_value(key="sequences")
    if not sequences:
        logger.warning("No sequences found in system, skipping conformer handling.")
        return ccd_map

    ligand_idx = -1
    for seq_idx, seq in enumerate(sequences):
        ligand_info = seq.get("ligand")
        if not ligand_info:
            logger.debug("No ligand found in sequence %d, skipping.", seq_idx)
            continue

        ligand_id = ligand_info.get("id")
        if isinstance(ligand_id, list):
            ligand_id = ligand_id[0]
        ligand_idx += 1

        smiles_value = ligand_info.get("smiles")
        if not smiles_value:
            logger.warning("Ligand '%s' has no SMILES, skipping.", ligand_id)
            continue

        # Determine SDF file path
        if sdf_file:
            sdf_file_path = Path(sdf_file)
        else:
            sdf_file_path = Path(wrk_dir) / f"_ccd_{ligand_id}.sdf"

        # Prepare ligand conformers
        if conformers == "sdf" and sdf_file:
            logger.debug("Using provided SDF file for ligand %s: %s", ligand_id, sdf_file_path)
        else:
            ligand.smiles_to_sdf(data=smiles_value, output_sdf_path=str(sdf_file_path))
            if conformers == "2D":
                ligand.generate_2d_conformers(str(sdf_file_path))
            elif conformers == "3D":
                ligand.generate_3d_conformers(str(sdf_file_path))

        # Read the molecule from SDF
        mols = read.read_sdf(str(sdf_file_path))
        if not mols:
            raise ValueError(f"No valid molecules found in SDF: {sdf_file_path}")

        if conformers == "sdf":
            mol = mols[ligand_idx]
        else:
            mol = mols[0]

        # CCD conversion
        boltz_cache = opt_obj.find_value(key='cache') or '~/.boltz'
        boltz_path = Path(boltz_cache).expanduser()
        resname = _allocate_ccd_resname(smiles_value, boltz_path)
        try:
            ligand.mol_to_ccd(resname, mol, boltz_path=boltz_path)
            logger.info("Saved CCD for %s to %s/mols/", resname, boltz_path)
        except Exception as e:
            logger.error("Failed to convert molecule '%s' to CCD: %s", resname, e)

        # Update system for each ligand ID
        sys_obj.update_system(resname, path=["sequences", seq_idx, "ligand", "ccd"])
        sys_obj.delete_system_key(path=["sequences", seq_idx, "ligand"], keys_to_delete=["smiles"])
 
        ccd_map[ligand_id] = resname

        # Clean up temporary SDF if no explicit file and not in debug mode
        if sdf_file is None and sdf_file_path.exists() and not logger.isEnabledFor(logging.DEBUG):
            try:
                os.remove(sdf_file_path)
                logger.debug("Temporary SDF removed: %s", sdf_file_path)
            except Exception as e:
                logger.warning("Failed to remove temporary SDF: %s", e)

        # Remove temporary SDF if not provided and not in debug
        if sdf_file is None and sdf_file_path.exists() and not logger.isEnabledFor(logging.DEBUG):
            try:
                os.remove(sdf_file_path)
                logger.debug("Temporary SDF removed: %s", sdf_file_path)
            except Exception as e:
                logger.warning("Failed to remove temporary SDF %s: %s", sdf_file_path, e)

    return ccd_map

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

def cache_mols_from_sdf(file_path: str, property_id: str, on_conflict: str = 'overwrite', cache: str = "~/.boltz/"):
    """Batch convert SDF molecules to CCD format and add to cache.

    Reads all molecules from an SDF file, converts each to CCD format,
    and stores them in the Boltz cache directory for use in predictions.

    Parameters
    ----------
    file_path : str
        Path to the input SDF file.
    property_id : str
        SDF property name to use as molecule identifier.
    on_conflict : {'use_cache', 'overwrite'}, default='overwrite'
        Behavior when a CCD file already exists:
        - 'use_cache': Skip processing, keep existing CCD
        - 'overwrite': Replace existing CCD with new conversion
    cache : str, default="~/.boltz/"
        Path to Boltz cache directory.

    Raises
    ------
    FileNotFoundError
        If the input SDF file does not exist.
    ValueError
        If duplicate molecule IDs are found in the SDF file, or if
        on_conflict has an invalid value.

    Notes
    -----
    - Automatically downloads default Boltz cache if it doesn't exist
    - Duplicate IDs in the SDF file cause an error before processing
    - Logs statistics on success, failure, and skipped molecules
    - Failed conversions are logged but don't stop the batch process
    """
    file_path = os.path.expanduser(file_path)
    cache = os.path.expanduser(cache)
    mols_dir = os.path.join(cache, "mols")
    os.makedirs(mols_dir, exist_ok=True)

    if not os.path.exists(cache):
        logger.info(f"Cache path {cache} does not exist. Downloading default cache...")
        command.download_cache(cache)

    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"Input file not found: {file_path}")

    # Pre-check for duplicate IDs
    ids = []
    suppl = read.read_sdf(file_path)
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
        logger.info(f"Found {len(overlap_ids)} overlapping IDs with cache: {sorted(list(overlap_ids))[:10]}...")

    # Re-iterate to process molecules
    suppl = read.read_sdf(file_path)
    n_success, n_fail, n_skipped = 0, 0, 0
    for mol in suppl:
        if mol is None:
            continue
        mol_id = mol.GetProp(property_id)
        mol_file = Path(mols_dir) / f"{mol_id}.pkl"

        if mol_file.exists():
            if on_conflict == "use_cache":
                logger.info(f"Using cached CCD for {mol_id}")
                n_skipped += 1
                continue
            elif on_conflict == "overwrite":
                logger.info(f"Overwriting cached CCD for {mol_id}")
            else:
                raise ValueError(f"Invalid on_conflict mode: {on_conflict}")

        try:
            mol_to_ccd(mol_id, mol)
            n_success += 1
        except Exception as e:
            logger.error(f"Failed to process ID {mol_id}: {e}")
            n_fail += 1

    logger.info(f"Caching complete — Success: {n_success}, Failed: {n_fail}, Skipped: {n_skipped}")

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
