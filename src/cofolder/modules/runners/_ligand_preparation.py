"""Boltz-facing ligand conformer and CCD-cache preparation."""

from __future__ import annotations

import fcntl
import hashlib
import json
import logging
import os
from pathlib import Path

from rdkit import Chem

from cofolder.modules.entities import ligand
from cofolder.modules.utils import read

# Preserve the established logger identity during the behavior-preserving move.
logger = logging.getLogger("cofolder.modules.entities.ligand")

_CCD_ALPHABET = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
_CCD_SPACE = len(_CCD_ALPHABET) ** 5


def _base36_fixed5(value: int) -> str:
    """Encode an integer as a stable five-character base36 token."""

    value = value % _CCD_SPACE
    chars = []
    for _ in range(5):
        value, remainder = divmod(value, len(_CCD_ALPHABET))
        chars.append(_CCD_ALPHABET[remainder])
    return "".join(reversed(chars))


def _canonical_smiles(smiles: str) -> str:
    """Return canonical SMILES when possible, otherwise the stripped input."""

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    except Exception:
        pass
    return str(smiles).strip()


def _allocate_ccd_resname(smiles: str, boltz_path: Path) -> str:
    """Allocate a stable five-character CCD ID in a shared cache."""

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
                with open(map_path, "r", encoding="utf-8") as handle:
                    loaded = json.load(handle)
                if isinstance(loaded, dict):
                    mapping["smiles_to_resname"] = dict(
                        loaded.get("smiles_to_resname", {})
                    )
                    mapping["resname_to_smiles"] = dict(
                        loaded.get("resname_to_smiles", {})
                    )
            except Exception:
                logger.warning(
                    "Failed reading CCD map at %s. Rebuilding mapping in memory.",
                    map_path,
                )

        smiles_to_resname = mapping["smiles_to_resname"]
        resname_to_smiles = mapping["resname_to_smiles"]
        existing = smiles_to_resname.get(canonical)
        if existing:
            return existing

        used_ids = set(resname_to_smiles)
        used_ids.update(path.stem for path in mols_dir.glob("*.pkl"))

        seed = int(hashlib.sha1(canonical.encode("utf-8")).hexdigest(), 16)
        seed %= _CCD_SPACE
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
        with open(tmp_path, "w", encoding="utf-8") as handle:
            json.dump(mapping, handle, indent=2, sort_keys=True)
        os.replace(tmp_path, map_path)

        return chosen


def prepare_ligand_conformers(
    sys_obj,
    *,
    cache_path: str | Path,
    wrk_dir: str | Path,
    conformers: str | None = None,
    sdf_file: str | Path | None = None,
    logger: logging.Logger,
) -> dict[str, str]:
    """Prepare Boltz ligand conformers and update the runner system in place."""

    logger.info("Creating conformer: %s", conformers)
    ccd_map: dict[str, str] = {}

    sequences = sys_obj.find_value(key="sequences")
    if not sequences:
        logger.warning("No sequences found in system, skipping conformer handling.")
        return ccd_map

    boltz_path = Path(cache_path).expanduser()
    ligand_idx = -1
    for seq_idx, sequence in enumerate(sequences):
        ligand_info = sequence.get("ligand")
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

        if sdf_file:
            sdf_file_path = Path(sdf_file)
        else:
            sdf_file_path = Path(wrk_dir) / f"_ccd_{ligand_id}.sdf"

        if conformers == "sdf" and sdf_file:
            logger.debug(
                "Using provided SDF file for ligand %s: %s",
                ligand_id,
                sdf_file_path,
            )
        else:
            ligand.smiles_to_sdf(
                data=smiles_value,
                output_sdf_path=str(sdf_file_path),
            )
            if conformers == "2D":
                ligand.generate_2d_conformers(str(sdf_file_path))
            elif conformers == "3D":
                ligand.generate_3d_conformers(str(sdf_file_path))

        mols = read.read_sdf(str(sdf_file_path))
        if not mols:
            raise ValueError(f"No valid molecules found in SDF: {sdf_file_path}")

        mol = mols[ligand_idx] if conformers == "sdf" else mols[0]
        resname = _allocate_ccd_resname(smiles_value, boltz_path)
        try:
            ligand.mol_to_ccd(resname, mol, boltz_path=boltz_path)
            logger.info("Saved CCD for %s to %s/mols/", resname, boltz_path)
        except Exception as exc:
            logger.error(
                "Failed to convert molecule '%s' to CCD: %s",
                resname,
                exc,
            )

        sys_obj.update_system(
            resname,
            path=["sequences", seq_idx, "ligand", "ccd"],
        )
        sys_obj.delete_system_key(
            path=["sequences", seq_idx, "ligand"],
            keys_to_delete=["smiles"],
        )
        ccd_map[ligand_id] = resname

        if (
            sdf_file is None
            and sdf_file_path.exists()
            and not logger.isEnabledFor(logging.DEBUG)
        ):
            try:
                os.remove(sdf_file_path)
                logger.debug("Temporary SDF removed: %s", sdf_file_path)
            except Exception as exc:
                logger.warning("Failed to remove temporary SDF: %s", exc)

    return ccd_map


__all__ = ["prepare_ligand_conformers"]
