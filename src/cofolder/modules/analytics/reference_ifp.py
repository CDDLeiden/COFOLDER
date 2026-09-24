"""Typed extraction and comparison of protein--ligand interaction fingerprints."""

from __future__ import annotations

import json
import math
import subprocess
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path

from Bio.Align import PairwiseAligner
from Bio.PDB import MMCIFParser, PDBParser
from Bio.SeqUtils import seq1

DEFAULT_PROLIF_INTERACTIONS = (
    "Hydrophobic",
    "HBAcceptor",
    "HBDonor",
    "PiStacking",
    "Anionic",
    "Cationic",
    "CationPi",
    "PiCation",
    "VdWContact",
)
_WATER_NAMES = frozenset({"HOH", "WAT", "DOD"})


class IFPTaxonomy(StrEnum):
    DISTANCE = "distance"
    PROLIF = "prolif"


class IFPSimilarityMetric(StrEnum):
    JACCARD = "jaccard"
    REFERENCE_COVERAGE = "reference_coverage"


class IFPComparisonStatus(StrEnum):
    COMPARABLE = "comparable"
    NOT_EVALUABLE = "not_evaluable"


class IFPMappingStatus(StrEnum):
    MAPPED = "mapped"
    UNMAPPABLE = "unmappable"


@dataclass(frozen=True, slots=True, order=True)
class ResidueIdentity:
    chain_id: str
    residue_number: int
    insertion_code: str = ""
    residue_name: str | None = field(default=None, compare=False)


@dataclass(frozen=True, slots=True, order=True)
class LigandIdentity:
    chain_id: str
    residue_number: int
    insertion_code: str = ""
    residue_name: str | None = field(default=None, compare=False)


@dataclass(frozen=True, slots=True, order=True)
class InteractionKey:
    receptor: ResidueIdentity
    interaction_type: str

    def __str__(self) -> str:
        suffix = f"{self.receptor.residue_number}{self.receptor.insertion_code}"
        return f"{self.receptor.chain_id}:{suffix}:{self.interaction_type}"

    @classmethod
    def parse(cls, value: str) -> InteractionKey:
        try:
            chain_id, residue_token, interaction_type = str(value).strip().split(":", 2)
        except ValueError as exc:
            raise ReferenceIFPInputError(
                f"Invalid interaction key {value!r}; expected CHAIN:RESNUM[ICODE]:TYPE."
            ) from exc
        index = 0
        if residue_token.startswith("-"):
            index = 1
        while index < len(residue_token) and residue_token[index].isdigit():
            index += 1
        number_token = residue_token[:index]
        insertion_code = residue_token[index:]
        if (
            not chain_id
            or not number_token
            or len(insertion_code) > 1
            or not interaction_type
        ):
            raise ReferenceIFPInputError(
                f"Invalid interaction key {value!r}; expected CHAIN:RESNUM[ICODE]:TYPE."
            )
        return cls(
            receptor=ResidueIdentity(chain_id, int(number_token), insertion_code),
            interaction_type=_normalize_interaction_type(interaction_type),
        )


@dataclass(frozen=True, slots=True)
class AtomIdentity:
    """Stable atom identity reported for one typed interaction occurrence."""

    chain_id: str
    residue_number: int
    insertion_code: str
    residue_name: str | None
    atom_name: str
    element: str | None = None
    atom_serial: int | None = None
    source_index: int | None = None

    def to_dict(self) -> dict[str, object]:
        return {
            "chain_id": self.chain_id,
            "residue_number": self.residue_number,
            "insertion_code": self.insertion_code,
            "residue_name": self.residue_name,
            "atom_name": self.atom_name,
            "element": self.element,
            "atom_serial": self.atom_serial,
            "source_index": self.source_index,
        }


@dataclass(frozen=True, slots=True)
class InteractionGeometry:
    """Named ProLIF geometry value with an explicit physical unit."""

    name: str
    value: float
    unit: str

    def __post_init__(self) -> None:
        if not self.name or not math.isfinite(self.value):
            raise ReferenceIFPInputError(
                "Interaction geometry must be named and finite."
            )

    def to_dict(self) -> dict[str, object]:
        return {"name": self.name, "value": self.value, "unit": self.unit}


@dataclass(frozen=True, slots=True)
class InteractionEvent:
    """One atom-level occurrence contributing to a deduplicated interaction key."""

    interaction: InteractionKey
    ligand_atoms: tuple[AtomIdentity, ...]
    protein_atoms: tuple[AtomIdentity, ...]
    ligand_role: str
    protein_role: str
    geometry: tuple[InteractionGeometry, ...] = ()

    def to_dict(self) -> dict[str, object]:
        return {
            "interaction_key": str(self.interaction),
            "interaction_type": self.interaction.interaction_type,
            "ligand_role": self.ligand_role,
            "protein_role": self.protein_role,
            "ligand_atoms": [atom.to_dict() for atom in self.ligand_atoms],
            "protein_atoms": [atom.to_dict() for atom in self.protein_atoms],
            "geometry": [measurement.to_dict() for measurement in self.geometry],
        }


@dataclass(frozen=True, slots=True)
class LigandSelector:
    chain_id: str | None = None
    residue_number: int | None = None
    insertion_code: str = ""


@dataclass(frozen=True, slots=True)
class IFPExtractionConfig:
    taxonomy: IFPTaxonomy = IFPTaxonomy.DISTANCE
    distance_cutoff_angstrom: float = 5.0
    prolif_interactions: tuple[str, ...] = DEFAULT_PROLIF_INTERACTIONS
    prolif_timeout_seconds: float = 120.0

    def __post_init__(self) -> None:
        object.__setattr__(self, "taxonomy", IFPTaxonomy(self.taxonomy))
        if (
            not math.isfinite(self.distance_cutoff_angstrom)
            or self.distance_cutoff_angstrom <= 0
        ):
            raise ReferenceIFPInputError(
                "distance_cutoff_angstrom must be positive and finite."
            )
        if (
            not math.isfinite(self.prolif_timeout_seconds)
            or self.prolif_timeout_seconds <= 0
        ):
            raise ReferenceIFPInputError(
                "prolif_timeout_seconds must be positive and finite."
            )
        unknown = set(self.prolif_interactions) - set(DEFAULT_PROLIF_INTERACTIONS)
        if unknown:
            raise ReferenceIFPInputError(
                f"Unsupported ProLIF interactions: {sorted(unknown)}"
            )


DEFAULT_IFP_EXTRACTION_CONFIG = IFPExtractionConfig()


@dataclass(frozen=True, slots=True)
class InteractionFingerprint:
    taxonomy: IFPTaxonomy
    ligand: LigandIdentity
    receptor_chains: tuple[str, ...]
    interactions: frozenset[InteractionKey]
    source_path: Path
    events: tuple[InteractionEvent, ...] = ()

    def serialized_interactions(self) -> list[str]:
        return [str(item) for item in sorted(self.interactions)]

    def serialized_events(self) -> list[dict[str, object]]:
        """Return occurrence-level events in a deterministic JSON-safe form."""

        return [event.to_dict() for event in self.events]


@dataclass(frozen=True, slots=True)
class ReferenceIdentityMapping:
    status: IFPMappingStatus
    chain_mapping: Mapping[str, str]
    residue_mapping: Mapping[ResidueIdentity, ResidueIdentity]
    reference_ligand: LigandIdentity
    predicted_ligand: LigandIdentity
    failures: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class IFPComparison:
    status: IFPComparisonStatus
    mapping: ReferenceIdentityMapping
    similarities: Mapping[IFPSimilarityMetric, float]
    matched_interactions: tuple[InteractionKey, ...]
    missing_interactions: tuple[InteractionKey, ...]
    extra_interactions: tuple[InteractionKey, ...]
    ignored_predicted_interactions: tuple[InteractionKey, ...]


class ReferenceIFPError(Exception):
    """Base class for reference-IFP failures."""


class ReferenceIFPInputError(ReferenceIFPError, ValueError):
    """Raised for invalid reference-IFP input or configuration."""


class ReferenceEntitySelectionError(ReferenceIFPInputError):
    """Raised when a receptor or ligand cannot be selected unambiguously."""


class ReferenceIFPExtractionError(ReferenceIFPError, RuntimeError):
    """Raised when fingerprint extraction fails."""


class ReferenceIFPMappingError(ReferenceIFPError):
    """Raised for invalid mapping API usage."""


class ProLIFWorkerError(ReferenceIFPExtractionError):
    """Base class for isolated ProLIF worker failures."""


class ProLIFWorkerTimeoutError(ProLIFWorkerError):
    """Raised when the isolated ProLIF worker times out."""


class ProLIFWorkerCrashError(ProLIFWorkerError):
    """Raised when the isolated ProLIF worker exits abnormally."""


def _load_structure(path: Path):
    path = Path(path)
    if not path.exists() or not path.is_file():
        raise ReferenceIFPInputError(f"Structure file does not exist: {path}")
    if path.suffix.lower() in {".cif", ".mmcif"}:
        parser = MMCIFParser(QUIET=True)
    elif path.suffix.lower() == ".pdb":
        parser = PDBParser(QUIET=True)
    else:
        raise ReferenceIFPInputError(f"Unsupported structure format: {path}")
    try:
        return parser.get_structure(path.stem, str(path))
    except Exception as exc:
        raise ReferenceIFPExtractionError(
            f"Unable to read structure {path}: {exc}"
        ) from exc


def _residue_identity(chain_id: str, residue) -> ResidueIdentity:
    return ResidueIdentity(
        chain_id=str(chain_id),
        residue_number=int(residue.id[1]),
        insertion_code=str(residue.id[2]).strip(),
        residue_name=str(residue.get_resname()).strip().upper() or None,
    )


def _eligible_ligands(
    model, selector: LigandSelector | None
) -> list[tuple[object, object]]:
    selected = []
    for chain in model:
        if (
            selector is not None
            and selector.chain_id is not None
            and chain.id != selector.chain_id
        ):
            continue
        for residue in chain:
            if str(residue.get_resname()).strip().upper() in _WATER_NAMES:
                continue
            is_explicit_chain = selector is not None and selector.chain_id is not None
            if residue.id[0] == " " and not is_explicit_chain:
                continue
            if selector is not None and selector.residue_number is not None:
                if int(residue.id[1]) != selector.residue_number:
                    continue
                if str(residue.id[2]).strip() != selector.insertion_code.strip():
                    continue
            selected.append((chain, residue))
    return selected


def _select_ligand(model, selector: LigandSelector | None):
    candidates = _eligible_ligands(model, selector)
    if not candidates:
        raise ReferenceEntitySelectionError("ligand_not_found")
    if len(candidates) != 1:
        reason = (
            "ambiguous_reference_ligand"
            if selector is None
            else "ambiguous_ligand_selector"
        )
        raise ReferenceEntitySelectionError(reason)
    return candidates[0]


def _protein_residues(model, receptor_chains: Sequence[str] | None, ligand_chain: str):
    requested = {str(item) for item in receptor_chains} if receptor_chains else None
    found: dict[str, list[object]] = {}
    for chain in model:
        if requested is not None and chain.id not in requested:
            continue
        if requested is None and chain.id == ligand_chain:
            continue
        residues = [res for res in chain if res.id[0] == " " and "CA" in res]
        if residues:
            found[str(chain.id)] = residues
    if requested is not None and set(found) != requested:
        raise ReferenceEntitySelectionError("receptor_chain_not_found")
    if not found:
        raise ReferenceEntitySelectionError("receptor_chain_not_found")
    return found


def extract_interaction_fingerprint(
    structure_path: Path,
    *,
    ligand: LigandSelector | None = None,
    receptor_chains: Sequence[str] | None = None,
    config: IFPExtractionConfig = DEFAULT_IFP_EXTRACTION_CONFIG,
) -> InteractionFingerprint:
    """Extract a normalized fingerprint from a PDB or mmCIF complex."""

    path = Path(structure_path)
    if config.taxonomy is IFPTaxonomy.PROLIF:
        return _extract_prolif(path, ligand, receptor_chains, config)

    structure = _load_structure(path)
    model = next(iter(structure), None)
    if model is None:
        raise ReferenceIFPExtractionError(f"Structure contains no model: {path}")
    ligand_chain, ligand_residue = _select_ligand(model, ligand)
    proteins = _protein_residues(model, receptor_chains, str(ligand_chain.id))
    ligand_atoms = [
        atom
        for atom in ligand_residue
        if str(getattr(atom, "element", "")).strip().upper() != "H"
    ]
    if not ligand_atoms:
        raise ReferenceIFPExtractionError("selected_ligand_has_no_heavy_atoms")
    cutoff_squared = config.distance_cutoff_angstrom**2
    interactions: set[InteractionKey] = set()
    for chain_id, residues in proteins.items():
        for residue in residues:
            contact = any(
                float(((protein_atom.coord - ligand_atom.coord) ** 2).sum())
                <= cutoff_squared
                for protein_atom in residue
                if str(getattr(protein_atom, "element", "")).strip().upper() != "H"
                for ligand_atom in ligand_atoms
            )
            if contact:
                interactions.add(
                    InteractionKey(
                        _residue_identity(chain_id, residue), "distance_contact"
                    )
                )
    return InteractionFingerprint(
        taxonomy=IFPTaxonomy.DISTANCE,
        ligand=LigandIdentity(
            str(ligand_chain.id),
            int(ligand_residue.id[1]),
            str(ligand_residue.id[2]).strip(),
            str(ligand_residue.get_resname()).strip().upper() or None,
        ),
        receptor_chains=tuple(sorted(proteins)),
        interactions=frozenset(interactions),
        source_path=path,
    )


def _extract_prolif(path, ligand, receptor_chains, config):
    from cofolder.modules.utils._optional_dependencies import (
        require_analysis_dependency,
    )

    require_analysis_dependency(
        "MDAnalysis",
        feature="ProLIF interaction fingerprints",
    )
    require_analysis_dependency(
        "prolif",
        feature="ProLIF interaction fingerprints",
    )
    # Resolve entities in the safe parent process before importing ProLIF in the
    # worker. This produces the same explicit selection diagnostics as distance IFPs.
    model = next(iter(_load_structure(path)), None)
    if model is None:
        raise ReferenceIFPExtractionError(f"Structure contains no model: {path}")
    selected_chain, selected_residue = _select_ligand(model, ligand)
    selected_receptors = _protein_residues(
        model, receptor_chains, str(selected_chain.id)
    )
    resolved_ligand = LigandSelector(
        chain_id=str(selected_chain.id),
        residue_number=int(selected_residue.id[1]),
        insertion_code=str(selected_residue.id[2]).strip(),
    )
    request = {
        "structure_path": str(path),
        "ligand": {
            "chain_id": resolved_ligand.chain_id,
            "residue_number": resolved_ligand.residue_number,
            "insertion_code": resolved_ligand.insertion_code,
        },
        "receptor_chains": sorted(selected_receptors),
        "interactions": list(config.prolif_interactions),
    }
    try:
        completed = subprocess.run(
            [sys.executable, "-m", "cofolder.modules.analytics._prolif_worker"],
            input=json.dumps(request),
            text=True,
            capture_output=True,
            timeout=config.prolif_timeout_seconds,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        raise ProLIFWorkerTimeoutError("prolif_worker_timeout") from exc
    if completed.returncode != 0:
        raise ProLIFWorkerCrashError(
            f"prolif_worker_crashed: exit={completed.returncode}; {completed.stderr[-500:]}"
        )
    try:
        payload = json.loads(completed.stdout)
        ligand_payload = payload["ligand"]
        interactions = frozenset(
            InteractionKey.parse(item) for item in payload["interactions"]
        )
        events = tuple(
            _interaction_event_from_payload(item) for item in payload.get("events", ())
        )
    except (
        AttributeError,
        KeyError,
        TypeError,
        ValueError,
        json.JSONDecodeError,
    ) as exc:
        raise ProLIFWorkerError("prolif_worker_malformed_output") from exc
    if events and frozenset(event.interaction for event in events) != interactions:
        raise ProLIFWorkerError("prolif_worker_malformed_output")
    return InteractionFingerprint(
        taxonomy=IFPTaxonomy.PROLIF,
        ligand=LigandIdentity(**ligand_payload),
        receptor_chains=tuple(payload["receptor_chains"]),
        interactions=interactions,
        source_path=path,
        events=events,
    )


def _interaction_event_from_payload(payload: Mapping[str, object]) -> InteractionEvent:
    def atom(value: Mapping[str, object]) -> AtomIdentity:
        return AtomIdentity(
            chain_id=str(value["chain_id"]),
            residue_number=int(value["residue_number"]),
            insertion_code=str(value.get("insertion_code") or ""),
            residue_name=(
                str(value["residue_name"]) if value.get("residue_name") else None
            ),
            atom_name=str(value["atom_name"]),
            element=str(value["element"]) if value.get("element") else None,
            atom_serial=(
                int(value["atom_serial"])
                if value.get("atom_serial") is not None
                else None
            ),
            source_index=(
                int(value["source_index"])
                if value.get("source_index") is not None
                else None
            ),
        )

    geometry = tuple(
        InteractionGeometry(
            name=str(item["name"]),
            value=float(item["value"]),
            unit=str(item["unit"]),
        )
        for item in payload.get("geometry", ())
    )
    ligand_atoms = tuple(atom(item) for item in payload["ligand_atoms"])
    protein_atoms = tuple(atom(item) for item in payload["protein_atoms"])
    if not ligand_atoms or not protein_atoms:
        raise ValueError("Interaction events require ligand and protein atoms.")
    return InteractionEvent(
        interaction=InteractionKey.parse(str(payload["interaction_key"])),
        ligand_atoms=ligand_atoms,
        protein_atoms=protein_atoms,
        ligand_role=str(payload["ligand_role"]),
        protein_role=str(payload["protein_role"]),
        geometry=geometry,
    )


def _chain_data(path: Path) -> dict[str, tuple[str, list[ResidueIdentity]]]:
    model = next(iter(_load_structure(path)))
    output = {}
    for chain in model:
        residues = [res for res in chain if res.id[0] == " " and "CA" in res]
        if not residues:
            continue
        identities = [_residue_identity(chain.id, residue) for residue in residues]
        sequence = "".join(
            seq1(item.residue_name or "UNK", undef_code="X") for item in identities
        )
        output[str(chain.id)] = (sequence, identities)
    return output


def _alignment_pairs(ref_sequence, pred_sequence):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(ref_sequence, pred_sequence)[0]
    pairs = []
    coordinates = alignment.coordinates
    for index in range(coordinates.shape[1] - 1):
        ref_start, ref_end = int(coordinates[0, index]), int(coordinates[0, index + 1])
        pred_start, pred_end = int(coordinates[1, index]), int(
            coordinates[1, index + 1]
        )
        span = min(ref_end - ref_start, pred_end - pred_start)
        pairs.extend(
            (ref_start + offset, pred_start + offset) for offset in range(max(0, span))
        )
    return float(alignment.score), pairs


def map_reference_identities(
    reference: InteractionFingerprint,
    prediction: InteractionFingerprint,
    *,
    reference_structure_path: Path,
    predicted_structure_path: Path,
    chain_hints: Mapping[str, str] | None = None,
) -> ReferenceIdentityMapping:
    """Map reference receptor residues into prediction identity space."""

    failures: list[str] = []
    if reference.taxonomy != prediction.taxonomy:
        failures.append("taxonomy_mismatch")
    ref_data = _chain_data(Path(reference_structure_path))
    pred_data = _chain_data(Path(predicted_structure_path))
    hints = dict(chain_hints or {})
    chain_mapping: dict[str, str] = {}
    used_prediction: set[str] = set()
    residue_mapping: dict[ResidueIdentity, ResidueIdentity] = {}

    for ref_chain in reference.receptor_chains:
        candidates = (
            [hints[ref_chain]]
            if ref_chain in hints
            else (
                [ref_chain]
                if ref_chain in prediction.receptor_chains
                and ref_chain not in used_prediction
                else [
                    chain
                    for chain in prediction.receptor_chains
                    if chain not in used_prediction
                ]
            )
        )
        candidates = [chain for chain in candidates if chain in pred_data]
        if ref_chain not in ref_data or not candidates:
            failures.append("receptor_chain_unmapped")
            continue
        scored = []
        for pred_chain in candidates:
            score, pairs = _alignment_pairs(
                ref_data[ref_chain][0], pred_data[pred_chain][0]
            )
            scored.append((score, pred_chain, pairs))
        best_score = max(item[0] for item in scored)
        best = [item for item in scored if math.isclose(item[0], best_score)]
        if len(best) != 1:
            failures.append("ambiguous_chain_mapping")
            continue
        _, pred_chain, pairs = best[0]
        if pred_chain in used_prediction:
            failures.append("prediction_chain_reused")
            continue
        used_prediction.add(pred_chain)
        chain_mapping[ref_chain] = pred_chain
        ref_residues, pred_residues = ref_data[ref_chain][1], pred_data[pred_chain][1]
        residue_mapping.update({ref_residues[i]: pred_residues[j] for i, j in pairs})

    for interaction in reference.interactions:
        if interaction.receptor not in residue_mapping:
            failures.append("interaction_residue_unmapped")
            break
    failures = list(dict.fromkeys(failures))
    return ReferenceIdentityMapping(
        status=IFPMappingStatus.UNMAPPABLE if failures else IFPMappingStatus.MAPPED,
        chain_mapping=chain_mapping,
        residue_mapping=residue_mapping,
        reference_ligand=reference.ligand,
        predicted_ligand=prediction.ligand,
        failures=tuple(failures),
    )


def compare_interaction_fingerprints(
    reference: InteractionFingerprint,
    prediction: InteractionFingerprint,
    mapping: ReferenceIdentityMapping,
) -> IFPComparison:
    """Compare two fingerprints after reference-to-prediction identity mapping."""

    if (
        reference.taxonomy != prediction.taxonomy
        or mapping.status is IFPMappingStatus.UNMAPPABLE
    ):
        return IFPComparison(
            IFPComparisonStatus.NOT_EVALUABLE,
            mapping,
            {},
            (),
            (),
            (),
            tuple(sorted(prediction.interactions)),
        )
    if not reference.interactions:
        return IFPComparison(
            IFPComparisonStatus.NOT_EVALUABLE,
            mapping,
            {},
            (),
            (),
            tuple(sorted(prediction.interactions)),
            (),
        )

    mapped: dict[InteractionKey, InteractionKey] = {
        item: InteractionKey(
            mapping.residue_mapping[item.receptor], item.interaction_type
        )
        for item in reference.interactions
    }
    mapped_residues = frozenset(mapping.residue_mapping.values())
    relevant_prediction = {
        item for item in prediction.interactions if item.receptor in mapped_residues
    }
    ignored = prediction.interactions - relevant_prediction
    mapped_reference = set(mapped.values())
    matched_reference = tuple(
        sorted(ref for ref, pred in mapped.items() if pred in relevant_prediction)
    )
    missing_reference = tuple(
        sorted(set(reference.interactions) - set(matched_reference))
    )
    extras = tuple(sorted(relevant_prediction - mapped_reference))
    intersection_size = len(mapped_reference & relevant_prediction)
    union_size = len(mapped_reference | relevant_prediction)
    similarities = {
        IFPSimilarityMetric.JACCARD: intersection_size / union_size,
        IFPSimilarityMetric.REFERENCE_COVERAGE: intersection_size
        / len(mapped_reference),
    }
    return IFPComparison(
        IFPComparisonStatus.COMPARABLE,
        mapping,
        similarities,
        matched_reference,
        missing_reference,
        extras,
        tuple(sorted(ignored)),
    )


def _normalize_interaction_type(value: str) -> str:
    token = "".join(
        character.lower() if character.isalnum() else "_" for character in str(value)
    )
    while "__" in token:
        token = token.replace("__", "_")
    canonical = (
        "hydrophobic",
        "hb_acceptor",
        "hb_donor",
        "pi_stacking",
        "anionic",
        "cationic",
        "cation_pi",
        "pi_cation",
        "vdw_contact",
    )
    aliases = {
        name.lower(): normalized
        for name, normalized in zip(DEFAULT_PROLIF_INTERACTIONS, canonical, strict=True)
    }
    aliases["distance_contact"] = "distance_contact"
    normalized = aliases.get(token.replace("_", ""), token.strip("_"))
    if not normalized:
        raise ReferenceIFPInputError("Interaction type must not be empty.")
    return normalized
