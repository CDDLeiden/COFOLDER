"""Tutorial helpers for the MAPK14 structure-gated Oracle walkthrough.

This module deliberately lives in ``tutorials``.  It demonstrates how to compose
the public Oracle callback and interaction-fingerprint APIs without making the
target-specific MAPK14 policy part of COFOLDER's general-purpose core API.
"""

from __future__ import annotations

import hashlib
import json
import math
import subprocess
import sys
import urllib.request
from collections.abc import Callable, Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any

import pandas as pd
import yaml
from Bio.PDB import MMCIFParser, PDBParser
from rdkit import Chem
from rdkit.Chem import Draw

ASSET_DIR = Path(__file__).resolve().parent / "assets" / "structure_gated_oracle"
CONFIG_PATH = ASSET_DIR / "config.yaml"
LIGANDS_PATH = ASSET_DIR / "ligands.csv"
ILLUSTRATIVE_RESULTS_PATH = ASSET_DIR / "illustrative_results.csv"

AUDIT_COLUMNS = (
    "candidate_id",
    "prediction_id",
    "expected_role",
    "structural_classification",
    "atp_site_contact",
    "back_pocket_contact",
    "dfg_out",
    "dfg_d1_angstrom",
    "dfg_d2_angstrom",
    "glu71_sidechain_hbond",
    "asp168_backbone_hbond",
    "met109_backbone_hbond",
    "matched_interactions",
    "matched_interaction_keys",
    "missing_interaction_keys",
    "native_affinity_score",
    "raw_score",
    "bounded_score",
    "scalar_reward",
    "score_rank",
    "structure_gated_rank",
    "evaluation_status",
    "evaluation_reason",
    "structure_path",
)


def load_config(path: str | Path = CONFIG_PATH) -> dict[str, Any]:
    """Load the versioned MAPK14 tutorial policy."""

    value = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    if not isinstance(value, dict) or value.get("version") != 2:
        raise ValueError("Expected structure-gated tutorial configuration version 2.")
    return value


def load_ligands(path: str | Path = LIGANDS_PATH) -> pd.DataFrame:
    """Load and chemically validate the six-member literature panel."""

    frame = pd.read_csv(path)
    if len(frame) != 6 or frame["candidate_id"].duplicated().any():
        raise ValueError(
            "The MAPK14 tutorial panel must contain six unique candidates."
        )
    canonical: list[str] = []
    for candidate_id, smiles in zip(frame["candidate_id"], frame["canonical_smiles"]):
        molecule = Chem.MolFromSmiles(str(smiles))
        if molecule is None:
            raise ValueError(f"Invalid SMILES for {candidate_id!r}.")
        canonical.append(Chem.MolToSmiles(molecule, isomericSmiles=True))
    frame = frame.copy()
    frame["canonical_smiles"] = canonical
    return frame


def molecule_grid_svg(ligands: pd.DataFrame) -> str:
    """Render the tutorial panel as an SVG suitable for a marimo Markdown cell."""

    molecules = [Chem.MolFromSmiles(value) for value in ligands["canonical_smiles"]]
    legends = [
        f"{row.name}\n{row.role}"
        for row in ligands[["name", "role"]].itertuples(index=False)
    ]
    return str(
        Draw.MolsToGridImage(
            molecules,
            molsPerRow=3,
            subImgSize=(360, 260),
            legends=legends,
            useSVG=True,
        )
    )


def bounded_score(value: float) -> float:
    """Map a finite score monotonically into the open interval ``(0, 1)``."""

    score = _finite_float(value, "raw_score")
    return 0.5 + math.atan(score) / math.pi


def structure_gated_reward(
    *,
    is_type_ii: bool,
    matched_interactions: int,
    raw_score: float,
    key_interaction_count: int,
) -> float:
    """Encode the tutorial's lexicographic ordering as one Oracle scalar."""

    if isinstance(matched_interactions, bool) or not isinstance(
        matched_interactions, int
    ):
        raise TypeError("matched_interactions must be an integer.")
    if not 0 <= matched_interactions <= key_interaction_count:
        raise ValueError("matched_interactions must be between zero and K.")
    if key_interaction_count < 1:
        raise ValueError("At least one key interaction is required.")
    return (
        (key_interaction_count + 1) * int(bool(is_type_ii))
        + matched_interactions
        + bounded_score(raw_score)
    )


def match_key_interactions(
    observed: Iterable[str], configured: Sequence[str]
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Return configured matches and misses, counting every key at most once."""

    observed_set = frozenset(str(item) for item in observed)
    configured_unique = tuple(dict.fromkeys(str(item) for item in configured))
    matched = tuple(item for item in configured_unique if item in observed_set)
    missing = tuple(item for item in configured_unique if item not in observed_set)
    return matched, missing


def match_interaction_features(
    events: Iterable[Mapping[str, Any]], features: Sequence[Mapping[str, Any]]
) -> tuple[tuple[str, ...], tuple[str, ...], dict[str, list[dict[str, Any]]]]:
    """Match atom-resolved events to the frozen literature feature definitions."""

    evidence: dict[str, list[dict[str, Any]]] = {}
    matched: list[str] = []
    missing: list[str] = []
    event_rows = [dict(event) for event in events]
    for feature in features:
        feature_id = str(feature["id"])
        allowed_atoms = {str(atom) for atom in feature["receptor_atoms"]}
        matches = [
            event
            for event in event_rows
            if event.get("interaction_key") == feature["key"]
            and event.get("ligand_role") == feature["ligand_role"]
            and event.get("protein_role") == feature["protein_role"]
            and any(
                atom.get("atom_name") in allowed_atoms
                for atom in event.get("protein_atoms", [])
            )
        ]
        evidence[feature_id] = matches
        (matched if matches else missing).append(str(feature["key"]))
    return tuple(matched), tuple(missing), evidence


def _prepared_interactions(
    structure_path: Path,
    *,
    smiles: str,
    ligand_chain: str,
    ligand_number: int,
    ligand_name: str,
    policy: Mapping[str, Any],
) -> dict[str, Any]:
    """Run chemistry preparation and ProLIF outside the notebook process."""

    request = {
        "structure_path": str(structure_path.resolve()),
        "smiles": smiles,
        "ligand": {
            "chain_id": ligand_chain,
            "residue_number": int(ligand_number),
            "residue_name": ligand_name,
        },
        "preparation": policy["interaction_policy"]["preparation"],
        "interactions": policy["interaction_policy"]["prolif_interactions"],
    }
    worker = Path(__file__).with_name("_structure_gated_prolif_worker.py")
    completed = subprocess.run(
        [sys.executable, str(worker)],
        input=json.dumps(request),
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode:
        detail = completed.stderr.strip().splitlines()[-1] if completed.stderr else ""
        raise ValueError(f"isolated preparation failed: {detail}")
    try:
        return json.loads(completed.stdout)
    except json.JSONDecodeError as exc:
        raise ValueError("isolated preparation returned invalid JSON") from exc


def rank_candidates(
    frame: pd.DataFrame, *, config: Mapping[str, Any] | None = None
) -> pd.DataFrame:
    """Add score and structure-gated ranks while retaining non-evaluable rows."""

    policy = dict(config or load_config())
    keys = tuple(policy["interaction_policy"]["features"])
    result = frame.copy()
    required = {
        "candidate_id",
        "raw_score",
        "matched_interactions",
        "structural_classification",
        "evaluation_status",
    }
    missing = sorted(required - set(result.columns))
    if missing:
        raise ValueError(f"Ranking input is missing columns: {missing}")

    result["bounded_score"] = math.nan
    result["scalar_reward"] = math.nan
    result["score_rank"] = pd.array([pd.NA] * len(result), dtype="Int64")
    result["structure_gated_rank"] = pd.array([pd.NA] * len(result), dtype="Int64")
    evaluable = result["evaluation_status"].eq("evaluable")
    for index in result.index[evaluable]:
        raw_score = _finite_float(result.at[index, "raw_score"], "raw_score")
        matched = int(result.at[index, "matched_interactions"])
        is_type_ii = result.at[index, "structural_classification"] == "type_II"
        result.at[index, "bounded_score"] = bounded_score(raw_score)
        result.at[index, "scalar_reward"] = structure_gated_reward(
            is_type_ii=is_type_ii,
            matched_interactions=matched,
            raw_score=raw_score,
            key_interaction_count=len(keys),
        )

    evaluable_indices = list(result.index[evaluable])
    score_order = sorted(
        evaluable_indices,
        key=lambda index: (
            -float(result.at[index, "raw_score"]),
            str(result.at[index, "candidate_id"]),
        ),
    )
    gated_order = sorted(
        evaluable_indices,
        key=lambda index: (
            -int(result.at[index, "structural_classification"] == "type_II"),
            -int(result.at[index, "matched_interactions"]),
            -float(result.at[index, "raw_score"]),
            str(result.at[index, "candidate_id"]),
        ),
    )
    for rank, index in enumerate(score_order, start=1):
        result.at[index, "score_rank"] = rank
    for rank, index in enumerate(gated_order, start=1):
        result.at[index, "structure_gated_rank"] = rank
    return result


def analyze_pose(
    structure_path: str | Path,
    *,
    raw_score: float,
    candidate_id: str,
    prediction_id: str,
    canonical_smiles: str | None = None,
    native_affinity_score: float | None = None,
    expected_role: str = "",
    config: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Evaluate one predicted pose and return one auditable ranking row."""

    policy = dict(config or load_config())
    path = Path(structure_path)
    base = _audit_base(candidate_id, prediction_id, expected_role, raw_score, path)
    try:
        score = _finite_float(raw_score, "raw_score")
        if not path.is_file():
            raise ValueError(f"structure file does not exist: {path}")
        target = policy["target"]
        predicate = policy["structural_predicate"]
        model = _load_model(path)
        chain = model[str(target["protein_chain"])]
        ligand_atoms = _ligand_heavy_atoms(model, str(target["ligand_chain"]))
        cutoff = float(predicate["contact_cutoff_angstrom"])
        atp_contact = _contacts_any(
            chain, ligand_atoms, predicate["atp_site_residues"], cutoff
        )
        back_contact = _contacts_any(
            chain, ligand_atoms, predicate["back_pocket_residues"], cutoff
        )
        remote_contact = _contacts_any(
            chain, ligand_atoms, predicate["remote_um101_region"], cutoff
        )
        dfg_d1, dfg_d2, dfg_out = _dfg_conformation(
            chain, predicate["dfg_conformation"]
        )
        is_type_ii = atp_contact and back_contact and dfg_out

        if canonical_smiles is None:
            ligand_table = load_ligands().set_index("candidate_id")
            if candidate_id not in ligand_table.index:
                raise ValueError("canonical_smiles is required for a new candidate")
            canonical_smiles = str(ligand_table.at[candidate_id, "canonical_smiles"])
        ligand_residue = _single_ligand_residue(model, str(target["ligand_chain"]))
        prepared = _prepared_interactions(
            path,
            smiles=canonical_smiles,
            ligand_chain=str(target["ligand_chain"]),
            ligand_number=int(ligand_residue.id[1]),
            ligand_name=str(ligand_residue.resname),
            policy=policy,
        )
        features = tuple(policy["interaction_policy"]["features"])
        matched, missing, evidence = match_interaction_features(
            prepared["events"], features
        )
        classification = (
            "type_II"
            if is_type_ii
            else "remote_site_candidate"
            if remote_contact and not atp_contact
            else "other_evaluable_mode"
        )
        reward = structure_gated_reward(
            is_type_ii=is_type_ii,
            matched_interactions=len(matched),
            raw_score=score,
            key_interaction_count=len(features),
        )
        feature_states = {
            str(feature["id"]): (
                "matched" if evidence[str(feature["id"])] else "geometrically_failed"
            )
            for feature in features
        }
        matched_events = [event for matches in evidence.values() for event in matches]
        return {
            **base,
            "structural_classification": classification,
            "atp_site_contact": atp_contact,
            "back_pocket_contact": back_contact,
            "dfg_out": dfg_out,
            "dfg_d1_angstrom": dfg_d1,
            "dfg_d2_angstrom": dfg_d2,
            "matched_interactions": len(matched),
            "matched_interaction_keys": ";".join(matched),
            "missing_interaction_keys": ";".join(missing),
            "native_affinity_score": (
                score if native_affinity_score is None else native_affinity_score
            ),
            "bounded_score": bounded_score(score),
            "scalar_reward": reward,
            **feature_states,
            "interaction_evidence": evidence,
            "diagnostic_interactions": [
                event for event in prepared["events"] if event not in matched_events
            ],
            "preparation_audit": prepared["preparation"],
            "evaluation_status": "evaluable",
            "evaluation_reason": "",
        }
    except Exception as exc:  # one failed diagnostic must remain visible in the audit
        return {
            **base,
            "evaluation_status": "not_evaluable",
            "evaluation_reason": f"{type(exc).__name__}: {exc}",
        }


def make_oracle_scoring_function(
    *,
    candidate_id: str,
    expected_role: str,
    canonical_smiles: str | None = None,
    audit_sink: list[dict[str, Any]],
    config: Mapping[str, Any] | None = None,
) -> Callable[[Any], float]:
    """Build an ``Oracle(scoring_function=...)`` callback for one candidate."""

    policy = dict(config or load_config())
    selector = str(policy["execution"]["score_metric"])
    score_multiplier = float(policy["execution"].get("score_multiplier", 1.0))

    def score(context: Any) -> float:
        raw_score = context.aggregated_metrics.get(selector)
        if raw_score is None:
            raise ValueError(f"Oracle output did not contain {selector!r}.")
        structure_path, prediction_id = _single_prediction_path(context)
        oriented_score = float(raw_score) * score_multiplier
        row = analyze_pose(
            structure_path,
            raw_score=oriented_score,
            candidate_id=candidate_id,
            prediction_id=prediction_id,
            canonical_smiles=canonical_smiles,
            native_affinity_score=float(raw_score),
            expected_role=expected_role,
            config=policy,
        )
        audit_sink.append(row)
        if row["evaluation_status"] != "evaluable":
            raise ValueError(row["evaluation_reason"])
        return float(row["scalar_reward"])

    return score


def fetch_references(
    destination: str | Path,
    *,
    config: Mapping[str, Any] | None = None,
) -> dict[str, Path]:
    """Download reference mmCIF files and verify their pinned SHA-256 hashes."""

    policy = dict(config or load_config())
    output = Path(destination)
    output.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    for name, reference in policy["references"].items():
        path = output / f"{reference['pdb_id']}.cif"
        if not path.exists():
            urllib.request.urlretrieve(str(reference["url"]), path)
        digest = sha256_file(path)
        if digest != reference["sha256"]:
            raise ValueError(
                f"Checksum mismatch for {reference['pdb_id']}: {digest}; "
                f"expected {reference['sha256']}"
            )
        paths[name] = path
    return paths


def reanalyze_saved_predictions(
    run_root: str | Path,
    *,
    config: Mapping[str, Any] | None = None,
) -> pd.DataFrame:
    """Reanalyse six completed Oracle runs without invoking the prediction backend."""

    policy = dict(config or load_config())
    ligand_table = load_ligands().set_index("candidate_id")
    rows: list[dict[str, Any]] = []
    for candidate_id, ligand in ligand_table.iterrows():
        result_dir = Path(run_root) / candidate_id / "oracle_run" / "results"
        metrics_path = result_dir / "metrics.csv"
        structures = sorted((result_dir / "structures").glob("*.cif"))
        if not metrics_path.is_file() or len(structures) != 1:
            rows.append(
                {
                    **_audit_base(
                        str(candidate_id),
                        "unavailable",
                        str(ligand["role"]),
                        math.nan,
                        structures[0] if structures else result_dir / "missing.cif",
                    ),
                    "evaluation_status": "not_evaluable",
                    "evaluation_reason": "saved run requires one structure and metrics.csv",
                }
            )
            continue
        metrics = pd.read_csv(metrics_path)
        selected = metrics[
            metrics["metric_name"].eq("affinity_pred_value")
            & metrics["status"].eq("computed")
            & metrics["chain_id"].eq(str(policy["target"]["ligand_chain"]))
        ]
        if len(selected) != 1 or pd.isna(selected.iloc[0]["value"]):
            rows.append(
                {
                    **_audit_base(
                        str(candidate_id),
                        structures[0].stem,
                        str(ligand["role"]),
                        math.nan,
                        structures[0],
                    ),
                    "evaluation_status": "not_evaluable",
                    "evaluation_reason": "exactly one computed ligand affinity is required",
                }
            )
            continue
        native_score = float(selected.iloc[0]["value"])
        rows.append(
            analyze_pose(
                structures[0],
                raw_score=native_score
                * float(policy["execution"].get("score_multiplier", 1.0)),
                native_affinity_score=native_score,
                candidate_id=str(candidate_id),
                prediction_id=structures[0].stem,
                expected_role=str(ligand["role"]),
                canonical_smiles=str(ligand["canonical_smiles"]),
                config=policy,
            )
        )
    return pd.DataFrame(rows)


def reference_geometry_report(
    references: Mapping[str, str | Path],
    *,
    config: Mapping[str, Any] | None = None,
) -> pd.DataFrame:
    """Check frozen geometry, pocket contacts, and keys against both references."""

    policy = dict(config or load_config())
    target = policy["target"]
    conformation = policy["structural_predicate"]["dfg_conformation"]
    rows = []
    for name, path in references.items():
        model = _load_model(Path(path))
        protein_chain = model[str(target["protein_chain"])]
        ligand_name = str(policy["references"][name]["ligand_residue"])
        ligand_residues = [
            residue for residue in protein_chain if residue.resname == ligand_name
        ]
        if len(ligand_residues) != 1:
            raise ValueError(
                f"Expected one {ligand_name} residue in {path}, found "
                f"{len(ligand_residues)}."
            )
        ligand_residue = ligand_residues[0]
        ligand_atoms = tuple(
            atom
            for atom in ligand_residue
            if str(getattr(atom, "element", "")).upper() != "H"
        )
        predicate = policy["structural_predicate"]
        cutoff = float(predicate["contact_cutoff_angstrom"])
        dfg_d1, dfg_d2, dfg_out = _dfg_conformation(protein_chain, conformation)
        candidate_id = str(policy["references"][name]["ligand_candidate_id"])
        ligand_table = load_ligands().set_index("candidate_id")
        prepared = _prepared_interactions(
            Path(path),
            smiles=str(ligand_table.at[candidate_id, "canonical_smiles"]),
            ligand_chain=str(target["protein_chain"]),
            ligand_number=int(ligand_residue.id[1]),
            ligand_name=str(ligand_residue.resname),
            policy=policy,
        )
        features = tuple(policy["interaction_policy"]["features"])
        matched, _, evidence = match_interaction_features(prepared["events"], features)
        rows.append(
            {
                "reference": name,
                "pdb_id": policy["references"][name]["pdb_id"],
                "atp_site_contact": _contacts_any(
                    protein_chain,
                    ligand_atoms,
                    predicate["atp_site_residues"],
                    cutoff,
                ),
                "back_pocket_contact": _contacts_any(
                    protein_chain,
                    ligand_atoms,
                    predicate["back_pocket_residues"],
                    cutoff,
                ),
                "dfg_d1_angstrom": dfg_d1,
                "dfg_d2_angstrom": dfg_d2,
                "dfg_out": dfg_out,
                "matched_key_interactions": ";".join(matched),
                **{
                    str(feature["id"]): (
                        "matched"
                        if evidence[str(feature["id"])]
                        else "geometrically_failed"
                    )
                    for feature in features
                },
                "interaction_evidence": evidence,
                "preparation_audit": prepared["preparation"],
            }
        )
    return pd.DataFrame(rows)


def write_audit_bundle(
    frame: pd.DataFrame,
    output_dir: str | Path,
    *,
    inputs: Iterable[str | Path] = (),
    metadata: Mapping[str, Any] | None = None,
) -> tuple[Path, Path, Path]:
    """Write summary, atom-level interaction evidence, and provenance."""

    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    ranked = rank_candidates(frame)
    for column in AUDIT_COLUMNS:
        if column not in ranked:
            ranked[column] = pd.NA
    csv_path = output / "structure_gated_oracle_audit.csv"
    ranked.loc[:, AUDIT_COLUMNS].to_csv(csv_path, index=False)
    evidence_path = output / "structure_gated_oracle_interactions.json"
    evidence_path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "candidates": [
                    {
                        "candidate_id": row.get("candidate_id"),
                        "prediction_id": row.get("prediction_id"),
                        "evaluation_status": row.get("evaluation_status"),
                        "features": row.get("interaction_evidence", {}),
                        "diagnostics": row.get("diagnostic_interactions", []),
                        "preparation": row.get("preparation_audit", {}),
                    }
                    for row in ranked.to_dict(orient="records")
                ],
            },
            indent=2,
            sort_keys=True,
            default=str,
        )
        + "\n",
        encoding="utf-8",
    )
    provenance = {
        "schema_version": 2,
        "policy_sha256": sha256_file(CONFIG_PATH),
        "input_sha256": {str(Path(path)): sha256_file(path) for path in inputs},
        **dict(metadata or {}),
    }
    json_path = output / "structure_gated_oracle_provenance.json"
    json_path.write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return csv_path, evidence_path, json_path


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _audit_base(
    candidate_id: str,
    prediction_id: str,
    expected_role: str,
    raw_score: Any,
    path: Path,
) -> dict[str, Any]:
    return {
        "candidate_id": candidate_id,
        "prediction_id": prediction_id,
        "expected_role": expected_role,
        "structural_classification": pd.NA,
        "atp_site_contact": pd.NA,
        "back_pocket_contact": pd.NA,
        "dfg_out": pd.NA,
        "dfg_d1_angstrom": math.nan,
        "dfg_d2_angstrom": math.nan,
        "glu71_sidechain_hbond": pd.NA,
        "asp168_backbone_hbond": pd.NA,
        "met109_backbone_hbond": pd.NA,
        "matched_interactions": pd.NA,
        "matched_interaction_keys": "",
        "missing_interaction_keys": "",
        "native_affinity_score": raw_score,
        "raw_score": raw_score,
        "bounded_score": math.nan,
        "scalar_reward": math.nan,
        "score_rank": pd.NA,
        "structure_gated_rank": pd.NA,
        "structure_path": str(path),
        "interaction_evidence": {},
        "diagnostic_interactions": [],
        "preparation_audit": {},
    }


def _finite_float(value: Any, name: str) -> float:
    if isinstance(value, bool):
        raise TypeError(f"{name} must be a finite number, not bool.")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite.")
    return result


def _load_model(path: Path):
    parser = (
        MMCIFParser(QUIET=True)
        if path.suffix.lower() in {".cif", ".mmcif"}
        else PDBParser(QUIET=True)
    )
    structure = parser.get_structure(path.stem, str(path))
    model = next(structure.get_models(), None)
    if model is None:
        raise ValueError("structure contains no model")
    return model


def _ligand_heavy_atoms(model: Any, chain_id: str) -> tuple[Any, ...]:
    if chain_id not in model:
        raise ValueError(f"ligand chain {chain_id!r} was not found")
    atoms = tuple(
        atom
        for residue in model[chain_id]
        if residue.resname not in {"HOH", "WAT", "DOD"}
        for atom in residue
        if str(getattr(atom, "element", "")).upper() != "H"
    )
    if not atoms:
        raise ValueError(f"ligand chain {chain_id!r} contains no heavy atoms")
    return atoms


def _single_ligand_residue(model: Any, chain_id: str):
    if chain_id not in model:
        raise ValueError(f"ligand chain {chain_id!r} was not found")
    residues = [
        residue
        for residue in model[chain_id]
        if residue.resname not in {"HOH", "WAT", "DOD"}
    ]
    if len(residues) != 1:
        raise ValueError(
            f"expected one ligand residue in chain {chain_id!r}, found {len(residues)}"
        )
    return residues[0]


def _protein_residue(chain: Any, number: int):
    matches = [
        residue
        for residue in chain
        if residue.id[0] == " " and residue.id[1] == int(number)
    ]
    if len(matches) != 1:
        raise ValueError(
            f"expected one protein residue numbered {number}, found {len(matches)}"
        )
    return matches[0]


def _contacts_any(
    chain: Any, ligand_atoms: Sequence[Any], residues: Sequence[int], cutoff: float
) -> bool:
    cutoff_squared = cutoff * cutoff
    for number in residues:
        residue = _protein_residue(chain, int(number))
        for protein_atom in residue:
            if str(getattr(protein_atom, "element", "")).upper() == "H":
                continue
            for ligand_atom in ligand_atoms:
                delta = protein_atom.coord - ligand_atom.coord
                if float((delta * delta).sum()) <= cutoff_squared:
                    return True
    return False


def _atom_distance(chain: Any, definition: Mapping[str, Any]) -> float:
    residue_1, atom_1 = definition["atom_1"]
    residue_2, atom_2 = definition["atom_2"]
    first = _protein_residue(chain, int(residue_1))[str(atom_1)]
    second = _protein_residue(chain, int(residue_2))[str(atom_2)]
    return round(float(first - second), 3)


def _dfg_conformation(
    chain: Any, definition: Mapping[str, Any]
) -> tuple[float, float, bool]:
    """Apply the manuscript SI classical DFG-out D1/D2 definition."""

    d1 = _atom_distance(chain, definition["d1"])
    d2 = _atom_distance(chain, definition["d2"])
    is_dfg_out = d1 <= float(definition["d1"]["threshold_angstrom"]) and d2 >= float(
        definition["d2"]["threshold_angstrom"]
    )
    return d1, d2, is_dfg_out


def _single_prediction_path(context: Any) -> tuple[Path, str]:
    candidates: list[tuple[str, str]] = []
    for frame in (context.system_metrics, context.chain_metrics):
        if frame is None or frame.empty:
            continue
        path_column = next(
            (name for name in ("cif_file", "structure_path") if name in frame), None
        )
        if path_column is None:
            continue
        for _, row in frame.iterrows():
            value = row.get(path_column)
            if pd.isna(value):
                continue
            prediction_id = str(row.get("model_name", Path(str(value)).stem))
            candidates.append((str(value), prediction_id))
    unique = list(dict.fromkeys(candidates))
    unique_paths = list(dict.fromkeys(value for value, _ in unique))
    if len(unique_paths) != 1:
        raise ValueError(
            "The tutorial callback requires exactly one predicted structure per ligand; "
            f"found {len(unique_paths)}. Configure one repeat and one diffusion sample."
        )
    raw_path = Path(unique_paths[0])
    alternatives = (
        raw_path,
        context.run_dir / raw_path,
        context.run_dir / "results" / raw_path,
        context.run_dir / "results" / "structures" / raw_path,
    )
    resolved = next((path for path in alternatives if path.is_file()), alternatives[0])
    prediction_id = next(
        identifier for value, identifier in unique if value == unique_paths[0]
    )
    return resolved, prediction_id
