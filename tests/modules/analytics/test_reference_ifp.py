"""Tests for reusable reference-complex interaction fingerprints."""

import json
import subprocess
from contextlib import ExitStack
from pathlib import Path

import gemmi
import pytest

from cofolder.modules.analytics.ifp_filtering import (
    IFPFilterMode,
    ReferenceIFPFilterPolicy,
    evaluate_reference_ifp_filter,
)
from cofolder.modules.analytics.reference_ifp import (
    IFPComparisonStatus,
    IFPExtractionConfig,
    IFPSimilarityMetric,
    IFPTaxonomy,
    InteractionKey,
    LigandSelector,
    ProLIFWorkerCrashError,
    ProLIFWorkerTimeoutError,
    ReferenceEntitySelectionError,
    compare_interaction_fingerprints,
    extract_interaction_fingerprint,
    map_reference_identities,
)
from cofolder.modules.contracts import (
    METRIC_CATALOG,
    EvidenceRegime,
    OptimizationDirection,
)
from tests.modules.analytics.test_reproduction import (
    _write_predicted_pdb,
    _write_reference_pdb,
    _write_reference_with_duplicate_ligands_pdb,
)


def test_distance_extraction_mapping_and_similarity(temp_dir):
    reference_path = temp_dir / "reference.pdb"
    prediction_path = temp_dir / "prediction.pdb"
    _write_reference_pdb(reference_path)
    _write_predicted_pdb(prediction_path, ((1.3, 1.2, 0.0), (2.7, 1.2, 0.0)))

    reference = extract_interaction_fingerprint(reference_path)
    prediction = extract_interaction_fingerprint(
        prediction_path, ligand=LigandSelector(chain_id="Z")
    )
    mapping = map_reference_identities(
        reference,
        prediction,
        reference_structure_path=reference_path,
        predicted_structure_path=prediction_path,
    )
    comparison = compare_interaction_fingerprints(reference, prediction, mapping)

    assert reference.serialized_interactions() == ["A:1:distance_contact"]
    assert comparison.status is IFPComparisonStatus.COMPARABLE
    assert comparison.similarities[IFPSimilarityMetric.JACCARD] == 1.0
    assert comparison.similarities[IFPSimilarityMetric.REFERENCE_COVERAGE] == 1.0
    assert comparison.missing_interactions == ()


def test_nonmatching_prediction_reports_missing_reference_contact(temp_dir):
    reference_path = temp_dir / "reference.pdb"
    prediction_path = temp_dir / "prediction.pdb"
    _write_reference_pdb(reference_path)
    _write_predicted_pdb(
        prediction_path, ((100.0, 100.0, 100.0), (102.0, 100.0, 100.0))
    )
    reference = extract_interaction_fingerprint(reference_path)
    prediction = extract_interaction_fingerprint(
        prediction_path, ligand=LigandSelector(chain_id="Z")
    )
    mapping = map_reference_identities(
        reference,
        prediction,
        reference_structure_path=reference_path,
        predicted_structure_path=prediction_path,
    )
    comparison = compare_interaction_fingerprints(reference, prediction, mapping)

    assert comparison.similarities[IFPSimilarityMetric.JACCARD] == 0.0
    assert [str(item) for item in comparison.missing_interactions] == [
        "A:1:distance_contact"
    ]
    outcome = evaluate_reference_ifp_filter(
        reference,
        comparison,
        ReferenceIFPFilterPolicy(
            IFPFilterMode.SIMILARITY,
            IFPSimilarityMetric.JACCARD,
            threshold=0.0,
        ),
    )
    assert outcome.passed is True  # inclusive boundary


def test_required_subset_policy_and_key_parser(temp_dir):
    reference_path = temp_dir / "reference.pdb"
    prediction_path = temp_dir / "prediction.pdb"
    _write_reference_pdb(reference_path)
    _write_predicted_pdb(prediction_path, ((1.3, 1.2, 0.0), (2.7, 1.2, 0.0)))
    reference = extract_interaction_fingerprint(reference_path)
    prediction = extract_interaction_fingerprint(
        prediction_path, ligand=LigandSelector(chain_id="Z")
    )
    mapping = map_reference_identities(
        reference,
        prediction,
        reference_structure_path=reference_path,
        predicted_structure_path=prediction_path,
    )
    comparison = compare_interaction_fingerprints(reference, prediction, mapping)
    required = frozenset({InteractionKey.parse("A:1:distance_contact")})
    outcome = evaluate_reference_ifp_filter(
        reference,
        comparison,
        ReferenceIFPFilterPolicy(
            IFPFilterMode.REQUIRED, required_interactions=required
        ),
    )
    assert outcome.status == "accepted"


def test_ambiguous_reference_ligand_requires_selector(temp_dir):
    path = temp_dir / "duplicate.pdb"
    _write_reference_with_duplicate_ligands_pdb(path)
    with pytest.raises(
        ReferenceEntitySelectionError, match="ambiguous_reference_ligand"
    ):
        extract_interaction_fingerprint(path)
    selected = extract_interaction_fingerprint(
        path, ligand=LigandSelector(chain_id="L")
    )
    assert selected.ligand.chain_id == "L"


def test_tied_sequence_mapping_is_explicitly_unmappable(temp_dir):
    reference_path = temp_dir / "reference.pdb"
    prediction_path = temp_dir / "prediction.pdb"
    base_prediction = temp_dir / "base_prediction.pdb"
    _write_reference_pdb(reference_path)
    _write_predicted_pdb(base_prediction, ((1.3, 1.2, 0.0), (2.7, 1.2, 0.0)))
    lines = base_prediction.read_text(encoding="utf-8").splitlines()
    protein = [line for line in lines if line.startswith("ATOM")]
    ligand = [line for line in lines if line.startswith("HETATM")]
    prediction_path.write_text(
        "\n".join(
            [line[:21] + "X" + line[22:] for line in protein]
            + [line[:21] + "Y" + line[22:] for line in protein]
            + ligand
            + ["TER", "END", ""]
        ),
        encoding="utf-8",
    )
    reference = extract_interaction_fingerprint(reference_path)
    prediction = extract_interaction_fingerprint(
        prediction_path,
        ligand=LigandSelector(chain_id="Z"),
        receptor_chains=("X", "Y"),
    )
    mapping = map_reference_identities(
        reference,
        prediction,
        reference_structure_path=reference_path,
        predicted_structure_path=prediction_path,
    )
    comparison = compare_interaction_fingerprints(reference, prediction, mapping)
    assert "ambiguous_chain_mapping" in mapping.failures
    assert comparison.status is IFPComparisonStatus.NOT_EVALUABLE


def test_prolif_worker_is_opt_in_and_failure_is_contained(temp_dir, monkeypatch):
    path = temp_dir / "reference.pdb"
    _write_reference_pdb(path)

    def must_not_run(*args, **kwargs):
        raise AssertionError("worker launched")

    monkeypatch.setattr(subprocess, "run", must_not_run)
    assert extract_interaction_fingerprint(path).taxonomy is IFPTaxonomy.DISTANCE

    payload = {
        "ligand": {
            "chain_id": "L",
            "residue_number": 1,
            "insertion_code": "",
            "residue_name": "LIG",
        },
        "receptor_chains": ["A"],
        "interactions": ["A:1:HBAcceptor"],
        "events": [
            {
                "interaction_key": "A:1:HBAcceptor",
                "ligand_role": "acceptor",
                "protein_role": "donor",
                "ligand_atoms": [
                    {
                        "chain_id": "L",
                        "residue_number": 1,
                        "insertion_code": "",
                        "residue_name": "LIG",
                        "atom_name": "O1",
                        "element": "O",
                        "atom_serial": 10,
                        "source_index": 9,
                    }
                ],
                "protein_atoms": [
                    {
                        "chain_id": "A",
                        "residue_number": 1,
                        "insertion_code": "",
                        "residue_name": "ASP",
                        "atom_name": "N",
                        "element": "N",
                        "atom_serial": 1,
                        "source_index": 0,
                    }
                ],
                "geometry": [
                    {"name": "distance", "value": 2.9, "unit": "angstrom"},
                    {"name": "DHA_angle", "value": 165.0, "unit": "degree"},
                ],
            }
        ],
    }
    monkeypatch.setattr(
        subprocess,
        "run",
        lambda *args, **kwargs: subprocess.CompletedProcess(
            args[0], 0, json.dumps(payload), ""
        ),
    )
    prolif = extract_interaction_fingerprint(
        path, config=IFPExtractionConfig(IFPTaxonomy.PROLIF)
    )
    assert prolif.serialized_interactions() == ["A:1:hb_acceptor"]
    assert len(prolif.events) == 1
    event = prolif.events[0]
    assert event.protein_atoms[0].residue_name == "ASP"
    assert event.protein_atoms[0].atom_name == "N"
    assert (event.ligand_role, event.protein_role) == ("acceptor", "donor")
    assert {item.name: item.value for item in event.geometry} == {
        "distance": 2.9,
        "DHA_angle": 165.0,
    }
    assert prolif.serialized_events()[0]["interaction_type"] == "hb_acceptor"

    def timeout(*args, **kwargs):
        raise subprocess.TimeoutExpired(args[0], 1)

    monkeypatch.setattr(subprocess, "run", timeout)
    with pytest.raises(ProLIFWorkerTimeoutError):
        extract_interaction_fingerprint(
            path, config=IFPExtractionConfig(IFPTaxonomy.PROLIF)
        )

    monkeypatch.setattr(
        subprocess,
        "run",
        lambda *args, **kwargs: subprocess.CompletedProcess(
            args[0], -11, "", "segfault"
        ),
    )
    with pytest.raises(ProLIFWorkerCrashError):
        extract_interaction_fingerprint(
            path, config=IFPExtractionConfig(IFPTaxonomy.PROLIF)
        )


def test_real_prolif_pdb_and_mmcif_extraction_are_equivalent(temp_dir):
    cif_path = Path("src/cofolder/resources/examples/4HJO.cif")
    pdb_path = temp_dir / "4HJO.pdb"
    gemmi.read_structure(str(cif_path)).write_pdb(str(pdb_path))
    options = {
        "ligand": LigandSelector(chain_id="A", residue_number=1001),
        "receptor_chains": ("A",),
        "config": IFPExtractionConfig(IFPTaxonomy.PROLIF),
    }

    from_cif = extract_interaction_fingerprint(cif_path, **options)
    from_pdb = extract_interaction_fingerprint(pdb_path, **options)

    assert from_cif.serialized_interactions() == from_pdb.serialized_interactions()
    assert from_cif.serialized_events() == from_pdb.serialized_events()
    assert from_cif.events
    first = from_cif.events[0]
    assert first.ligand_atoms[0].residue_name == "AQ4"
    assert first.protein_atoms[0].atom_name
    assert all(item.value >= 0 for event in from_cif.events for item in event.geometry)


def test_mmcif_worker_aliases_long_chain_ids_and_cleans_temporary_pdb(temp_dir):
    from cofolder.modules.analytics._prolif_worker import (
        _prepare_topology,
        _selection,
    )

    structure = gemmi.read_structure("src/cofolder/resources/examples/4HJO.cif")
    for model in structure:
        for chain in model:
            chain.name = "protein_chain"
    cif_path = temp_dir / "long-chain.cif"
    structure.make_mmcif_document().write_file(str(cif_path))

    with ExitStack() as stack:
        pdb_path, original_to_alias, alias_to_original = _prepare_topology(
            cif_path, stack
        )
        assert pdb_path.is_file()
        assert original_to_alias == {"protein_chain": "A"}
        assert alias_to_original == {"A": "protein_chain"}
        assert gemmi.read_structure(str(pdb_path))[0][0].name == "A"
    assert not pdb_path.exists()
    assert _selection("A", 168, "B").endswith("resid 168 and icode B")


def test_reference_filter_metrics_are_registered_for_both_reference_regimes():
    similarity = METRIC_CATALOG["ifp_filter_similarity"]
    assert similarity.direction is OptimizationDirection.MAXIMIZE
    assert similarity.allowed_evidence_regimes == {
        EvidenceRegime.REFERENCE_STRUCTURE,
        EvidenceRegime.CUSTOM_POCKET,
    }
    for name in (
        "ifp_filter_similarity_metric",
        "ifp_filter_policy",
        "ifp_filter_taxonomy",
        "ifp_filter_required_interactions",
        "ifp_filter_missing_interactions",
        "ifp_filter_mapping_status",
        "ifp_filter_mapping_failures",
    ):
        assert name in METRIC_CATALOG
