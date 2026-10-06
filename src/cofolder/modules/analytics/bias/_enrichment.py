"""Cached RCSB and MMseqs enrichment for bias datasets."""

from __future__ import annotations

from collections.abc import Callable
from rdkit import DataStructs
from pandas.errors import EmptyDataError
from pathlib import Path
from cofolder.modules.analytics.bias_training import _resolve_mmseqs_bin
import json
import logging
from functools import lru_cache
import pandas as pd
import subprocess
import tempfile
from urllib.request import urlopen

from ._common import (
    BIAS_LIGAND_VIEW_MIN_SIMILARITY,
    BIAS_PAIRING_COMPONENT_1_ONLY,
    BIAS_PAIRING_COMPONENT_2_ONLY,
    BIAS_PAIRING_LIGAND_ONLY,
    BIAS_PAIRING_PAIRED,
    BIAS_PAIRING_PROTEIN_ONLY,
    BIAS_PROGRESS_LOG_EVERY,
    BIAS_PROGRESS_WRITE_EVERY,
    BIAS_PROTEIN_VIEW_MIN_SIMILARITY,
    CORE_ENTRY_URL,
    CORE_NONPOLY_URL,
    FASTA_URL,
    _norm_id,
    _normalize_smiles,
    _path_cache_token,
)

from ._datasets import (
    _build_protein_lookup_index,
    _collapse_pair_value,
    _finalize_bias_training_dataframe,
    _finalize_same_type_pair_dataframe,
    _reference_ligand_similarity_rows,
    _reference_pdb_id_from_bias_row,
    _reference_pdb_id_from_same_type_row,
    _reference_protein_similarity_rows,
)

from ._references import (
    _components_cif_smiles_index,
    _smiles_from_ccd_cache,
)

from ._similarity import (
    _morgan_fp_from_smiles,
)


def _fetch_json(url: str, timeout: int = 30) -> dict | list:
    with urlopen(url, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8"))


def _fetch_text(url: str, timeout: int = 30) -> str:
    with urlopen(url, timeout=timeout) as response:
        return response.read().decode("utf-8", errors="replace")


def _extract_nonpoly_comp_id(nonpoly_payload: dict) -> str | None:
    if not isinstance(nonpoly_payload, dict):
        return None
    container = nonpoly_payload.get("rcsb_nonpolymer_entity_container_identifiers", {})
    if isinstance(container, dict):
        comp = container.get("nonpolymer_comp_id")
        if comp:
            return str(comp).upper()
    ent = nonpoly_payload.get("pdbx_entity_nonpoly", {})
    if isinstance(ent, dict):
        comp = ent.get("comp_id")
        if comp:
            return str(comp).upper()
    return None


@lru_cache(maxsize=4096)
def _entry_ligand_ids(pdb_id: str, timeout: int = 30) -> tuple[str, ...]:
    try:
        payload = _fetch_json(CORE_ENTRY_URL.format(pdb_id=pdb_id), timeout=timeout)
    except Exception:
        return tuple()
    ids_obj = (
        payload.get("rcsb_entry_container_identifiers", {})
        if isinstance(payload, dict)
        else {}
    )
    entity_ids = (
        ids_obj.get("non_polymer_entity_ids", []) if isinstance(ids_obj, dict) else []
    )
    out: set[str] = set()
    for ent_id in entity_ids:
        try:
            nonpoly = _fetch_json(
                CORE_NONPOLY_URL.format(pdb_id=pdb_id, entity_id=ent_id),
                timeout=timeout,
            )
        except Exception:
            continue
        comp = _extract_nonpoly_comp_id(nonpoly)
        if comp:
            out.add(comp)
    return tuple(sorted(out))


@lru_cache(maxsize=4096)
def _entry_fasta_sequences(pdb_id: str, timeout: int = 30) -> tuple[str, ...]:
    try:
        txt = _fetch_text(FASTA_URL.format(pdb_id=pdb_id), timeout=timeout)
    except Exception:
        return tuple()
    sequences: list[str] = []
    buffer: list[str] = []
    for raw_line in txt.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if buffer:
                sequences.append("".join(buffer))
                buffer = []
            continue
        buffer.append(line)
    if buffer:
        sequences.append("".join(buffer))
    return tuple(sequence for sequence in sequences if sequence)


@lru_cache(maxsize=65536)
def _resolve_pdb_ligand_smiles(
    ccd_id: str,
    boltz_cache_token: str,
    components_cif_token: str,
) -> str:
    normalized_ccd_id = _norm_id(ccd_id)
    smiles = _components_cif_smiles_index(components_cif_token).get(normalized_ccd_id)
    if smiles:
        return smiles
    if boltz_cache_token:
        cached = _smiles_from_ccd_cache(normalized_ccd_id, Path(boltz_cache_token))
        if cached:
            return _normalize_smiles(cached)
    return ""


@lru_cache(maxsize=4096)
def _cached_pdb_ligand_entries(
    pdb_id: str,
    boltz_cache_token: str,
    components_cif_token: str,
    timeout: int = 30,
) -> tuple[tuple[str, str], ...]:
    entries: list[tuple[str, str]] = []
    for ligand_id in _entry_ligand_ids(pdb_id, timeout=timeout):
        smiles = _resolve_pdb_ligand_smiles(
            ligand_id,
            boltz_cache_token,
            components_cif_token,
        )
        if not smiles:
            continue
        entries.append((str(ligand_id).strip(), smiles))
    return tuple(entries)


@lru_cache(maxsize=8192)
def _cached_pdb_protein_similarity_rows(
    pdb_id: str,
    query_sequence: str,
    timeout: int = 30,
) -> tuple[tuple[str, float], ...]:
    if not query_sequence:
        return tuple()
    target_sequences = tuple(
        sequence
        for sequence in _entry_fasta_sequences(pdb_id, timeout=timeout)
        if sequence
    )
    if not target_sequences:
        return tuple()

    resolved_mmseqs = _resolved_mmseqs_bin_token()
    if not resolved_mmseqs:
        return tuple()

    with tempfile.TemporaryDirectory(prefix="cofolder-bias-mmseqs-") as tmpdir:
        tmp_root = Path(tmpdir)
        query_fasta = tmp_root / "query.fasta"
        target_fasta = tmp_root / "target.fasta"
        _write_fasta_records(query_fasta, [("query", str(query_sequence))])
        _write_fasta_records(
            target_fasta,
            [
                (f"{str(pdb_id).strip().upper()}_{index}", sequence)
                for index, sequence in enumerate(target_sequences, start=1)
            ],
        )

        qdb = tmp_root / "query_db"
        tdb = tmp_root / "target_db"
        rdb = tmp_root / "result_db"
        out_tsv = tmp_root / "result.tsv"
        search_tmp = tmp_root / "search_tmp"
        search_tmp.mkdir(parents=True, exist_ok=True)

        commands = [
            [resolved_mmseqs, "createdb", str(query_fasta), str(qdb)],
            [resolved_mmseqs, "createdb", str(target_fasta), str(tdb)],
            [
                resolved_mmseqs,
                "search",
                str(qdb),
                str(tdb),
                str(rdb),
                str(search_tmp),
                "--threads",
                "1",
                "--max-seqs",
                str(max(1, len(target_sequences))),
            ],
            [
                resolved_mmseqs,
                "convertalis",
                str(qdb),
                str(tdb),
                str(rdb),
                str(out_tsv),
                "--format-output",
                "target,pident,tseq",
            ],
        ]

        try:
            for command in commands:
                subprocess.run(
                    command,
                    check=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    timeout=max(30, int(timeout)),
                )
        except (
            OSError,
            subprocess.CalledProcessError,
            subprocess.TimeoutExpired,
        ) as exc:
            logging.getLogger("cofolder.modules.analytics.bias").warning(
                "MMseqs protein fallback failed for pdb_id=%s: %s",
                str(pdb_id).strip().upper(),
                exc,
            )
            return tuple()

        similarity_by_sequence: dict[str, float] = {}
        if out_tsv.exists():
            try:
                hits = pd.read_csv(
                    out_tsv,
                    sep="\t",
                    header=None,
                    names=["target", "pident", "tseq"],
                )
            except EmptyDataError:
                hits = pd.DataFrame(columns=["target", "pident", "tseq"])
            if not hits.empty:
                hits["pident"] = pd.to_numeric(hits["pident"], errors="coerce").fillna(
                    0.0
                )
                hits["tseq"] = hits["tseq"].astype(str)
                for _, hit in hits.iterrows():
                    sequence = str(hit["tseq"])
                    similarity = float(hit["pident"])
                    current = similarity_by_sequence.get(sequence)
                    if current is None or similarity > current:
                        similarity_by_sequence[sequence] = similarity

    return tuple(
        (sequence, float(similarity_by_sequence.get(sequence, 0.0)))
        for sequence in target_sequences
    )


@lru_cache(maxsize=1)
def _resolved_mmseqs_bin_token() -> str:
    return _resolve_mmseqs_bin() or ""


def _write_fasta_records(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for record_id, sequence in records:
            handle.write(f">{record_id}\n{sequence}\n")


def _pdb_ligand_similarity_rows(
    *,
    pdb_id: str,
    query_smiles: str,
    boltz_cache_path: Path,
    components_cif_path: Path | None,
    exclude_ligand_ids: set[str] | None = None,
    timeout: int = 30,
) -> list[dict[str, object]]:
    query_fp = _morgan_fp_from_smiles(query_smiles)
    if query_fp is None:
        return []
    excluded_ids = {
        _norm_id(ligand_id)
        for ligand_id in (exclude_ligand_ids or set())
        if str(ligand_id).strip()
    }
    boltz_cache_token = _path_cache_token(boltz_cache_path)
    components_cif_token = _path_cache_token(components_cif_path)
    rows: list[dict[str, object]] = []
    for ligand_id, smiles in _cached_pdb_ligand_entries(
        str(pdb_id).strip().upper(),
        boltz_cache_token,
        components_cif_token,
        timeout=timeout,
    ):
        if _norm_id(ligand_id) in excluded_ids:
            continue
        ligand_fp = _morgan_fp_from_smiles(smiles)
        if ligand_fp is None:
            continue
        rows.append(
            {
                "ligand_id": ligand_id,
                "smiles": smiles,
                "ecfp_similarity": float(
                    DataStructs.TanimotoSimilarity(query_fp, ligand_fp)
                ),
            }
        )
    return rows


def _pdb_protein_similarity_rows(
    *,
    pdb_id: str,
    query_sequence: str,
    timeout: int = 30,
) -> list[dict[str, object]]:
    if not query_sequence:
        return []
    return [
        {
            "sequence": sequence,
            "sequence_similarity": similarity,
            "sequence_similarity_pairwise": pd.NA,
            "sequence_similarity_method": "mmseqs_pident",
        }
        for sequence, similarity in _cached_pdb_protein_similarity_rows(
            str(pdb_id).strip().upper(),
            str(query_sequence),
            timeout=timeout,
        )
    ]


def _write_progress_checkpoint(
    *,
    rows: list[dict[str, object]],
    pending_rows: list[dict[str, object]],
    output_path: Path | None,
    finalize: Callable[[pd.DataFrame], pd.DataFrame],
    progress_label: str,
    completed_lookups: int,
    total_lookups: int,
) -> None:
    logger = logging.getLogger("cofolder.modules.analytics.bias")
    if output_path is None:
        logger.info(
            "Bias lookup progress for %s: %d/%d lookups completed",
            progress_label,
            completed_lookups,
            total_lookups,
        )
        return

    checkpoint_df = finalize(pd.DataFrame([*rows, *pending_rows]))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    checkpoint_df.to_csv(output_path, index=False)
    logger.info(
        "Bias lookup progress for %s: %d/%d lookups completed; wrote checkpoint to %s",
        progress_label,
        completed_lookups,
        total_lookups,
        output_path,
    )


def _log_bias_lookup_progress(
    *,
    progress_label: str,
    completed_lookups: int,
    total_lookups: int,
) -> None:
    logging.getLogger("cofolder.modules.analytics.bias").info(
        "Bias lookup progress for %s: %d/%d lookups completed",
        progress_label,
        completed_lookups,
        total_lookups,
    )


def _flatten_resolved_bias_rows(
    resolved_rows: dict[int, list[dict[str, object]]],
) -> list[dict[str, object]]:
    flattened: list[dict[str, object]] = []
    for index in sorted(resolved_rows):
        flattened.extend(resolved_rows[index])
    return flattened


def _pending_bias_rows(
    *,
    resolved_rows: dict[int, list[dict[str, object]]],
    input_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    return [
        row_dict
        for index, row_dict in enumerate(input_rows)
        if index not in resolved_rows
    ]


def _enrich_mixed_bias_row_with_pdb_backfill(
    *,
    row_dict: dict[str, object],
    protein_queries: dict[str, str | None],
    ligand_queries: dict[str, str | None],
    protein_lookup_index: dict[str, list[dict[str, object]]] | None = None,
    boltz_cache_path: Path,
    components_cif_path: Path | None,
    ligand_reference_df: pd.DataFrame | None = None,
    ligand_reference_index: dict[str, list[dict[str, str]]] | None = None,
    timeout: int = 30,
) -> list[dict[str, object]]:
    row = pd.Series(row_dict)
    pdb_id = _reference_pdb_id_from_bias_row(row)
    if not pdb_id:
        return [row_dict]

    generated_rows: list[dict[str, object]] = []
    if row.get("pairing_status") == BIAS_PAIRING_PROTEIN_ONLY and pd.isna(
        row.get("ligand_pdb_id")
    ):
        ligand_query_id = str(row.get("query_ligand_chain_id")).strip()
        ligand_query_smiles = ligand_queries.get(ligand_query_id)
        if ligand_query_smiles:
            ligand_matches = _reference_ligand_similarity_rows(
                pdb_id=pdb_id,
                query_smiles=ligand_query_smiles,
                ligand_reference_df=ligand_reference_df,
                ligand_reference_index=ligand_reference_index,
            )
            if not ligand_matches:
                ligand_matches = _pdb_ligand_similarity_rows(
                    pdb_id=pdb_id,
                    query_smiles=ligand_query_smiles,
                    boltz_cache_path=boltz_cache_path,
                    components_cif_path=components_cif_path,
                    timeout=timeout,
                )
            for ligand_match in ligand_matches:
                updated_row = dict(row_dict)
                updated_row["pairing_status"] = BIAS_PAIRING_PAIRED
                updated_row["pdb_id"] = pdb_id
                updated_row["ligand_pdb_id"] = pdb_id
                updated_row["ligand_release_date"] = updated_row.get(
                    "protein_release_date", pd.NA
                )
                updated_row["ligand_source"] = updated_row.get("ligand_source", pd.NA)
                if (
                    pd.isna(updated_row["ligand_source"])
                    or not str(updated_row["ligand_source"]).strip()
                ):
                    updated_row["ligand_source"] = "public"
                updated_row["ligand_dataset_name"] = updated_row.get(
                    "ligand_dataset_name", pd.NA
                )
                if (
                    pd.isna(updated_row["ligand_dataset_name"])
                    or not str(updated_row["ligand_dataset_name"]).strip()
                ):
                    updated_row["ligand_dataset_name"] = "public"
                updated_row["ligand_id"] = ligand_match["ligand_id"]
                updated_row["smiles"] = ligand_match["smiles"]
                updated_row["ecfp_similarity"] = ligand_match["ecfp_similarity"]
                updated_row["plot_ecfp_similarity"] = float(
                    ligand_match["ecfp_similarity"]
                )
                updated_row["source"] = _collapse_pair_value(
                    [
                        updated_row.get("protein_source"),
                        updated_row.get("ligand_source"),
                    ],
                    mixed_label="mixed",
                )
                updated_row["dataset_name"] = _collapse_pair_value(
                    [
                        updated_row.get("protein_dataset_name"),
                        updated_row.get("ligand_dataset_name"),
                    ]
                )
                generated_rows.append(updated_row)

    if row.get("pairing_status") == BIAS_PAIRING_LIGAND_ONLY and pd.isna(
        row.get("protein_pdb_id")
    ):
        protein_query_id = str(row.get("query_protein_chain_id")).strip()
        protein_query_sequence = protein_queries.get(protein_query_id)
        if protein_query_sequence:
            protein_matches = _reference_protein_similarity_rows(
                pdb_id=pdb_id,
                protein_lookup_index=protein_lookup_index,
            )
            used_protein_lookup = bool(protein_matches)
            if not protein_matches:
                if not _resolved_mmseqs_bin_token():
                    logging.getLogger("cofolder.modules.analytics.bias").warning(
                        "MMseqs protein fallback unavailable for mixed bias row %s "
                        "(query_protein_chain_id=%s, pdb_id=%s); leaving row ligand_only",
                        str(row.get("query_pair_id")).strip(),
                        protein_query_id,
                        pdb_id,
                    )
                    return [row_dict]
                protein_matches = _pdb_protein_similarity_rows(
                    pdb_id=pdb_id,
                    query_sequence=protein_query_sequence,
                    timeout=timeout,
                )

            for protein_match in protein_matches:
                sequence_similarity = pd.to_numeric(
                    pd.Series([protein_match.get("sequence_similarity")]),
                    errors="coerce",
                ).iloc[0]
                if pd.isna(sequence_similarity):
                    continue
                updated_row = dict(row_dict)
                updated_row["pairing_status"] = BIAS_PAIRING_PAIRED
                updated_row["pdb_id"] = pdb_id
                updated_row["protein_pdb_id"] = pdb_id
                updated_row["protein_release_date"] = protein_match.get(
                    "release_date",
                    updated_row.get("ligand_release_date", pd.NA),
                )
                updated_row["protein_source"] = protein_match.get(
                    "source",
                    updated_row.get("protein_source", pd.NA),
                )
                if (
                    pd.isna(updated_row["protein_source"])
                    or not str(updated_row["protein_source"]).strip()
                ):
                    updated_row["protein_source"] = "public"
                updated_row["protein_dataset_name"] = protein_match.get(
                    "dataset_name",
                    updated_row.get("protein_dataset_name", pd.NA),
                )
                if (
                    pd.isna(updated_row["protein_dataset_name"])
                    or not str(updated_row["protein_dataset_name"]).strip()
                ):
                    updated_row["protein_dataset_name"] = "public"
                updated_row["protein_source_structure_path"] = protein_match.get(
                    "source_structure_path",
                    updated_row.get("protein_source_structure_path", pd.NA),
                )
                updated_row["protein_source_reference_path"] = protein_match.get(
                    "source_reference_path",
                    updated_row.get("protein_source_reference_path", pd.NA),
                )
                updated_row["sequence"] = protein_match["sequence"]
                updated_row["sequence_similarity"] = float(sequence_similarity)
                updated_row["sequence_similarity_pairwise"] = pd.NA
                updated_row["sequence_similarity_method"] = "mmseqs_pident"
                updated_row["plot_sequence_similarity"] = (
                    float(sequence_similarity) / 100.0
                )
                updated_row["plot_sequence_similarity_pairwise"] = pd.NA
                updated_row["source"] = _collapse_pair_value(
                    [
                        updated_row.get("protein_source"),
                        updated_row.get("ligand_source"),
                    ],
                    mixed_label="mixed",
                )
                updated_row["dataset_name"] = _collapse_pair_value(
                    [
                        updated_row.get("protein_dataset_name"),
                        updated_row.get("ligand_dataset_name"),
                    ]
                )
                if (
                    not used_protein_lookup
                    and float(sequence_similarity) > BIAS_PROTEIN_VIEW_MIN_SIMILARITY
                ):
                    logging.getLogger("cofolder.modules.analytics.bias").warning(
                        "Protein MMseqs fallback exceeded reference threshold for mixed bias row %s "
                        "(query_protein_chain_id=%s, pdb_id=%s, sequence_similarity=%.3f); "
                        "PDB absent from protein_training_data reference view",
                        str(row.get("query_pair_id")).strip(),
                        protein_query_id,
                        pdb_id,
                        float(sequence_similarity),
                    )
                generated_rows.append(updated_row)

    return generated_rows or [row_dict]


def _enrich_same_type_ligand_pair_row_with_pdb_backfill(
    *,
    row_dict: dict[str, object],
    ligand_queries: dict[str, str | None],
    boltz_cache_path: Path,
    components_cif_path: Path | None,
    ligand_reference_df: pd.DataFrame | None = None,
    ligand_reference_index: dict[str, list[dict[str, str]]] | None = None,
    timeout: int = 30,
) -> list[dict[str, object]]:
    row = pd.Series(row_dict)
    pdb_id = _reference_pdb_id_from_same_type_row(row)
    if not pdb_id:
        return [row_dict]

    generated_rows: list[dict[str, object]] = []
    if row.get("pairing_status") == BIAS_PAIRING_COMPONENT_1_ONLY:
        component_1_similarity = pd.to_numeric(
            pd.Series([row.get("component_1_similarity")]),
            errors="coerce",
        ).iloc[0]
        component_2_query_id = str(row.get("component_2_id")).strip()
        component_2_query_smiles = ligand_queries.get(component_2_query_id)
        if (
            pd.notna(component_1_similarity)
            and float(component_1_similarity) > BIAS_LIGAND_VIEW_MIN_SIMILARITY
            and component_2_query_smiles
        ):
            exclude_ligand_ids = (
                {str(row.get("_component_1_ligand_id")).strip()}
                if pd.notna(row.get("_component_1_ligand_id"))
                and str(row.get("_component_1_ligand_id")).strip()
                else set()
            )
            ligand_matches = _reference_ligand_similarity_rows(
                pdb_id=pdb_id,
                query_smiles=component_2_query_smiles,
                ligand_reference_df=ligand_reference_df,
                ligand_reference_index=ligand_reference_index,
                exclude_ligand_ids=exclude_ligand_ids,
            )
            if not ligand_matches:
                ligand_matches = _pdb_ligand_similarity_rows(
                    pdb_id=str(pdb_id),
                    query_smiles=str(component_2_query_smiles),
                    boltz_cache_path=boltz_cache_path,
                    components_cif_path=components_cif_path,
                    exclude_ligand_ids=set(exclude_ligand_ids),
                    timeout=timeout,
                )
            for ligand_match in ligand_matches:
                updated_row = dict(row_dict)
                updated_row["pairing_status"] = BIAS_PAIRING_PAIRED
                updated_row["reference_pdb_id"] = pdb_id
                if (
                    pd.isna(updated_row.get("component_2_reference_key"))
                    or not str(updated_row.get("component_2_reference_key")).strip()
                ):
                    updated_row["component_2_reference_key"] = updated_row.get(
                        "reference_key"
                    )
                updated_row["component_2_pdb_id"] = pdb_id
                updated_row["component_2_release_date"] = updated_row.get(
                    "component_1_release_date", pd.NA
                )
                if (
                    pd.isna(updated_row.get("component_2_source"))
                    or not str(updated_row.get("component_2_source")).strip()
                ):
                    updated_row["component_2_source"] = "public"
                if (
                    pd.isna(updated_row.get("component_2_dataset_name"))
                    or not str(updated_row.get("component_2_dataset_name")).strip()
                ):
                    updated_row["component_2_dataset_name"] = "public"
                updated_row["component_2_similarity"] = ligand_match["ecfp_similarity"]
                updated_row["plot_component_2_similarity"] = float(
                    ligand_match["ecfp_similarity"]
                )
                updated_row["_component_2_ligand_id"] = ligand_match["ligand_id"]
                updated_row["_component_2_smiles"] = ligand_match["smiles"]
                updated_row["source"] = _collapse_pair_value(
                    [
                        updated_row.get("component_1_source"),
                        updated_row.get("component_2_source"),
                    ],
                    mixed_label="mixed",
                )
                updated_row["dataset_name"] = _collapse_pair_value(
                    [
                        updated_row.get("component_1_dataset_name"),
                        updated_row.get("component_2_dataset_name"),
                    ]
                )
                generated_rows.append(updated_row)

    if row.get("pairing_status") == BIAS_PAIRING_COMPONENT_2_ONLY:
        component_2_similarity = pd.to_numeric(
            pd.Series([row.get("component_2_similarity")]),
            errors="coerce",
        ).iloc[0]
        component_1_query_id = str(row.get("component_1_id")).strip()
        component_1_query_smiles = ligand_queries.get(component_1_query_id)
        if (
            pd.notna(component_2_similarity)
            and float(component_2_similarity) > BIAS_LIGAND_VIEW_MIN_SIMILARITY
            and component_1_query_smiles
        ):
            exclude_ligand_ids = (
                {str(row.get("_component_2_ligand_id")).strip()}
                if pd.notna(row.get("_component_2_ligand_id"))
                and str(row.get("_component_2_ligand_id")).strip()
                else set()
            )
            ligand_matches = _reference_ligand_similarity_rows(
                pdb_id=pdb_id,
                query_smiles=component_1_query_smiles,
                ligand_reference_df=ligand_reference_df,
                ligand_reference_index=ligand_reference_index,
                exclude_ligand_ids=exclude_ligand_ids,
            )
            if not ligand_matches:
                ligand_matches = _pdb_ligand_similarity_rows(
                    pdb_id=str(pdb_id),
                    query_smiles=str(component_1_query_smiles),
                    boltz_cache_path=boltz_cache_path,
                    components_cif_path=components_cif_path,
                    exclude_ligand_ids=set(exclude_ligand_ids),
                    timeout=timeout,
                )
            for ligand_match in ligand_matches:
                updated_row = dict(row_dict)
                updated_row["pairing_status"] = BIAS_PAIRING_PAIRED
                updated_row["reference_pdb_id"] = pdb_id
                if (
                    pd.isna(updated_row.get("component_1_reference_key"))
                    or not str(updated_row.get("component_1_reference_key")).strip()
                ):
                    updated_row["component_1_reference_key"] = updated_row.get(
                        "reference_key"
                    )
                updated_row["component_1_pdb_id"] = pdb_id
                updated_row["component_1_release_date"] = updated_row.get(
                    "component_2_release_date", pd.NA
                )
                if (
                    pd.isna(updated_row.get("component_1_source"))
                    or not str(updated_row.get("component_1_source")).strip()
                ):
                    updated_row["component_1_source"] = "public"
                if (
                    pd.isna(updated_row.get("component_1_dataset_name"))
                    or not str(updated_row.get("component_1_dataset_name")).strip()
                ):
                    updated_row["component_1_dataset_name"] = "public"
                updated_row["component_1_similarity"] = ligand_match["ecfp_similarity"]
                updated_row["plot_component_1_similarity"] = float(
                    ligand_match["ecfp_similarity"]
                )
                updated_row["_component_1_ligand_id"] = ligand_match["ligand_id"]
                updated_row["_component_1_smiles"] = ligand_match["smiles"]
                updated_row["source"] = _collapse_pair_value(
                    [
                        updated_row.get("component_1_source"),
                        updated_row.get("component_2_source"),
                    ],
                    mixed_label="mixed",
                )
                updated_row["dataset_name"] = _collapse_pair_value(
                    [
                        updated_row.get("component_1_dataset_name"),
                        updated_row.get("component_2_dataset_name"),
                    ]
                )
                generated_rows.append(updated_row)

    return generated_rows or [row_dict]


def _enrich_mixed_bias_dataset_with_pdb_backfill(
    df: pd.DataFrame,
    *,
    protein_queries: dict[str, str | None],
    ligand_queries: dict[str, str | None],
    protein_lookup_view: pd.DataFrame | None = None,
    boltz_cache_path: Path,
    components_cif_path: Path | None,
    ligand_reference_df: pd.DataFrame | None = None,
    ligand_reference_index: dict[str, list[dict[str, str]]] | None = None,
    progress_output_path: Path | None = None,
    progress_label: str | None = None,
    timeout: int = 30,
) -> pd.DataFrame:
    if df.empty:
        return _finalize_bias_training_dataframe(df)

    input_rows = [row.to_dict() for _, row in df.iterrows()]
    total_lookups = 0
    for row_dict in input_rows:
        row = pd.Series(row_dict)
        pdb_id = _reference_pdb_id_from_bias_row(row)
        if (
            row.get("pairing_status") == BIAS_PAIRING_PROTEIN_ONLY
            and pd.isna(row.get("ligand_pdb_id"))
            and pdb_id
            and ligand_queries.get(str(row.get("query_ligand_chain_id")).strip())
        ):
            total_lookups += 1
        elif (
            row.get("pairing_status") == BIAS_PAIRING_LIGAND_ONLY
            and pd.isna(row.get("protein_pdb_id"))
            and pdb_id
            and protein_queries.get(str(row.get("query_protein_chain_id")).strip())
        ):
            total_lookups += 1

    progress_name = progress_label or "mixed_bias_backfill"
    protein_lookup_index = _build_protein_lookup_index(protein_lookup_view)
    _log_bias_lookup_progress(
        progress_label=progress_name,
        completed_lookups=0,
        total_lookups=total_lookups,
    )

    enriched_rows: list[dict[str, object]] = []
    completed_lookups = 0
    for index, row_dict in enumerate(input_rows):
        row = pd.Series(row_dict)
        pdb_id = _reference_pdb_id_from_bias_row(row)
        lookup_attempted = False

        if (
            row.get("pairing_status") == BIAS_PAIRING_PROTEIN_ONLY
            and pd.isna(row.get("ligand_pdb_id"))
            and pdb_id
            and ligand_queries.get(str(row.get("query_ligand_chain_id")).strip())
        ):
            lookup_attempted = True
        elif (
            row.get("pairing_status") == BIAS_PAIRING_LIGAND_ONLY
            and pd.isna(row.get("protein_pdb_id"))
            and pdb_id
            and protein_queries.get(str(row.get("query_protein_chain_id")).strip())
        ):
            lookup_attempted = True

        enriched_rows.extend(
            _enrich_mixed_bias_row_with_pdb_backfill(
                row_dict=row_dict,
                protein_queries=protein_queries,
                ligand_queries=ligand_queries,
                protein_lookup_index=protein_lookup_index,
                boltz_cache_path=boltz_cache_path,
                components_cif_path=components_cif_path,
                ligand_reference_df=ligand_reference_df,
                ligand_reference_index=ligand_reference_index,
                timeout=timeout,
            )
        )

        if lookup_attempted:
            completed_lookups += 1
            if (
                completed_lookups % BIAS_PROGRESS_LOG_EVERY == 0
                or completed_lookups == total_lookups
            ):
                _log_bias_lookup_progress(
                    progress_label=progress_name,
                    completed_lookups=completed_lookups,
                    total_lookups=total_lookups,
                )
            if (
                completed_lookups % BIAS_PROGRESS_WRITE_EVERY == 0
                or completed_lookups == total_lookups
            ):
                _write_progress_checkpoint(
                    rows=enriched_rows,
                    pending_rows=input_rows[index + 1 :],
                    output_path=progress_output_path,
                    finalize=_finalize_bias_training_dataframe,
                    progress_label=progress_name,
                    completed_lookups=completed_lookups,
                    total_lookups=total_lookups,
                )

    return _finalize_bias_training_dataframe(pd.DataFrame(enriched_rows))


def _enrich_same_type_ligand_pair_dataset_with_pdb_backfill(
    df: pd.DataFrame,
    *,
    ligand_queries: dict[str, str | None],
    boltz_cache_path: Path,
    components_cif_path: Path | None,
    ligand_reference_df: pd.DataFrame | None = None,
    ligand_reference_index: dict[str, list[dict[str, str]]] | None = None,
    progress_output_path: Path | None = None,
    progress_label: str | None = None,
    timeout: int = 30,
) -> pd.DataFrame:
    if df.empty:
        return df.copy()

    input_rows = [row.to_dict() for _, row in df.iterrows()]
    total_lookups = 0
    for row_dict in input_rows:
        row = pd.Series(row_dict)
        pdb_id = _reference_pdb_id_from_same_type_row(row)
        if row.get("pairing_status") == BIAS_PAIRING_COMPONENT_1_ONLY and pdb_id:
            component_1_similarity = pd.to_numeric(
                pd.Series([row.get("component_1_similarity")]),
                errors="coerce",
            ).iloc[0]
            component_2_query_id = str(row.get("component_2_id")).strip()
            if (
                pd.notna(component_1_similarity)
                and float(component_1_similarity) > BIAS_LIGAND_VIEW_MIN_SIMILARITY
                and ligand_queries.get(component_2_query_id)
            ):
                total_lookups += 1
        elif row.get("pairing_status") == BIAS_PAIRING_COMPONENT_2_ONLY and pdb_id:
            component_2_similarity = pd.to_numeric(
                pd.Series([row.get("component_2_similarity")]),
                errors="coerce",
            ).iloc[0]
            component_1_query_id = str(row.get("component_1_id")).strip()
            if (
                pd.notna(component_2_similarity)
                and float(component_2_similarity) > BIAS_LIGAND_VIEW_MIN_SIMILARITY
                and ligand_queries.get(component_1_query_id)
            ):
                total_lookups += 1

    progress_name = progress_label or "ligand_pair_backfill"
    _log_bias_lookup_progress(
        progress_label=progress_name,
        completed_lookups=0,
        total_lookups=total_lookups,
    )

    enriched_rows: list[dict[str, object]] = []
    completed_lookups = 0
    for index, row_dict in enumerate(input_rows):
        row = pd.Series(row_dict)
        pdb_id = _reference_pdb_id_from_same_type_row(row)
        lookup_attempted = False

        if row.get("pairing_status") == BIAS_PAIRING_COMPONENT_1_ONLY and pdb_id:
            component_1_similarity = pd.to_numeric(
                pd.Series([row.get("component_1_similarity")]),
                errors="coerce",
            ).iloc[0]
            component_2_query_id = str(row.get("component_2_id")).strip()
            if (
                pd.notna(component_1_similarity)
                and float(component_1_similarity) > BIAS_LIGAND_VIEW_MIN_SIMILARITY
                and ligand_queries.get(component_2_query_id)
            ):
                lookup_attempted = True
        elif row.get("pairing_status") == BIAS_PAIRING_COMPONENT_2_ONLY and pdb_id:
            component_2_similarity = pd.to_numeric(
                pd.Series([row.get("component_2_similarity")]),
                errors="coerce",
            ).iloc[0]
            component_1_query_id = str(row.get("component_1_id")).strip()
            if (
                pd.notna(component_2_similarity)
                and float(component_2_similarity) > BIAS_LIGAND_VIEW_MIN_SIMILARITY
                and ligand_queries.get(component_1_query_id)
            ):
                lookup_attempted = True

        enriched_rows.extend(
            _enrich_same_type_ligand_pair_row_with_pdb_backfill(
                row_dict=row_dict,
                ligand_queries=ligand_queries,
                boltz_cache_path=boltz_cache_path,
                components_cif_path=components_cif_path,
                ligand_reference_df=ligand_reference_df,
                ligand_reference_index=ligand_reference_index,
                timeout=timeout,
            )
        )

        if lookup_attempted:
            completed_lookups += 1
            if (
                completed_lookups % BIAS_PROGRESS_LOG_EVERY == 0
                or completed_lookups == total_lookups
            ):
                _log_bias_lookup_progress(
                    progress_label=progress_name,
                    completed_lookups=completed_lookups,
                    total_lookups=total_lookups,
                )
            if (
                completed_lookups % BIAS_PROGRESS_WRITE_EVERY == 0
                or completed_lookups == total_lookups
            ):
                _write_progress_checkpoint(
                    rows=enriched_rows,
                    pending_rows=input_rows[index + 1 :],
                    output_path=progress_output_path,
                    finalize=_finalize_same_type_pair_dataframe,
                    progress_label=progress_name,
                    completed_lookups=completed_lookups,
                    total_lookups=total_lookups,
                )

    return pd.DataFrame(enriched_rows)


__all__ = [
    "_cached_pdb_ligand_entries",
    "_cached_pdb_protein_similarity_rows",
    "_enrich_mixed_bias_dataset_with_pdb_backfill",
    "_enrich_mixed_bias_row_with_pdb_backfill",
    "_enrich_same_type_ligand_pair_dataset_with_pdb_backfill",
    "_enrich_same_type_ligand_pair_row_with_pdb_backfill",
    "_entry_fasta_sequences",
    "_entry_ligand_ids",
    "_extract_nonpoly_comp_id",
    "_fetch_json",
    "_fetch_text",
    "_flatten_resolved_bias_rows",
    "_log_bias_lookup_progress",
    "_pdb_ligand_similarity_rows",
    "_pdb_protein_similarity_rows",
    "_pending_bias_rows",
    "_resolve_pdb_ligand_smiles",
    "_resolved_mmseqs_bin_token",
    "_write_fasta_records",
    "_write_progress_checkpoint",
]
