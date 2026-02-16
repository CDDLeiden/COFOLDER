#!/usr/bin/env python3
"""Build bias training CSVs from CCD-first hit expansion.

Workflow:
1) Read ligand queries from the target `system.yaml`.
2) Compare query ligands to CCD (components.cif) by ECFP/Tanimoto.
3) Keep CCD hits above a threshold.
4) Fetch only PDB entries that contain those CCD hits.
5) Filter entries by release cutoff and write:
   - protein_training_data.csv (pdb_id, release_date, sequence)
   - ligand_training_data.csv (pdb_id, release_date, ligand_id, smiles)

This implements the intended bias bootstrap strategy:
- Start from CCD ligand space first.
- Expand to PDB entries only for identified ligand hits.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import date
import json
from pathlib import Path
import re
import shutil
import subprocess
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import pandas as pd
from Bio.Align import PairwiseAligner
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit import RDLogger
import yaml
import gemmi

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
CORE_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
CORE_NONPOLY_URL = "https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{pdb_id}/{entity_id}"
FASTA_URL = "https://www.rcsb.org/fasta/entry/{pdb_id}/download"
MORGAN_FP = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
# Silence verbose RDKit parser noise for invalid CCD descriptors we intentionally skip.
RDLogger.DisableLog("rdApp.error")
ALIGNER = PairwiseAligner(mode="global")
ALIGNER.match_score = 1.0
ALIGNER.mismatch_score = 0.0
ALIGNER.open_gap_score = 0.0
ALIGNER.extend_gap_score = 0.0


@dataclass(frozen=True)
class LigandQuery:
    chain_ids: tuple[str, ...]
    smiles: str


def _fetch_json(url: str, timeout: int) -> dict | list:
    with urlopen(url, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


def _post_json(url: str, payload: dict, timeout: int) -> dict:
    req = Request(
        url=url,
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    with urlopen(req, timeout=timeout) as resp:
        body = resp.read().decode("utf-8", errors="replace")
        try:
            parsed = json.loads(body)
        except json.JSONDecodeError as e:
            snippet = body[:200].replace("\n", " ")
            raise ValueError(
                f"Non-JSON response from {url} (first 200 chars: {snippet!r})"
            ) from e
        if not isinstance(parsed, dict):
            raise ValueError(f"Unexpected JSON response type from {url}: {type(parsed).__name__}")
        return parsed


def _fetch_text(url: str, timeout: int) -> str:
    with urlopen(url, timeout=timeout) as resp:
        return resp.read().decode("utf-8", errors="replace")


def _parse_iso_date(value: str | None) -> date | None:
    if not value:
        return None
    try:
        return date.fromisoformat(str(value)[:10])
    except Exception:
        return None


def _as_list(value) -> list[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return [str(x) for x in value if str(x)]
    return [str(value)]


def _extract_single_ligand_query(
    system_path: Path,
    ccd_smiles_by_id: dict[str, str],
) -> LigandQuery | None:
    data = yaml.safe_load(system_path.read_text(encoding="utf-8"))
    seqs = data.get("sequences") if isinstance(data, dict) else None
    if not isinstance(seqs, list):
        return None

    # Single-query mode: first ligand definition with usable smiles (direct or CCD-mapped).
    for item in seqs:
        if not isinstance(item, dict):
            continue
        lig = item.get("ligand")
        if not isinstance(lig, dict):
            continue
        chain_ids = _as_list(lig.get("id"))
        smiles = _normalize_smiles(lig.get("smiles"))
        if not smiles:
            # Also support ligand definitions that only provide CCD IDs.
            ccd_ids = [x.upper() for x in _as_list(lig.get("ccd"))]
            for ccd_id in ccd_ids:
                mapped = _normalize_smiles(ccd_smiles_by_id.get(ccd_id))
                if not mapped:
                    continue
                return LigandQuery(
                    chain_ids=tuple(chain_ids),
                    smiles=mapped,
                )
            continue
        return LigandQuery(
            chain_ids=tuple(chain_ids),
            smiles=smiles,
        )

    return None


def _load_ccd_smiles(components_cif_path: Path) -> pd.DataFrame:
    """Load CCD smiles from components.cif using gemmi.

    Preference order per CCD ID:
    1) SMILES_CANONICAL
    2) SMILES
    """
    doc = gemmi.cif.read_file(str(components_cif_path))
    rows: list[dict[str, str]] = []

    for block in doc:
        ccd_id = (block.find_value("_chem_comp.id") or "").strip().upper()
        if not ccd_id:
            continue

        tab = block.find(
            "_pdbx_chem_comp_descriptor.",
            ["type", "descriptor"],
        )

        smiles_canonical = None
        smiles_generic = None
        if tab:
            for row in tab:
                d_type = str(row[0]).strip().upper()
                desc = _normalize_smiles(str(row[1]))
                if not desc or desc in {"?", "."}:
                    continue
                if d_type == "SMILES_CANONICAL":
                    smiles_canonical = desc
                    break
                if d_type == "SMILES" and smiles_generic is None:
                    smiles_generic = desc

        smiles = smiles_canonical or smiles_generic
        if not smiles:
            continue

        rows.append({"ligand_id": ccd_id, "smiles": smiles})

    out = pd.DataFrame(rows).drop_duplicates(subset=["ligand_id"], keep="first")
    return out


def _normalize_smiles(smiles: str | None) -> str:
    if smiles is None:
        return ""
    s = str(smiles).strip()
    # Some CCD descriptors are quoted; strip one or more matching wrappers.
    while len(s) >= 2 and (
        (s[0] == '"' and s[-1] == '"') or (s[0] == "'" and s[-1] == "'")
    ):
        s = s[1:-1].strip()
    # Some CCD descriptors contain trailing semicolon wrappers.
    s = s.strip(";").strip()
    return s


def _looks_like_supported_smiles(smiles: str) -> bool:
    if not smiles:
        return False
    # Exclude common non-RDKit descriptor patterns seen in CCD.
    if "|" in smiles:
        return False
    if "\n" in smiles or "\r" in smiles:
        return False
    if smiles.startswith(";") or smiles.endswith(";"):
        return False
    return True


def _morgan_fp(smiles: str):
    clean = _normalize_smiles(smiles)
    if not _looks_like_supported_smiles(clean):
        return None
    mol = Chem.MolFromSmiles(clean)
    if mol is None:
        return None
    return MORGAN_FP.GetFingerprint(mol)


def _find_ccd_hits(
    query_smiles: str,
    ccd_df: pd.DataFrame,
    threshold: float,
    top_k: int,
) -> list[tuple[str, str, float]]:
    qfp = _morgan_fp(query_smiles)
    if qfp is None:
        return []

    scored: list[tuple[str, str, float]] = []
    for _, row in ccd_df.iterrows():
        ccd_id = str(row["ligand_id"])
        smiles = _normalize_smiles(str(row["smiles"]))
        fp = row.get("fp")
        if fp is None:
            fp = _morgan_fp(smiles)
        if fp is None:
            continue
        sim = float(DataStructs.TanimotoSimilarity(qfp, fp))
        if sim >= threshold:
            scored.append((ccd_id, smiles, sim))

    scored.sort(key=lambda x: x[2], reverse=True)
    if top_k > 0:
        scored = scored[:top_k]
    return scored


def _search_entries_for_ccd(ccd_id: str, timeout: int) -> list[str]:
    payload = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id",
                "operator": "exact_match",
                "value": ccd_id,
            },
        },
        "return_type": "entry",
        "request_options": {
            "results_verbosity": "compact",
            "paginate": {"start": 0, "rows": 10000},
        },
    }
    res = _post_json(SEARCH_URL, payload, timeout=timeout)
    out: list[str] = []
    result_set = res.get("result_set", [])
    if not isinstance(result_set, list):
        return []
    for item in result_set:
        if isinstance(item, dict):
            ident = str(item.get("identifier", "")).upper()
        else:
            ident = str(item).upper()
        if re.fullmatch(r"[A-Z0-9]{4}", ident):
            out.append(ident)
    return sorted(set(out))


def _entry_release_date(pdb_id: str, timeout: int) -> date | None:
    payload = _fetch_json(CORE_ENTRY_URL.format(pdb_id=pdb_id), timeout=timeout)
    accession = payload.get("rcsb_accession_info", {}) if isinstance(payload, dict) else {}
    return _parse_iso_date(accession.get("initial_release_date"))


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


def _entry_ligand_ids(pdb_id: str, timeout: int) -> list[str]:
    try:
        payload = _fetch_json(CORE_ENTRY_URL.format(pdb_id=pdb_id), timeout=timeout)
    except Exception:
        return []
    ids_obj = payload.get("rcsb_entry_container_identifiers", {}) if isinstance(payload, dict) else {}
    entity_ids = ids_obj.get("non_polymer_entity_ids", []) if isinstance(ids_obj, dict) else []
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
    return sorted(out)


def _entry_fasta_sequences(pdb_id: str, timeout: int) -> list[str]:
    txt = _fetch_text(FASTA_URL.format(pdb_id=pdb_id), timeout=timeout)
    seqs: list[str] = []
    buf: list[str] = []
    for line in txt.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if buf:
                seqs.append("".join(buf))
                buf = []
            continue
        buf.append(line)
    if buf:
        seqs.append("".join(buf))
    return [s for s in seqs if s]


def _sequence_identity_percent(query: str, target: str) -> float:
    if not query or not target:
        return 0.0
    score = ALIGNER.score(query, target)
    denom = max(len(query), len(target))
    if denom == 0:
        return 0.0
    return float(100.0 * float(score) / float(denom))


def _read_single_fasta_sequence(path: Path) -> str:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    seq_parts: list[str] = []
    for line in lines:
        line = line.strip()
        if not line or line.startswith(">"):
            continue
        seq_parts.append(line)
    return "".join(seq_parts)


def _resolve_query_sequence(seq_arg: str | None, fasta_arg: Path | None) -> str | None:
    if seq_arg:
        return str(seq_arg).strip()
    if fasta_arg is not None:
        seq = _read_single_fasta_sequence(fasta_arg).strip()
        return seq or None
    return None


def _extract_protein_query_from_system(system_path: Path) -> str | None:
    data = yaml.safe_load(system_path.read_text(encoding="utf-8"))
    seqs = data.get("sequences") if isinstance(data, dict) else None
    if not isinstance(seqs, list):
        return None

    for item in seqs:
        if not isinstance(item, dict):
            continue
        prot = item.get("protein")
        if not isinstance(prot, dict):
            continue
        seq = (prot.get("fasta") or prot.get("sequence") or "").strip()
        if seq:
            return seq
    return None


def _guess_mmseqs_db_from_components(components_cif: Path) -> Path:
    # expected fetch layout:
    # <output_root>/ccd/components.cif
    # <output_root>/mmseqs/pdb/db
    output_root = components_cif.parent.parent
    return output_root / "mmseqs" / "pdb" / "db"


def _run_mmseqs(
    mmseqs_bin: str,
    query_fasta: Path,
    target_db: Path,
    workers: int,
    tmp_root: Path,
    max_seqs: int,
) -> pd.DataFrame:
    if shutil.which(mmseqs_bin) is None and not Path(mmseqs_bin).exists():
        raise RuntimeError(
            f"mmseqs binary not found: {mmseqs_bin}. Provide --mmseqs_bin with full path."
        )

    tmp_root.mkdir(parents=True, exist_ok=True)
    qdb = tmp_root / "query_db"
    rdb = tmp_root / "result_db"
    out_tsv = tmp_root / "result.tsv"
    search_tmp = tmp_root / "search_tmp"
    # Clean stale MMseqs artifacts so reruns don't fail on existing db files.
    for stem in ("query_db", "result_db"):
        for p in tmp_root.glob(f"{stem}*"):
            if p.is_dir():
                shutil.rmtree(p, ignore_errors=True)
            else:
                p.unlink(missing_ok=True)
    out_tsv.unlink(missing_ok=True)
    shutil.rmtree(search_tmp, ignore_errors=True)
    search_tmp.mkdir(parents=True, exist_ok=True)

    cmds = [
        [mmseqs_bin, "createdb", str(query_fasta), str(qdb)],
        [
            mmseqs_bin,
            "search",
            str(qdb),
            str(target_db),
            str(rdb),
            str(search_tmp),
            "--threads",
            str(max(1, workers)),
            "--max-seqs",
            str(max(1, int(max_seqs))),
        ],
        [
            mmseqs_bin,
            "convertalis",
            str(qdb),
            str(target_db),
            str(rdb),
            str(out_tsv),
            "--format-output",
            "target,pident,tseq",
        ],
    ]
    for cmd in cmds:
        subprocess.run(cmd, check=True)

    if not out_tsv.exists():
        return pd.DataFrame(columns=["target", "pident", "tseq"])

    return pd.read_csv(
        out_tsv,
        sep="\t",
        header=None,
        names=["target", "pident", "tseq"],
    )


def _collect_ccd_hits(
    query: LigandQuery,
    ccd_df: pd.DataFrame,
    threshold: float,
    top_k: int,
) -> pd.DataFrame:
    if "fp" not in ccd_df.columns:
        ccd_df = ccd_df.copy()
        ccd_df["fp"] = ccd_df["smiles"].apply(_morgan_fp)
        ccd_df = ccd_df[ccd_df["fp"].notna()].copy()

    rows: list[dict[str, str | float]] = []
    query_chain_id = query.chain_ids[0] if query.chain_ids else ""
    hits = _find_ccd_hits(
        query_smiles=query.smiles,
        ccd_df=ccd_df,
        threshold=threshold,
        top_k=top_k,
    )
    for ccd_id, hit_smiles, sim in hits:
        rows.append(
            {
                "query_chain_id": query_chain_id,
                "ligand_pdbid": ccd_id,
                "hit_smiles": hit_smiles,
                "ecfp_similarity": sim,
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--system_path", type=Path, required=True)
    parser.add_argument("--components_cif", type=Path, required=True)
    parser.add_argument("--output_protein_csv", type=Path, required=True)
    parser.add_argument("--output_ligand_csv", type=Path, required=True)
    parser.add_argument("--output_bias_csv", type=Path, default=None)
    parser.add_argument("--output_hits_csv", type=Path, default=None)
    parser.add_argument("--release_cutoff", type=str, default="2023-06-01")
    parser.add_argument("--ligand_similarity_threshold", type=float, default=0.35)
    parser.add_argument("--top_k_per_query", type=int, default=200)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--timeout", type=int, default=30)
    parser.add_argument(
        "--query_sequence",
        type=str,
        default=None,
        help="Optional protein query sequence. If provided, protein output is filtered by similarity threshold.",
    )
    parser.add_argument(
        "--query_sequence_fasta",
        type=Path,
        default=None,
        help="Optional FASTA file containing the protein query sequence.",
    )
    parser.add_argument(
        "--protein_similarity_threshold",
        type=float,
        default=25.0,
        help="Minimum sequence identity percentage for protein output.",
    )
    parser.add_argument(
        "--mmseqs_db_path",
        type=Path,
        default=None,
        help="Path to MMseqs target DB (default: inferred from --components_cif as <training_data>/mmseqs/pdb/db).",
    )
    parser.add_argument(
        "--mmseqs_bin",
        type=str,
        default="mmseqs",
        help="MMseqs executable name or full path.",
    )
    parser.add_argument(
        "--mmseqs_max_seqs",
        type=int,
        default=10000,
        help="Maximum number of MMseqs hits to keep per query (default: 10000).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output CSVs (default: false).",
    )
    args = parser.parse_args()

    cutoff = _parse_iso_date(args.release_cutoff)
    if cutoff is None:
        raise ValueError("--release_cutoff must be YYYY-MM-DD")
    if args.query_sequence and args.query_sequence_fasta:
        raise ValueError("Use either --query_sequence or --query_sequence_fasta, not both.")

    ccd_df = _load_ccd_smiles(args.components_cif)
    if ccd_df.empty:
        raise ValueError("No CCD SMILES could be read from components.cif")
    ccd_smiles_by_id = {
        str(row["ligand_id"]).upper(): _normalize_smiles(str(row["smiles"]))
        for _, row in ccd_df.iterrows()
    }

    query = _extract_single_ligand_query(
        args.system_path,
        ccd_smiles_by_id=ccd_smiles_by_id,
    )
    if query is None:
        raise ValueError(
            "No ligand queries found. Provide ligand.smiles and/or ligand.ccd values in system.yaml."
        )

    print("[info] query_mode=single_smiles")
    print(f"[info] query_chain_ids={','.join(query.chain_ids)}")
    print(f"[info] ccd_with_smiles={len(ccd_df)}")

    output_bias_csv = args.output_bias_csv or args.output_ligand_csv.with_name("bias_training_data.csv")

    hits_df = _collect_ccd_hits(
        query=query,
        ccd_df=ccd_df,
        threshold=args.ligand_similarity_threshold,
        top_k=args.top_k_per_query,
    )
    if hits_df.empty:
        raise ValueError(
            "No CCD hits above threshold. Lower --ligand_similarity_threshold or inspect query ligands."
        )

    # Keep only best score per CCD ID across all query ligands.
    hits_df = hits_df.sort_values("ecfp_similarity", ascending=False)
    hits_best = hits_df.drop_duplicates(subset=["ligand_pdbid"], keep="first")
    hit_ccd_ids = sorted(hits_best["ligand_pdbid"].astype(str).unique())

    print(f"[info] ccd_hits={len(hit_ccd_ids)}")

    # For each hit CCD, fetch associated PDB entries.
    ccd_to_entries: dict[str, list[str]] = {}
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as ex:
        fut_map = {
            ex.submit(_search_entries_for_ccd, ccd_id, args.timeout): ccd_id
            for ccd_id in hit_ccd_ids
        }
        done = 0
        for fut in as_completed(fut_map):
            ccd_id = fut_map[fut]
            done += 1
            if done % 50 == 0:
                print(f"[map] ccd_done={done}/{len(fut_map)}", flush=True)
            try:
                ccd_to_entries[ccd_id] = fut.result()
            except Exception as e:
                print(f"[warn] failed ccd->pdb mapping for {ccd_id}: {e}")
                ccd_to_entries[ccd_id] = []

    candidate_pdb_ids = sorted({p for vals in ccd_to_entries.values() for p in vals})
    print(f"[info] candidate_pdb_entries={len(candidate_pdb_ids)}")

    # Fetch release dates and keep pre-cutoff entries only.
    pdb_release: dict[str, date] = {}
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as ex:
        fut_map = {
            ex.submit(_entry_release_date, pdb_id, args.timeout): pdb_id
            for pdb_id in candidate_pdb_ids
        }
        done = 0
        for fut in as_completed(fut_map):
            pdb_id = fut_map[fut]
            done += 1
            if done % 200 == 0:
                print(f"[release] checked={done}/{len(fut_map)}", flush=True)
            try:
                rel = fut.result()
            except (HTTPError, URLError):
                continue
            if rel is not None and rel < cutoff:
                pdb_release[pdb_id] = rel

    selected_pdb_ids = sorted(pdb_release.keys())
    print(f"[info] selected_pdb_entries_pre_cutoff={len(selected_pdb_ids)}")

    # Build ligand training rows (ligand branch only).
    hit_smiles_by_ccd = {
        str(row["ligand_pdbid"]): str(row["hit_smiles"])
        for _, row in hits_best.iterrows()
    }
    ccd_similarity = {
        str(row["ligand_pdbid"]): float(row["ecfp_similarity"])
        for _, row in hits_best.iterrows()
    }
    ligand_rows: list[dict[str, str]] = []
    for ccd_id, pdb_ids in ccd_to_entries.items():
        smiles = hit_smiles_by_ccd.get(ccd_id)
        if not smiles:
            continue
        for pdb_id in pdb_ids:
            rel = pdb_release.get(pdb_id)
            if rel is None:
                continue
            ligand_rows.append(
                {
                    "pdb_id": pdb_id,
                    "release_date": rel.isoformat(),
                    "ligand_id": ccd_id,
                    "smiles": smiles,
                    "ecfp_similarity": ccd_similarity.get(ccd_id),
                }
            )

    ligand_df = pd.DataFrame(ligand_rows).drop_duplicates(
        subset=["pdb_id", "ligand_id", "smiles"]
    )
    if not ligand_df.empty:
        ligand_df = ligand_df.sort_values(
            by=["ecfp_similarity", "pdb_id", "ligand_id"],
            ascending=[False, True, True],
        )
    if args.output_ligand_csv.exists() and not args.overwrite:
        try:
            ligand_df = pd.read_csv(args.output_ligand_csv)
            print(f"[info] reusing existing ligand CSV: {args.output_ligand_csv}")
        except Exception as e:
            print(f"[warn] failed to read existing ligand CSV, using newly computed data: {e}")

    # Build protein training rows independently from MMseqs DB (no ligand filtering).
    query_sequence = _resolve_query_sequence(args.query_sequence, args.query_sequence_fasta)
    if not query_sequence:
        query_sequence = _extract_protein_query_from_system(args.system_path)
        if query_sequence:
            print("[info] using protein query sequence from system.yaml")
    if not query_sequence:
        raise ValueError(
            "Protein query sequence not found. Provide --query_sequence / --query_sequence_fasta "
            "or include protein.fasta/protein.sequence in system.yaml."
        )

    if args.output_protein_csv.exists() and not args.overwrite:
        try:
            protein_df = pd.read_csv(args.output_protein_csv)
            if "sequence_similarity_query" in protein_df.columns and "sequence_similarity" not in protein_df.columns:
                protein_df = protein_df.rename(columns={"sequence_similarity_query": "sequence_similarity"})
            print(f"[info] reusing existing protein CSV: {args.output_protein_csv}")
        except Exception as e:
            print(f"[warn] failed to read existing protein CSV, recomputing with MMseqs: {e}")
            protein_df = pd.DataFrame()
    else:
        protein_df = pd.DataFrame()

    if protein_df.empty:
        mmseqs_db_path = args.mmseqs_db_path or _guess_mmseqs_db_from_components(args.components_cif)
        if not mmseqs_db_path.exists():
            raise ValueError(
                f"MMseqs DB path not found: {mmseqs_db_path}. "
                "Pass --mmseqs_db_path explicitly or run fetch_bias_training_data first."
            )

        protein_tmp = args.output_protein_csv.parent / "_mmseqs_bias_tmp"
        query_fasta_path = protein_tmp / "query.fasta"
        query_fasta_path.parent.mkdir(parents=True, exist_ok=True)
        query_fasta_path.write_text(f">query\n{query_sequence}\n", encoding="utf-8")

        mmseqs_hits = _run_mmseqs(
            mmseqs_bin=args.mmseqs_bin,
            query_fasta=query_fasta_path,
            target_db=mmseqs_db_path,
            workers=args.workers,
            tmp_root=protein_tmp,
            max_seqs=args.mmseqs_max_seqs,
        )
        if len(mmseqs_hits) >= int(args.mmseqs_max_seqs):
            print(
                "[warn] MMseqs hit count reached --mmseqs_max_seqs=%d; results may be truncated. "
                "Increase --mmseqs_max_seqs for broader coverage."
                % int(args.mmseqs_max_seqs)
            )
        if mmseqs_hits.empty:
            protein_df = pd.DataFrame(columns=["pdb_id", "release_date", "sequence", "sequence_similarity"])
        else:
            mmseqs_hits["target"] = mmseqs_hits["target"].astype(str)
            mmseqs_hits["sequence"] = mmseqs_hits["tseq"].astype(str)
            mmseqs_hits["sequence_similarity"] = pd.to_numeric(mmseqs_hits["pident"], errors="coerce")
            mmseqs_hits["pdb_id"] = mmseqs_hits["target"].str.split("_").str[0].str.upper()
            mmseqs_hits = mmseqs_hits.dropna(subset=["sequence_similarity", "pdb_id", "sequence"])
            mmseqs_hits = mmseqs_hits[
                mmseqs_hits["sequence_similarity"] >= float(args.protein_similarity_threshold)
            ].copy()

            unique_pdb = sorted(mmseqs_hits["pdb_id"].unique())
            pdb_release_all: dict[str, date] = {}
            with ThreadPoolExecutor(max_workers=max(1, args.workers)) as ex:
                fut_map = {
                    ex.submit(_entry_release_date, pdb_id, args.timeout): pdb_id
                    for pdb_id in unique_pdb
                }
                for fut in as_completed(fut_map):
                    pdb_id = fut_map[fut]
                    try:
                        rel = fut.result()
                    except Exception:
                        continue
                    if rel is not None and rel < cutoff:
                        pdb_release_all[pdb_id] = rel

            mmseqs_hits = mmseqs_hits[mmseqs_hits["pdb_id"].isin(pdb_release_all.keys())].copy()
            mmseqs_hits["release_date"] = mmseqs_hits["pdb_id"].map(
                lambda x: pdb_release_all[x].isoformat()
            )
            protein_df = mmseqs_hits[["pdb_id", "release_date", "sequence", "sequence_similarity"]].copy()
            protein_df = protein_df.drop_duplicates(subset=["pdb_id", "sequence"])
            protein_df = protein_df.sort_values(
                by=["sequence_similarity", "pdb_id"],
                ascending=[False, True],
            )
            print(f"[info] protein MMseqs hits kept: {len(protein_df)}")

    # Build combined bias training data per pdb_id with protein/ligand combinations.
    # Required output columns: pdb_id, sequence_similarity, ecfp_similarity
    query_ligand_fp = _morgan_fp(query.smiles)
    ligand_by_pdb: dict[str, list[float]] = {}
    for _, r in ligand_df.iterrows():
        pdb_id = str(r["pdb_id"]).upper()
        sim = pd.to_numeric(r.get("ecfp_similarity"), errors="coerce")
        if pd.isna(sim):
            continue
        ligand_by_pdb.setdefault(pdb_id, []).append(float(sim))

    protein_by_pdb: dict[str, list[float]] = {}
    for _, r in protein_df.iterrows():
        pdb_id = str(r["pdb_id"]).upper()
        sim = pd.to_numeric(r.get("sequence_similarity"), errors="coerce")
        if pd.isna(sim):
            continue
        protein_by_pdb.setdefault(pdb_id, []).append(float(sim))

    def _rows_for_pdb(
        pdb_id: str,
        seq_sims: list[float],
        lig_sims: list[float],
    ) -> list[dict[str, float | str | None]]:
        # Keep one-side-known entries with NaN counterpart until backfilled.
        seq_vals = seq_sims if seq_sims else [None]
        lig_vals = lig_sims if lig_sims else [None]
        out_rows: list[dict[str, float | str | None]] = []
        for ss in seq_vals:
            for ls in lig_vals:
                out_rows.append(
                    {
                        "pdb_id": pdb_id,
                        "sequence_similarity": (float(ss) if ss is not None else None),
                        "ecfp_similarity": (float(ls) if ls is not None else None),
                    }
                )
        return out_rows

    # Step 1: write initial combined table with already-known similarities.
    initial_rows: list[dict[str, float | str | None]] = []
    all_pdb_ids = sorted(set(ligand_by_pdb) | set(protein_by_pdb))
    for pdb_id in all_pdb_ids:
        seq_sims = list(protein_by_pdb.get(pdb_id, []))
        lig_sims = list(ligand_by_pdb.get(pdb_id, []))
        initial_rows.extend(_rows_for_pdb(pdb_id, seq_sims, lig_sims))

    bias_df = pd.DataFrame(initial_rows)

    # Ensure output parent dirs exist.
    args.output_protein_csv.parent.mkdir(parents=True, exist_ok=True)
    args.output_ligand_csv.parent.mkdir(parents=True, exist_ok=True)
    output_bias_csv.parent.mkdir(parents=True, exist_ok=True)

    if args.output_protein_csv.exists() and not args.overwrite:
        print(f"[info] keeping existing protein CSV (overwrite=false): {args.output_protein_csv}")
    else:
        protein_df.to_csv(args.output_protein_csv, index=False)

    if args.output_ligand_csv.exists() and not args.overwrite:
        print(f"[info] keeping existing ligand CSV (overwrite=false): {args.output_ligand_csv}")
    else:
        ligand_df.to_csv(args.output_ligand_csv, index=False)

    if output_bias_csv.exists() and not args.overwrite:
        try:
            existing_bias = pd.read_csv(output_bias_csv)
            if {"pdb_id", "sequence_similarity", "ecfp_similarity"} <= set(existing_bias.columns):
                bias_df = existing_bias.copy()
                print(f"[info] reusing existing bias CSV for backfill: {output_bias_csv}")
            else:
                print("[warn] existing bias CSV missing expected columns; rebuilding initial bias table")
                if not bias_df.empty:
                    bias_df = bias_df.sort_values(
                        by=["sequence_similarity", "ecfp_similarity", "pdb_id"],
                        ascending=[False, False, True],
                        na_position="last",
                    )
                bias_df.to_csv(output_bias_csv, index=False)
        except Exception as e:
            print(f"[warn] failed reading existing bias CSV, rebuilding: {e}")
            if not bias_df.empty:
                bias_df = bias_df.sort_values(
                    by=["sequence_similarity", "ecfp_similarity", "pdb_id"],
                    ascending=[False, False, True],
                    na_position="last",
                )
            bias_df.to_csv(output_bias_csv, index=False)
    else:
        if not bias_df.empty:
            bias_df = bias_df.sort_values(
                by=["sequence_similarity", "ecfp_similarity", "pdb_id"],
                ascending=[False, False, True],
                na_position="last",
            )
        bias_df.to_csv(output_bias_csv, index=False)
        print(f"[info] initial bias_csv written (pre-backfill): {output_bias_csv} rows={len(bias_df)}")

    # Step 2: backfill missing sides in chunks of 100 PDBs and save after each chunk.
    if not bias_df.empty and {"pdb_id", "sequence_similarity", "ecfp_similarity"} <= set(bias_df.columns):
        missing_mask = bias_df["sequence_similarity"].isna() | bias_df["ecfp_similarity"].isna()
        missing_pdb_ids = sorted(bias_df.loc[missing_mask, "pdb_id"].astype(str).str.upper().unique())
    else:
        missing_pdb_ids = [
            pdb_id
            for pdb_id in all_pdb_ids
            if (not protein_by_pdb.get(pdb_id)) or (not ligand_by_pdb.get(pdb_id))
        ]
    chunk_size = 100
    for start in range(0, len(missing_pdb_ids), chunk_size):
        chunk = missing_pdb_ids[start:start + chunk_size]
        replacement_rows: list[dict[str, float | str | None]] = []
        for pdb_id in chunk:
            seq_sims = list(protein_by_pdb.get(pdb_id, []))
            lig_sims = list(ligand_by_pdb.get(pdb_id, []))

            if not lig_sims and query_ligand_fp is not None:
                for ccd_id in _entry_ligand_ids(pdb_id, args.timeout):
                    smiles = ccd_smiles_by_id.get(ccd_id)
                    if not smiles:
                        continue
                    fp = _morgan_fp(smiles)
                    if fp is None:
                        continue
                    lig_sims.append(float(DataStructs.TanimotoSimilarity(query_ligand_fp, fp)))

            if not seq_sims and query_sequence:
                try:
                    seqs = _entry_fasta_sequences(pdb_id, args.timeout)
                except Exception:
                    seqs = []
                for s in seqs:
                    seq_sims.append(_sequence_identity_percent(query_sequence, str(s)))

            replacement_rows.extend(_rows_for_pdb(pdb_id, seq_sims, lig_sims))

        # Replace chunk rows and save intermediate result.
        keep_df = bias_df[~bias_df["pdb_id"].isin(chunk)].copy() if not bias_df.empty else pd.DataFrame()
        repl_df = pd.DataFrame(replacement_rows)
        bias_df = pd.concat([keep_df, repl_df], ignore_index=True)
        if not bias_df.empty:
            bias_df = bias_df.sort_values(
                by=["sequence_similarity", "ecfp_similarity", "pdb_id"],
                ascending=[False, False, True],
                na_position="last",
            )
        bias_df.to_csv(output_bias_csv, index=False)
        print(
            f"[info] bias backfill progress: {min(start + chunk_size, len(missing_pdb_ids))}/"
            f"{len(missing_pdb_ids)} pdbs, rows={len(bias_df)}"
        )

    # Finalize: unresolved missing values mean no ligand/protein was found.
    # Set these to 0.0 as explicit "no similarity evidence" values.
    if not bias_df.empty:
        missing_before = int(
            bias_df["sequence_similarity"].isna().sum()
            + bias_df["ecfp_similarity"].isna().sum()
        )
        if missing_before > 0:
            bias_df["sequence_similarity"] = pd.to_numeric(
                bias_df["sequence_similarity"], errors="coerce"
            ).fillna(0.0)
            bias_df["ecfp_similarity"] = pd.to_numeric(
                bias_df["ecfp_similarity"], errors="coerce"
            ).fillna(0.0)
            bias_df = bias_df.sort_values(
                by=["sequence_similarity", "ecfp_similarity", "pdb_id"],
                ascending=[False, False, True],
            )
            bias_df.to_csv(output_bias_csv, index=False)
            print(
                f"[info] finalized bias_csv: filled unresolved missing similarities with 0.0 "
                f"(values_filled={missing_before})"
            )

    # Save ligand-hit provenance for plotting/debugging.
    hits_out = args.output_hits_csv or args.output_ligand_csv.with_name("ccd_ligand_hits.csv")
    hits_best_sorted = hits_best.sort_values("ecfp_similarity", ascending=False)
    hits_best_sorted.to_csv(hits_out, index=False)

    print(f"[done] protein_csv={args.output_protein_csv} rows={len(protein_df)}")
    print(f"[done] ligand_csv={args.output_ligand_csv} rows={len(ligand_df)}")
    print(f"[done] bias_csv={output_bias_csv} rows={len(bias_df)}")
    print(f"[done] hits_csv={hits_out} rows={len(hits_best_sorted)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
