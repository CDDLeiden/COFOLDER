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
import os
from pathlib import Path
import re
import shutil
import subprocess
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit import RDLogger
import yaml
import gemmi

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
CORE_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
CORE_NONPOLY_URL = "https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{pdb_id}/{entity_id}"
MORGAN_FP = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
MMSEQS_PIDENT_METHOD = "mmseqs_pident"
UNAVAILABLE_SIMILARITY_METHOD = "unavailable"
COMBINED_BIAS_COLUMNS = [
    "pdb_id",
    "sequence_similarity",
    "sequence_similarity_method",
    "ecfp_similarity",
]
# Silence verbose RDKit parser noise for invalid CCD descriptors we intentionally skip.
RDLogger.DisableLog("rdApp.error")


@dataclass(frozen=True)
class LigandQuery:
    chain_id: str
    smiles: str
    source: str
    source_ccd: str | None = None


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


def _extract_ligand_queries(
    system_path: Path,
    ccd_smiles_by_id: dict[str, str],
    selected_chains: set[str] | None = None,
) -> list[LigandQuery]:
    data = yaml.safe_load(system_path.read_text(encoding="utf-8"))
    seqs = data.get("sequences") if isinstance(data, dict) else None
    if not isinstance(seqs, list):
        return []

    out: list[LigandQuery] = []

    for item in seqs:
        if not isinstance(item, dict):
            continue
        lig = item.get("ligand")
        if not isinstance(lig, dict):
            continue
        chain_ids = [c.strip().upper() for c in _as_list(lig.get("id")) if c and c.strip()]
        if selected_chains:
            chain_ids = [c for c in chain_ids if c in selected_chains]
            if not chain_ids:
                continue
        smiles = _normalize_smiles(lig.get("smiles"))
        if smiles:
            for chain_id in chain_ids:
                out.append(
                    LigandQuery(
                        chain_id=chain_id,
                        smiles=smiles,
                        source="smiles",
                    )
                )
            continue

        # Also support ligand definitions that only provide CCD IDs.
        ccd_ids = [x.upper() for x in _as_list(lig.get("ccd"))]
        mapped_smiles = None
        mapped_ccd = None
        for ccd_id in ccd_ids:
            mapped = _normalize_smiles(ccd_smiles_by_id.get(ccd_id))
            if not mapped:
                continue
            mapped_smiles = mapped
            mapped_ccd = ccd_id
            break

        if not mapped_smiles:
            joined_ccd = ",".join(ccd_ids) if ccd_ids else "<none>"
            joined_chains = ",".join(chain_ids) if chain_ids else "<none>"
            print(f"[warn] no CCD->SMILES mapping found for chains={joined_chains} ccd={joined_ccd}")
            continue

        for chain_id in chain_ids:
            print(
                f"[info] ccd_to_smiles chain={chain_id} ccd={mapped_ccd} smiles={mapped_smiles}"
            )
            out.append(
                LigandQuery(
                    chain_id=chain_id,
                    smiles=mapped_smiles,
                    source="ccd",
                    source_ccd=mapped_ccd,
                )
            )

    return out


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
        if sim > threshold:
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


def _extract_protein_query_chain_id_from_system(system_path: Path) -> str | None:
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
        ids = _as_list(prot.get("id"))
        for cid in ids:
            token = str(cid).strip().upper()
            if token:
                return token
    return None


def _order_ligand_training_columns(df: pd.DataFrame) -> pd.DataFrame:
    preferred = [
        "query_chain_id",
        "pdb_id",
        "release_date",
        "ligand_id",
        "ecfp_similarity",
        "smiles",
    ]
    cols = [c for c in preferred if c in df.columns] + [c for c in df.columns if c not in preferred]
    return df.loc[:, cols].copy()


def _order_protein_training_columns(df: pd.DataFrame) -> pd.DataFrame:
    preferred = [
        "query_chain_id",
        "pdb_id",
        "release_date",
        "sequence_similarity",
        "sequence",
    ]
    cols = [c for c in preferred if c in df.columns] + [c for c in df.columns if c not in preferred]
    return df.loc[:, cols].copy()


def _normalize_mmseqs_hits(mmseqs_hits: pd.DataFrame) -> pd.DataFrame:
    """Normalize MMseqs output without applying the protein similarity threshold."""
    columns = ["pdb_id", "sequence_similarity", "sequence"]
    if mmseqs_hits.empty:
        return pd.DataFrame(columns=columns)

    normalized = mmseqs_hits.copy()
    normalized["target"] = normalized["target"].astype("string")
    normalized["sequence"] = normalized["tseq"].astype("string")
    normalized["sequence_similarity"] = pd.to_numeric(normalized["pident"], errors="coerce")
    normalized["pdb_id"] = normalized["target"].str.split("_").str[0].str.upper()
    normalized = normalized.dropna(subset=columns)
    normalized = normalized[
        normalized["pdb_id"].astype(str).str.strip().ne("")
        & normalized["sequence"].astype(str).str.strip().ne("")
    ]
    return normalized.loc[:, columns].drop_duplicates(subset=["pdb_id", "sequence"])


def _protein_similarity_lookup(protein_hits: pd.DataFrame) -> dict[str, list[float]]:
    lookup: dict[str, list[float]] = {}
    for _, row in protein_hits.iterrows():
        pdb_id = str(row.get("pdb_id", "")).strip().upper()
        similarity = pd.to_numeric(row.get("sequence_similarity"), errors="coerce")
        if not pdb_id or pd.isna(similarity):
            continue
        lookup.setdefault(pdb_id, []).append(float(similarity))
    return lookup


def _combined_rows_for_pdb(
    pdb_id: str,
    sequence_similarities: list[float],
    ligand_similarities: list[float],
) -> list[dict[str, float | str | None]]:
    """Create combined rows while keeping protein similarity provenance explicit."""
    sequence_values = sequence_similarities if sequence_similarities else [None]
    ligand_values = ligand_similarities if ligand_similarities else [None]
    rows: list[dict[str, float | str | None]] = []
    for sequence_similarity in sequence_values:
        for ligand_similarity in ligand_values:
            rows.append(
                {
                    "pdb_id": pdb_id,
                    "sequence_similarity": (
                        float(sequence_similarity) if sequence_similarity is not None else None
                    ),
                    "sequence_similarity_method": (
                        MMSEQS_PIDENT_METHOD
                        if sequence_similarity is not None
                        else UNAVAILABLE_SIMILARITY_METHOD
                    ),
                    "ecfp_similarity": (
                        float(ligand_similarity) if ligand_similarity is not None else None
                    ),
                }
            )
    return rows


def _read_reusable_bias_csv(path: Path) -> pd.DataFrame | None:
    """Return a current-schema bias CSV, or None when it must be rebuilt."""
    try:
        existing = pd.read_csv(path)
    except Exception as exc:
        print(f"[warn] failed reading existing bias CSV, rebuilding: {exc}")
        return None

    if not set(COMBINED_BIAS_COLUMNS) <= set(existing.columns):
        print(
            "[warn] existing bias CSV missing expected columns or similarity provenance; "
            "rebuilding initial bias table"
        )
        return None
    methods = existing["sequence_similarity_method"].astype(str)
    valid_methods = {MMSEQS_PIDENT_METHOD, UNAVAILABLE_SIMILARITY_METHOD}
    if not methods.isin(valid_methods).all():
        print("[warn] existing bias CSV has invalid similarity provenance; rebuilding")
        return None
    similarities = pd.to_numeric(existing["sequence_similarity"], errors="coerce")
    inconsistent = (
        (similarities.notna() & methods.ne(MMSEQS_PIDENT_METHOD))
        | (similarities.isna() & methods.ne(UNAVAILABLE_SIMILARITY_METHOD))
    )
    if inconsistent.any():
        print("[warn] existing bias CSV has inconsistent similarity provenance; rebuilding")
        return None
    return existing


def _guess_mmseqs_db_from_components(components_cif: Path) -> Path:
    # expected fetch layout:
    # <output_root>/ccd/components.cif
    # <output_root>/mmseqs/pdb/db
    output_root = components_cif.parent.parent
    return output_root / "mmseqs" / "pdb" / "db"


def _resolve_mmseqs_bin(mmseqs_bin: str | None = None) -> str | None:
    candidates: list[str] = []
    if mmseqs_bin:
        candidates.append(str(mmseqs_bin))
    explicit_env = os.environ.get("COFOLDER_MMSEQS_BIN")
    if explicit_env:
        candidates.insert(0, explicit_env)
    candidates.extend(
        [
            str(Path.home() / ".cofolder/vendor/mmseqs/bin/mmseqs"),
            str(Path(__file__).resolve().parents[1] / "vendor/mmseqs/bin/mmseqs"),
            str(Path(__file__).resolve().parents[1] / "vendor/mmseqs/mmseqs"),
        ]
    )
    from_path = shutil.which("mmseqs")
    if from_path:
        candidates.append(from_path)

    for cand in candidates:
        if not cand:
            continue
        p = Path(cand).expanduser()
        if p.exists() and p.is_file():
            return str(p)
        found = shutil.which(cand)
        if found:
            return found
    return None


def _run_mmseqs(
    mmseqs_bin: str,
    query_fasta: Path,
    target_db: Path,
    workers: int,
    tmp_root: Path,
    max_seqs: int,
) -> pd.DataFrame:
    resolved_mmseqs = _resolve_mmseqs_bin(mmseqs_bin)
    if resolved_mmseqs is None:
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
        [resolved_mmseqs, "createdb", str(query_fasta), str(qdb)],
        [
            resolved_mmseqs,
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
            resolved_mmseqs,
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
    hits = _find_ccd_hits(
        query_smiles=query.smiles,
        ccd_df=ccd_df,
        threshold=threshold,
        top_k=top_k,
    )
    for ccd_id, hit_smiles, sim in hits:
        rows.append(
            {
                "query_chain_id": query.chain_id,
                "ligand_pdbid": ccd_id,
                "hit_smiles": hit_smiles,
                "ecfp_similarity": sim,
            }
        )
    return pd.DataFrame(rows)


def _write_chain_ligand_csvs(output_ligand_csv: Path, ligand_df: pd.DataFrame) -> None:
    if "query_chain_id" not in ligand_df.columns:
        return

    chain_ids = sorted(
        ligand_df["query_chain_id"].dropna().astype(str).str.strip().unique().tolist()
    )
    stem = output_ligand_csv.stem
    suffix = output_ligand_csv.suffix or ".csv"
    for chain_id in chain_ids:
        if not chain_id:
            continue
        chain_path = output_ligand_csv.with_name(f"{stem}_{chain_id}{suffix}")
        sub = ligand_df[
            ligand_df["query_chain_id"].astype(str).str.strip() == chain_id
        ].copy()
        sub.to_csv(chain_path, index=False)
        print(f"[done] chain_ligand_csv[{chain_id}]={chain_path} rows={len(sub)}")


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
    parser.add_argument(
        "--skip_protein_mmseqs",
        action="store_true",
        help=(
            "Skip MMseqs protein expansion and emit an empty protein_training_data.csv "
            "with required columns."
        ),
    )
    parser.add_argument(
        "--skip_bias_csv",
        action="store_true",
        help="Do not write bias_training_data.csv.",
    )
    parser.add_argument(
        "--skip_ligand_ecfp",
        action="store_true",
        help="Skip ligand ECFP/CCD expansion and emit an empty ligand_training_data.csv.",
    )
    parser.add_argument(
        "--ligand_chains",
        nargs="+",
        default=None,
        help="Optional ligand chain IDs to process (e.g. --ligand_chains E F).",
    )
    args = parser.parse_args()

    cutoff = _parse_iso_date(args.release_cutoff)
    if cutoff is None:
        raise ValueError("--release_cutoff must be YYYY-MM-DD")
    if args.query_sequence and args.query_sequence_fasta:
        raise ValueError("Use either --query_sequence or --query_sequence_fasta, not both.")

    ccd_smiles_by_id: dict[str, str] = {}
    ccd_to_entries: dict[str, list[str]] = {}
    queries: list[LigandQuery] = []
    hits_best_sorted = pd.DataFrame(
        columns=["query_chain_id", "ligand_pdbid", "hit_smiles", "ecfp_similarity"]
    )
    output_bias_csv = args.output_bias_csv or args.output_ligand_csv.with_name("bias_training_data.csv")

    if args.skip_ligand_ecfp:
        print("[info] skipping ligand ECFP/CCD expansion (--skip_ligand_ecfp)")
        ligand_df = pd.DataFrame(
            columns=["query_chain_id", "pdb_id", "release_date", "ligand_id", "ecfp_similarity", "smiles"]
        )
    else:
        ccd_df = _load_ccd_smiles(args.components_cif)
        if ccd_df.empty:
            raise ValueError("No CCD SMILES could be read from components.cif")
        ccd_smiles_by_id = {
            str(row["ligand_id"]).upper(): _normalize_smiles(str(row["smiles"]))
            for _, row in ccd_df.iterrows()
        }

        selected_ligand_chains: set[str] | None = None
        if args.ligand_chains:
            selected_ligand_chains = set()
            for value in args.ligand_chains:
                for token in str(value).split(","):
                    token = token.strip().upper()
                    if token:
                        selected_ligand_chains.add(token)
            print(f"[info] selected_ligand_chains={','.join(sorted(selected_ligand_chains))}")

        queries = _extract_ligand_queries(
            args.system_path,
            ccd_smiles_by_id=ccd_smiles_by_id,
            selected_chains=selected_ligand_chains,
        )
        if not queries:
            raise ValueError(
                "No ligand queries found. Provide ligand.smiles and/or ligand.ccd values in system.yaml."
            )

        print("[info] query_mode=multi_chain")
        print(f"[info] query_chain_ids={','.join(q.chain_id for q in queries)}")
        print(f"[info] ccd_with_smiles={len(ccd_df)}")

        hit_frames: list[pd.DataFrame] = []
        for query in queries:
            print(f"[info] collecting_ccd_hits chain={query.chain_id} source={query.source}")
            hit_frames.append(
                _collect_ccd_hits(
                    query=query,
                    ccd_df=ccd_df,
                    threshold=args.ligand_similarity_threshold,
                    top_k=args.top_k_per_query,
                )
            )
        hits_df = pd.concat(hit_frames, ignore_index=True) if hit_frames else pd.DataFrame()
        if hits_df.empty:
            raise ValueError(
                "No CCD hits above threshold. Lower --ligand_similarity_threshold or inspect query ligands."
            )

        # Keep only best score per CCD ID for each query chain independently.
        hits_df = hits_df.sort_values("ecfp_similarity", ascending=False)
        hits_best = hits_df.drop_duplicates(subset=["query_chain_id", "ligand_pdbid"], keep="first")
        hit_ccd_ids = sorted(hits_best["ligand_pdbid"].astype(str).unique())
        hits_best_sorted = hits_best.sort_values("ecfp_similarity", ascending=False).copy()

        print(f"[info] ccd_hits={len(hit_ccd_ids)}")

        # For each hit CCD, fetch associated PDB entries.
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
        ligand_rows: list[dict[str, str]] = []
        for _, hit in hits_best.iterrows():
            chain_id = str(hit["query_chain_id"]).strip()
            ccd_id = str(hit["ligand_pdbid"]).strip().upper()
            smiles = str(hit["hit_smiles"])
            sim = float(hit["ecfp_similarity"])
            pdb_ids = ccd_to_entries.get(ccd_id, [])
            for pdb_id in pdb_ids:
                rel = pdb_release.get(pdb_id)
                if rel is None:
                    continue
                ligand_rows.append(
                    {
                        "query_chain_id": chain_id,
                        "pdb_id": pdb_id,
                        "release_date": rel.isoformat(),
                        "ligand_id": ccd_id,
                        "smiles": smiles,
                        "ecfp_similarity": sim,
                    }
                )

        ligand_columns = [
            "query_chain_id",
            "pdb_id",
            "release_date",
            "ligand_id",
            "ecfp_similarity",
            "smiles",
        ]
        ligand_df = pd.DataFrame(ligand_rows, columns=ligand_columns)
        if not ligand_df.empty:
            ligand_df = ligand_df.drop_duplicates(
                subset=["query_chain_id", "pdb_id", "ligand_id", "smiles"]
            )
        if not ligand_df.empty:
            ligand_df = ligand_df.sort_values(
                by=["query_chain_id", "ecfp_similarity", "pdb_id", "ligand_id"],
                ascending=[True, False, True, True],
            )
        if args.output_ligand_csv.exists() and not args.overwrite:
            try:
                ligand_df = pd.read_csv(args.output_ligand_csv)
                print(f"[info] reusing existing ligand CSV: {args.output_ligand_csv}")
            except Exception as e:
                print(f"[warn] failed to read existing ligand CSV, using newly computed data: {e}")
    ligand_df = _order_ligand_training_columns(ligand_df)

    # Persist ligand outputs immediately after compute/reuse.
    args.output_ligand_csv.parent.mkdir(parents=True, exist_ok=True)
    if args.output_ligand_csv.exists() and not args.overwrite:
        print(f"[info] keeping existing ligand CSV (overwrite=false): {args.output_ligand_csv}")
    else:
        ligand_df.to_csv(args.output_ligand_csv, index=False)
    _write_chain_ligand_csvs(args.output_ligand_csv, ligand_df)

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

    protein_query_chain_id = _extract_protein_query_chain_id_from_system(args.system_path)
    raw_mmseqs_hits = pd.DataFrame(columns=["pdb_id", "sequence_similarity", "sequence"])
    needs_complete_mmseqs_lookup = not args.skip_bias_csv
    should_run_mmseqs = (
        not args.skip_protein_mmseqs
        and (protein_df.empty or needs_complete_mmseqs_lookup)
    )

    if should_run_mmseqs:
        mmseqs_db_path = args.mmseqs_db_path or _guess_mmseqs_db_from_components(args.components_cif)
        if not mmseqs_db_path.exists():
            raise ValueError(
                f"MMseqs DB path not found: {mmseqs_db_path}. "
                "Pass --mmseqs_db_path explicitly, run fetch_bias_training_data first, "
                "or use --skip_protein_mmseqs."
            )

        protein_tmp = args.output_protein_csv.parent / "_mmseqs_bias_tmp"
        query_fasta_path = protein_tmp / "query.fasta"
        query_fasta_path.parent.mkdir(parents=True, exist_ok=True)
        query_fasta_path.write_text(f">query\n{query_sequence}\n", encoding="utf-8")

        mmseqs_output = _run_mmseqs(
            mmseqs_bin=args.mmseqs_bin,
            query_fasta=query_fasta_path,
            target_db=mmseqs_db_path,
            workers=args.workers,
            tmp_root=protein_tmp,
            max_seqs=args.mmseqs_max_seqs,
        )
        if len(mmseqs_output) >= int(args.mmseqs_max_seqs):
            print(
                "[warn] MMseqs hit count reached --mmseqs_max_seqs=%d; results may be truncated. "
                "Increase --mmseqs_max_seqs for broader coverage."
                % int(args.mmseqs_max_seqs)
            )
        raw_mmseqs_hits = _normalize_mmseqs_hits(mmseqs_output)

    if protein_df.empty:
        if args.skip_protein_mmseqs:
            print(
                "[info] skipping MMseqs protein expansion (--skip_protein_mmseqs); "
                "protein CSV will be empty"
            )
        if raw_mmseqs_hits.empty:
            protein_df = pd.DataFrame(
                columns=["query_chain_id", "pdb_id", "release_date", "sequence_similarity", "sequence"]
            )
        else:
            thresholded_hits = raw_mmseqs_hits[
                raw_mmseqs_hits["sequence_similarity"] > float(args.protein_similarity_threshold)
            ].copy()

            unique_pdb = sorted(thresholded_hits["pdb_id"].unique())
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

            thresholded_hits = thresholded_hits[
                thresholded_hits["pdb_id"].isin(pdb_release_all)
            ].copy()
            thresholded_hits["release_date"] = thresholded_hits["pdb_id"].map(
                lambda pdb_id: pdb_release_all[pdb_id].isoformat()
            )
            protein_df = thresholded_hits[
                ["pdb_id", "release_date", "sequence_similarity", "sequence"]
            ].copy()
            protein_df = protein_df.sort_values(
                by=["sequence_similarity", "pdb_id"],
                ascending=[False, True],
            )
            print(f"[info] protein MMseqs hits kept: {len(protein_df)}")
    if "query_chain_id" not in protein_df.columns:
        protein_df["query_chain_id"] = protein_query_chain_id or ""
    else:
        protein_df["query_chain_id"] = protein_df["query_chain_id"].fillna(protein_query_chain_id or "")
    protein_df = _order_protein_training_columns(protein_df)

    # The standalone protein table is thresholded, but the combined table must use
    # every MMseqs pident available for ligand-side PDB lookup.
    protein_lookup_source = raw_mmseqs_hits if not raw_mmseqs_hits.empty else protein_df

    # Build combined bias training data per pdb_id with protein/ligand combinations.
    query_fps_by_chain = {
        q.chain_id: _morgan_fp(q.smiles)
        for q in queries
    } if queries else {}
    ligand_by_pdb: dict[str, list[float]] = {}
    for _, r in ligand_df.iterrows():
        pdb_id = str(r["pdb_id"]).upper()
        sim = pd.to_numeric(r.get("ecfp_similarity"), errors="coerce")
        if pd.isna(sim):
            continue
        ligand_by_pdb.setdefault(pdb_id, []).append(float(sim))

    protein_by_pdb = _protein_similarity_lookup(protein_lookup_source)
    standalone_protein_pdb_ids = set(_protein_similarity_lookup(protein_df))

    # Step 1: write initial combined table with already-known similarities.
    initial_rows: list[dict[str, float | str | None]] = []
    all_pdb_ids = sorted(set(ligand_by_pdb) | standalone_protein_pdb_ids)
    for pdb_id in all_pdb_ids:
        seq_sims = list(protein_by_pdb.get(pdb_id, []))
        lig_sims = list(ligand_by_pdb.get(pdb_id, []))
        initial_rows.extend(_combined_rows_for_pdb(pdb_id, seq_sims, lig_sims))

    bias_df = pd.DataFrame(initial_rows, columns=COMBINED_BIAS_COLUMNS)

    # Ensure output parent dirs exist.
    args.output_protein_csv.parent.mkdir(parents=True, exist_ok=True)
    if not args.skip_bias_csv:
        output_bias_csv.parent.mkdir(parents=True, exist_ok=True)

    if args.output_protein_csv.exists() and not args.overwrite:
        print(f"[info] keeping existing protein CSV (overwrite=false): {args.output_protein_csv}")
    else:
        protein_df.to_csv(args.output_protein_csv, index=False)

    # Always mirror training references into results/bias_train.
    common_root = Path(
        os.path.commonpath(
            [str(args.output_protein_csv.parent.resolve()), str(args.output_ligand_csv.parent.resolve())]
        )
    )
    bias_train_dir = common_root / "results" / "bias_train"
    bias_train_dir.mkdir(parents=True, exist_ok=True)
    protein_mirror = bias_train_dir / "protein_training_data.csv"
    ligand_mirror = bias_train_dir / "ligand_training_data.csv"

    if protein_mirror.exists() and not args.overwrite:
        print(f"[info] keeping existing mirrored protein CSV (overwrite=false): {protein_mirror}")
    else:
        protein_df.to_csv(protein_mirror, index=False)

    if ligand_mirror.exists() and not args.overwrite:
        print(f"[info] keeping existing mirrored ligand CSV (overwrite=false): {ligand_mirror}")
    else:
        ligand_df.to_csv(ligand_mirror, index=False)

    if args.skip_bias_csv:
        print("[info] skipping bias CSV generation (--skip_bias_csv)")
    else:
        if output_bias_csv.exists() and not args.overwrite:
            existing_bias = _read_reusable_bias_csv(output_bias_csv)
            if existing_bias is not None:
                bias_df = existing_bias.copy()
                print(f"[info] reusing existing bias CSV for backfill: {output_bias_csv}")
            else:
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

        # Step 2: backfill only missing ligand similarities. Protein similarity is
        # exclusively sourced from the complete MMseqs lookup above.
        if not bias_df.empty and set(COMBINED_BIAS_COLUMNS) <= set(bias_df.columns):
            missing_mask = bias_df["ecfp_similarity"].isna()
            missing_pdb_ids = sorted(bias_df.loc[missing_mask, "pdb_id"].astype(str).str.upper().unique())
        else:
            missing_pdb_ids = [
                pdb_id
                for pdb_id in all_pdb_ids
                if not ligand_by_pdb.get(pdb_id)
            ]
        chunk_size = 100
        for start in range(0, len(missing_pdb_ids), chunk_size):
            chunk = missing_pdb_ids[start:start + chunk_size]
            replacement_rows: list[dict[str, float | str | None]] = []
            for pdb_id in chunk:
                seq_sims = list(protein_by_pdb.get(pdb_id, []))
                lig_sims = list(ligand_by_pdb.get(pdb_id, []))

                if not lig_sims:
                    # Backfill against all ligand query chains independently.
                    for ccd_id in _entry_ligand_ids(pdb_id, args.timeout):
                        smiles = ccd_smiles_by_id.get(ccd_id)
                        if not smiles:
                            continue
                        fp = _morgan_fp(smiles)
                        if fp is None:
                            continue
                        for qfp in query_fps_by_chain.values():
                            if qfp is None:
                                continue
                            lig_sims.append(float(DataStructs.TanimotoSimilarity(qfp, fp)))

                replacement_rows.extend(
                    _combined_rows_for_pdb(pdb_id, seq_sims, lig_sims)
                )

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

        # Finalize unresolved ligand values as 0.0. Protein values remain null
        # with explicit provenance when MMseqs returned no result.
        if not bias_df.empty:
            missing_ligand_before = int(bias_df["ecfp_similarity"].isna().sum())
            if missing_ligand_before > 0:
                bias_df["ecfp_similarity"] = pd.to_numeric(
                    bias_df["ecfp_similarity"], errors="coerce"
                ).fillna(0.0)
                bias_df = bias_df.sort_values(
                    by=["sequence_similarity", "ecfp_similarity", "pdb_id"],
                    ascending=[False, False, True],
                )
                bias_df.to_csv(output_bias_csv, index=False)
                print(
                    f"[info] finalized bias_csv: filled unresolved ligand similarities with 0.0 "
                    f"(values_filled={missing_ligand_before})"
                )

    # Save ligand-hit provenance for plotting/debugging.
    hits_out = args.output_hits_csv or args.output_ligand_csv.with_name("ccd_ligand_hits.csv")
    hits_best_sorted.to_csv(hits_out, index=False)

    print(f"[done] protein_csv={args.output_protein_csv} rows={len(protein_df)}")
    print(f"[done] ligand_csv={args.output_ligand_csv} rows={len(ligand_df)}")
    print(f"[done] mirrored_protein_csv={protein_mirror}")
    print(f"[done] mirrored_ligand_csv={ligand_mirror}")
    if args.skip_bias_csv:
        print("[done] bias_csv=skipped")
    else:
        print(f"[done] bias_csv={output_bias_csv} rows={len(bias_df)}")
    print(f"[done] hits_csv={hits_out} rows={len(hits_best_sorted)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
