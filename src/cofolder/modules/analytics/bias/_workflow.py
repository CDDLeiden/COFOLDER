"""Stable orchestration entry point for bias metrics."""

from __future__ import annotations

from cofolder.modules.utils.timing import DebugTimingCollector
from pathlib import Path
import logging
from contextlib import nullcontext
from cofolder.modules.analytics.bias_database import parse_bias_release_policy
import pandas as pd

from ._common import (
    BIAS_PLOT_FILE_STEM,
    BIAS_PLOT_LIGAND_THRESHOLD,
    BIAS_PLOT_SEQUENCE_THRESHOLD,
    _norm_id,
)

from ._similarity import (
    _best_ligand_hit,
    _best_pairwise_protein_hit,
    _best_protein_hit,
    _reference_rows_for_query,
)

from ._datasets import (
    _build_ligand_reference_index,
    _combined_mixed_pair_bias_training_dataframe,
)

from ._artifacts import (
    _clear_pair_specific_artifacts,
    _materialize_bias_training_views,
    _materialize_reference_landscape_summary,
    _pair_specific_bias_training_datasets,
    _write_bias_plot_artifacts,
    _write_same_type_plot_artifacts,
)

from ._references import (
    _extract_system_queries,
    _load_ligand_training,
    _load_protein_training,
)


def apply_bias_metrics(
    system_df: pd.DataFrame,
    chain_df: pd.DataFrame,
    sys_obj,
    protein_training_data_path: Path | str | None,
    ligand_training_data_path: Path | str | None,
    custom_protein_reference_path: Path | str | None = None,
    custom_ligand_reference_path: Path | str | None = None,
    release_cutoff: str = "2023-06-01",
    bias_chains: set[str] | list[str] | None = None,
    protein_top_n: int = 100,
    protein_similarity_threshold: float = BIAS_PLOT_SEQUENCE_THRESHOLD,
    ligand_similarity_threshold: float = BIAS_PLOT_LIGAND_THRESHOLD,
    boltz_cache_path: Path | str | None = None,
    output_dir: Path | str | None = None,
    logger: logging.Logger | None = None,
    timings: DebugTimingCollector | None = None,
    components_cif_path: Path | str | None = None,
    strict_query_resolution: bool = False,
):
    """Compute lightweight chain-level bias metrics only."""
    if logger is None:
        logger = logging.getLogger("cofolder.modules.analytics.bias")

    release_policy = parse_bias_release_policy(release_cutoff)
    cutoff = release_policy.cutoff
    protein_training_path = (
        Path(protein_training_data_path)
        if protein_training_data_path is not None
        else None
    )
    ligand_training_path = (
        Path(ligand_training_data_path)
        if ligand_training_data_path is not None
        else None
    )
    custom_protein_path = (
        Path(custom_protein_reference_path)
        if custom_protein_reference_path is not None
        else None
    )
    custom_ligand_path = (
        Path(custom_ligand_reference_path)
        if custom_ligand_reference_path is not None
        else None
    )
    proteins_df = _load_protein_training(
        protein_training_path,
        cutoff,
        top_n=protein_top_n,
        custom_reference_path=custom_protein_path,
    )
    ligands_df = _load_ligand_training(
        ligand_training_path,
        cutoff,
        custom_reference_path=custom_ligand_path,
    )
    cache_path = (
        Path(boltz_cache_path).expanduser()
        if boltz_cache_path
        else Path("~/.boltz").expanduser()
    )
    components_path = (
        Path(components_cif_path).expanduser()
        if components_cif_path is not None
        else None
    )
    logger.info(
        "Bias release policy applied (%s): protein_rows=%d ligand_rows=%d",
        release_policy.label,
        len(proteins_df),
        len(ligands_df),
    )
    selected_chains = {
        str(x).strip().upper() for x in (bias_chains or []) if str(x).strip()
    }
    if selected_chains:
        logger.info("Bias restricted to selected chains: %s", sorted(selected_chains))

    protein_queries, ligand_queries, unresolved_ccd_by_chain = _extract_system_queries(
        sys_obj,
        ligands_df,
        boltz_cache_path=cache_path,
        components_cif_path=components_path,
    )
    if selected_chains:
        protein_queries = {
            cid: seq
            for cid, seq in protein_queries.items()
            if str(cid).strip().upper() in selected_chains
        }
        ligand_queries = {
            cid: smi
            for cid, smi in ligand_queries.items()
            if str(cid).strip().upper() in selected_chains
        }
        unresolved_ccd_by_chain = {
            cid: ccd_ids
            for cid, ccd_ids in unresolved_ccd_by_chain.items()
            if str(cid).strip().upper() in selected_chains
        }
    if strict_query_resolution and unresolved_ccd_by_chain:
        details = ", ".join(
            f"{chain_id} ({'/'.join(ccd_ids)})"
            for chain_id, ccd_ids in sorted(unresolved_ccd_by_chain.items())
        )
        raise ValueError(
            "Could not resolve SMILES for CCD-backed ligand chain(s): "
            f"{details}. Provide --bias_training_components_cif, use ligand training data "
            "with matching ligand_id values, or prepare the local CCD cache."
        )

    protein_best: dict[str, float | None] = {}
    protein_pairwise_best: dict[str, float | None] = {}
    protein_timer = (
        timings.measure("scores.bias_metrics.protein_similarity", logger=logger)
        if timings is not None
        else nullcontext()
    )
    with protein_timer:
        for chain_id, query_seq in protein_queries.items():
            normalized_chain_id = str(chain_id).strip().upper()
            query_references = _reference_rows_for_query(proteins_df, chain_id)
            protein_best[normalized_chain_id] = (
                _best_protein_hit(query_seq, query_references) if query_seq else None
            )
            protein_pairwise_best[normalized_chain_id] = (
                _best_pairwise_protein_hit(query_seq, query_references)
                if query_seq
                else None
            )

    ligand_best: dict[str, float | None] = {}
    ligand_timer = (
        timings.measure("scores.bias_metrics.ligand_similarity", logger=logger)
        if timings is not None
        else nullcontext()
    )
    with ligand_timer:
        for chain_id, query_smiles in ligand_queries.items():
            query_references = _reference_rows_for_query(ligands_df, chain_id)
            ligand_best[str(chain_id).strip().upper()] = (
                _best_ligand_hit(query_smiles, query_references)
                if query_smiles
                else None
            )

    if "bias_prot_sim_train" not in chain_df.columns:
        chain_df["bias_prot_sim_train"] = pd.NA
    if "bias_prot_sim_train_pairwise" not in chain_df.columns:
        chain_df["bias_prot_sim_train_pairwise"] = pd.NA
    if "bias_lig_sim_train" not in chain_df.columns:
        chain_df["bias_lig_sim_train"] = pd.NA

    for idx, row in chain_df.iterrows():
        entity = str(row.get("ENTITY_TYPE"))
        chain_id = str(row.get("CHAIN_ID")).strip().upper()
        if selected_chains and chain_id not in selected_chains:
            continue
        if entity == "protein":
            chain_df.at[idx, "bias_prot_sim_train"] = protein_best.get(chain_id)
            chain_df.at[idx, "bias_prot_sim_train_pairwise"] = (
                protein_pairwise_best.get(chain_id)
            )
        elif entity == "ligand":
            sim = ligand_best.get(chain_id)
            mol_id = row.get("ligand_molecule_id")
            if pd.notna(mol_id) and str(mol_id).strip():
                ref = ligands_df[
                    ligands_df.get("ligand_id").astype(str) == _norm_id(mol_id)
                ]
                if not ref.empty:
                    query_smiles = ligand_queries.get(chain_id)
                    if (
                        query_smiles
                        and ref["smiles"].astype(str).eq(str(query_smiles)).any()
                    ):
                        sim = 1.0
                    elif "ecfp_similarity" in ref.columns:
                        ecfp = pd.to_numeric(
                            ref["ecfp_similarity"], errors="coerce"
                        ).dropna()
                        if len(ecfp):
                            sim = float(ecfp.max())
                    if sim is None:
                        query_smiles = str(ref.iloc[0]["smiles"])
                        sim = _best_ligand_hit(query_smiles, ligands_df)
            chain_df.at[idx, "bias_lig_sim_train"] = sim

    prot_vals = pd.to_numeric(
        chain_df.get("bias_prot_sim_train"), errors="coerce"
    ).dropna()
    prot_pairwise_vals = pd.to_numeric(
        chain_df.get("bias_prot_sim_train_pairwise"), errors="coerce"
    ).dropna()
    lig_vals = pd.to_numeric(
        chain_df.get("bias_lig_sim_train"), errors="coerce"
    ).dropna()

    if len(prot_vals):
        system_df["bias_prot_sim_train_max"] = float(prot_vals.max())
    if len(prot_pairwise_vals):
        system_df["bias_prot_sim_train_pairwise_max"] = float(prot_pairwise_vals.max())
    if len(lig_vals):
        system_df["bias_lig_sim_train_max"] = float(lig_vals.max())

    if output_dir is not None:
        ligand_reference_index = _build_ligand_reference_index(ligands_df)
        output_chain_df = chain_df
        if selected_chains:
            output_chain_df = chain_df[
                chain_df["CHAIN_ID"]
                .astype(str)
                .str.strip()
                .str.upper()
                .isin(selected_chains)
            ].copy()
        output_root = Path(output_dir)
        protein_view, ligand_views, protein_lookup_view, ligand_lookup_views = (
            _materialize_bias_training_views(
                chain_df=output_chain_df,
                ligands_df=ligands_df,
                proteins_df=_load_protein_training(
                    protein_training_path,
                    cutoff,
                    top_n=0,
                    custom_reference_path=custom_protein_path,
                ),
                protein_queries=protein_queries,
                ligand_queries=ligand_queries,
                boltz_cache_path=cache_path,
                output_dir=output_root,
                components_cif_path=components_path,
            )
        )
        _clear_pair_specific_artifacts(output_root)
        pair_artifacts = _pair_specific_bias_training_datasets(
            chain_df=output_chain_df,
            protein_queries=protein_queries,
            ligand_queries=ligand_queries,
            ligand_reference_df=ligands_df,
            ligand_reference_index=ligand_reference_index,
            protein_view=protein_view,
            ligand_views=ligand_views,
            protein_lookup_view=protein_lookup_view,
            ligand_lookup_views=ligand_lookup_views,
            boltz_cache_path=cache_path,
            components_cif_path=components_path,
            output_dir=output_root,
        )
        bias_training_df = _combined_mixed_pair_bias_training_dataframe(pair_artifacts)
        bias_training_df.to_csv(output_root / "bias_training_data.csv", index=False)
        _write_bias_plot_artifacts(
            bias_training_df=bias_training_df,
            output_dir=output_root,
            file_stem=BIAS_PLOT_FILE_STEM,
            sequence_threshold=protein_similarity_threshold,
            ligand_threshold=ligand_similarity_threshold,
        )
        for pair_artifact in pair_artifacts:
            pair_df = pair_artifact["dataframe"]
            pair_df.to_csv(
                output_root / str(pair_artifact["data_file"]),
                index=False,
            )
            if pair_artifact["pair_type"] == "mixed":
                _write_bias_plot_artifacts(
                    bias_training_df=pair_df,
                    output_dir=output_root,
                    file_stem=str(pair_artifact["plot_file_stem"]),
                    sequence_threshold=protein_similarity_threshold,
                    ligand_threshold=ligand_similarity_threshold,
                )
            else:
                _write_same_type_plot_artifacts(
                    pair_df=pair_df,
                    output_dir=output_root,
                    file_stem=str(pair_artifact["plot_file_stem"]),
                    component_type=(
                        "protein"
                        if pair_artifact["pair_type"] == "protein_pair"
                        else "ligand"
                    ),
                    threshold=(
                        protein_similarity_threshold
                        if pair_artifact["pair_type"] == "protein_pair"
                        else ligand_similarity_threshold
                    ),
                )
        _materialize_reference_landscape_summary(
            chain_df=output_chain_df,
            proteins_df=proteins_df,
            ligands_df=ligands_df,
            protein_queries=protein_queries,
            ligand_queries=ligand_queries,
            output_dir=output_root,
        )
        logger.info("Bias training views written to %s", output_dir)

    logger.info(
        "Lightweight bias metrics computed: protein_max=%s ligand_max=%s (protein_top_n=%d)",
        (float(prot_vals.max()) if len(prot_vals) else None),
        (float(lig_vals.max()) if len(lig_vals) else None),
        protein_top_n,
    )

    return system_df, chain_df


__all__ = [
    "apply_bias_metrics",
]
