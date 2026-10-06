"""Public, non-mutating workflow preflight reports."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from cofolder.modules.input import load_yaml_document
from cofolder.modules.input.system import System, iter_system_chains
from cofolder.modules.input.validation import WorkflowInputRequirements
from cofolder.modules.runners import get_runner


@dataclass(frozen=True, slots=True)
class PreflightReport:
    workflow: str
    ready: bool
    runner: str | None = None
    runner_available: bool | None = None
    capabilities: tuple[str, ...] = ()
    selected_chains: tuple[str, ...] = ()
    database_sources: tuple[str, ...] = ()
    release_cutoff: str | None = None
    release_policy: str | None = None
    query_hashes: tuple[str, ...] = ()
    bias_query_cache: str | None = None
    bias_cache_ready: bool | None = None
    expected_bias_backfill: bool | None = None
    protein_similarity_threshold: float | None = None
    ligand_similarity_threshold: float | None = None
    repeats: int = 1
    diffusion_samples: int = 1
    planned_executions: int = 0
    external_msa_required: bool = False
    output_dir: str = ""
    messages: tuple[str, ...] = field(default_factory=tuple)

    def format_text(self) -> str:
        lines = [
            f"preflight={'ready' if self.ready else 'not-ready'}",
            f"workflow={self.workflow}",
        ]
        for name, value in (
            ("runner", self.runner),
            ("runner_available", self.runner_available),
            ("capabilities", ",".join(self.capabilities)),
            ("selected_chains", ",".join(self.selected_chains)),
            ("database_sources", ",".join(self.database_sources)),
            ("release_cutoff", self.release_cutoff),
            ("release_policy", self.release_policy),
            ("query_hashes", ",".join(self.query_hashes)),
            ("bias_query_cache", self.bias_query_cache),
            ("bias_cache_ready", self.bias_cache_ready),
            ("expected_bias_backfill", self.expected_bias_backfill),
            ("protein_similarity_threshold", self.protein_similarity_threshold),
            ("ligand_similarity_threshold", self.ligand_similarity_threshold),
            ("repeats", self.repeats),
            ("diffusion_samples", self.diffusion_samples),
            ("planned_executions", self.planned_executions),
            ("external_msa_required", self.external_msa_required),
            ("output_dir", self.output_dir),
        ):
            if value not in (None, "", ()):
                lines.append(f"{name}={value}")
        lines.extend(f"message={message}" for message in self.messages)
        return "\n".join(lines)


def _system(path: str | Path) -> System:
    document = load_yaml_document(Path(path))
    if not isinstance(document.value, dict):
        raise ValueError("System YAML root must be a mapping.")
    return System(system=document.value)


def _chain_ids(system_obj: System) -> tuple[str, ...]:
    return tuple(sorted({chain.chain_id for chain in iter_system_chains(system_obj)}))


def _missing_msa(system_obj: System) -> bool:
    sequences = system_obj.find_value(key="sequences") or []
    return any(
        isinstance(entry, dict)
        and isinstance(entry.get("protein"), dict)
        and not str(entry["protein"].get("msa") or "").strip()
        for entry in sequences
    )


def prediction_preflight(
    *,
    workflow: str,
    system_path: str | Path,
    options_path: str | Path,
    runner_name: str,
    repeats: int,
    output_dir: str | Path,
    require_ligand: bool = False,
    assess_bias: bool = False,
    use_bias_databases: bool = True,
    protein_database_path: str | Path | None = None,
    ligand_database_path: str | Path | None = None,
    release_cutoff: str | None = None,
    bias_chains: set[str] | list[str] | tuple[str, ...] | None = None,
    bias_query_cache_path: str | Path | None = None,
    custom_bias_reference_path: str | Path | None = None,
    protein_similarity_threshold: float = 0.25,
    ligand_similarity_threshold: float = 0.35,
) -> tuple[PreflightReport, System, Any]:
    try:
        system_obj = _system(system_path)
        runner = get_runner(runner_name)
        options = runner.load_options(Path(options_path))
        runner.validate_system(
            system_obj,
            options,
            check_atom_names=False,
            source_path=Path(system_path),
            requirements=(
                WorkflowInputRequirements(require_protein=True, require_ligand=True)
                if require_ligand
                else None
            ),
        )
        models = runner.execution_models(options)
        available, availability_message = runner.check_availability()
        messages = [availability_message] if availability_message else []
        databases: list[str] = []
        release_policy = None
        query_hashes: list[str] = []
        cache_path_text = None
        cache_ready = None
        expected_backfill = None
        if assess_bias:
            from cofolder.modules.analytics.bias_database import (
                DEFAULT_LIGAND_DATABASE,
                DEFAULT_PROTEIN_DATABASE,
                LIGAND_TABLE_NAME,
                _selected_queries,
                bias_query_cache_readiness,
                parse_bias_release_policy,
                validate_bias_database_bundle,
            )

            policy = parse_bias_release_policy(release_cutoff or "2023-06-01")
            for name, value in (
                ("protein", protein_similarity_threshold),
                ("ligand", ligand_similarity_threshold),
            ):
                if not 0.0 <= float(value) <= 1.0:
                    raise ValueError(
                        f"Bias {name} similarity threshold must be between 0 and 1."
                    )
            release_policy = policy.mode
            selected = {
                str(value).strip().upper()
                for value in (bias_chains or ())
                if str(value).strip()
            }
            entities = {
                chain.entity_type
                for chain in iter_system_chains(system_obj)
                if not selected or chain.chain_id.strip().upper() in selected
            }
            protein_bundle = None
            ligand_bundle = None
            ligand_smiles_by_id: dict[str, str] = {}
            if use_bias_databases:
                protein_path = protein_database_path or DEFAULT_PROTEIN_DATABASE
                ligand_path = ligand_database_path or DEFAULT_LIGAND_DATABASE
                if "protein" in entities:
                    protein_bundle = validate_bias_database_bundle(protein_path, "protein")
                    databases.append(
                        f"protein:{protein_bundle.root}:{protein_bundle.fingerprint}"
                    )
                if "ligand" in entities:
                    ligand_bundle = validate_bias_database_bundle(ligand_path, "ligand")
                    databases.append(
                        f"ligand:{ligand_bundle.root}:{ligand_bundle.fingerprint}"
                    )
                    import pandas as pd

                    source = pd.read_csv(ligand_bundle.root / LIGAND_TABLE_NAME)
                    if {"ligand_id", "smiles"} <= set(source.columns):
                        ligand_smiles_by_id = {
                            str(row["ligand_id"]).upper(): str(row["smiles"])
                            for _, row in source.dropna(
                                subset=["ligand_id", "smiles"]
                            ).iterrows()
                        }
            proteins, ligands = _selected_queries(
                system_obj, selected, ligand_smiles_by_id
            )
            import hashlib
            from rdkit import Chem

            query_hashes.extend(
                f"protein:{chain}:{hashlib.sha256(''.join(sequence.split()).upper().encode()).hexdigest()}"
                for chain, sequence in sorted(proteins.items())
            )
            for chain, smiles in sorted(ligands.items()):
                molecule = Chem.MolFromSmiles(smiles)
                canonical = Chem.MolToSmiles(molecule, isomericSmiles=True) if molecule else smiles
                query_hashes.append(
                    f"ligand:{chain}:{hashlib.sha256(canonical.encode()).hexdigest()}"
                )
            cache_root, cache_ready, cache_message = bias_query_cache_readiness(
                bias_query_cache_path
            )
            cache_path_text = f"{cache_root} ({cache_message})"
            expected_backfill = bool(
                protein_bundle is not None
                and protein_bundle.sequence_index_path is not None
                and proteins
                and ligands
            )
        if assess_bias and custom_bias_reference_path:
            from cofolder.modules.analytics.bias_database import validate_custom_bias_reference_bundle
            custom = validate_custom_bias_reference_bundle(custom_bias_reference_path)
            databases.append(f"custom:{custom.root}:{custom.fingerprint}")
        report = PreflightReport(
            workflow=workflow,
            ready=bool(available),
            runner=runner_name,
            runner_available=available,
            capabilities=tuple(sorted(runner.capabilities)),
            selected_chains=_chain_ids(system_obj),
            database_sources=tuple(databases),
            release_cutoff=release_cutoff if assess_bias else None,
            release_policy=release_policy,
            query_hashes=tuple(query_hashes),
            bias_query_cache=cache_path_text,
            bias_cache_ready=cache_ready,
            expected_bias_backfill=expected_backfill,
            protein_similarity_threshold=(
                float(protein_similarity_threshold) if assess_bias else None
            ),
            ligand_similarity_threshold=(
                float(ligand_similarity_threshold) if assess_bias else None
            ),
            repeats=int(repeats),
            diffusion_samples=len(models),
            planned_executions=int(repeats) * len(models),
            external_msa_required=_missing_msa(system_obj),
            output_dir=str(Path(output_dir) / "results"),
            messages=tuple(messages),
        )
        return report, system_obj, options
    except Exception as exc:
        return (
            PreflightReport(
                workflow=workflow,
                ready=False,
                runner=runner_name,
                release_cutoff=release_cutoff if assess_bias else None,
                output_dir=str(Path(output_dir) / "results"),
                messages=(str(exc),),
            ),
            System(system={"sequences": []}),
            None,
        )


__all__ = ["PreflightReport", "prediction_preflight"]
