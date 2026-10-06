"""Strict, source-aware loading for runner options YAML."""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import Any, Literal

import yaml

from cofolder.modules.contracts.models import JSONValue

FieldToken = str | int


class InputValidationError(ValueError):
    """Base class for actionable errors raised before backend execution."""

    error_code = "input_validation_failed"

    def __init__(
        self,
        message: str,
        *,
        source_path: Path | str | None = None,
        field_path: tuple[FieldToken, ...] = (),
        line: int | None = None,
        column: int | None = None,
        entity_id: str | None = None,
        chain_id: str | None = None,
        source_record_id: str | None = None,
        residue_position: int | None = None,
        character: str | None = None,
    ) -> None:
        super().__init__(message)
        self.source_path = Path(source_path) if source_path is not None else None
        self.field_path = tuple(field_path)
        self.line = line
        self.column = column
        self.entity_id = entity_id
        self.chain_id = chain_id
        self.source_record_id = source_record_id
        self.residue_position = residue_position
        self.character = character

    @property
    def details(self) -> dict[str, JSONValue]:
        details: dict[str, JSONValue] = {}
        if self.source_path is not None:
            details["source_path"] = str(self.source_path)
        if self.field_path:
            details["field_path"] = list(self.field_path)
        if self.line is not None:
            details["line"] = self.line
        if self.column is not None:
            details["column"] = self.column
        if self.entity_id is not None:
            details["entity_id"] = self.entity_id
        if self.chain_id is not None:
            details["chain_id"] = self.chain_id
        if self.source_record_id is not None:
            details["source_record_id"] = self.source_record_id
        if self.residue_position is not None:
            details["residue_position"] = self.residue_position
        if self.character is not None:
            details["character"] = self.character
        return details


class YamlLoadError(InputValidationError):
    error_code = "yaml_load_failed"


class OptionsValidationError(InputValidationError):
    error_code = "options_validation_failed"


class SystemInputValidationError(InputValidationError):
    error_code = "system_input_validation_failed"


class SequenceValidationError(SystemInputValidationError):
    error_code = "sequence_validation_failed"


class MsaValidationError(SystemInputValidationError):
    error_code = "msa_validation_failed"


class LigandValidationError(SystemInputValidationError):
    error_code = "ligand_validation_failed"


class LigandSelectionError(SystemInputValidationError):
    error_code = "ligand_selection_failed"


@dataclass(frozen=True, slots=True)
class YamlDocument:
    value: Any
    source_path: Path
    locations: Mapping[tuple[FieldToken, ...], tuple[int, int]] = field(
        default_factory=dict
    )

    def location(self, path: tuple[FieldToken, ...]) -> tuple[int | None, int | None]:
        line, column = self.locations.get(path, (None, None))
        return line, column


def _node_locations(
    node: yaml.Node,
    *,
    path: tuple[FieldToken, ...] = (),
    output: dict[tuple[FieldToken, ...], tuple[int, int]],
) -> None:
    output[path] = (node.start_mark.line + 1, node.start_mark.column + 1)
    if isinstance(node, yaml.MappingNode):
        for key_node, value_node in node.value:
            key = key_node.value
            child_path = (*path, key)
            output[child_path] = (
                value_node.start_mark.line + 1,
                value_node.start_mark.column + 1,
            )
            _node_locations(value_node, path=child_path, output=output)
    elif isinstance(node, yaml.SequenceNode):
        for index, value_node in enumerate(node.value):
            _node_locations(value_node, path=(*path, index), output=output)


def load_yaml_document(path: Path | str) -> YamlDocument:
    source_path = Path(path)
    try:
        text = source_path.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as exc:
        raise YamlLoadError(
            f"Unable to read YAML file {source_path}: {exc}", source_path=source_path
        ) from exc
    try:
        value = yaml.safe_load(text)
        node = yaml.compose(text, Loader=yaml.SafeLoader)
    except yaml.YAMLError as exc:
        mark = getattr(exc, "problem_mark", None)
        line = mark.line + 1 if mark is not None else None
        column = mark.column + 1 if mark is not None else None
        location = f" at line {line}, column {column}" if line is not None else ""
        problem = getattr(exc, "problem", None) or str(exc)
        raise YamlLoadError(
            f"Malformed YAML in {source_path}{location}: {problem}",
            source_path=source_path,
            line=line,
            column=column,
        ) from exc
    locations: dict[tuple[FieldToken, ...], tuple[int, int]] = {}
    if node is not None:
        _node_locations(node, output=locations)
    return YamlDocument(value=value, source_path=source_path, locations=locations)


@dataclass(frozen=True, slots=True)
class RunnerRuntimeOptions:
    cache_path: Path | None = None
    diffusion_samples: int = 1
    executable: str | None = None
    subcommand: str | None = None
    extra_args: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class RunnerOptions:
    version: Literal[1]
    runtime: RunnerRuntimeOptions
    runner: Mapping[str, JSONValue]

    def find_value(self, key: str | None = None, path: list[str | int] | None = None) -> Any:
        combined: dict[str, Any] = {
            "runtime": {
                "cache_path": str(self.runtime.cache_path) if self.runtime.cache_path else None,
                "diffusion_samples": self.runtime.diffusion_samples,
                "executable": self.runtime.executable,
                "subcommand": self.runtime.subcommand,
                "extra_args": list(self.runtime.extra_args),
            },
            "runner": dict(self.runner),
        }
        if path:
            current: Any = combined
            for token in path:
                current = current[token] if isinstance(current, dict) else current[int(token)]
            return current
        aliases = {"cache": "cache_path", "samples_per_seed": "diffusion_samples"}
        key = aliases.get(str(key), key)
        runtime_value = getattr(self.runtime, str(key), None)
        if runtime_value is not None:
            return str(runtime_value) if isinstance(runtime_value, Path) else runtime_value
        return self.runner.get(str(key))


@dataclass(frozen=True, slots=True)
class OptionField:
    value_types: tuple[type, ...]
    required: bool = False
    choices: frozenset[JSONValue] | None = None
    minimum: float | None = None
    fields: Mapping[str, OptionField] | None = None
    allow_unknown_fields: bool = False
    item_types: tuple[type, ...] | None = None


@dataclass(frozen=True, slots=True)
class RunnerOptionsSchema:
    runtime_fields: Mapping[str, OptionField]
    runner_fields: Mapping[str, OptionField]

    def _error(
        self, document: YamlDocument, path: tuple[FieldToken, ...], message: str
    ) -> OptionsValidationError:
        line, column = document.location(path)
        dotted = ".".join(str(token) for token in path)
        return OptionsValidationError(
            f"Invalid options in {document.source_path} at '{dotted}': {message}",
            source_path=document.source_path,
            field_path=path,
            line=line,
            column=column,
        )

    def _validate_mapping(
        self,
        document: YamlDocument,
        section: str,
        value: Any,
        specs: Mapping[str, OptionField],
    ) -> dict[str, Any]:
        if not isinstance(value, dict):
            raise self._error(document, (section,), "must be a mapping")
        non_string_keys = [key for key in value if not isinstance(key, str)]
        if non_string_keys:
            raise self._error(document, (section,), "keys must be strings")
        unknown = sorted(set(value) - set(specs))
        if unknown:
            raise self._error(
                document, (section, unknown[0]), f"unknown key {unknown[0]!r}"
            )
        for name, spec in specs.items():
            if spec.required and name not in value:
                raise self._error(document, (section, name), "required key is missing")
            if name not in value:
                continue
            self._validate_field(document, (section, name), value[name], spec)
        return dict(value)

    def _validate_field(
        self,
        document: YamlDocument,
        path: tuple[FieldToken, ...],
        item: Any,
        spec: OptionField,
    ) -> None:
        if type(item) not in spec.value_types:
            expected = ", ".join(item_type.__name__ for item_type in spec.value_types)
            raise self._error(document, path, f"must be of type {expected}")
        if spec.choices is not None and item not in spec.choices:
            raise self._error(
                document, path, f"must be one of {sorted(spec.choices)!r}"
            )
        if spec.minimum is not None and item is not None:
            if isinstance(item, float) and not math.isfinite(item):
                raise self._error(document, path, "must be finite")
            if item < spec.minimum:
                raise self._error(
                    document, path, f"must be at least {spec.minimum:g}"
                )
        if isinstance(item, list) and spec.item_types is not None:
            for index, child in enumerate(item):
                if type(child) not in spec.item_types:
                    expected = ", ".join(
                        item_type.__name__ for item_type in spec.item_types
                    )
                    raise self._error(
                        document, (*path, index), f"must be of type {expected}"
                    )
        if isinstance(item, dict):
            if any(not isinstance(key, str) for key in item):
                raise self._error(document, path, "keys must be strings")
            fields = spec.fields or {}
            unknown = sorted(set(item) - set(fields))
            if unknown and not spec.allow_unknown_fields:
                raise self._error(
                    document, (*path, unknown[0]), f"unknown key {unknown[0]!r}"
                )
            for name, child_spec in fields.items():
                if child_spec.required and name not in item:
                    raise self._error(
                        document, (*path, name), "required key is missing"
                    )
                if name in item:
                    self._validate_field(
                        document, (*path, name), item[name], child_spec
                    )

    def validate(self, document: YamlDocument) -> RunnerOptions:
        root = document.value
        if not isinstance(root, dict):
            raise self._error(document, (), "document root must be a mapping")
        if any(not isinstance(key, str) for key in root):
            raise self._error(document, (), "root keys must be strings")
        unknown = sorted(set(root) - {"version", "runtime", "runner"})
        if unknown:
            raise self._error(document, (unknown[0],), f"unknown key {unknown[0]!r}")
        if root.get("version") != 1 or type(root.get("version")) is not int:
            raise self._error(document, ("version",), "must be the integer 1")
        if "runner" not in root:
            raise self._error(document, ("runner",), "required key is missing")
        runtime = self._validate_mapping(
            document, "runtime", root.get("runtime", {}), self.runtime_fields
        )
        runner = self._validate_mapping(
            document, "runner", root["runner"], self.runner_fields
        )
        extra_args = runtime.get("extra_args", [])
        if any(type(item) is not str for item in extra_args):
            raise self._error(document, ("runtime", "extra_args"), "must contain only strings")
        return RunnerOptions(
            version=1,
            runtime=RunnerRuntimeOptions(
                cache_path=Path(runtime["cache_path"]).expanduser()
                if runtime.get("cache_path")
                else None,
                diffusion_samples=runtime.get("diffusion_samples", 1),
                executable=runtime.get("executable"),
                subcommand=runtime.get("subcommand"),
                extra_args=tuple(extra_args),
            ),
            runner=MappingProxyType(runner),
        )


def load_runner_options(
    path: Path | str, *, schema: RunnerOptionsSchema
) -> RunnerOptions:
    return schema.validate(load_yaml_document(path))


_POSITIVE_INT = OptionField((int,), minimum=1)
_NONNEGATIVE_INT = OptionField((int,), minimum=0)
_BOOL = OptionField((bool,))
_STRING = OptionField((str,))
_NUMBER = OptionField((int, float), minimum=0)
_NULLABLE_STRING = OptionField((str, type(None)))
_NULLABLE_INT = OptionField((int, type(None)), minimum=0)
_STRING_LIST = OptionField((list,), item_types=(str,))
_INT_LIST = OptionField((list,), item_types=(int,))

COMMON_RUNTIME_FIELDS: Mapping[str, OptionField] = MappingProxyType(
    {"cache_path": _STRING, "diffusion_samples": _POSITIVE_INT}
)

BOLTZ_RUNNER_FIELDS: Mapping[str, OptionField] = MappingProxyType(
    {
        "checkpoint": _STRING,
        "devices": _POSITIVE_INT,
        "accelerator": OptionField((str,), choices=frozenset({"gpu", "cpu", "tpu"})),
        "recycling_steps": _NONNEGATIVE_INT,
        "sampling_steps": _POSITIVE_INT,
        "max_parallel_samples": _POSITIVE_INT,
        "step_scale": _NUMBER,
        "write_full_pae": _BOOL,
        "write_full_pde": _BOOL,
        "output_format": OptionField((str,), choices=frozenset({"pdb", "mmcif"})),
        "num_workers": _NONNEGATIVE_INT,
        "override": _BOOL,
        "msa_server_url": _STRING,
        "msa_pairing_strategy": OptionField((str,), choices=frozenset({"greedy", "complete"})),
        "msa_server_username": _STRING,
        "msa_server_password": _STRING,
        "api_key_header": _STRING,
        "api_key_value": _STRING,
        "use_potentials": _BOOL,
        "method": _STRING,
        "preprocessing_threads": _POSITIVE_INT,
        "affinity_mw_correction": _BOOL,
        "sampling_steps_affinity": _POSITIVE_INT,
        "diffusion_samples_affinity": _POSITIVE_INT,
        "affinity_checkpoint": _STRING,
        "max_msa_seqs": _POSITIVE_INT,
        "subsample_msa": _BOOL,
        "num_subsampled_msa": _POSITIVE_INT,
        "no_kernels": _BOOL,
        "write_embeddings": _BOOL,
    }
)

BOLTZ_OPTIONS_SCHEMA = RunnerOptionsSchema(COMMON_RUNTIME_FIELDS, BOLTZ_RUNNER_FIELDS)

_OPENFOLD_MSA_MAX_COUNTS = OptionField(
    (dict,),
    fields=MappingProxyType(
        {
            name: _NONNEGATIVE_INT
            for name in (
                "uniref90_hits",
                "uniprot_hits",
                "bfd_uniclust_hits",
                "bfd_uniref_hits",
                "cfdb_uniref30",
                "mgnify_hits",
                "rfam_hits",
                "rnacentral_hits",
                "nt_hits",
                "nucleotide_collection_hits",
                "concat_cfdb_uniref100_filtered",
                "mmseqs_colabfold",
                "colabfold_main",
                "colabfold_paired",
            )
        }
    ),
)

_OPENFOLD_DATASET_FIELDS = MappingProxyType(
    {
        "ccd_file_path": _NULLABLE_STRING,
        "msa": OptionField(
            (dict,),
            fields=MappingProxyType(
                {
                    "max_rows_paired": _NONNEGATIVE_INT,
                    "max_rows": _NONNEGATIVE_INT,
                    "subsample_with_bands": _BOOL,
                    "min_chains_paired_partial": _NONNEGATIVE_INT,
                    "pairing_mask_keys": _STRING_LIST,
                    "moltypes": _INT_LIST,
                    "max_seq_counts": _OPENFOLD_MSA_MAX_COUNTS,
                    "msas_to_pair": _STRING_LIST,
                    "aln_order": _STRING_LIST,
                    "paired_msa_order": _STRING_LIST,
                }
            ),
        ),
        "template": OptionField(
            (dict,),
            fields=MappingProxyType(
                {"n_templates": _NONNEGATIVE_INT, "take_top_k": _BOOL}
            ),
        ),
        "distogram": OptionField(
            (dict,),
            fields=MappingProxyType(
                {"min_bin": _NUMBER, "max_bin": _NUMBER, "n_bins": _POSITIVE_INT}
            ),
        ),
        "pocket_sampling": OptionField(
            (dict,),
            fields=MappingProxyType(
                {
                    "enabled": _BOOL,
                    "num_parents": _POSITIVE_INT,
                    "candidates": _POSITIVE_INT,
                    "noise_frac": _NUMBER,
                    "ligand_jitter": _NUMBER,
                    "center_jitter": _NUMBER,
                    "surface_jitter": _NUMBER,
                    "vdw_buffer": _NUMBER,
                    "diversity_rmsd": _NUMBER,
                    "rdkit_num_conformers": _POSITIVE_INT,
                    "rdkit_conformer_rng": _NONNEGATIVE_INT,
                    "rdkit_conformer_prune_rmsd": _NUMBER,
                    "rdkit_conformer_max_iters": _POSITIVE_INT,
                }
            ),
        ),
    }
)

_OPENFOLD_TEMPLATE_FIELDS = MappingProxyType(
    {
        "mode": OptionField((str,), choices=frozenset({"predict"})),
        "moltypes": _INT_LIST,
        "max_sequences_parse": _POSITIVE_INT,
        "max_seq_id": OptionField((int, float, type(None)), minimum=0),
        "min_align": OptionField((int, float, type(None)), minimum=0),
        "min_len": OptionField((int, float, type(None)), minimum=0),
        "max_release_date": _NULLABLE_STRING,
        "min_release_date_diff": _NULLABLE_INT,
        "max_templates": _NONNEGATIVE_INT,
        "fetch_missing_structures": _BOOL,
        "create_precache": _BOOL,
        "preparse_structures": _BOOL,
        "create_logs": _BOOL,
        "n_processes": _POSITIVE_INT,
        "chunksize": _POSITIVE_INT,
        "preprocess_timeout": _POSITIVE_INT,
        "structure_directory": _NULLABLE_STRING,
        "structure_file_format": OptionField(
            (str,), choices=frozenset({"cif", "pdb", "cif.gz"})
        ),
        "output_directory": _NULLABLE_STRING,
        "precache_directory": _NULLABLE_STRING,
        "structure_array_directory": _NULLABLE_STRING,
        "cache_directory": _NULLABLE_STRING,
        "log_directory": _NULLABLE_STRING,
        "ccd_file_path": _NULLABLE_STRING,
    }
)

OPENFOLD3_OPTIONS_SCHEMA = RunnerOptionsSchema(
    MappingProxyType(
        {
            **COMMON_RUNTIME_FIELDS,
            "executable": _STRING,
            "subcommand": _STRING,
            "extra_args": OptionField((list,)),
        }
    ),
    MappingProxyType(
        {
            "experiment_settings": OptionField(
                (dict,), fields=MappingProxyType({"use_templates": _BOOL})
            ),
            "pl_trainer_args": OptionField(
                (dict,),
                fields=MappingProxyType(
                    {
                        "max_epochs": _POSITIVE_INT,
                        "accelerator": OptionField(
                            (str,), choices=frozenset({"gpu", "cpu", "tpu", "auto"})
                        ),
                        "precision": OptionField((str, int)),
                        "num_nodes": _POSITIVE_INT,
                        "devices": _POSITIVE_INT,
                        "profiler": _NULLABLE_STRING,
                        "log_every_n_steps": _POSITIVE_INT,
                        "enable_checkpointing": _BOOL,
                        "enable_model_summary": _BOOL,
                        "deepspeed_config_path": _NULLABLE_STRING,
                        "distributed_timeout": _STRING,
                        "mpi_plugin": _BOOL,
                    }
                ),
            ),
            "model_update": OptionField(
                (dict,),
                fields=MappingProxyType(
                    {
                        "presets": _STRING_LIST,
                        # OpenFold intentionally defines this node as a model
                        # configuration update tree. Its children are dynamic,
                        # while every surrounding native section is explicit.
                        "custom": OptionField((dict,), allow_unknown_fields=True),
                    }
                ),
            ),
            "data_module_args": OptionField(
                (dict,),
                fields=MappingProxyType(
                    {
                        "batch_size": _POSITIVE_INT,
                        "num_workers": _NONNEGATIVE_INT,
                        "num_workers_validation": _NONNEGATIVE_INT,
                        "epoch_len": _POSITIVE_INT,
                    }
                ),
            ),
            "dataset_config_kwargs": OptionField(
                (dict,), fields=_OPENFOLD_DATASET_FIELDS
            ),
            "output_writer_settings": OptionField(
                (dict,),
                fields=MappingProxyType(
                    {
                        "structure_format": OptionField(
                            (str,), choices=frozenset({"cif", "pdb", "cif.gz"})
                        ),
                        "full_confidence_output_format": OptionField(
                            (str,), choices=frozenset({"json", "npz"})
                        ),
                        "full_confidence_output_dtype": OptionField(
                            (str,), choices=frozenset({"float32", "float16"})
                        ),
                        "write_features": _BOOL,
                        "write_latent_outputs": _BOOL,
                        "write_full_confidence_scores": _BOOL,
                    }
                ),
            ),
            "msa_computation_settings": OptionField(
                (dict,),
                fields=MappingProxyType(
                    {
                        "msa_file_format": OptionField(
                            (str,), choices=frozenset({"npz", "a3m", "csv"})
                        ),
                        "server_user_agent": _STRING,
                        "server_url": _STRING,
                        "save_mappings": _BOOL,
                        "msa_output_directory": _NULLABLE_STRING,
                        "cleanup_msa_dir": _BOOL,
                        "save_openfold_outputs": _BOOL,
                        "save_colabfold_outputs": _BOOL,
                        "colabfold_output_dir": _NULLABLE_STRING,
                    }
                ),
            ),
            "template_preprocessor_settings": OptionField(
                (dict,), fields=_OPENFOLD_TEMPLATE_FIELDS
            ),
        }
    ),
)
