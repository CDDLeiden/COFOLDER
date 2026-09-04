"""Single-ligand Oracle scoring backed by the Validate workflow."""

from __future__ import annotations

import copy
import json
import logging
import math
import operator
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Any, Callable, Mapping, Sequence

import pandas as pd

from cofolder.modules.input import system
from cofolder.modules.utils import read, write
from cofolder.recipes._metrics import (
    collect_qualified_metric_values,
    read_metric_frames,
)
from cofolder.recipes.validate import DEFAULT_SCORING_FUNCTIONS, Validate

logger = logging.getLogger(__name__)


class _MetricUnavailableError(ValueError):
    """A selected metric is absent or has no finite numeric values."""


@dataclass(frozen=True, slots=True)
class OracleGate:
    """A threshold comparison applied to one aggregated Oracle metric."""

    metric: str
    operator: str
    threshold: float

    def __post_init__(self) -> None:
        metric = str(self.metric).strip()
        comparison = str(self.operator).strip().lower()
        try:
            threshold = float(self.threshold)
        except (TypeError, ValueError) as exc:
            raise ValueError("OracleGate.threshold must be a finite number.") from exc
        if not metric:
            raise ValueError("OracleGate.metric must not be empty.")
        if comparison not in {"gt", "ge", "lt", "le"}:
            raise ValueError("OracleGate.operator must be one of: gt, ge, lt, le.")
        if not math.isfinite(threshold):
            raise ValueError("OracleGate.threshold must be a finite number.")
        object.__setattr__(self, "metric", metric)
        object.__setattr__(self, "operator", comparison)
        object.__setattr__(self, "threshold", threshold)


@dataclass(frozen=True, slots=True)
class OracleGatePolicy:
    """The score transformation used when one or more Oracle gates fail."""

    mode: str
    value: float

    def __post_init__(self) -> None:
        mode = str(self.mode).strip().lower()
        try:
            value = float(self.value)
        except (TypeError, ValueError) as exc:
            raise ValueError("OracleGatePolicy.value must be a finite number.") from exc
        if mode not in {"downweight", "fixed_penalty", "non_binder"}:
            raise ValueError(
                "OracleGatePolicy.mode must be one of: downweight, fixed_penalty, non_binder."
            )
        if not math.isfinite(value):
            raise ValueError("OracleGatePolicy.value must be a finite number.")
        if mode == "downweight" and not 0 <= value <= 1:
            raise ValueError("A downweight policy value must be in [0, 1].")
        object.__setattr__(self, "mode", mode)
        object.__setattr__(self, "value", value)


@dataclass(frozen=True, slots=True)
class OracleScoreContext:
    """Inputs made available to an arbitrary Python Oracle scoring function."""

    query_smiles: str
    run_dir: Path
    system_metrics: pd.DataFrame
    chain_metrics: pd.DataFrame
    aggregated_metrics: Mapping[str, float]


class Oracle:
    """Run Validate once for one ligand and return a finite scalar score."""

    _AGGREGATORS = {"first", "mean", "max", "min", "median"}
    _COMPARISONS = {
        "gt": operator.gt,
        "ge": operator.ge,
        "lt": operator.lt,
        "le": operator.le,
    }

    def __init__(
        self,
        wrk_dir: str,
        system_path: str,
        options_path: str,
        runner: str = "boltz2",
        input_smiles: str | None = None,
        input_mol_file: str | None = None,
        output_metric: str | None = None,
        aggregate: str = "first",
        repeats: int = 1,
        seed: int | None = None,
        scoring_functions: list[str] | None = None,
        assess_robustness: bool = True,
        assess_bias: bool = False,
        protein_training_data_path: str | None = None,
        ligand_training_data_path: str | None = None,
        bias_release_cutoff: str = "2023-06-01",
        bias_ligand_similarity_threshold: float = 0.35,
        bias_chains: list[str] | None = None,
        build_bias_training_data: bool = False,
        bias_training_components_cif: str | None = None,
        conformers: str | None = None,
        sdf_file: str | None = None,
        reference_path: str | None = None,
        pocket_coverage_reference: str | None = None,
        reproduction_metrics: list[str] | None = None,
        score_components: Mapping[str, float] | None = None,
        scoring_function: Callable[[OracleScoreContext], float] | None = None,
        score_gates: Sequence[OracleGate] | None = None,
        gate_policy: OracleGatePolicy | None = None,
    ):
        self.wrk_dir = Path(wrk_dir)
        self.system_path = Path(system_path)
        self.options_path = Path(options_path)
        self.runner = str(runner)
        self.input_smiles = input_smiles.strip() if input_smiles else None
        self.input_mol_file = Path(input_mol_file) if input_mol_file else None
        self.output_metric = str(output_metric).strip() if output_metric else None
        self.aggregate = str(aggregate).strip().lower()
        self.score_components = self._normalize_components(score_components)
        self.scoring_function = scoring_function
        self.score_gates = tuple(score_gates or ())
        self.gate_policy = gate_policy

        effective_scoring_functions = (
            scoring_functions
            if scoring_functions is not None
            else sorted(DEFAULT_SCORING_FUNCTIONS)
        )
        self.validate_kwargs: dict[str, Any] = {
            "repeats": repeats,
            "seed": seed,
            "scoring_functions": effective_scoring_functions,
            "assess_robustness": assess_robustness,
            "assess_bias": assess_bias,
            "protein_training_data_path": protein_training_data_path,
            "ligand_training_data_path": ligand_training_data_path,
            "bias_release_cutoff": bias_release_cutoff,
            "bias_ligand_similarity_threshold": bias_ligand_similarity_threshold,
            "bias_chains": bias_chains,
            "build_bias_training_data": build_bias_training_data,
            "bias_training_components_cif": bias_training_components_cif,
            "conformers": conformers,
            "sdf_file": sdf_file,
            "reference_path": reference_path,
            "pocket_coverage_reference": pocket_coverage_reference,
            "reproduction_metrics": reproduction_metrics,
        }
        self.base_system = read.read_yaml(path=self.system_path)
        self.query_ligand_chain = self._find_query_ligand_chain(self.base_system)
        self.logger = logging.getLogger("cofolder.oracle")
        self._validate_config()

    @staticmethod
    def _normalize_components(
        components: Mapping[str, float] | None,
    ) -> dict[str, float]:
        normalized: dict[str, float] = {}
        for metric, raw_weight in (components or {}).items():
            name = str(metric).strip()
            if not name:
                raise ValueError("score_components metric names must not be empty.")
            try:
                weight = float(raw_weight)
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"Weight for score component '{name}' must be a finite number."
                ) from exc
            if not math.isfinite(weight):
                raise ValueError(
                    f"Weight for score component '{name}' must be a finite number."
                )
            normalized[name] = weight
        return normalized

    def _validate_config(self) -> None:
        if bool(self.input_smiles) == bool(self.input_mol_file):
            raise ValueError(
                "Provide exactly one input: --input_smiles or --input_mol_file."
            )
        if self.input_mol_file is not None:
            if not self.input_mol_file.exists():
                raise ValueError(
                    f"--input_mol_file does not exist: {self.input_mol_file}"
                )
            if not self.input_mol_file.is_file():
                raise ValueError(
                    f"--input_mol_file is not a file: {self.input_mol_file}"
                )
        sources = sum(
            (
                bool(self.output_metric),
                bool(self.score_components),
                self.scoring_function is not None,
            )
        )
        if sources != 1:
            raise ValueError(
                "Provide exactly one score source: output_metric, score_components, or scoring_function."
            )
        if self.scoring_function is not None and not callable(self.scoring_function):
            raise ValueError("scoring_function must be callable.")
        if self.aggregate not in self._AGGREGATORS:
            raise ValueError(
                f"Invalid --aggregate '{self.aggregate}'. Choose from: {sorted(self._AGGREGATORS)}"
            )
        if any(not isinstance(gate, OracleGate) for gate in self.score_gates):
            raise ValueError("score_gates must contain OracleGate instances.")
        if self.score_gates and not isinstance(self.gate_policy, OracleGatePolicy):
            raise ValueError("gate_policy is required when score_gates are configured.")
        if not self.score_gates and self.gate_policy is not None:
            raise ValueError("gate_policy requires at least one score gate.")

        selectors = [
            *self.score_components,
            *(gate.metric for gate in self.score_gates),
        ]
        if self.output_metric:
            selectors.append(self.output_metric)
        for selector in selectors:
            self._validate_metric_dependency(selector)

    def _validate_metric_dependency(self, selector: str) -> None:
        metric = str(selector).rsplit("__", 1)[-1]
        scoring = set(self.validate_kwargs.get("scoring_functions") or [])
        if metric == "ifp_distance":
            if "ifp_distance" not in scoring:
                raise ValueError(
                    f"Metric '{selector}' requires one of scoring functions "
                    "['ifp_distance']."
                )
            raise ValueError(
                "Metric 'ifp_distance' is vector-valued and cannot be an Oracle scalar; "
                "use pocket_coverage_custom for an IFP-reference overlap score."
            )
        if metric.startswith("bias_"):
            if not self.validate_kwargs.get("assess_bias"):
                raise ValueError(f"Metric '{selector}' requires assess_bias=True.")
            return

        required_group: set[str] | None = None
        if metric in {"sasa", "sasa_norm_heavy"}:
            required_group = {"sasa", "sasa_normalized"}
        elif metric in {
            "affinity_pred_value",
            "affinity_probability_binary",
            "pIC50",
            "IC50_M",
            "pIC50_kcal_per_mol",
        }:
            required_group = {"affinity_metrics", "affinity_metrics_ext"}
        elif metric.startswith("pair_chains_iptm_") or metric in {
            "chains_ptm",
            "ptm",
            "iptm",
            "confidence_score",
        }:
            required_group = {"confidence_metrics"}
        elif metric == "pocket_coverage_custom":
            required_group = {"ifp_distance"}
            if not self.validate_kwargs.get("pocket_coverage_reference"):
                raise ValueError(
                    f"Metric '{selector}' requires pocket_coverage_reference."
                )
            reproduction = self.validate_kwargs.get("reproduction_metrics")
            if reproduction is not None and "pocket_coverage" not in reproduction:
                raise ValueError(
                    f"Metric '{selector}' requires the pocket_coverage reproduction metric."
                )
        elif metric == "pocket_coverage_ref":
            if not self.validate_kwargs.get("reference_path"):
                raise ValueError(f"Metric '{selector}' requires reference_path.")
            reproduction = self.validate_kwargs.get("reproduction_metrics")
            if reproduction is not None and "pocket_coverage" not in reproduction:
                raise ValueError(
                    f"Metric '{selector}' requires the pocket_coverage reproduction metric."
                )
        if required_group is not None and not (required_group & scoring):
            raise ValueError(
                f"Metric '{selector}' requires one of scoring functions {sorted(required_group)}, "
                f"but enabled are {sorted(scoring)}."
            )

    def run(self) -> float:
        self.wrk_dir.mkdir(parents=True, exist_ok=True)
        run_dir = self.wrk_dir / "oracle_run"
        run_dir.mkdir(parents=True, exist_ok=True)
        smiles = self.input_smiles or self._smiles_from_mol_file(self.input_mol_file)
        if not smiles:
            raise ValueError("Failed to resolve valid input SMILES for oracle run.")

        sys_obj = system.System(system=copy.deepcopy(self.base_system))
        ligand_path = self._find_first_ligand_smiles_path(sys_obj.system)
        if ligand_path is None:
            raise ValueError("No ligand SMILES field found in system YAML to update.")
        sys_obj.update_system(value=smiles, path=ligand_path)
        row_system_path = run_dir / "oracle_system.yaml"
        write.write_yaml(sys_obj, path=row_system_path)
        Validate(
            wrk_dir=str(run_dir),
            system_path=str(row_system_path),
            options_path=str(self.options_path),
            runner=self.runner,
            **self.validate_kwargs,
        ).run()

        system_df, chain_df = read_metric_frames(run_dir)
        raw_metrics = collect_qualified_metric_values(system_df, chain_df)
        context = OracleScoreContext(
            query_smiles=smiles,
            run_dir=run_dir,
            system_metrics=system_df.copy(),
            chain_metrics=chain_df.copy(),
            aggregated_metrics=MappingProxyType(
                self._aggregate_all_metrics(raw_metrics)
            ),
        )
        base_value, component_values, score_mode = self._calculate_base_score(
            context, raw_metrics
        )
        final_value, gate_values, failed_gates, action = self._apply_gates(
            base_value, raw_metrics
        )
        result = {
            "output_metric": self.output_metric,
            "aggregate": self.aggregate,
            "value": final_value,
            "score_mode": score_mode,
            "base_value": base_value,
            "component_values": json.dumps(component_values, sort_keys=True),
            "gate_pass": not failed_gates,
            "gate_values": json.dumps(gate_values, sort_keys=True),
            "failed_gates": json.dumps(failed_gates),
            "gate_action": action,
        }
        pd.DataFrame([result]).to_csv(self.wrk_dir / "oracle_result.csv", index=False)
        return final_value

    def _calculate_base_score(self, context, raw_metrics):
        if self.output_metric:
            value = self._metric_value(self.output_metric, raw_metrics)
            return value, {self.output_metric: value}, "metric"
        if self.score_components:
            values = {
                name: self._metric_value(name, raw_metrics)
                for name in self.score_components
            }
            score = sum(
                values[name] * weight for name, weight in self.score_components.items()
            )
            return self._finite_score(score, "Composite score"), values, "composite"
        assert self.scoring_function is not None
        return (
            self._finite_score(
                self.scoring_function(context), "Custom scoring function"
            ),
            {},
            "custom",
        )

    def _apply_gates(self, base_value, raw_metrics):
        gate_values: dict[str, float | None] = {}
        failed: list[str] = []
        for gate in self.score_gates:
            try:
                value = self._metric_value(gate.metric, raw_metrics)
            except _MetricUnavailableError:
                value = None
            gate_values[gate.metric] = value
            if value is None or not self._COMPARISONS[gate.operator](
                value, gate.threshold
            ):
                failed.append(gate.metric)
        if not failed:
            return base_value, gate_values, failed, "none"
        assert self.gate_policy is not None
        final = (
            base_value * self.gate_policy.value
            if self.gate_policy.mode == "downweight"
            else self.gate_policy.value
        )
        return (
            self._finite_score(final, "Gate-adjusted score"),
            gate_values,
            failed,
            self.gate_policy.mode,
        )

    def _aggregate_all_metrics(self, raw_metrics):
        aggregated: dict[str, float] = {}
        for selector, values in raw_metrics.items():
            try:
                aggregated[selector] = self._aggregate_values(selector, values)
            except ValueError:
                continue
        return aggregated

    def _metric_value(self, selector, raw_metrics):
        qualified = self._resolve_selector(selector, raw_metrics)
        return self._aggregate_values(qualified, raw_metrics[qualified])

    def _resolve_selector(self, selector, raw_metrics):
        name = str(selector).strip()
        if "__" in name:
            if name not in raw_metrics:
                raise _MetricUnavailableError(
                    f"Metric '{name}' was not found in Oracle outputs."
                )
            return name
        if self.query_ligand_chain:
            ligand_name = f"ligand_{self.query_ligand_chain}__{name}"
            if ligand_name in raw_metrics:
                return ligand_name
        ligand_matches = sorted(
            metric
            for metric in raw_metrics
            if metric.startswith("ligand_") and metric.endswith(f"__{name}")
        )
        if len(ligand_matches) == 1:
            return ligand_matches[0]
        system_name = f"system__{name}"
        if system_name in raw_metrics:
            return system_name
        matches = sorted(
            metric for metric in raw_metrics if metric.endswith(f"__{name}")
        )
        if len(matches) == 1:
            return matches[0]
        if len(matches) > 1:
            raise ValueError(
                f"Metric '{name}' is ambiguous; use a qualified selector such as '{matches[0]}'."
            )
        raise _MetricUnavailableError(
            f"Metric '{name}' was not found in Oracle outputs."
        )

    def _aggregate_values(self, selector, values):
        numeric = pd.to_numeric(pd.Series(list(values)), errors="coerce").dropna()
        finite = [float(value) for value in numeric if math.isfinite(float(value))]
        if not finite:
            raise _MetricUnavailableError(
                f"Metric '{selector}' was not found with finite numeric values in Oracle outputs."
            )
        if self.aggregate == "first":
            return finite[0]
        if self.aggregate == "mean":
            return float(sum(finite) / len(finite))
        if self.aggregate == "max":
            return max(finite)
        if self.aggregate == "min":
            return min(finite)
        if self.aggregate == "median":
            ordered = sorted(finite)
            middle = len(ordered) // 2
            return (
                ordered[middle]
                if len(ordered) % 2
                else float((ordered[middle - 1] + ordered[middle]) / 2)
            )
        raise ValueError(f"Unsupported aggregate: {self.aggregate}")

    @staticmethod
    def _finite_score(value, source):
        if isinstance(value, bool):
            raise ValueError(f"{source} must return a finite numeric scalar, not bool.")
        try:
            score = float(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{source} must return a finite numeric scalar.") from exc
        if not math.isfinite(score):
            raise ValueError(f"{source} must return a finite numeric scalar.")
        return score

    @staticmethod
    def _find_query_ligand_chain(system_dict):
        path = Oracle._find_first_ligand_smiles_path(system_dict)
        if path is None:
            return None
        ligand = system_dict["sequences"][path[1]]["ligand"]
        identifier = ligand.get("id") if isinstance(ligand, dict) else None
        if isinstance(identifier, list):
            return str(identifier[0]) if len(identifier) == 1 else None
        return str(identifier) if identifier is not None else None

    @staticmethod
    def _find_first_ligand_smiles_path(system_dict):
        sequences = system_dict.get("sequences", [])
        if not isinstance(sequences, list):
            return None
        for index, entry in enumerate(sequences):
            ligand = entry.get("ligand") if isinstance(entry, dict) else None
            if isinstance(ligand, dict) and "smiles" in ligand:
                return ["sequences", index, "ligand", "smiles"]
        return None

    @staticmethod
    def _smiles_from_mol_file(path):
        if path is None:
            return None
        try:
            from rdkit import Chem
        except Exception as exc:
            raise ValueError("RDKit is required to parse --input_mol_file.") from exc
        suffix = path.suffix.lower()
        mol = None
        if suffix == ".sdf":
            mol = next(
                (
                    item
                    for item in Chem.SDMolSupplier(str(path), removeHs=False)
                    if item is not None
                ),
                None,
            )
        elif suffix in {".mol", ".mol2"}:
            mol = Chem.MolFromMolFile(str(path), removeHs=False)
        else:
            mol = next(
                (
                    item
                    for item in Chem.SDMolSupplier(str(path), removeHs=False)
                    if item is not None
                ),
                None,
            )
            if mol is None:
                mol = Chem.MolFromMolFile(str(path), removeHs=False)
        return None if mol is None else Chem.MolToSmiles(mol)
