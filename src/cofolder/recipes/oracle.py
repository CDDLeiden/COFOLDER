"""Oracle recipe implemented as a single-run Validate wrapper."""

from __future__ import annotations

import copy
import logging
from pathlib import Path
from typing import Any

import pandas as pd

from cofolder.modules.input import system
from cofolder.modules.utils import read, write
from cofolder.recipes.validate import Validate

logger = logging.getLogger(__name__)


class Oracle:
    """Run validate once for a single input ligand and return a single score."""

    _AGGREGATORS = {"first", "mean", "max", "min", "median"}

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
    ):
        self.wrk_dir = Path(wrk_dir)
        self.system_path = Path(system_path)
        self.options_path = Path(options_path)
        self.runner = str(runner)

        self.input_smiles = input_smiles.strip() if input_smiles else None
        self.input_mol_file = Path(input_mol_file) if input_mol_file else None
        self.output_metric = str(output_metric).strip() if output_metric else None
        self.aggregate = str(aggregate).strip().lower()

        self.validate_kwargs: dict[str, Any] = {
            "repeats": repeats,
            "seed": seed,
            "scoring_functions": scoring_functions,
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
        self.logger = logging.getLogger("cofolder.oracle")
        self._validate_config()

    def _validate_config(self) -> None:
        if bool(self.input_smiles) == bool(self.input_mol_file):
            raise ValueError("Provide exactly one input: --input_smiles or --input_mol_file.")
        if self.input_mol_file is not None:
            if not self.input_mol_file.exists():
                raise ValueError(f"--input_mol_file does not exist: {self.input_mol_file}")
            if not self.input_mol_file.is_file():
                raise ValueError(f"--input_mol_file is not a file: {self.input_mol_file}")
        if not self.output_metric:
            raise ValueError("--output_metric is required.")
        if self.aggregate not in self._AGGREGATORS:
            raise ValueError(
                f"Invalid --aggregate '{self.aggregate}'. "
                f"Choose from: {sorted(self._AGGREGATORS)}"
            )
        self._validate_metric_supported_by_scoring_functions()

    def _validate_metric_supported_by_scoring_functions(self) -> None:
        metric = str(self.output_metric)
        scoring = set(self.validate_kwargs.get("scoring_functions") or [])
        assess_bias = bool(self.validate_kwargs.get("assess_bias"))

        if metric.startswith("bias_"):
            if not assess_bias:
                raise ValueError(
                    f"Metric '{metric}' requires --assess_bias/--no-assess_bias to be enabled."
                )
            return

        required_group = None
        if metric == "ifp_distance":
            required_group = {"ifp_distance"}
        elif metric in {"sasa", "sasa_norm_heavy"}:
            required_group = {"sasa", "sasa_normalized"}
        elif metric in {"affinity_pred_value", "affinity_probability_binary", "pIC50", "IC50_M", "pIC50_kcal_per_mol"}:
            required_group = {"affinity_metrics", "affinity_metrics_ext"}
        elif metric.startswith("pair_chains_iptm_") or metric in {"chains_ptm", "ptm", "iptm", "confidence_score"}:
            required_group = {"confidence_metrics"}

        if required_group is None:
            # Unknown metric: skip strict check and defer to runtime extraction.
            return
        if not (required_group & scoring):
            raise ValueError(
                f"Metric '{metric}' requires one of scoring functions {sorted(required_group)}, "
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

        validator = Validate(
            wrk_dir=str(run_dir),
            system_path=str(row_system_path),
            options_path=str(self.options_path),
            runner=self.runner,
            **self.validate_kwargs,
        )
        validator.run()

        value = self._extract_single_value(run_dir=run_dir)
        out = pd.DataFrame(
            [
                {
                    "output_metric": self.output_metric,
                    "aggregate": self.aggregate,
                    "value": value,
                }
            ]
        )
        out.to_csv(self.wrk_dir / "oracle_result.csv", index=False)
        return float(value)

    def _extract_single_value(self, run_dir: Path) -> float:
        metric = self.output_metric
        values: list[float] = []
        for csv_path in (run_dir / "results" / "system_metrics.csv", run_dir / "results" / "chain_metrics.csv"):
            if not csv_path.exists():
                continue
            df = pd.read_csv(csv_path)
            if metric not in df.columns:
                continue
            series = pd.to_numeric(df[metric], errors="coerce").dropna()
            values.extend([float(x) for x in series.tolist()])

        if not values:
            raise ValueError(
                f"Metric '{metric}' was not found with numeric values in oracle outputs. "
                "Ensure the metric exists and corresponding scoring function is enabled."
            )

        if self.aggregate == "first":
            return values[0]
        if self.aggregate == "mean":
            return float(sum(values) / len(values))
        if self.aggregate == "max":
            return float(max(values))
        if self.aggregate == "min":
            return float(min(values))
        if self.aggregate == "median":
            sorted_vals = sorted(values)
            n = len(sorted_vals)
            mid = n // 2
            if n % 2 == 1:
                return float(sorted_vals[mid])
            return float((sorted_vals[mid - 1] + sorted_vals[mid]) / 2.0)

        raise ValueError(f"Unsupported aggregate: {self.aggregate}")

    @staticmethod
    def _find_first_ligand_smiles_path(system_dict: dict[str, Any]) -> list[Any] | None:
        sequences = system_dict.get("sequences", [])
        if not isinstance(sequences, list):
            return None
        for i, entry in enumerate(sequences):
            if not isinstance(entry, dict):
                continue
            lig = entry.get("ligand")
            if isinstance(lig, dict) and "smiles" in lig:
                return ["sequences", i, "ligand", "smiles"]
        return None

    @staticmethod
    def _smiles_from_mol_file(path: Path | None) -> str | None:
        if path is None:
            return None
        try:
            from rdkit import Chem
        except Exception as exc:  # pragma: no cover - env-dependent
            raise ValueError("RDKit is required to parse --input_mol_file.") from exc

        suffix = path.suffix.lower()
        mol = None
        if suffix == ".sdf":
            supplier = Chem.SDMolSupplier(str(path), removeHs=False)
            mol = next((m for m in supplier if m is not None), None)
        elif suffix in {".mol", ".mol2"}:
            mol = Chem.MolFromMolFile(str(path), removeHs=False)
        else:
            # Try SDF loader first, then MOL loader.
            supplier = Chem.SDMolSupplier(str(path), removeHs=False)
            mol = next((m for m in supplier if m is not None), None)
            if mol is None:
                mol = Chem.MolFromMolFile(str(path), removeHs=False)

        if mol is None:
            return None
        return Chem.MolToSmiles(mol)
