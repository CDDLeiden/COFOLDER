# Adding New Runners

This tutorial shows how to add a new prediction backend to COFOLDER without changing the shared `validate`, `screen`, or `oracle` workflows.

This is the authoritative contributor-facing guide for runner authoring. Other docs may summarize the runner system or preserve compatibility notes, but new runners should follow this document's contract and examples.

The key idea is:

- COFOLDER handles workflow orchestration, input mutation, repeats, result gathering, and analytics.
- A runner handles backend-specific dependency checks, option loading, optional backend-specific system preparation, execution, and normalization into COFOLDER's canonical bundle.

## Where Runners Live

Runner modules are discovered from `src/cofolder/modules/runners/`.

Each runner file must expose a module-level `RUNNER` object. Discovery is handled by `cofolder.modules.runners.discover_runners()`, which scans the package and registers `RUNNER.name`.

Examples in the current codebase:

- `src/cofolder/modules/runners/boltz1_runner.py`
- `src/cofolder/modules/runners/boltz_runner.py`
- `src/cofolder/modules/runners/boltz2_runner.py`
- `src/cofolder/modules/runners/boltz_community_runner.py`

For the current contract, `Boltz2Runner` backed by `BoltzRunner` is the in-tree reference implementation. It already exercises the typed runtime surface, canonical normalized bundle layout, explicit `metric_outcomes`, and early boundary validation path used by the shared workflows.

`Boltz1Runner` is also on the supported contract path, but it is intentionally the narrower-capability example: it supports confidence metrics only, and affinity-related metric groups stay explicit `unsupported` outcomes rather than simulated payloads.

`BoltzCommunityRunner` follows the same full-capability supported path as `Boltz2Runner`, but it keeps a distinct coexistence rule: it must fail clearly if a PyPI `boltz` package-line install is present in the same environment.

## Runner Lifecycle

At a high level, COFOLDER calls a runner in four steps:

1. `check_availability()` / `ensure_available()`
2. `load_options(options_path)`
3. `prepare_system(system_obj, options_obj, wrk_dir, conformers, sdf_file, logger)` if the backend needs preparation
4. `run(request)`

The supported contract types live in `src/cofolder/modules/runners/contracts.py`, and `src/cofolder/modules/runners/base.py` provides the default no-op preparation path for simple runners.

Compatibility note:

- Import and author against `RunnerExecutionRequest`, `RunnerExecutionResult`, and `RunnerPreparationResult`.
- Legacy alias names such as `RunnerRequest`, `RunnerResult`, and `RunnerPreparation` may remain importable for compatibility, but they are not the recommended authoring surface for new code or docs.

## Input Provided By COFOLDER

COFOLDER gives the runner two kinds of input:

### 1. Preparation-time input

`prepare_system(...)` receives:

- `system_obj`: parsed system YAML object, already loaded by COFOLDER
- `options_obj`: runner-specific options object returned by `load_options()`
- `wrk_dir`: the shared raw working directory for the run
- `conformers`: optional conformer mode requested by the user
- `sdf_file`: optional SDF path provided by the user
- `logger`: workflow logger

The runner returns a `RunnerPreparationResult` object:

```python
@dataclass(slots=True)
class RunnerRuntime:
    diffusion_samples: int = 1
    cache_path: str | None = None
    model_name: str | None = None


@dataclass(slots=True)
class RunnerPreparationResult:
    system_obj: Any
    options_obj: Any
    warnings: list[str] = field(default_factory=list)
    runtime: RunnerRuntime = field(default_factory=RunnerRuntime)
```

Use this stage for backend-specific preparation, for example:

- forcing a backend-specific option such as `model: boltz2`
- converting ligands into backend cache formats
- enriching the system object before COFOLDER writes the final YAML used for execution
- recording explicit runtime metadata such as cache location or diffusion sample count

If your backend does not need custom preparation, inherit `BaseRunner.prepare_system(...)` as-is. The default implementation is a no-op and already returns the original system, options, and empty typed runtime metadata.

### 2. Execution-time input

`run(...)` receives a `RunnerExecutionRequest`:

```python
@dataclass(slots=True)
class RunnerExecutionRequest:
    runner_name: str
    system_name: str
    system_path: Path
    system_obj: Any
    options_path: Path
    options_obj: Any
    repeat: int
    seed: int
    repeat_dir: Path
    raw_dir: Path
    logger: logging.Logger
    timings: Any | None = None
    label_prefix: str | None = None
    runtime: RunnerRuntime = field(default_factory=RunnerRuntime)
```

Important fields:

- `system_path`: the YAML path COFOLDER wrote after shared input handling and runner preparation
- `repeat_dir`: runner-owned directory for one repeat, typically `raw/repeat_<n>/`
- `seed`: per-repeat seed resolved by COFOLDER
- `timings`: optional timing collector used by the shared debug summary
- `label_prefix`: per-repeat label for timing metrics

## Output Required By COFOLDER

The runner must normalize backend output into the canonical structure below inside `request.repeat_dir / "normalized"`.

Required files:

- `normalized/system_metrics.csv`
- `normalized/chain_metrics.csv`
- `normalized/manifest.json`
- `normalized/structures/` containing copied `.cif`, `.mmcif`, or `.pdb` files

The runner then returns a `RunnerExecutionResult`:

```python
@dataclass(slots=True)
class RunnerExecutionResult:
    runner_name: str
    raw_output_dir: Path
    normalized_dir: Path
    structures_dir: Path
    system_metrics_path: Path
    chain_metrics_path: Path
    manifest_path: Path
    diffusion_samples: int = 1
    capabilities: set[str] = field(default_factory=set)
    warnings: list[str] = field(default_factory=list)
    runtime: RunnerRuntime = field(default_factory=RunnerRuntime)
    sample_records: list[dict[str, Any]] = field(default_factory=list)
```

## Metric Outcome States

The canonical normalized bundle now distinguishes runner metric-group outcomes with four exact states:

- `computed`: the required normalized payload is present and schema-valid
- `unsupported`: the selected runner does not support the requested metric group
- `missing`: the runner declares support, but the required normalized payload is absent
- `failed`: the runner attempted production, but the metric group did not complete successfully

Unsupported requested groups are the only runner-boundary case that should warn and continue. `missing`, `failed`, malformed bundles, and contradictory states stop the workflow after runner execution and before shared gather or analytics proceed.

Shared validation no longer derives requested metric-group outcomes from CSV shape or runner capabilities. If a requested runner metric group matters to workflow behavior, the runner must emit an explicit `metric_outcomes` entry for it.

`Boltz2Runner` is the reference example of the supported authoring path: it emits explicit `metric_outcomes` directly from runner-owned normalization code, and shared validation only verifies that authored state instead of inventing backend-specific logic in recipes.

## Boundary Validation

COFOLDER validates the normalized bundle immediately after `run(...)` returns and before shared gather steps begin.

That validation checks:

- required bundle files and directories exist under `normalized/`
- canonical CSV columns are present
- structure files and `sample_records` agree
- explicit metric outcomes, if supplied, are consistent with runner capabilities and payload shape

If validation fails, the workflow stops after the attempted runner execution, and partial artifacts remain in `raw/repeat_<n>/` for debugging.

After that, COFOLDER takes over again:

- normalized structures are gathered into `results/structures/`
- per-repeat CSVs are merged into `results/system_metrics.csv` and `results/chain_metrics.csv`
- shared analytics such as structure metrics, reproduction metrics, bias, and robustness run on the normalized bundle

The merge step is implemented in `src/cofolder/modules/utils/gather.py`.

`manifest.json` may still include `runtime_context` for debugging and provenance, but shared recipe logic should rely on the typed `RunnerRuntime` contract rather than re-reading loose runtime dictionaries.

## Canonical CSV Expectations

At minimum, the normalized CSVs should include these common columns:

### `system_metrics.csv`

- `cif_file`
- `model_name`
- `repeat`
- `diffusion_sample`

Optional runner-native metrics can then be added, for example:

- `ptm`
- `iptm`
- `confidence_score`

### `chain_metrics.csv`

- `conf_chain_id`
- `cif_file`
- `model_name`
- `repeat`
- `diffusion_sample`

COFOLDER later adds shared chain metadata such as:

- `CHAIN_ID`
- `ENTITY_TYPE`
- `ligand_molecule_id`

If your runner supports affinity outputs, add the normalized chain-level columns COFOLDER already expects:

- `affinity_pred_value`
- `affinity_probability_binary`
- `pIC50`
- `IC50_M`
- `pIC50_kcal_per_mol`

## Capabilities

Each runner declares the metric groups it can natively provide through `capabilities`.

Current shared runner metric groups are:

- `confidence_metrics`
- `affinity_metrics`
- `affinity_metrics_ext`

If a user requests a metric group the selected runner does not support, COFOLDER warns and continues. The corresponding output columns are added as empty values so downstream schemas remain stable.

This means you should declare only the metrics your backend can truly normalize.

## Minimal Runner Skeleton

The following skeleton shows the expected shape:

```python
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from cofolder.modules.runners import (
    BaseRunner,
    RunnerExecutionRequest,
    RunnerExecutionResult,
    RunnerMetricOutcome,
    RunnerRuntime,
)


class MyRunner(BaseRunner):
    name = "my-runner"
    capabilities = {"confidence_metrics"}

    def check_availability(self) -> tuple[bool, str | None]:
        return self.check_distribution_available(
            distribution_name="my-backend",
            missing_message=(
                "The selected 'my-runner' backend is not installed. "
                'Install it with `pip install "cofolder[my-runner]"`.'
            ),
        )

    def load_options(self, options_path: Path):
        # Return any runner-specific options object you want.
        return {"options_path": str(options_path)}

    def run(self, request: RunnerExecutionRequest) -> RunnerExecutionResult:
        raw_output_dir = request.repeat_dir / "my_backend_output"
        normalized_dir = request.repeat_dir / "normalized"
        structures_dir = normalized_dir / "structures"
        structures_dir.mkdir(parents=True, exist_ok=True)
        raw_output_dir.mkdir(parents=True, exist_ok=True)

        # 1. Run your backend here.
        # 2. Copy normalized structure files into normalized/structures/.
        # 3. Build canonical CSVs.
        structure_name = f"{request.repeat}_{request.system_name}_model_0.cif"
        (structures_dir / structure_name).write_text("data_demo", encoding="utf-8")

        system_df = pd.DataFrame([
            {
                "cif_file": structure_name,
                "model_name": request.system_name,
                "repeat": request.repeat,
                "diffusion_sample": 0,
                "ptm": 0.82,
                "iptm": 0.76,
                "confidence_score": 0.91,
            }
        ])
        chain_df = pd.DataFrame([
            {
                "conf_chain_id": 0,
                "cif_file": structure_name,
                "model_name": request.system_name,
                "repeat": request.repeat,
                "diffusion_sample": 0,
                "chains_ptm": 0.88,
            }
        ])

        system_metrics_path = normalized_dir / "system_metrics.csv"
        chain_metrics_path = normalized_dir / "chain_metrics.csv"
        system_df.to_csv(system_metrics_path, index=False)
        chain_df.to_csv(chain_metrics_path, index=False)

        sample_records = [
            {
                "repeat": request.repeat,
                "diffusion_sample": 0,
                "cif_file": structure_name,
            }
        ]

        manifest_path = normalized_dir / "manifest.json"
        manifest = {
            "runner": self.name,
            "capabilities": sorted(self.capabilities),
            "repeat": request.repeat,
            "raw_output_dir": str(raw_output_dir),
            "normalized_dir": str(normalized_dir),
            "system_metrics_path": str(system_metrics_path),
            "chain_metrics_path": str(chain_metrics_path),
            "structures_dir": str(structures_dir),
            "diffusion_samples": 1,
            "runtime_context": {"diffusion_samples": 1},
            "warnings": [],
            "sample_records": sample_records,
        }
        manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

        return RunnerExecutionResult(
            runner_name=self.name,
            raw_output_dir=raw_output_dir,
            normalized_dir=normalized_dir,
            structures_dir=structures_dir,
            system_metrics_path=system_metrics_path,
            chain_metrics_path=chain_metrics_path,
            manifest_path=manifest_path,
            diffusion_samples=1,
            capabilities=set(self.capabilities),
            runtime=RunnerRuntime(diffusion_samples=1),
            sample_records=sample_records,
            metric_outcomes={
                "confidence_metrics": RunnerMetricOutcome(state="computed"),
            },
        )


RUNNER = MyRunner()
```

## Recommended Development Steps

When adding a new runner, this is the safest sequence:

1. Add the runner module under `src/cofolder/modules/runners/`.
2. Implement `check_availability()` first, so CLI selection fails early with a useful install message.
3. Implement `load_options()` and only override `prepare_system()` if you truly need backend-specific mutation.
4. Implement `run()` to write a complete normalized bundle even before all metrics are available.
5. Declare only the metric groups you can truly support in `capabilities`.
6. Add targeted tests for:
   - runner discovery
   - availability checking
   - normalized bundle creation
   - any version- or backend-line-specific behavior

## Testing a New Runner

At minimum, add:

- one discovery test to confirm the runner appears in `list_runner_names()`
- one availability test to confirm the install hint is actionable
- one bundle test to confirm `system_metrics.csv`, `chain_metrics.csv`, `manifest.json`, and normalized structures are produced

The existing runner tests are good templates:

- `tests/modules/runners/test_boltz1_runner.py`
- `tests/modules/runners/test_boltz_runner.py`
- `tests/modules/runners/test_boltz_community_runner.py`
- `tests/modules/runners/test_contracts.py`
- `tests/modules/runners/test_validators.py`
- `tests/modules/runners/test_registry.py`

## Practical Notes

- Prefer normalizing backend-specific output inside the runner, not in `validate`.
- Keep raw backend artifacts inside the runner-owned repeat directory and expose them through `manifest.json`.
- If your backend does not support a metric group, do not fake support. Let COFOLDER handle the warning and empty-column behavior.
- If your backend line has version-specific behavior, encode that directly in `check_availability()`, as the current Boltz runners do.

## Summary

A new COFOLDER runner only needs to do three things well:

1. Check that the correct backend is installed.
2. Execute the backend from the standardized `RunnerExecutionRequest`.
3. Write the normalized bundle that the shared analytics pipeline expects.

If you keep that boundary clean, the existing `validate`, `screen`, `oracle`, gather, and analytics code can stay unchanged.
