# Validate recipe: Python API and file flow

For the database-backed `assess_bias` inputs, see
[Providing bias training data](../../user-guide/bias-training-data.md).
These inputs include `whole` snapshot selection, the shared query cache, and a
prepared custom-complex supplement bundle.

This is a code-oriented reference for reconstructing the current `Validate` recipe
without going through the CLI. It describes the implementation as it exists in
`src/cofolder/recipes/validate.py`; it is not a promise that private helper names
will remain stable.

## What Validate does

`Validate` is COFOLDER's foundational prediction workflow. It loads one molecular
system and one runner-options document, executes the selected backend for one or
more repeats, normalizes all successful backend outputs, applies shared analytics,
and writes the versioned public result bundle.

`Screen` and `Oracle` both delegate their actual prediction work to `Validate`.
Consequently, a replacement runner or shared analysis added at the Validate
boundary is normally inherited by those recipes as well.

## Smallest Python call

```python
from cofolder.recipes.validate import Validate

workflow = Validate(
    wrk_dir="runs/example",
    system_path="system.yaml",
    options_path="options.yaml",
    runner="boltz2",
    repeats=1,
    seed=1234,
)
workflow.run()
```

`run()` communicates through files and normally returns `None`. Read the canonical
outputs from `<wrk_dir>/results/`; do not expect the DataFrames to be returned.

## The journey through the code and files

1. **Construction selects the runner.**

   `Validate.__init__()` stores all paths and requested analyses, assigns a run ID,
   resolves the runner with `cofolder.modules.runners.get_runner()`, and creates the
   initial public `OutputIdentity`. Runner discovery scans modules in
   `cofolder.modules.runners` for a module-level `RUNNER` object.

2. **The working area is created.**

   `run()` enters `_run_impl()`. The recipe creates `<wrk_dir>/raw/`. Raw and
   backend-private files live here; public files are written later under
   `<wrk_dir>/results/`.

3. **Repeats and seeds become an explicit plan.**

   Without an injected `RunnerExecutionPlan`, `helpers.get_seeds()` resolves the
   requested or generated base seed and one derived seed per repeat. The recipe
   records these as `SeedPlan` and `RepeatSeedProvenance` objects. A runner may
   adjust a backend-facing seed through `resolve_effective_seed()`, but it must
   preserve the original provenance and state why it changed the seed.

4. **Runner options are loaded by the runner.**

   `runner.load_options(options_path)` parses and validates the options YAML using
   that runner's schema. The supported current shape is a typed `RunnerOptions` or
   a runner-owned wrapper around it. This is the first backend-specific seam.

5. **The system YAML is loaded and relative MSA paths are resolved.**

   `load_yaml_document()` reads the YAML and retains source locations for useful
   validation errors. The mapping is wrapped in `modules.input.system.System`.
   `resolve_declared_msa_paths()` resolves declared MSA files relative to the input
   YAML's directory. A deep copy is used for the run so the source document is not
   modified.

6. **The unprepared system is checked against runner capabilities.**

   `runner.validate_system(..., check_atom_names=False)` applies shared input
   validation using the runner's declared entity and constraint capabilities. This
   catches unsupported DNA, RNA, ligand, bond, pocket, or contact inputs before
   backend preparation.

7. **Backend availability and provenance are established.**

   `runner.ensure_available()` checks that the backend can run.
   `runner.detect_backend_identity()` records the runner name, backend distribution,
   and installed version status in the public identity.

8. **The runner prepares the system.**

   `runner.prepare_system()` may convert conformers, populate backend-specific
   chemistry data, or alter backend options. It returns `RunnerPreparationResult`:
   the prepared system, prepared options, warnings, and typed `RunnerRuntime`
   metadata such as cache path, sample count, and model name. Simple runners can use
   the no-op implementation inherited from `BaseRunner`.

9. **The prepared system is validated and written.**

   The recipe validates again with `check_atom_names=True`, then writes the effective
   system to `<wrk_dir>/raw/<original-system-filename>`. This generated YAML—not the
   source YAML—is passed to the backend.

10. **One typed runner request is executed per repeat.**

    For each seed, the recipe creates `<wrk_dir>/raw/repeat_<n>/` and a
    `RunnerExecutionRequest`. The request contains the generated YAML, prepared
    options, repeat and seed provenance, backend identity, chain-identity mapping,
    runtime metadata, logger, and timing collector. `runner.run(request)` owns the
    backend invocation and normalization for that repeat.

    A runner writes its private normalized boundary under
    `raw/repeat_<n>/normalized/` and returns `RunnerExecutionResult`. If supported,
    generated protein MSAs are captured after a run and injected into later repeats.
    In a Screen parent workflow, the same cache is also shared between compounds.

11. **Every runner result is validated before analytics.**

    `validate_runner_bundle()` checks the normalized directory layout, canonical CSV
    columns, structures, sample records, chain identities, backend/seed provenance,
    companion artifacts, typed records, and explicit metric-group outcomes. A bad
    repeat becomes a typed output-validation failure. Other valid repeats may still
    continue; the recipe stops if none remain usable.

12. **Valid repeat outputs are gathered.**

    `modules.utils.gather.gather_structures()` copies normalized structures into
    `<wrk_dir>/results/structures/`. `merge_runner_results()` merges repeat-level
    `system_metrics.csv` and `chain_metrics.csv` into DataFrames. `add_chain_info()`
    attaches semantic chain and entity identities. Unsupported optional runner
    metric groups receive explicit empty columns rather than invented values.

13. **Shared analytics run on the normalized DataFrames and structures.**

    Depending on configuration, `Structure` adds distance or ProLIF interaction
    fingerprints and SASA metrics. `scaffold_reproduction_metrics()` handles the
    optional reference structure and custom-pocket comparisons. Optional bias uses
    `BiasAssessmentWorkflow.apply()`. Robustness aggregation runs when more than one
    repeat or diffusion sample exists. An analytics failure is recorded without
    necessarily discarding otherwise usable backend results.

14. **The public bundle is written.**

    `_write_public_results()` passes the final DataFrames, evidence, failures,
    artifacts, backend identity, and seed plan to `bundle_from_frames()`, followed by
    `write_public_bundle()`. The canonical public location is
    `<wrk_dir>/results/`, including the manifest and `records.jsonl`, derived metric
    tables, and any structure, alignment, matrix, IFP, or bias directories that were
    produced.

15. **Failures remain machine-readable.**

    A failure before a usable result is available is converted with
    `failure_from_exception()`, written as a failed public bundle, and raised as
    `WorkflowExecutionError`. Inspect `exc.failures` and `exc.output_dir` when using
    the API programmatically.

## Where to replace one component

- **Prediction backend:** add or select another runner. This is the designed stable
  extension seam and does not require rewriting `Validate`.
- **Options interpretation or backend command:** implement it inside the runner's
  `load_options()` and `run()` boundary.
- **Ligand/backend preparation:** implement `runner.prepare_system()` and declare
  `ligand_preparation_capabilities`.
- **MSA reuse:** implement the three optional runner MSA methods. The recipe owns
  when they are called.
- **Bundle validation:** `validate_runner_bundle()` is the boundary shared by all
  runners. Replacing it changes the contract, not merely one backend.
- **Structure, reproduction, bias, or robustness analytics:** these are currently
  direct calls in `Validate`, not injected plugins. To substitute only one while
  preserving behavior, keep its DataFrame inputs/outputs and typed failure handling
  equivalent. The proposed recipe refactor intends to make these boundaries smaller.
- **Public serialization:** preserve `PublicOutputBundle`, identities, provenance,
  metric states, and `write_public_bundle()` if downstream consumers must continue
  to read the same outputs.

## Reconstruction checklist

A Python-only reconstruction is behaviorally equivalent when it preserves: typed
option and system validation; runner capability checks; seed and backend provenance;
preparation before final atom-name validation; one request per repeat; immediate
normalized-bundle validation; partial-repeat behavior; chain identity; optional MSA
reuse; requested-versus-unavailable metric states; analytics failure records; and
the versioned public result bundle.
