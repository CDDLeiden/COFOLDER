# Oracle recipe: Python API and file flow

For database-backed bias inputs, see
[Providing bias training data](../../user-guide/bias-training-data.md).
The cache path is forwarded to candidate Validate runs; protein searches are shared
and ligand entries use canonical molecular identity.

This reference reconstructs the current `Oracle` recipe in
`src/cofolder/recipes/oracle.py`. Oracle turns one candidate ligand and the output of
one child `Validate` workflow into one finite scalar.

## Smallest Python call

```python
from cofolder.recipes.oracle import Oracle

score = Oracle(
    wrk_dir="runs/oracle",
    system_path="system.yaml",
    options_path="options.yaml",
    runner="boltz2",
    input_smiles="CCO",
    output_metric="affinity_pred_value",
    aggregate="mean",
    repeats=1,
    seed=1234,
).run()
```

Exactly one of `input_smiles` and `input_mol_file` is required. Exactly one scoring
mode is also required: `output_metric`, `score_components`, or a Python
`scoring_function` callback. `run()` returns a `float` and also writes a public
bundle under `<wrk_dir>/results/`.

## The journey through the code and files

1. **Construction records the candidate and scoring policy.**

   Oracle stores the base system, options, runner, ligand selector, and the arguments
   forwarded to Validate. It normalizes optional weighted score components. Gates
   are represented by `OracleGate`; their failure transformation is represented by
   `OracleGatePolicy`.

2. **Configuration is validated without invoking the backend.**

   `_validate_config()` checks the one-of input rule, one-of scoring rule, aggregate
   name, callable scoring function, gate policy, and dependencies of every selected
   metric. For example, a bias metric requires `assess_bias=True`, reference metrics
   require a reference, and a vector-valued interaction fingerprint cannot be used
   directly as a scalar.

3. **The candidate becomes validated SMILES.**

   A direct SMILES string is used as-is initially. An SDF, MOL, or MOL2 input is read
   with RDKit and converted to SMILES. Oracle creates `<wrk_dir>/oracle_run/` for the
   child workflow.

4. **The base system is loaded and checked with the selected runner.**

   `load_yaml_document()` reads the system, `resolve_declared_msa_paths()` resolves
   relative MSA inputs, and `get_runner()` obtains the adapter. The runner loads the
   options and validates that the system contains both protein and ligand inputs
   that it can consume.

5. **The ligand entity is selected and replaced.**

   If `ligand_chain` is omitted, Oracle infers it only when exactly one ligand entity
   exists. `resolve_ligand_target()` identifies the entity, `validate_smiles()`
   creates a source-aware normalized ligand, and `replace_ligand_smiles()` changes
   that entity in a copy of the validated system.

6. **The child system is written and Validate is called.**

   Oracle checks backend availability, writes
   `<wrk_dir>/oracle_run/oracle_system.yaml`, and constructs `Validate` with that
   file, the original options path, runner name, repeats, seed, requested analyses,
   and optional evidence. The entire prediction and analysis path is therefore the
   same one used by a direct Validate call.

7. **Metrics are read and qualified.**

   `read_metric_frames()` reads the child's versioned public records.
   `collect_qualified_metric_values()` exposes values with names such as
   `system__confidence_score` and `ligand_B__affinity_pred_value`. A bare metric name
   first prefers the selected ligand, then the system, and succeeds elsewhere only
   when the match is unambiguous.

8. **One base scalar is calculated.**

   - `output_metric` selects one metric and applies `first`, `mean`, `max`, `min`, or
     `median` across its finite values.
   - `score_components` aggregates each selected metric and returns their weighted
     sum.
   - `scoring_function` receives `OracleScoreContext`, containing query SMILES, child
     run directory, copies of both metric DataFrames, and all aggregated qualified
     metrics. It must return a finite numeric scalar.

9. **Optional gates transform the score.**

   Each gate aggregates its selected metric and performs `gt`, `ge`, `lt`, or `le`.
   If a gate fails, the policy can downweight the score, apply a fixed penalty, or
   return a configured non-binder value. There are no universal scientific default
   thresholds; callers choose them.

10. **The scalar and provenance are published.**

    Oracle returns the final float. It also writes a success record,
    `oracle_raw_score`, `oracle_score`, and—when adjusted—
    `oracle_gate_adjusted_score` into `<wrk_dir>/results/`. The candidate's public ID
    is a SHA-256 digest of canonicalized SMILES where possible, avoiding raw ligand
    disclosure in the identity.

11. **A failure is both written and raised.**

    Input, child Validate, unavailable-metric, callback, or gate failures become a
    typed failed Oracle bundle and a `WorkflowExecutionError`.

## Where to replace one component

- **Backend and prediction:** select or add a runner; Oracle delegates prediction to
  Validate.
- **Candidate source conversion:** replace the MOL/SDF-to-SMILES step while retaining
  source-aware ligand validation and unambiguous target selection.
- **Metric source:** implement a custom `scoring_function`. This is already a public
  Python injection seam and is the least invasive way to replace scalar scoring.
- **Metric aggregation:** the five reducers are internal. Replacing one requires
  preserving finite-value filtering and the qualified metric-selection rules.
- **Gating:** construct different `OracleGate` and `OracleGatePolicy` values without
  modifying the recipe. A new policy mode would require a recipe change and public
  contract tests.
- **Child workflow:** a replacement for Validate must provide equivalent public
  metric frames and typed failures; merely writing an arbitrary CSV is insufficient.

## Reconstruction checklist

Preserve exclusive input/scoring modes, runner-aware base-system validation,
unambiguous ligand selection, source-aware ligand replacement, child Validate
execution, qualified metric names, finite scalar enforcement, all five reducers,
custom callbacks, gates and penalties, hashed candidate identity, and public failure
serialization.
