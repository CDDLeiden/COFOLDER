# Screen recipe: Python API and file flow

For database-backed bias inputs, see
[Providing bias training data](../../user-guide/bias-training-data.md).
The cache path is forwarded to every child Validate run so invariant protein
searches are reused across compounds.

This reference reconstructs the current `Screen` recipe in
`src/cofolder/recipes/screen.py`. `Screen` is a library-expansion and aggregation
layer around `Validate`; it is not a separate prediction implementation.

## Smallest Python call

```python
from cofolder.recipes.screen import Screen

results = Screen(
    wrk_dir="runs/screen",
    system_path="system.yaml",
    options_path="options.yaml",
    runner="boltz2",
    library="compounds.csv",
    library_format="csv",
    ligand_chain="B",
    smiles_column="smiles",
    col_id="compound_id",
    repeats=1,
    seed=1234,
).run()
```

`run()` returns the detailed screening `pandas.DataFrame`. The versioned files under
`<wrk_dir>/results/` remain the canonical persisted result.

For a multi-parameter CSV, pass mappings instead of `smiles_column`:

```python
results = Screen(
    wrk_dir="runs/variants",
    system_path="system.yaml",
    options_path="options.yaml",
    library="variants.csv",
    col_id="experiment",
    mappings=[
        "protein_sequence=sequences.0.protein.sequence",
        "ligand_smiles=sequences.1.ligand.smiles",
    ],
).run()
```

Each mapping targets an existing template field, and every CSV row supplies one
complete combination. The row system is validated before it reaches `Validate`.
If a protein sequence varies while the template supplies an MSA, the row must also
map a matching `msa` value. Screen validates all supplied alignment queries before
starting any child workflow. Without supplied alignments, Boltz-family runners cache
and reuse generated MSAs by normalized sequence and generation settings.

## The journey through the code and files

1. **Construction separates Screen concerns from Validate concerns.**

   Screen stores library parsing, ligand selection, optional IFP filtering and
   clustering settings. Prediction and shared scoring arguments are collected in
   `validate_kwargs`. It resolves the runner immediately with `get_runner()`. For a
   runner with MSA reuse, it reserves `<wrk_dir>/shared/msa/<runner>/`.

2. **Screen-specific configuration is validated before prediction.**

   `_validate_config()` requires a library, infers CSV/SDF/MOL from
   the suffix when necessary, checks either ID/SMILES inputs or general CSV mappings, and
   checks IFP filter and cluster thresholds. Reference-complex or custom-pocket IFP
   policy is parsed before any backend call.

3. **The base system and options are loaded once.**

   `load_yaml_document()` reads the base system; declared MSA paths are resolved
   relative to it. `runner.load_options()` loads typed backend options. Screen then
   constructs one immutable `RunnerExecutionPlan` containing backend identity, seed
   plan, and the Cartesian product of repeats and runner-declared model/sample slots.
   This same plan is passed to every child Validate run, so all compounds are
   compared on the same execution axes.

4. **The base system is validated and any ligand target is resolved.**

   `runner.validate_system()` applies the shared requirement that the system contain
   a protein and, in single-ligand mode, a ligand. `resolve_ligand_target()` locates
   the ligand entity named by `ligand_chain`. Mapped mode can run protein-only
   systems unless ligand-specific analysis is requested.

5. **The input library becomes canonical members or failures.**

   `load_compound_library()` reads CSV, SDF, or MOL input and returns ordered
   outcomes. Each outcome retains source identity, a stable execution ID, metadata,
   normalized ligand data, coordinate mode, and a safe execution-directory name.
   Duplicate-ID policy is applied here. Invalid members remain explicit failures so
   the rest of the library can continue.

   In mapped mode, `load_mapped_system_library()` parses raw CSV text, applies every
   mapped value to a fresh template copy, and retains the same source identity and
   execution-directory conventions. Completed row systems are runner-validated
   before prediction.

6. **Each compound gets its own child workspace.**

   For every outcome, Screen creates
   `<wrk_dir>/<outcome.execution_directory>/` and plans
   `screen_system.yaml`. A valid member replaces only the selected ligand in a copy
   of the base `System` using `replace_ligand_smiles()`.

   An SDF/MOL member may preserve its molblock as `source_ligand.sdf`. Based on the
   selected runner's conformer capabilities, Screen either passes those coordinates,
   asks Validate to generate coordinates, or falls back to native SMILES handling.

7. **Screen delegates prediction to Validate.**

   The per-compound system is written to `screen_system.yaml`, then Screen constructs
   a `Validate` object with the common options path, runner name, shared execution
   plan, shared MSA directory, and `validate_kwargs`. `validator.run()` performs all
   backend preparation, repeat execution, normalized-bundle validation, analytics,
   and child public serialization described in the Validate reference.

8. **Failures are isolated by compound and execution slot.**

   A bad library member or failed child workflow is converted into Screen failure
   information without immediately stopping other compounds. Child Validate failure
   records are rebased onto Screen identities. Each planned repeat/model/sample slot
   becomes `success`, `failed`, or `unavailable`; this prevents a backend from
   silently returning fewer samples than were planned.

9. **Child metrics are read through the public contract.**

   `recipes._metrics.read_metric_frames()` prefers each child
   `<run_dir>/results/records.jsonl` and reconstructs system and chain DataFrames.
   It falls back to older CSV files only for compatibility. Screen selects and
   qualifies metrics as `system__<metric>` or
   `<entity-type>_<chain>__<metric>` so identically named chain metrics remain
   distinguishable.

10. **Optional IFP decisions are applied.**

    Each successful prediction can be compared with a reference complex or custom
    pocket. Screen records similarity, required-interaction, pass/fail, and
    not-evaluable states. Optional clustering groups compatible interaction
    fingerprints after all rows have been processed and writes the cluster summary,
    native SciPy linkage matrix, and deterministic leaf order. ProLIF filtering or
    clustering also writes occurrence-level atom identities, roles, and geometry to
    `results/ifp_interaction_events.jsonl`.

11. **Screen-level outputs are assembled.**

    The source-to-execution mapping is written atomically to
    `<wrk_dir>/results/compound_members.csv`. `_write_public_results()` creates one
    typed `ExecutionRecord` per compound/repeat/model/sample, rebases child metric and
    failure records to Screen identities, adds Screen-derived metrics, and writes the
    public bundle. The manifest points to `executions.csv`,
    `compound_members.csv`, mapped-mode `system_mappings.json`, and the optional IFP
    cluster and event artifacts.

12. **At least one success is required.**

    Partial screens return their DataFrame and publish both successes and failures.
    If no compound produces a usable result, `run()` raises
    `WorkflowExecutionError` after writing the failed public output.

## Where to replace one component

- **Backend:** use another discovered runner name. Validate and Screen share the
  runner contract.
- **Library reader or duplicate policy:** preserve the `CompoundMember`, source
  identity, failure, execution-ID, and execution-directory semantics returned by
  `load_compound_library()`.
- **Ligand substitution:** preserve `resolve_ligand_target()` and
  `replace_ligand_smiles()` behavior so only the intended entity changes.
- **Per-compound prediction:** `Validate` is the deliberate seam. A replacement must
  still produce the same public metric frames, execution cardinality, provenance,
  and typed failures expected by Screen.
- **Metric selection:** preserve the qualified naming scheme used by
  `recipes._metrics`; Oracle uses the same convention.
- **IFP filtering or clustering:** these are currently internal methods rather than
  injected services. A substitute must retain explicit pass, fail, and
  not-evaluable states and identity-bearing fingerprints.
- **Aggregation/serialization:** preserve one terminal execution record for every
  planned slot and the stable member mapping.

## Reconstruction checklist

Preserve early base-system validation, one shared execution plan, stable compound
identity, isolated child workspaces, shared MSA reuse, child Validate output
contracts, partial-failure continuation, exact planned cardinality, qualified metric
names, optional IFP decision states, and Screen-level public serialization.
