# Bias recipe: Python API and file flow

This reference reconstructs the standalone `Bias` recipe and its reusable
`BiasAssessmentWorkflow` from `src/cofolder/recipes/bias.py`.

Bias is different from the other recipes: it does **not** select a runner or perform
structure prediction. It assesses how closely query proteins and ligands resemble
configured training/reference data. Validate can call the same assessment workflow
after prediction to add bias metrics to its result frames.

## Smallest Python call

```python
from cofolder.recipes.bias import Bias

system_metrics, chain_metrics = Bias(
    wrk_dir="runs/bias",
    system_path="system.yaml",
    protein_training_data_path="references/proteins.csv",
    ligand_training_data_path="references/ligands.csv",
    bias_release_cutoff="2023-06-01",
).run()
```

`run()` returns two `pandas.DataFrame` objects and writes a versioned public bundle
under `<wrk_dir>/results/`.

## The journey through the code and files

1. **Construction loads the query system.**

   `Bias.__init__()` stores public and optional custom reference paths, date and
   ligand-similarity policy, selected chains, and training-data build settings. It
   reads the system YAML through `modules.utils.read.read_yaml()` and wraps the result
   in `modules.input.system.System`.

2. **The standalone workspace is prepared.**

   `_run_impl()` creates `<wrk_dir>/results/bias_train/`. This directory holds
   supporting bias artifacts. It removes two obsolete nested metric CSV names so a
   rerun cannot mix old and current output layouts.

3. **The molecular system becomes minimal metric frames.**

   `build_bias_dataframes()` creates one system row and one chain row per explicitly
   identified protein or ligand chain. Each chain receives `CHAIN_ID`, `ENTITY_ID`,
   `ENTITY_TYPE`, and a ligand/protein molecule identity. These frames provide the
   same shape that the reusable workflow receives from Validate, without pretending
   that a prediction occurred.

4. **The reusable assessment workflow resolves its sources.**

   `BiasAssessmentWorkflow.apply()` validates the ISO release cutoff and restricts
   source requirements to `bias_chains` when configured. Protein queries require a
   public protein-training file or custom protein reference; ligand queries require
   the corresponding ligand source. Standalone Bias is strict and raises when a
   needed source is missing. Validate uses non-strict mode and records a warning
   while retaining its original metric frames.

5. **Optional training data is built first.**

   With `build_bias_training_data=True`, `_build_training_data()` resolves
   `components.cif`, decides independently whether protein MMseqs and ligand ECFP
   protocols are needed for the selected chains, and calls
   `run_build_bias_training_data()`. Invalid selected ligand SMILES disable only the
   ligand ECFP portion. Generated supporting data is placed under `bias_train` or at
   the explicitly requested paths.

6. **The analytics implementation adds bias metrics.**

   `modules.analytics.bias.apply_bias_metrics()` receives the two DataFrames, query
   `System`, public/custom reference sources, release cutoff, chain selection, Boltz
   cache/CCD paths, output directory, and timing/logger context. It performs the
   protein and ligand similarity/reference work and returns enriched system and
   chain DataFrames. Supporting enrichment, dataset, checkpoint, and plot artifacts
   remain under `<wrk_dir>/results/bias_train/`.

7. **The public result bundle is written.**

   The standalone recipe describes configured datasets as evidence, identifies the
   supporting directory as an artifact, converts the enriched frames with
   `bundle_from_frames()`, and calls `write_public_bundle()` for
   `<wrk_dir>/results/`. It then returns the two enriched DataFrames.

8. **Failures remain public.**

   An error is converted into a typed `bias_analytics_failed` record, written to the
   failed bundle, and raised as `WorkflowExecutionError`.

## Reusing bias inside another Python workflow

If you already have canonical system and chain metric frames, instantiate
`BiasAssessmentWorkflow` directly and call `apply(system_df=..., chain_df=...,
boltz_cache_path=...)`. The required chain columns are `CHAIN_ID`, `ENTITY_ID`, and
`ENTITY_TYPE`, plus the execution identity columns used by the public metric
contract. This is the component boundary used internally by Validate.

## Where to replace one component

- **Query-frame construction:** replace `build_bias_dataframes()` only for standalone
  Bias. Validate supplies real normalized frames instead.
- **Reference-data construction:** replace or wrap
  `run_build_bias_training_data()` while retaining its explicit output paths,
  protein/ligand protocol selection, cutoff, and failure behavior.
- **Bias science:** `apply_bias_metrics()` is the stable orchestration seam named by
  the cleanup design. A substitute should accept and return the same frames, preserve
  chain identities and missing-value meanings, and write equivalent provenance and
  supporting artifacts.
- **Reference sources:** public training CSVs and custom reference paths are already
  selectable through the Python constructor; no recipe change is required.
- **Serialization:** preserve the `bias_metrics` request classification, evidence
  records, artifact reference, and public typed failures.

## Reconstruction checklist

Preserve the no-prediction nature of standalone Bias, explicit chain identity,
chain-scoped source requirements, strict standalone versus non-strict Validate
behavior, independent protein and ligand build protocols, public and custom sources,
cutoff and similarity policies, supporting artifact families, returned DataFrames,
and versioned public records.
