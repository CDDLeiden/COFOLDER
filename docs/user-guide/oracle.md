# Oracle Command and Python API

The Oracle workflow runs `validate` for one query ligand in a fixed system and
returns one finite scalar. The command line exposes standard single-metric calls;
the Python API additionally supports weighted objectives, structure gates, and
arbitrary scoring functions.

## Single-metric command

```bash
cofolder oracle \
  -s examples/system.yaml -o examples/options.yaml \
  --input_smiles "CCO" \
  --output_metric ligand_B__affinity_pred_value \
  --aggregate first \
  -w oracle_affinity
```

Exactly one of `--input_smiles` and `--input_mol_file` is required. `--aggregate`
can be `first`, `mean`, `max`, `min`, or `median` and reduces values across repeats
and diffusion samples.

Metric selectors use the same qualification as Screen: `system__<metric>` or
`<entity>_<chain ID>__<metric>`. Bare names remain supported and prefer the query
ligand, then the system. Qualify a selector whenever a system has multiple matching
chains.

## Ligand training-set similarity

Bias is ligand similarity/proximity to a supplied training reference, not a binding
score or a separately inferred error estimate.

```bash
cofolder oracle \
  -s examples/system.yaml -o examples/options.yaml \
  --input_smiles "CCO" \
  --output_metric ligand_B__bias_lig_sim_train \
  --assess_bias --bias_chains B \
  --ligand_training_data_path \
    src/cofolder/acceptance/data/ligand_training_data.csv \
  -w oracle_bias
```

## Custom pocket/IFP-reference coverage

Raw `ifp_distance` is a vector and cannot be returned as a scalar. Supply a pocket
reference and request its numerical coverage instead:

```bash
cofolder oracle \
  -s examples/system.yaml -o examples/options.yaml \
  --input_smiles "CCO" \
  --scoring_functions ifp_distance \
  --reproduction_metrics pocket_coverage \
  --pocket_coverage_reference "A25 G48 Y51" \
  --output_metric ligand_B__pocket_coverage_custom \
  -w oracle_pocket
```

The pocket residue list is illustrative. Derive target-specific definitions and
thresholds during system validation.

## Composite and structure-gated scoring

Composite objectives and gates are deliberately Python-only. A weighted composite is
the unnormalized sum of each aggregated metric times its weight.

```python
from cofolder.recipes.oracle import Oracle, OracleGate, OracleGatePolicy

oracle = Oracle(
    wrk_dir="oracle_gated",
    system_path="examples/system.yaml",
    options_path="examples/options.yaml",
    input_smiles="CCO",
    aggregate="mean",
    scoring_functions=[
        "confidence_metrics",
        "affinity_metrics_ext",
        "ifp_distance",
        "sasa_normalized",
    ],
    pocket_coverage_reference="A25 G48 Y51",
    reproduction_metrics=["pocket_coverage"],
    score_components={
        "ligand_B__pIC50": 1.0,
        "system__confidence_score": 1.0,
    },
    score_gates=[
        OracleGate("ligand_B__pocket_coverage_custom", "ge", 0.60),
        OracleGate("ligand_B__sasa_norm_heavy", "le", 2.0),
    ],
    gate_policy=OracleGatePolicy("downweight", 0.25),
)
score = oracle.run()
```

All gates must pass. If any metric is missing or fails its comparison, the policy is
applied once. The available policies are:

- `OracleGatePolicy("downweight", factor)` multiplies the base score by a factor in
  `[0, 1]`.
- `OracleGatePolicy("fixed_penalty", value)` returns the supplied penalty.
- `OracleGatePolicy("non_binder", value)` assigns the supplied non-binder score and
  records that interpretation in the audit output.

## Arbitrary Python scoring

A scoring function receives an `OracleScoreContext` with the query SMILES, run path,
raw system and chain DataFrames, and a read-only mapping of aggregated qualified
metrics:

```python
from cofolder.recipes.oracle import Oracle


def custom_score(context):
    affinity = context.aggregated_metrics["ligand_B__pIC50"]
    confidence = context.aggregated_metrics["system__confidence_score"]
    return affinity * confidence


score = Oracle(
    wrk_dir="oracle_custom",
    system_path="examples/system.yaml",
    options_path="examples/options.yaml",
    input_smiles="CCO",
    scoring_functions=["confidence_metrics", "affinity_metrics_ext"],
    scoring_function=custom_score,
).run()
```

Exactly one base-score source is allowed: `output_metric`, `score_components`, or
`scoring_function`.

## Output and audit trail

`Oracle.run()` returns the final `float` and writes the versioned public contract
under `<wrk_dir>/results/`. The base, final, and gate-adjusted values are long-form
Oracle metric records. Raw Validate outputs remain under `<wrk_dir>/oracle_run/results/`.

## Metric prerequisites

- Affinity and confidence selectors require their corresponding scoring groups.
- `sasa` and `sasa_norm_heavy` require SASA scoring.
- `bias_*` requires `assess_bias=True` and appropriate training references.
- Custom pocket coverage requires `ifp_distance`, a pocket reference, and the
  `pocket_coverage` reproduction metric.
- Reference-derived pocket coverage requires `reference_path`.

No confidence, affinity, structural, or bias metric is a universal pass/fail test.
Calibrate objectives and gates against the validated target system.

## Related

- [Bias Command](bias.md)
- [Screen Command](screen.md)
- [Oracle API Reference](../api/recipes/oracle.md)
