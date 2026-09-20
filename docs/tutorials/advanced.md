# Advanced Features

This guide uses only supported COFOLDER 1.0 interfaces. Start with the
[Basic Usage](basic.md) and [Virtual Screening](screening.md) tutorials first.

## Multi-Chain Systems

Protein, DNA, RNA, and ligand entities use the same versioned system schema:

```yaml
version: 1
sequences:
  - protein:
      id: A
      sequence: ACDEFGHIK
  - protein:
      id: B
      sequence: LMNPQRSTV
  - ligand:
      id: L
      smiles: CCO
```

Use a list of IDs only when multiple chains share one entity. Separate ligand entries
represent separate chemistry and require `--ligand_chain` for Screen or Oracle when
the replacement target would otherwise be ambiguous.

## Repeats and Sampling

Use `--repeats` for workflow-level repeats and configure backend diffusion samples in
the options YAML:

```bash
cofolder validate -s examples/system.yaml -o examples/options.yaml \
  --repeats 3 --seed 2026 -w ./repeat-output
```

```yaml
version: 1
runtime:
  cache_path: ./cache/.boltz
  diffusion_samples: 5
runner:
  recycling_steps: 3
```

The manifest records requested, resolved, derived, and effective seeds. Metric rows
retain repeat, model, and sample identity.

## Covalent Binding

Use a CCD ligand when atom names are needed by a bond constraint:

```yaml
version: 1
sequences:
  - protein:
      id: A
      sequence: ACDEFGHIK
  - ligand:
      id: B
      ccd: COV
constraints:
  - bond:
      atom1: [A, 2, SG]
      atom2: [B, 1, C12]
```

Residue positions are one-based. Both atom names must exist in the selected protein
residue or preprocessed ligand CCD. See [Ligand Handling](ligands.md) for custom CCDs.

## Supported Python Workflows

The recipe classes mirror the CLI and write the same public records:

```python
from cofolder.recipes.validate import Validate

Validate(
    wrk_dir="validate-output",
    system_path="examples/system.yaml",
    options_path="examples/options.yaml",
    runner="boltz2",
    repeats=2,
    seed=2026,
).run()
```

```python
from cofolder.recipes.screen import Screen

screen_results = Screen(
    wrk_dir="screen-output",
    system_path="examples/system_screen.yaml",
    options_path="examples/options.yaml",
    library="examples/ligand_screen.csv",
    ligand_chain="B",
    col_id="Name",
    smiles_column="SMILES",
    merge_data="pIC50",
).run()
```

`Screen.run()` returns a convenience `DataFrame`; the files under `results/` retain
the authoritative versioned record contract.

## Custom Oracle Scores and Gates

```python
from cofolder.recipes.oracle import Oracle, OracleGate, OracleGatePolicy


def custom_score(context):
    affinity = context.aggregated_metrics["ligand_B__pIC50"]
    confidence = context.aggregated_metrics["system__confidence_score"]
    return affinity * confidence


score = Oracle(
    wrk_dir="oracle-output",
    system_path="examples/system.yaml",
    options_path="examples/options.yaml",
    input_smiles="CCO",
    aggregate="mean",
    scoring_functions=["confidence_metrics", "affinity_metrics_ext", "sasa_normalized"],
    scoring_function=custom_score,
    score_gates=[OracleGate("ligand_B__sasa_norm_heavy", "le", 2.0)],
    gate_policy=OracleGatePolicy("downweight", 0.25),
).run()
```

The threshold is illustrative and must be calibrated for the target system. Supported
reducers are `first`, `mean`, `max`, `min`, and `median`; gate policies are
`downweight`, `fixed_penalty`, and `non_binder`.

## Analyze Public Metrics

Read the long-form output directly instead of relying on runner-owned files:

```python
import pandas as pd

metrics = pd.read_csv("screen-output/results/metrics.csv")
computed = metrics[metrics["status"] == "computed"]

confidence = computed[computed["metric_name"] == "confidence_score"]
affinity = computed[computed["metric_name"] == "affinity_pred_value"]
distance_ifp = computed[computed["metric_name"] == "ifp_distance"]
```

For reference-free contact-pattern grouping, enable `--cluster_ifps` on Screen and
read `results/ifp_cluster_summary.csv`. For reference-overlap decisions, configure the
documented IFP filter options instead of calling analytics implementation details.

## Parallel Jobs and Recovery

COFOLDER 1.0 does not expose an in-process multi-GPU scheduler or checkpoint API.
For independent library shards, create each CSV explicitly and submit one supported
Screen command per GPU through the site scheduler:

```bash
CUDA_VISIBLE_DEVICES=0 cofolder screen -s examples/system_screen.yaml \
  -o examples/options.yaml -c library-part-1.csv --col_id Name \
  --smiles_column SMILES --ligand_chain B -w ./screen-part-1

CUDA_VISIBLE_DEVICES=1 cofolder screen -s examples/system_screen.yaml \
  -o examples/options.yaml -c library-part-2.csv --col_id Name \
  --smiles_column SMILES --ligand_chain B -w ./screen-part-2
```

Each work directory has its own manifest, execution records, and failures. Recovery
means rerunning the failed source members in a new, auditable Screen invocation; no
public checkpoint-resume interface is claimed.

## Post-v1 Scope

Protein-iteration/selectivity screening, ligand soaking, and a standalone molecule-
generation workflow are not implemented COFOLDER 1.0 capabilities. Molecule generators
may integrate externally through the Oracle interface, but this documentation does not
present those three ambitions as available release functionality.

## Related

- [Oracle Guide](../user-guide/oracle.md)
- [Public Output Contract](../user-guide/public-output-contract.md)
- [Configuration Guide](../getting-started/configuration.md)
