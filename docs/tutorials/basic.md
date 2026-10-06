# Basic Usage Tutorial

This tutorial walks through basic COFOLDER usage for protein-ligand co-folding.

## Prerequisites

- COFOLDER installed
- CUDA-compatible GPU
- Basic understanding of protein-ligand systems

## Tutorial Overview

We'll predict the structure of ibuprofen bound to a protein.

## Step 1: Prepare System File

Create `system.yaml`:

```yaml
version: 1
sequences:
  - protein:
      id: "protein_1"
      fasta: "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL"
  - ligand:
      id: B
      smiles: "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
```

## Step 2: Configure Runner Options

Create `options.yaml`:

```yaml
version: 1
runtime:
  cache_path: ./cache/.boltz
  diffusion_samples: 1
runner:
  devices: 1
  recycling_steps: 3
```

## Step 3: Run Prediction

Check the complete input contract without starting inference:

```bash
cofolder validate -s system.yaml -o options.yaml \
  -w ./tutorial_output --preflight_only
```

After preflight reports `preflight=ready`, run the prediction:

```bash
cofolder validate -s system.yaml -o options.yaml -w ./tutorial_output
```

## Step 4: Examine Output

Check the output directory:

```bash
ls -l tutorial_output/
```

You should see:
- `raw/` - runner-owned raw execution artifacts
- `results/records.jsonl` - authoritative versioned records
- `results/{executions,successes,metrics,failures}.csv` - long-form views
- `results/manifest.json` - invocation, backend, evidence, and completion metadata
- `results/structures/` - gathered output structures
- log files in the working directory

## Step 5: Visualize Results

Load one of the gathered structure files in PyMOL, ChimeraX, or your favorite molecular viewer:

```bash
pymol tutorial_output/results/structures/*.cif
```

## Step 6: Check Confidence

Examine the confidence scores:

```python
import pandas as pd

metrics_df = pd.read_csv("tutorial_output/results/metrics.csv")
confidence = metrics_df[
    (metrics_df["metric_name"] == "confidence_score")
    & (metrics_df["status"] == "computed")
]

print(confidence[["repeat_id", "model_id", "sample_id", "value"]].head())
```

`metrics.csv` is deliberately long-form: metric names live in `metric_name` and
their numeric result lives in `value`. Do not expect the retired wide
`system_metrics.csv` or `chain_metrics.csv` files in `results/`.

## Next Steps

### Try Different Options

Experiment with different runner parameters:

```yaml
version: 1
runtime:
  cache_path: ./cache/.boltz
  diffusion_samples: 5
runner:
  recycling_steps: 5
```

### Add 3D Conformer Generation

```bash
cofolder validate \
  -s system.yaml \
  -o options.yaml \
  --conformers 3D
```

### Enable Debug Logging

```bash
cofolder validate -s system.yaml -o options.yaml -d
```

## Troubleshooting

### Common Issues

**GPU Out of Memory:**
- Reduce `diffusion_samples`
- Use smaller protein sequences
- Check GPU memory with `nvidia-smi`

**Invalid SMILES:**
- Validate SMILES with RDKit
- Check for special characters
- Use canonical SMILES

**CCD Errors:**
- Keep CCD identifiers ≤ 5 characters
- Avoid special characters
- Use unique identifiers

## Further Reading

- [Validate Command Guide](../user-guide/validate.md)
- [Configuration Reference](../getting-started/configuration.md)
- [Virtual Screening Tutorial](screening.md)
