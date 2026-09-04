# Virtual Screening Tutorial

Learn how to screen compound libraries with COFOLDER.

## Overview

This tutorial demonstrates high-throughput virtual screening of a ligand library against a protein target.

## Prerequisites

- Completed [Basic Tutorial](basic.md)
- Compound library (CSV or SDF format)
- Target protein structure or sequence

## Scenario

Screen a library of 100 drug-like compounds against a protein kinase.

## Step 1: Prepare Library

Create `library.csv`:

```csv
compound_id,smiles,mw,logp
CMPD001,CC(C)Cc1ccc(cc1)C(C)C(=O)O,206.28,3.50
CMPD002,CC(=O)OC1=CC=CC=C1C(=O)O,180.16,1.19
CMPD003,CN1C=NC2=C1C(=O)N(C(=O)N2C)C,194.19,-0.07
```

## Step 2: Setup System Template

Create `system_template.yaml`:

```yaml
version: 1
sequences:
  - protein:
      id: "kinase"
      fasta: "MENLNMDLLYMAAAVMMGLAAIGAAIGIGILGGKFLEGAARQPDLIPLLRTQFFIVKGGLPFMMVGMAG..."
  - ligand:
      smiles: "PLACEHOLDER"  # Will be replaced
      ccd: "LIG"
```

## Step 3: Configure Screening

Create `screening_options.yaml`:

```yaml
out_dir: screening_results
devices: [0]
num_models: 1
recycling_steps: 3
diffusion_samples: 1
```

## Step 4: Run Screening

```bash
cofolder screen \
  -s system_template.yaml \
  -o screening_options.yaml \
  -c library.csv \
  --col_id compound_id \
  --variable sequences,1,ligand,smiles --col_variable smiles \
  --merge_data "mw,logp" \
  -w ./screening_output
```

## Step 5: Monitor Progress

The screen command will process each compound sequentially. Monitor with:

```bash
tail -f screening_output/log.log
```

For Boltz-family runners, the first ligand also resolves any missing protein MSA.
COFOLDER stages it in `screening_output/shared/msa/<runner>/` and injects that path into all
later ligand systems and repeats. If `system_template.yaml` already contains an `msa`
path, COFOLDER preserves it and skips generation. Relative paths are interpreted from
the directory containing the original system YAML.

## Step 6: Analyze Results

Results are saved in a CSV file:

```python
import pandas as pd

results = pd.read_csv("screening_output/screen_results_with_scores.csv")

# Sort by a system-level confidence metric
top_hits = results.sort_values("system__confidence_score", ascending=False).head(10)

print("Top 10 Compounds:")
print(top_hits[["compound_id", "system__confidence_score", "mw", "logp"]])
```

The same table is returned directly when calling `Screen.run()` from Python. Its
stable columns include model affinity/pIC50 and binding likelihood when supported,
confidence scores, ligand SASA, distance IFPs, bias similarities, and configured
pocket/reference metrics. Metrics that are unavailable for a runner remain empty.

### Discover contact-pattern families without a reference

Add clustering when no experimental pose or predefined pocket is available:

```bash
cofolder screen \
  -s system_template.yaml -o screening_options.yaml -c library.csv \
  --col_id compound_id \
  --variable sequences,1,ligand,smiles --col_variable smiles \
  --scoring_functions confidence_metrics affinity_metrics affinity_metrics_ext \
    sasa sasa_normalized ifp_distance \
  --cluster_ifps \
  --ifp_cluster_similarity_threshold 0.5 \
  -w ./screening_output
```

Read `ifp_cluster_summary.csv` to inspect cluster sizes, member IDs, medoids,
consensus fingerprints, and within-cluster Jaccard similarity. The row-level
`ifp_cluster_id` is also present in both consolidated screening CSVs. This workflow
does not need `--reference_path` or `--pocket_coverage_reference`.

Reference-overlap filtering is a separate, non-destructive decision layer. It marks
every result as accepted, rejected, not evaluable, or not applied and never removes
predictions. You may enable filtering and clustering together when both reference
agreement and reference-free contact-pattern diversity matter.

## Step 7: Visualize Top Hits

```python
import matplotlib.pyplot as plt
import seaborn as sns

# Plot confidence distribution
plt.figure(figsize=(10, 6))
sns.histplot(results["system__confidence_score"].dropna(), bins=30)
plt.xlabel('Confidence Score')
plt.ylabel('Count')
plt.title('Screening Confidence Distribution')
plt.savefig('confidence_dist.png')
```

## Advanced Options

### Multi-GPU Screening

Split library and run parallel jobs:

```bash
# GPU 0
cofolder screen -s system.yaml -o options.yaml \
  -c library_part1.csv --col_id compound_id \
  --variable sequences,1,ligand,smiles --col_variable smiles

# GPU 1 (in parallel)
cofolder screen -s system.yaml -o options.yaml \
  -c library_part2.csv --col_id compound_id \
  --variable sequences,1,ligand,smiles --col_variable smiles
```

### Multi-Variable Mapping

Update multiple fields in one run:

```bash
cofolder screen \
  -s system.yaml \
  -o options.yaml \
  -c library.csv \
  --col_id compound_id \
  --variable sequences,1,ligand,smiles --col_variable smiles \
  --variable sequences,1,ligand,ccd --col_variable ccd
```

## Best Practices

### Library Preparation

1. **Validate SMILES**: Ensure all SMILES are valid
   ```python
   from rdkit import Chem
   valid = [Chem.MolFromSmiles(s) is not None for s in smiles_list]
   ```

2. **Remove duplicates**: Check for duplicate structures
3. **Filter properties**: Apply drug-likeness filters
4. **Standardize**: Use consistent tautomers and protonation states

### Performance Optimization

1. **Batch processing**: Split large libraries
2. **Parallel execution**: Use multiple GPUs
3. **2D conformers**: Faster for initial screening
4. **Disk space**: Ensure sufficient space for outputs

### Result Analysis

1. **Confidence thresholds**: Filter low-confidence predictions
2. **Visual inspection**: Manually check top hits
3. **Diversity**: Consider structural diversity in hits
4. **Follow-up**: Validate top candidates experimentally

## Troubleshooting

**Slow processing:**
- Reduce `recycling_steps` or `diffusion_samples`
- Use 2D instead of 3D conformers
- Split library across multiple GPUs

**High failure rate:**
- Validate SMILES before screening
- Check for incompatible functional groups
- Ensure CCD identifiers are valid

**Disk space issues:**
- Clean up intermediate files
- Use compression for large outputs
- Monitor disk usage during screening

## Next Steps

- [Validate top hits](../user-guide/validate.md) with the single-system workflow
- Try [advanced features](advanced.md) like custom scoring
- Explore [ligand handling](ligands.md) in detail

## Related

- [Screen Command Reference](../user-guide/screen.md)
- [Screen API Documentation](../api/recipes/screen.md)
