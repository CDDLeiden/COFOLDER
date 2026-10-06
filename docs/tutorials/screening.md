# Virtual Screening Tutorial

## Screen complete system variants

The regular workflow below replaces one ligand per row. To vary several system
parameters together, use the bundled `parameter_screen.csv`; each row maps to one
completed system:

```bash
cofolder screen \
  -s examples/system_screen.yaml \
  -o examples/options.yaml \
  -c examples/parameter_screen.csv \
  --col_id experiment \
  --map 'protein_sequence=sequences.0.protein.sequence' \
  --map 'ligand_smiles=sequences.1.ligand.smiles' \
  -w screen_variants
```

Run the same command with `--preflight_only` first to validate every row and inspect
the planned execution count without inference.

If the system template contains a fixed protein `msa`, map a matching MSA column
alongside every varying sequence. Otherwise preflight rejects the screen because the
template alignment belongs to only its original sequence. Removing the template MSA
enables automatic generation: each unique sequence is searched once and reused when
it appears again.

Learn how to screen compound libraries with COFOLDER.

## Overview

This tutorial demonstrates high-throughput virtual screening of a ligand library against a protein target.

## Prerequisites

- Completed [Basic Tutorial](basic.md)
- Compound library (CSV or SDF format)
- Target protein structure or sequence
- The `analysis` extra when using IFP clustering or the plotting examples below

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
      id: B
      smiles: "CCO"  # Valid template chemistry; replaced for each CSV row
```

## Step 3: Configure Screening

Create `screening_options.yaml`:

```yaml
version: 1
runtime:
  cache_path: ./cache/.boltz
  diffusion_samples: 1
runner:
  devices: 1
  recycling_steps: 3
```

## Step 4: Run Screening

```bash
cofolder screen \
  -s system_template.yaml \
  -o screening_options.yaml \
  -c library.csv \
  --col_id compound_id \
  --ligand_chain B --smiles_column smiles \
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

The public metric view is long-form. Select the computed system confidence rows,
rename `value`, and join source metadata explicitly:

```python
import pandas as pd

results = pd.read_csv("screening_output/results/metrics.csv")
library = pd.read_csv("library.csv")

confidence = results[
    (results["metric_name"] == "confidence_score")
    & (results["status"] == "computed")
    & (results["chain_id"].isna())
][["compound_id", "value"]].rename(columns={"value": "confidence_score"})

top_hits = (
    confidence.merge(library, left_on="compound_id", right_on="compound_id")
    .sort_values("confidence_score", ascending=False)
    .head(10)
)

print("Top 10 Compounds:")
print(top_hits[["compound_id", "confidence_score", "mw", "logp"]])
```

`records.jsonl` is authoritative; `metrics.csv` is its tabular metric view. Each row
retains metric status, scope, provenance, unit, and direction. Unsupported, missing,
failed, and unrequested metrics therefore remain distinguishable from numerical zero.
When calling `Screen.run()` from Python, its convenience return value is a wide
`DataFrame`; that return value is not a replacement public file format.

### SDF and MOL libraries

The packaged ethanol files exercise structural library ingestion without converting
them to CSV first:

```bash
cofolder screen -s system_template.yaml -o screening_options.yaml \
  -c examples/ethanol.sdf --ligand_chain B --id_property ID \
  -w ./sdf_screen_output

cofolder screen -s system_template.yaml -o screening_options.yaml \
  -c examples/ethanol.mol --ligand_chain B --id_property ID \
  -w ./mol_screen_output
```

For SDF, every `$$$$` record becomes one outcome. A missing `--id_property` falls
back to `record_000001`, `record_000002`, and so on. MOL input contains one record.
Duplicate identifiers are rejected by default; select `--duplicate_id_policy suffix`
or `source_index` when deterministic rewriting is intended. Malformed structural
records are retained as source-mapped failures while valid records continue.

### Discover contact-pattern families without a reference

Add clustering when no experimental pose or predefined pocket is available:

```bash
cofolder screen \
  -s system_template.yaml -o screening_options.yaml -c library.csv \
  --col_id compound_id \
  --ligand_chain B --smiles_column smiles \
  --scoring_functions confidence_metrics affinity_metrics affinity_metrics_ext \
    sasa sasa_normalized ifp_distance \
  --cluster_ifps \
  --ifp_cluster_similarity_threshold 0.5 \
  -w ./screening_output
```

Read `results/ifp_cluster_summary.csv` to inspect cluster sizes, member IDs, medoids,
consensus fingerprints, and within-cluster Jaccard similarity. The row-level
`ifp_cluster_id` is also present in both consolidated screening CSVs. This workflow
does not need `--reference_path` or `--pocket_coverage_reference`.

For publication plots, read `results/ifp_cluster_linkage.csv` as the four-column
SciPy linkage matrix and use `results/ifp_cluster_leaf_order.csv` for the recorded
leaf positions and stable execution-key labels. This reproduces the tree generated
by Screen instead of reclustering the exported fingerprints. With
`--ifp_taxonomy prolif`, `results/ifp_interaction_events.jsonl` also provides the
atom identities, roles, distances, and angles needed for atom-specific interaction
tables.

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
sns.histplot(confidence["confidence_score"].dropna(), bins=30)
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
  --ligand_chain B --smiles_column smiles

# GPU 1 (in parallel)
cofolder screen -s system.yaml -o options.yaml \
  -c library_part2.csv --col_id compound_id \
  --ligand_chain B --smiles_column smiles
```

### Fixed-System Ligand Replacement

Select one ligand entity and provide its SMILES column:

```bash
cofolder screen \
  -s system.yaml \
  -o options.yaml \
  -c library.csv \
  --col_id compound_id \
  --ligand_chain B --smiles_column smiles
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

## Post-v1 Scope

COFOLDER 1.0 screens a ligand library against a fixed system. Protein-iteration or
selectivity screening is not an implemented 1.0 workflow and is not implied by this
tutorial.

## Related

- [Screen Command Reference](../user-guide/screen.md)
- [Screen API Documentation](../api/recipes/screen.md)
