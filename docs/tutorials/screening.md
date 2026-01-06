# Virtual Screening Tutorial

Learn how to screen compound libraries with Boltz-Lab.

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
boltz-lab screen \
  -s system_template.yaml \
  -b screening_options.yaml \
  -v "sequences,1,ligand,smiles" \
  -c library.csv \
  --col_variable smiles \
  --col_id compound_id \
  --merge_data "mw,logp" \
  -w ./screening_output
```

## Step 5: Monitor Progress

The screen command will process each compound sequentially. Monitor with:

```bash
tail -f screening_output/screening.log
```

## Step 6: Analyze Results

Results are saved in a CSV file:

```python
import pandas as pd

results = pd.read_csv('screening_output/screening_results.csv')

# Sort by confidence score
top_hits = results.sort_values('confidence', ascending=False).head(10)

print("Top 10 Compounds:")
print(top_hits[['compound_id', 'confidence', 'mw', 'logp']])
```

## Step 7: Visualize Top Hits

```python
import matplotlib.pyplot as plt
import seaborn as sns

# Plot confidence distribution
plt.figure(figsize=(10, 6))
sns.histplot(results['confidence'], bins=30)
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
boltz-lab screen -s system.yaml -b options.yaml \
  -c library_part1.csv --col_variable smiles --col_id compound_id \
  -v "sequences,1,ligand,smiles"

# GPU 1 (in parallel)
boltz-lab screen -s system.yaml -b options.yaml \
  -c library_part2.csv --col_variable smiles --col_id compound_id \
  -v "sequences,1,ligand,smiles"
```

### SDF Input

Screen from SDF file:

```bash
boltz-lab screen \
  -s system.yaml \
  -b options.yaml \
  -v "sequences,1,ligand,smiles" \
  --variable_sdf library.sdf \
  --property_id ID \
  --generate_conformers 3D
```

### 3D Conformer Generation

Generate 3D conformers for better accuracy:

```bash
boltz-lab screen \
  -s system.yaml \
  -b options.yaml \
  -v "sequences,1,ligand,smiles" \
  -c library.csv \
  --col_variable smiles \
  --col_id compound_id \
  --generate_conformers 3D
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

- [Evaluate top hits](../user-guide/evaluate.md) against reference structures
- Try [advanced features](advanced.md) like custom scoring
- Explore [ligand handling](ligands.md) in detail

## Related

- [Screen Command Reference](../user-guide/screen.md)
- [Screen API Documentation](../api/recipes/screen.md)
