# Basic Usage Tutorial

This tutorial walks through basic Boltz-Lab usage for protein-ligand co-folding.

## Prerequisites

- Boltz-Lab installed
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
      smiles: "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
      ccd: "IBP"
```

## Step 2: Configure Boltz Options

Create `options.yaml`:

```yaml
out_dir: output
devices: [0]
num_models: 1
recycling_steps: 3
diffusion_samples: 1
```

## Step 3: Run Prediction

```bash
boltz-lab predict -s system.yaml -b options.yaml -w ./tutorial_output
```

## Step 4: Examine Output

Check the output directory:

```bash
ls -l tutorial_output/
```

You should see:
- `predictions.cif` - Predicted structure
- `confidence_model_0.json` - Confidence scores
- Log files

## Step 5: Visualize Results

Load `predictions.cif` in PyMOL, ChimeraX, or your favorite molecular viewer:

```bash
pymol tutorial_output/predictions.cif
```

## Step 6: Check Confidence

Examine the confidence scores:

```python
import json

with open('tutorial_output/confidence_model_0.json') as f:
    confidence = json.load(f)

print(f"Overall confidence: {confidence['overall']}")
```

## Next Steps

### Try Different Options

Experiment with different Boltz parameters:

```yaml
# Higher quality prediction
recycling_steps: 5
diffusion_samples: 5

# Ensemble prediction
num_models: 5
```

### Add 3D Conformer Generation

```bash
boltz-lab predict \
  -s system.yaml \
  -b options.yaml \
  --generate_conformers 3D
```

### Enable Debug Logging

```bash
boltz-lab predict -s system.yaml -b options.yaml -d
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

- [Predict Command Guide](../user-guide/predict.md)
- [Configuration Reference](../getting-started/configuration.md)
- [Virtual Screening Tutorial](screening.md)
