# Advanced Features

Explore advanced capabilities and customization options in COFOLDER.

## Multi-Chain Systems

### Protein-Protein-Ligand Complexes

```yaml
version: 1
sequences:
  - protein:
      id: "chain_A"
      fasta: "SEQUENCE_A..."
  - protein:
      id: "chain_B"
      fasta: "SEQUENCE_B..."
  - ligand:
      smiles: "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
      ccd: "IBP"
```

### Multiple Ligands

```yaml
version: 1
sequences:
  - protein:
      id: "protein"
      fasta: "SEQUENCE..."
  - ligand:
      smiles: "SMILES_1"
      ccd: "LIG1"
  - ligand:
      smiles: "SMILES_2"
      ccd: "LIG2"
```

## Ensemble Predictions

Generate multiple predictions for uncertainty estimation:

```yaml
# options.yaml
out_dir: ensemble_output
devices: [0]
num_models: 5           # Generate 5 models
recycling_steps: 5
diffusion_samples: 10   # 10 samples per model
```

## Custom Boltz Parameters

### High-Quality Predictions

```yaml
out_dir: high_quality
devices: [0]
num_models: 5
recycling_steps: 10
diffusion_samples: 20
sampling_steps: 500
diffusion_temperature: 0.8
```

### Fast Screening

```yaml
out_dir: fast_screen
devices: [0]
num_models: 1
recycling_steps: 1
diffusion_samples: 1
sampling_steps: 100
```

## Covalent Binding

For covalent inhibitors, specify the covalent bond:

```yaml
version: 1
sequences:
  - protein:
      id: "protease"
      fasta: "SEQUENCE..."
  - ligand:
      smiles: "COVALENT_SMILES"
      ccd: "COV"
      covalent:
        protein_residue: "CYS145"
        ligand_atom: 12
```

## Custom Scoring Functions

### Implement Custom Oracle

```python
from cofolder.recipes.oracle import Oracle
import numpy as np


class CustomOracle(Oracle):
    def score(self, prediction):
        # Custom scoring logic
        confidence = prediction['confidence']
        rmsd = self.calculate_rmsd(prediction)

        # Combined score
        score = confidence * np.exp(-rmsd)
        return score

    def calculate_rmsd(self, prediction):
        # RMSD calculation
        pass
```

## Batch Processing

### Parallel Screening

```python
import multiprocessing as mp
from functools import partial

def screen_batch(compounds, gpu_id):
    # Configure to use specific GPU
    options = load_options()
    options['devices'] = [gpu_id]

    # Run screening
    results = []
    for compound in compounds:
        result = run_prediction(compound, options)
        results.append(result)
    return results

# Split compounds across GPUs
n_gpus = 4
compound_batches = split_list(compounds, n_gpus)

with mp.Pool(n_gpus) as pool:
    func = partial(screen_batch)
    all_results = pool.starmap(func,
        [(batch, i) for i, batch in enumerate(compound_batches)])
```

## Result Post-Processing

### Aggregate Multiple Predictions

```python
import pandas as pd
import glob

# Collect all results
result_files = glob.glob("output/*/results.csv")
dfs = [pd.read_csv(f) for f in result_files]
combined = pd.concat(dfs, ignore_index=True)

# Calculate statistics
stats = combined.groupby('compound_id').agg({
    'confidence': ['mean', 'std', 'min', 'max'],
    'rmsd': ['mean', 'std']
})
```

### Extract Interaction Fingerprints

```python
from cofolder.modules.analytics import calculate_ifp

predictions = load_predictions("output/")
reference = load_reference("reference.pdb")

for pred in predictions:
    ifp = calculate_ifp(pred, reference)
    overlap = ifp.overlap()
    print(f"IFP overlap: {overlap:.2f}")
```

## Integration with Other Tools

### PyMOL Automation

```python
import pymol
from pymol import cmd

def visualize_results(predictions, reference=None):
    cmd.load(predictions[0], "pred1")

    if reference:
        cmd.load(reference, "ref")
        cmd.align("pred1", "ref")

    cmd.show("cartoon", "pred1")
    cmd.show("sticks", "organic")
    cmd.png("visualization.png", dpi=300)
```

### RDKit Integration

```python
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

def analyze_ligand(smiles):
    mol = Chem.MolFromSmiles(smiles)

    # Calculate properties
    props = {
        'MW': Descriptors.MolWt(mol),
        'LogP': Descriptors.MolLogP(mol),
        'TPSA': Descriptors.TPSA(mol),
        'HBD': Descriptors.NumHDonors(mol),
        'HBA': Descriptors.NumHAcceptors(mol),
    }

    return props
```

## Performance Optimization

### GPU Memory Management

```python
import torch

def clear_gpu_cache():
    torch.cuda.empty_cache()
    torch.cuda.synchronize()

# Clear between predictions
for compound in compounds:
    result = predict(compound)
    clear_gpu_cache()
```

### Disk I/O Optimization

```python
import tempfile
import shutil

def predict_with_tmpdir(compound):
    with tempfile.TemporaryDirectory() as tmpdir:
        # Run prediction in temp directory
        result = run_prediction(compound, work_dir=tmpdir)

        # Copy only essential results
        shutil.copy(f"{tmpdir}/predictions.cif", "output/")

    return result
```

## Debugging

### Detailed Logging

```python
import logging

logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('debug.log'),
        logging.StreamHandler()
    ]
)
```

### Checkpoint Recovery

```python
import pickle

def save_checkpoint(state, filename):
    with open(filename, 'wb') as f:
        pickle.dump(state, f)

def load_checkpoint(filename):
    with open(filename, 'rb') as f:
        return pickle.load(f)

# In screening loop
for i, compound in enumerate(compounds):
    if i % 100 == 0:
        save_checkpoint({'index': i, 'results': results},
                       f'checkpoint_{i}.pkl')
```

## Related

- [API Reference](../api/recipes/validate.md)
- [Configuration Guide](../getting-started/configuration.md)
- [Contributing](../contributing.md)
