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
      id: L
      smiles: "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
```

### Multiple Ligands

```yaml
version: 1
sequences:
  - protein:
      id: "protein"
      fasta: "SEQUENCE..."
  - ligand:
      id: L1
      smiles: "CCO"
  - ligand:
      id: L2
      smiles: "CCN"
```

## Ensemble Predictions

Generate multiple predictions for uncertainty estimation:

```yaml
# options.yaml
version: 1
runtime:
  cache_path: ./cache/.boltz
  diffusion_samples: 10
runner:
  devices: 1
  recycling_steps: 5
```

## Custom Boltz Parameters

### High-Quality Predictions

```yaml
version: 1
runtime:
  cache_path: ./cache/.boltz
  diffusion_samples: 20
runner:
  devices: 1
  recycling_steps: 10
  sampling_steps: 500
  step_scale: 1.5
```

### Fast Screening

```yaml
version: 1
runtime:
  cache_path: ./cache/.boltz
  diffusion_samples: 1
runner:
  devices: 1
  recycling_steps: 1
  sampling_steps: 100
```

## Covalent Binding

For covalent inhibitors, specify the covalent bond:

```yaml
version: 1
sequences:
  - protein:
      id: A
      sequence: "SEQUENCE..."
  - ligand:
      id: B
      ccd: COV
constraints:
  - bond:
      atom1: [A, 145, SG]
      atom2: [B, 1, C12]
```

Residues are 1-based. Bond endpoints use `[chain_id, residue_id, atom_name]`;
the atom name must exist in the selected residue or ligand CCD.

## Custom Scoring Functions

### Implement a custom Oracle function

```python
from cofolder.recipes.oracle import Oracle


def custom_score(context):
    """Combine real, aggregated Validate outputs into one finite scalar."""
    affinity = context.aggregated_metrics["ligand_B__pIC50"]
    confidence = context.aggregated_metrics["system__confidence_score"]
    pocket_coverage = context.aggregated_metrics[
        "ligand_B__pocket_coverage_custom"
    ]
    return affinity * confidence * pocket_coverage


score = Oracle(
    wrk_dir="oracle_custom",
    system_path="examples/system.yaml",
    options_path="examples/options.yaml",
    input_smiles="CCO",
    scoring_functions=[
        "confidence_metrics",
        "affinity_metrics_ext",
        "ifp_distance",
    ],
    pocket_coverage_reference="A25 G48 Y51",
    reproduction_metrics=["pocket_coverage"],
    scoring_function=custom_score,
).run()
```

The callable receives raw metric DataFrames as well as qualified, aggregated numeric
metrics. See the [Oracle guide](../user-guide/oracle.md) for weighted composites,
SASA/IFP gates, and explicit down-weight, penalty, and non-binder policies.

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
        shutil.copytree(
            f"{tmpdir}/results/structures",
            "output/structures",
            dirs_exist_ok=True,
        )

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
