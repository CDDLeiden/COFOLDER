# Oracle Command

The `oracle` command uses Boltz as a scoring function for molecular design workflows.

## Basic Usage

```bash
boltz-lab oracle -s system.yaml -b options.yaml
```

## Arguments

### Required Arguments

- `-s, --system_path`: Path to system YAML file
- `-b, --boltz_options_path`: Path to Boltz options YAML file

### Optional Arguments

- `-w, --wrk_dir`: Working directory
- `-d, --debug`: Enable debug logging

## Use Cases

### Molecular Optimization

Use Boltz as a fitness function in optimization:

```python
from boltz_lab.recipes.oracle import Oracle

def score_molecule(smiles: str) -> float:
    # Update system with SMILES
    # Run oracle
    oracle = Oracle(
        wrk_dir="./tmp",
        system_path="system.yaml",
        options_path="options.yaml"
    )
    result = oracle.run()
    return result.score
```

### Design-Make-Test Cycles

Integrate with design workflows:

1. Design: Generate candidate molecules
2. Score: Use oracle to predict binding
3. Select: Choose top candidates
4. Make: Synthesize selected compounds
5. Test: Validate experimentally

### Active Learning

Use in active learning loops:

- Score unlabeled molecules
- Select informative samples
- Update model with new data
- Repeat

## Integration Examples

### With RDKit

```python
from rdkit import Chem
from rdkit.Chem import AllChem

# Generate analogs
mol = Chem.MolFromSmiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O")
analogs = generate_analogs(mol)

# Score with oracle
scores = []
for analog in analogs:
    smiles = Chem.MolToSmiles(analog)
    score = oracle_score(smiles)
    scores.append(score)
```

### With Optimization Algorithms

```python
from scipy.optimize import differential_evolution

def objective(params):
    smiles = params_to_smiles(params)
    return -oracle_score(smiles)  # Minimize negative score

result = differential_evolution(objective, bounds)
```

## Output

Returns prediction results that can be used as scores in optimization workflows.

## Tips

- Cache results to avoid redundant calculations
- Use batch processing for efficiency
- Consider uncertainty estimates
- Validate oracle scores experimentally

## Related

- [Oracle API Reference](../api/recipes/oracle.md)
- [Oracle Tutorial](../tutorials/basic.md)
