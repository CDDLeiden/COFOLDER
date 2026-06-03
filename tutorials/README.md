# Tutorials

This folder contains interactive notebook tutorials for users who want a hands-on path through COFOLDER.

The notebooks complement the written tutorials in `docs/tutorials/`. A good rule of thumb is:

- read `docs/` when you want the shortest supported explanation
- use these notebooks when you want to experiment step by step

## Recommended Order

1. `base.ipynb`
2. `system_validation.ipynb`
3. `virtual_screening.ipynb`
4. `oracle.ipynb`
5. `ligand_handling.ipynb`
6. `polymer_handling.ipynb`
7. `runners.ipynb`

## Notebook Roles

- `base.ipynb`: first hands-on walkthrough
- `system_validation.ipynb`: single-system validation workflow
- `virtual_screening.ipynb`: multi-ligand screening workflow
- `oracle.ipynb`: single-metric oracle workflow
- `ligand_handling.ipynb`: ligand preparation and utility patterns
- `polymer_handling.ipynb`: polymer and system-input handling
- `runners.ipynb`: backend and runner-facing concepts for advanced users and contributors

## Before You Open A Notebook

Install the package and any backend you need first:

```bash
python -m pip install -e .
python -m pip install -e ".[boltz2]"
```

Then open the notebooks in Jupyter, VS Code, or another notebook environment that uses the same Python environment.
