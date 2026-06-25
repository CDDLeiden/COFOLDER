# Tutorials

This folder contains the interactive COFOLDER tutorial notebooks in marimo format.

The notebooks complement the written tutorials in `docs/tutorials/`. A good rule of thumb is:

- read `docs/` when you want the shortest supported explanation
- use these marimo notebooks when you want a guided, executable walkthrough

## Launching The Tutorials

Install the notebook runtime and the backend you want to use:

```bash
python -m pip install -e ".[tutorials]"
python -m pip install -e ".[boltz2]"
marimo edit tutorials/bias.py
```

The `bias` parts of the tutorial set can run without a prediction backend, but the
`validate`, `screen`, `oracle`, and some ligand-helper demonstrations assume a backend
environment such as `.[boltz2]`.

## Recommended Order

1. `bias.py`
2. `validate.py`
3. `screen.py`
4. `oracle.py`
5. `ligand_handling.py`
6. `runners.py`

## Notebook Roles

- `bias.py`: standalone reference-overlap diagnostics workflow
- `validate.py`: single-system validation workflow
- `screen.py`: multi-ligand screening workflow
- `oracle.py`: single-metric oracle workflow
- `ligand_handling.py`: specific ligand preparation and conformer-handling tutorial
- `runners.py`: contributor notebook for authoring new runners plus backend-acceptance launch guidance
