# Tutorials

This folder contains the interactive COFOLDER tutorial notebooks in marimo format.

The notebooks complement the written tutorials in `docs/tutorials/`. A good rule of thumb is:

- read `docs/` when you want the shortest supported explanation
- use these marimo notebooks when you want a guided, executable walkthrough

## Launching The Tutorials

Install the notebook runtime and the backend you want to use:

```bash
python -m pip install -e ".[analysis,tutorials,boltz2]"
marimo edit tutorials/bias.py
```

The `bias` parts of the tutorial set can run without a prediction backend, but the
`validate`, `screen`, `oracle`, and some ligand-helper demonstrations assume a backend
environment such as `.[boltz2]`.

The backend system-input tutorials require a dedicated backend environment and a
CUDA-capable machine. For example:

```bash
python -m pip install -e ".[analysis,tutorials,boltz2]"
marimo edit tutorials/boltz_system_inputs.py
```

For OpenFold3, prepare its cache before launching the notebook:

```bash
python -m pip install -e ".[analysis,tutorials,openfold3]"
export OPENFOLD_CACHE="$PWD/cache/.openfold3-cache"
cofolder-tools setup-openfold3
marimo edit tutorials/openfold3_system_inputs.py
```

## Recommended Order

1. `bias.py`
2. `validate.py`
3. `screen.py`
4. `oracle.py`
5. `structure_gated_oracle.py`
6. `ligand_handling.py`
7. `runners.py`
8. `boltz_system_inputs.py` or `openfold3_system_inputs.py` in the matching backend environment

## Notebook Roles

- `bias.py`: standalone reference-overlap diagnostics workflow
- `validate.py`: single-system validation workflow
- `screen.py`: multi-ligand screening workflow
- `oracle.py`: single-metric oracle workflow
- `structure_gated_oracle.py`: MAPK14 pose-aware lexicographic Oracle with an
  offline walkthrough and optional Boltz2 execution
- `ligand_handling.py`: specific ligand preparation and conformer-handling tutorial
- `runners.py`: contributor notebook for authoring new runners plus backend-acceptance launch guidance
- `boltz_system_inputs.py`: real-backend walkthrough for protein, DNA, RNA, ligand,
  constraints, and runner-specific rejection across the Boltz family
- `openfold3_system_inputs.py`: real-backend walkthrough for nucleic acids, OpenFold3
  pocket translation, and rejection of unsupported bonds

The system-input notebooks keep inference behind explicit checkboxes, show the actual
backend input, and verify normalized chain metadata. The input-contract notebook in
`src/cofolder/acceptance/` remains the stricter maintainer lane for recording release
acceptance.
