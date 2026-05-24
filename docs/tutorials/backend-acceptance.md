# Backend Acceptance Tutorial

This guide documents the manual backend acceptance lane for COFOLDER maintainers.

Use it when you need fresh-environment validation of the real CLI workflows against a real backend, especially before promoting backend-, runner-, or CLI-adjacent changes toward `main`.

This lane is intentionally expensive:

- it is not part of routine `pytest`
- it is not part of CI-default checks
- it should be run deliberately when backend behavior or runner integration has changed

## When To Run This Lane

Run the backend acceptance notebooks before promoting changes that affect:

- runner implementations under `src/cofolder/modules/runners/`
- shared CLI behavior in `src/cofolder/cli.py`
- workflow behavior in `src/cofolder/recipes/validate.py`, `screen.py`, or `oracle.py`
- packaging or dependency behavior for `boltz1`, `boltz2`, `boltz-community`, or the `acceptance` extra

If you are adding or integrating a new runner, pair this guide with [Adding New Runners](runners.md).

## Python Version

Use a fresh conda environment with:

- Python `>=3.11`
- Python `3.11` or `3.12` strongly recommended for this acceptance lane
- do not use Python `3.13` for `boltz1` or `boltz2` acceptance environments in this project, because those extras are currently constrained to `<3.13` in `pyproject.toml`

## Clean Install Rules

Treat each backend as a separate environment:

- create one new conda environment per backend
- install exactly one backend extra per environment
- do not reuse an existing development environment
- do not install `boltz1`, `boltz2`, and `boltz-community` into the same environment
- if you accidentally mixed backend installs, remove that conda environment and start over with a new one

The commands below assume you are in the repository root and want the current checkout installed in editable mode.

## Full Clean Installs

### `boltz1`

```bash
conda create -n cofolder-acceptance-boltz1 python=3.12
conda activate cofolder-acceptance-boltz1
python -m pip install -e ".[acceptance,boltz1]"
python --version
cofolder --help
```

If Python `3.12` is not available in your conda setup, use Python `3.11` instead.

### `boltz2`

```bash
conda create -n cofolder-acceptance-boltz2 python=3.12
conda activate cofolder-acceptance-boltz2
python -m pip install -e ".[acceptance,boltz2]"
python --version
cofolder --help
```

### `boltz-community`

```bash
conda create -n cofolder-acceptance-boltz-community python=3.12
conda activate cofolder-acceptance-boltz-community
python -m pip install -e ".[acceptance,boltz-community]"
python --version
cofolder --help
```

If you prefer a non-editable install from the current checkout, replace:

- `python -m pip install -e ".[acceptance,boltz1]"`
- `python -m pip install -e ".[acceptance,boltz2]"`
- `python -m pip install -e ".[acceptance,boltz-community]"`

with:

- `python -m pip install ".[acceptance,boltz1]"`
- `python -m pip install ".[acceptance,boltz2]"`
- `python -m pip install ".[acceptance,boltz-community]"`

To discard an acceptance environment completely:

```bash
conda deactivate
conda env remove -n cofolder-acceptance-boltz1
```

Swap the environment name as needed for `boltz2` or `boltz-community`.

## Launching The Notebooks

Recommended launch pattern from inside the matching clean environment:

```bash
marimo edit "$(python -c 'import cofolder.acceptance.boltz2_backend_acceptance as nb; print(nb.__file__)')"
```

Swap the module name for:

- `cofolder.acceptance.boltz1_backend_acceptance`
- `cofolder.acceptance.boltz_community_backend_acceptance`

Each notebook opens safely by default:

- the expensive `validate`, `screen`, and `oracle` steps are gated behind explicit checkboxes
- each subprocess run is bound to the notebook workspace so relative cache paths stay isolated
- the packaged options file is rewritten per run so the backend cache points into that run's temporary workspace
- command output streams live to the notebook console while the command runs

## What The Notebooks Exercise

Each backend notebook runs the command-line workflows that matter for shared behavior:

- `cofolder validate`
- `cofolder screen`
- `cofolder oracle`

The notebooks call the real CLI in a clean environment instead of importing recipes directly, which makes them suitable as a pre-merge acceptance lane for runner and workflow changes.

## Expected Pass/Fail Semantics

- `boltz2` and `boltz-community`: `validate`, `screen`, and `oracle` should produce populated affinity-related outputs.
- `boltz1`: `validate` and `screen` should warn that affinity scoring groups are unsupported while leaving downstream affinity-related columns present-but-empty; `oracle` should use `confidence_score` as the supported positive-path metric.

Keep scientific wording disciplined when reviewing these outputs:

- validation metrics, confidence metrics, and affinity metrics are distinct concepts
- structural confidence should not be described as a proxy for binding affinity

## Optional Bias Setup

Bias-related setup is intentionally outside the first-pass backend acceptance lane.

If you need `--assess_bias`, install `mmseqs2` in the same environment and follow the bias setup documented in [Installation](../getting-started/installation.md). Keep that work separate from the base backend acceptance pass unless the change under review specifically affects bias behavior.
