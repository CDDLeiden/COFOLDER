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
- OpenFold3 installation guidance, runner wiring, or confidence-only workflow behavior

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

### `openfold3`

OpenFold3 now installs through a COFOLDER optional extra. Use a fresh environment, install the acceptance tooling plus the backend extra, then run the repo setup script before launching the expensive acceptance cells:

```bash
conda create -n cofolder-acceptance-openfold3 python=3.12
conda activate cofolder-acceptance-openfold3
python -m pip install -e ".[acceptance,openfold3]"
export OPENFOLD_CACHE="$PWD/.openfold3-cache"
scripts/setup_openfold3.sh
python --version
cofolder --help
```

Notes:

- OpenFold3 upstream currently recommends `pixi` for reproducible environments, but the upstream `openfold3` pip package is also documented and is what the COFOLDER `openfold3` extra installs.
- `scripts/setup_openfold3.sh` wraps upstream `setup_openfold`, exports `OPENFOLD_CACHE` to a standard location by default (`~/.openfold3`), and prepares the cache, model parameters, and CCD before you run the acceptance workflows.
- The wrapper also removes the current upstream checkpoint-choice ambiguity by answering the setup prompts explicitly: it uses `OPENFOLD_CACHE` for both path questions and selects parameter download choice `1` by default. Use `OPENFOLD3_PARAMETER_CHOICE=2 scripts/setup_openfold3.sh` if you want all published checkpoints instead.
- The wrapper answers the upstream integration-test prompt with `no` by default. Use `OPENFOLD3_RUN_INTEGRATION_TESTS=yes scripts/setup_openfold3.sh` only when you intentionally want those upstream tests to run during setup.
- Upstream docs note that first inference can also download default model parameters to `$HOME/.openfold3`, but for manual acceptance this project prefers the explicit setup script so environment readiness is checked before the expensive lane starts.
- This first-pass acceptance lane is confidence-only. It does not claim affinity support for OpenFold3.

If you prefer a non-editable install from the current checkout, replace:

- `python -m pip install -e ".[acceptance,boltz1]"`
- `python -m pip install -e ".[acceptance,boltz2]"`
- `python -m pip install -e ".[acceptance,boltz-community]"`
- `python -m pip install -e ".[acceptance,openfold3]"`

with:

- `python -m pip install ".[acceptance,boltz1]"`
- `python -m pip install ".[acceptance,boltz2]"`
- `python -m pip install ".[acceptance,boltz-community]"`
- `python -m pip install ".[acceptance,openfold3]"`

To discard an acceptance environment completely:

```bash
conda deactivate
conda env remove -n cofolder-acceptance-boltz1
```

Swap the environment name as needed for `boltz2`, `boltz-community`, or `openfold3`.

## Launching The Notebooks

Recommended launch pattern from inside the matching clean environment:

```bash
marimo edit "$(python -c 'import cofolder.acceptance.boltz2_backend_acceptance as nb; print(nb.__file__)')"
```

Swap the module name for:

- `cofolder.acceptance.boltz1_backend_acceptance`
- `cofolder.acceptance.boltz_community_backend_acceptance`
- `cofolder.acceptance.openfold3_backend_acceptance`

Each notebook opens safely by default:

- the expensive `validate`, `screen`, and `oracle` steps are gated behind explicit checkboxes
- each subprocess run is bound to the notebook workspace so relative cache paths stay isolated
- the packaged options file is rewritten per run so the backend cache points into that run's temporary workspace
- command output streams live to the notebook console while the command runs

For `openfold3`, the notebook uses a dedicated OpenFold3 options fixture and points `OPENFOLD_CACHE` at a notebook-local path, but you still need to run `scripts/setup_openfold3.sh` against that path before triggering the expensive workflow cells.

## What The Notebooks Exercise

Each backend notebook runs the command-line workflows that matter for shared behavior:

- `cofolder validate`
- `cofolder screen`
- `cofolder oracle`

The notebooks call the real CLI in a clean environment instead of importing recipes directly, which makes them suitable as a pre-merge acceptance lane for runner and workflow changes.

For `openfold3`, the positive path is:

- `validate` with `confidence_metrics`
- `screen` with `confidence_metrics`
- `oracle` with a real OpenFold3-produced confidence metric such as `sample_ranking_score`

## Expected Pass/Fail Semantics

- `boltz2` and `boltz-community`: `validate`, `screen`, and `oracle` should produce populated affinity-related outputs.
- `boltz1`: `validate` and `screen` should warn that affinity scoring groups are unsupported while leaving downstream affinity-related columns present-but-empty; `oracle` should use `confidence_score` as the supported positive-path metric.
- `openfold3`: `validate`, `screen`, and `oracle` should produce populated confidence outputs, including OpenFold3-native confidence fields. This lane does not request affinity groups, and it must not imply that structural-confidence outputs are affinity proxies.

Keep scientific wording disciplined when reviewing these outputs:

- validation metrics, confidence metrics, and affinity metrics are distinct concepts
- structural confidence should not be described as a proxy for binding affinity

## Optional Bias Setup

Bias-related setup is intentionally outside the first-pass backend acceptance lane.

If you need `--assess_bias`, install `mmseqs2` in the same environment and follow the bias setup documented in [Installation](../getting-started/installation.md). Keep that work separate from the base backend acceptance pass unless the change under review specifically affects bias behavior.
