# COFOLDER Documentation

COFOLDER provides command-line workflows for protein-ligand co-folding and pre-cofolding diagnostics.

The supported user-facing workflows are:

- `bias`: standalone reference-overlap diagnostics from `system.yaml`
- `validate`: single runner-backed co-folding workflow
- `screen`: validate-style workflow across a ligand library
- `oracle`: single-metric runner-backed scoring workflow

## Start Here

If you are new to COFOLDER, follow this order:

1. [Installation](getting-started/installation.md)
2. [Quick Start](getting-started/quickstart.md)
3. [Workflow Overview](user-guide/overview.md)
4. [Tutorial Map](tutorials/index.md)

## Choose Your Path

### I want the shortest path to a first run

- Read [Quick Start](getting-started/quickstart.md)
- Copy files from `examples/`
- Run `cofolder validate`

### I want to understand which command to use

- Read [User Guide Overview](user-guide/overview.md)
- Then jump to `bias`, `validate`, `screen`, or `oracle`

### I want step-by-step learning

- Start with [Basic Usage](tutorials/basic.md)
- Continue through the [Tutorial Map](tutorials/index.md)

### I want contributor-facing backend context

- Read [Adding New Runners](tutorials/runners.md)
- Use [Backend Acceptance](tutorials/backend-acceptance.md) before promoting backend or CLI-adjacent changes

## Repository Orientation

The public repo is organized around a few distinct user paths:

- `docs/`: curated reading path
- `examples/`: copy-and-run input files
- repository `tutorials/`: interactive notebooks
- `scripts/`: optional setup and bias-data helpers
- `legacy/`: archived historical material

For a concise map of those roles, see [Repository Tour](getting-started/repository-tour.md).

## Scientific Scope

COFOLDER keeps validation metrics, model-derived confidence metrics, and structure-derived diagnostics distinct. Structural confidence should not be described as a proxy for binding affinity.
