# Tutorial Map

COFOLDER tutorials come in two forms:

- written tutorials in `docs/tutorials/` for the shortest supported reading path
- marimo notebook tutorials in the repository `tutorials/` folder for interactive exploration

## Suggested Learning Path

1. Read [Quick Start](../getting-started/quickstart.md).
2. Read [User Guide Overview](../user-guide/overview.md).
3. Work through [Basic Usage](basic.md).
4. Continue with the workflow-specific or topic-specific tutorials below.

## Written Tutorials

- [Basic Usage](basic.md): first runner-backed workflow
- [Virtual Screening](screening.md): screening across a ligand library
- [Ligand Handling](ligands.md): ligand formats, conformers, and utilities
- [Advanced Features](advanced.md): advanced configuration and power-user features
- [Adding New Runners](runners.md): contributor-facing runner authoring guide
- [Backend Acceptance](backend-acceptance.md): clean-install validation for backend and CLI changes

## Notebook Tutorials

These live in the repository root under `tutorials/`:

- `bias.py`
- `validate.py`
- `screen.py`
- `oracle.py`
- `ligand_handling.py`
- `runners.py`

Launch them with marimo:

```bash
python -m pip install -e ".[tutorials]"
marimo edit tutorials/bias.py
```

Use the notebooks when you want to inspect intermediate outputs, execute cells step by step, or adapt the tutorial flow to your own systems. The workflow-first set centers on `bias`, `validate`, `screen`, and `oracle`, then adds cross-cutting ligand handling plus one contributor notebook for runners and backend acceptance.
