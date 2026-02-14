# AGENTS.md

Guidance for Codex-style coding agents working in this repository.

## Project Scope

- Repository: `cofolder`
- Default integration branch: `dev`
- Main package root: `src/cofolder`
- Tests root: `tests`
- Docs root: `docs`

## Primary Goals

- Implement features with minimal behavioral regressions.
- Keep CLI behavior stable unless change is explicitly requested.
- Add or update tests for changed behavior.
- Update docs for user-facing changes.

## Environment Setup

Use Python 3.10-3.12.

Install for development:

```bash
pip install -e ".[docs,test]"
```

If UI changes are involved, ensure UI dependencies are available:

```bash
pip install -e ".[ui]"
```

## Common Commands

Run tests:

```bash
pytest -q
```

Run a focused test file:

```bash
pytest -q tests/path/to/test_file.py
```

Show CLI help:

```bash
cofolder -h
cofolder validate -h
cofolder screen -h
cofolder oracle -h
```

Run docs locally:

```bash
mkdocs serve
```

Run UI locally:

```bash
./run_ui.sh
```

## Code Organization

- CLI entrypoint: `src/cofolder/cli.py`
- Recipes: `src/cofolder/recipes/`
- Domain modules: `src/cofolder/modules/`
- Streamlit app: `src/cofolder/ui/app.py`

## Architecture

COFOLDER consists of two software parts:

- A fully independent API in `modules/`
- CLI recipes that run common API protocols

Protocol hierarchy for recipes:

- `validate` is the base/core protocol
- `screen` and `oracle` are built on top of this protocol

When adding functionality:

- Put command workflows in `recipes/`.
- Put reusable logic in `modules/`.
- Keep CLI parsing/dispatching in `cli.py` thin.

## Development Conventions

- Follow existing style and PEP 8.
- Add type hints for new/changed function signatures when practical.
- Add or update NumPy-style docstrings for new public functions/classes.
- Keep functions focused; prefer small composable helpers.
- Avoid unrelated refactors in feature PRs.

## Testing Expectations

For behavior changes, add tests in the closest relevant location, for example:

- `tests/recipes/` for recipe-level behavior
- `tests/modules/` for module logic
- `tests/test_cli.py` for CLI-level changes
- `tests/test_ui.py` for UI behavior

At minimum before finishing work:

1. Run targeted tests for changed components.
2. Run `pytest -q` if feasible.
3. Report any tests not run and why.

## Documentation Expectations

Update docs when behavior or interfaces change:

- User guides under `docs/user-guide/`
- Getting-started material under `docs/getting-started/`
- API docs under `docs/api/`
- Tutorials when workflows change significantly

## Git Workflow

- Branch from `dev`.
- Use branch naming like `feature/<short-description>` or `fix/<short-description>`.
- Keep commits scoped and messages clear.
- Target PRs to `dev`.

## Agent Operating Rules

- Do not remove or rewrite unrelated user changes.
- Do not run destructive git commands unless explicitly requested.
- Prefer minimal diffs that solve the requested task.
- If assumptions are required, state them clearly in the final summary.
- If blocked by missing inputs or environment constraints, explain the blocker and propose the shortest next step.

## Completion Checklist

Before handing work back, confirm:

- Code changes compile/import successfully (as applicable).
- Relevant tests were added/updated and executed.
- Relevant docs were updated.
- Final summary includes changed files and validation performed.
