# Examples

This folder contains the smallest supported input files for learning and smoke-testing COFOLDER workflows.

Use these files when you want to copy a working starting point before editing it for your own system.

## Files

- `system.yaml`: basic single-system example for `validate` and `oracle`
- `options.yaml`: standard runner options example
- `system_screen.yaml`: screening template used with a ligand library
- `ligand_screen.csv`: small screening CSV used by the screen tutorial
- `system_covalent.yaml`: example of a more specialized system definition
- `4HJO.pdb` and `4HJO.cif`: structure fixtures that support examples and manual inspection

## Recommended Use

Start with:

1. `system.yaml`
2. `options.yaml`
3. one of the quick commands from the root `README.md` or `docs/getting-started/quickstart.md`

Then move to `system_screen.yaml` and `ligand_screen.csv` when you want to explore `screen`.

## Relationship To The Docs

- `docs/getting-started/quickstart.md` explains the first command to run
- `docs/tutorials/index.md` maps the written tutorials
- `tutorials/README.md` maps the notebook tutorials

The examples here are intended to stay small, readable, and copyable.
