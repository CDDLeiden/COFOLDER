# Examples

This folder contains the smallest supported input files for learning and smoke-testing COFOLDER workflows.

Use these files when you want to copy a working starting point before editing it for your own system.

## Files

- `system.yaml`: basic single-system example for `validate` and `oracle`
- `options.yaml`: standard runner options example
- `system_screen.yaml`: screening template used with a ligand library
- `ligand_screen.csv`: small screening CSV used by the screen tutorial
- `ifp_clustering_demo.py`: CPU-only deterministic IFP clustering demonstration
- `system_covalent.yaml`: example of a more specialized system definition
- `system_nucleic_acid.yaml`: protein/DNA/RNA/ligand system with a pocket constraint
- `4HJO.pdb` and `4HJO.cif`: structure fixtures that support examples and manual inspection

## Recommended Use

Start with:

1. `system.yaml`
2. `options.yaml`
3. one of the quick commands from the root `README.md` or `docs/getting-started/quickstart.md`

Then move to `system_screen.yaml` and `ligand_screen.csv` when you want to explore `screen`.

The copyable no-reference IFP clustering workflow is:

```bash
cofolder screen \
  -s examples/system_screen.yaml \
  -o examples/options.yaml \
  -c examples/ligand_screen.csv \
  --col_id Name \
  --variable sequences,1,ligand,smiles --col_variable SMILES \
  --scoring_functions ifp_distance \
  --cluster_ifps \
  -w ./screen_out
```

For Boltz-family runners the first row generates a missing fixed-protein MSA and all
later rows reuse it. To skip the MSA server entirely, add an `msa` path to the protein
entry in `system_screen.yaml`; relative paths are resolved from that YAML's directory.

## Relationship To The Docs

- `docs/getting-started/quickstart.md` explains the first command to run
- `docs/tutorials/index.md` maps the written tutorials
- `tutorials/README.md` maps the notebook tutorials

The examples here are intended to stay small, readable, and copyable.
