# Examples

This folder contains the smallest supported input files for learning and smoke-testing COFOLDER workflows.

Use these files when you want to copy a working starting point before editing it for your own system.

## Files

- `system.yaml`: basic single-system example for `validate` and `oracle`
- `options.yaml`: standard runner options example
- `system_screen.yaml`: screening template used with a ligand library
- `ligand_screen.csv`: small screening CSV used by the screen tutorial
- `parameter_screen.csv`: paired protein-sequence and ligand parameter screen
- `ethanol.sdf` and `ethanol.mol`: structural ligand inputs for conformer reuse and Screen
- `system_custom_ccd.yaml`: custom `ET5` CCD input populated from `ethanol.sdf`
- `ifp_clustering_demo.py`: CPU-only deterministic IFP clustering demonstration
- `system_covalent.yaml`: example of a more specialized system definition
- `system_nucleic_acid.yaml`: protein/DNA/RNA/ligand system with a pocket constraint
- `4HJO.pdb` and `4HJO.cif`: structure fixtures that support examples and manual inspection
- `custom_bias_complexes.yaml`: explicit selectors for preparing a custom supplement
- `bias_matrix.yaml`: the pre-2023-06-01, whole-snapshot,
  pre-2023-06-01-plus-custom, and custom-only comparison runs

## Recommended Use

Start with:

1. `system.yaml`
2. `options.yaml`
3. one of the quick commands from the root `README.md` or `docs/getting-started/quickstart.md`

Then move to `system_screen.yaml` and `ligand_screen.csv` when you want to explore `screen`.
Use `parameter_screen.csv` with repeatable `--map COLUMN=YAML_PATH` arguments when
each row should change several system fields together.
When varying a protein sequence, remove any fixed template `msa` to generate and
reuse one alignment per unique sequence, or add a matching MSA column and map it to
the same protein's `msa` field.

The ligand tutorial uses `ethanol.sdf`, `ethanol.mol`, and
`system_custom_ccd.yaml`. The custom CCD system becomes runnable after `ET5` is
explicitly populated into the cache; it is not expected to work before that setup.

Prepare and compare bias references with:

```bash
cofolder-tools prepare-bias-custom-complexes custom_bias_complexes.yaml --output_root custom-bias
cofolder-tools run-bias-matrix bias_matrix.yaml
```

The copyable no-reference IFP clustering workflow is:

```bash
cofolder screen \
  -s examples/system_screen.yaml \
  -o examples/options.yaml \
  -c examples/ligand_screen.csv \
  --col_id Name \
  --ligand_chain B --smiles_column SMILES \
  --scoring_functions ifp_distance \
  --cluster_ifps \
  -w ./screen_out
```

`ifp_clustering_demo.py` shows the corresponding Python API, including the native
linkage matrix and stable dendrogram leaf labels.

For Boltz-family runners the first row generates a missing fixed-protein MSA and all
later rows reuse it. To skip the MSA server entirely, add an `msa` path to the protein
entry in `system_screen.yaml`; relative paths are resolved from that YAML's directory.

## Relationship To The Docs

- `docs/getting-started/quickstart.md` explains the first command to run
- `docs/tutorials/index.md` maps the written tutorials
- `tutorials/README.md` maps the notebook tutorials

The examples here are intended to stay small, readable, and copyable.
