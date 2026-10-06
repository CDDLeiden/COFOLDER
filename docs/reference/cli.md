# Command-line reference

The installed `cofolder` command and `python -m cofolder` use the same parser,
defaults, validation, output contracts, and exit statuses. The console form is used
below for brevity.

```text
cofolder [-h] [-v] {bias,validate,screen,oracle} ...
```

Use `cofolder COMMAND --help` for the complete parser-generated option list for the
installed version. Unknown commands and invalid arguments exit with status `2`.
Input-validation failures also exit `2`; execution failures exit `1`; successful
execution and successful `--preflight_only` validation exit `0`.

## Shared runner-backed arguments

`validate`, `screen`, and `oracle` share these option groups:

| Group | Arguments and defaults |
| --- | --- |
| Required inputs | `-s/--system_path`, `-o/--options_path` |
| Execution | `-w/--wrk_dir` (current directory), `--runner` (`boltz2`), `--repeats` (`1`), `--seed` (generated when omitted), `--preflight_only` |
| Scoring | `--scoring_functions` (all supported groups except `ifp_prolif`), `--assess_robustness/--no-assess_robustness` (enabled), `--assess_bias/--no-assess_bias` (disabled) |
| Structure evidence | `--reference_path`, `--pocket_coverage_reference`, `--reproduction_metrics` |
| Ligand preparation | `--conformers {2D,3D,sdf}`, `--sdf_file` |
| Diagnostics | bias reference/source/cache controls and robustness thresholds; inspect `--help` for the installed names and defaults |
| Logging | `-d/--debug` |

Workflow-owned values such as runner selection, repeats, seeds, scoring groups,
evidence, and output paths belong on the command line. Backend runtime values belong
in the separate options YAML.

## `bias`

Runs reference-overlap diagnostics without a co-folding backend or options YAML.
`--system_path` is required. Supply one supported public source-bundle route, custom
reference route, or build route. `--preflight_only` validates the plan without
searches or materialization. Results are written below
`<wrk_dir>/results/bias_train/`.

```bash
cofolder bias --system_path system.yaml --wrk_dir bias-output \
  --custom_protein_reference_path custom-proteins.csv \
  --custom_ligand_reference_path custom-ligands.csv
```

See [Bias](../user-guide/bias.md) and
[Bias training data](../user-guide/bias-training-data.md) for the mutually exclusive
source modes and provenance controls.

## `validate`

Runs one system for one or more repeats. In addition to the shared arguments it
accepts no workflow-specific required argument. Public results are written below
`<wrk_dir>/results/`, with runner-owned artifacts below `<wrk_dir>/raw/`.

```bash
cofolder validate -s examples/system.yaml -o examples/options.yaml \
  --runner boltz2 --repeats 1 --seed 20260919 \
  --scoring_functions confidence_metrics -w validate-output
```

## `screen`

Adds `-c/--library`, `--library_format {csv,sdf,mol}`, `--ligand_chain`,
`--smiles_column`, `--col_id`, `--id_property`, `--duplicate_id_policy`, and
`--merge_data`. CSV requires ID and SMILES columns. SDF/MOL uses an identifier
property and deterministic source-record fallbacks. IFP filter and clustering flags
are documented in [Screen](../user-guide/screen.md).

```bash
cofolder screen -s examples/system_screen.yaml -o examples/options.yaml \
  -c examples/ligand_screen.csv --ligand_chain B \
  --col_id Name --smiles_column SMILES --preflight_only
```

## `oracle`

Requires exactly one of `--input_smiles` and `--input_mol_file`, plus
`--output_metric`. `--aggregate` is one of `first`, `mean`, `max`, `min`, or
`median` and defaults to `first`. `--ligand_chain` is optional only when the target
ligand is unambiguous.

```bash
cofolder oracle -s examples/system.yaml -o examples/options.yaml \
  --input_smiles CCO --output_metric system__confidence_score \
  --scoring_functions confidence_metrics --aggregate first --preflight_only
```

The requested scalar must be produced by an enabled scoring group and supported by
the selected runner. Missing, nonnumeric, empty, or vector-valued selectors fail
explicitly. See [Oracle](../user-guide/oracle.md) for composite Python scores and
gates.
