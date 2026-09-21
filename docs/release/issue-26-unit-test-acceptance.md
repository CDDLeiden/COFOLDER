# Issue 26 unit-test and backend acceptance

Issue [#26](https://github.com/CDDLeiden/COFOLDER/issues/26) was verified on
2026-09-21 against commit `1395cb50e62757772dc8ce0e71990f1baada0e26`.
The latest consolidated publication checklist in the issue is the controlling
scope; obsolete `boltz_lab`, Predict/Evaluate, global 60% coverage, and automatic
backend-CI criteria are not release gates.

## Automated release environment

A new Python 3.12 virtual environment was created under `/tmp`, and the current
checkout was installed with the declared test, analysis, tutorial, and development
extras. The final environment used Python 3.12.4, COFOLDER 1.0.0, pytest 9.1.1,
pytest-cov 7.1.0, RDKit 2026.3.6, and pandas 3.0.6. `pip check` passed.

The initial `.[test,analysis]` run produced 704 passes and one collection failure
because `tests/test_tutorials.py` intentionally imports Marimo from the separate
`tutorials` extra. After installing that declared extra, the complete supported
suite passed; no test was skipped or weakened.

Commands and results:

```text
python -m pytest tests/ -v
705 passed, 1 warning

python -m pytest tests/ --cov=cofolder --cov-report=term-missing --cov-report=xml
705 passed, 1 warning; TOTAL 78%

python scripts/run_test_lane.py all
705 passed, 1 warning

python -m ruff check src tests scripts tutorials examples
All checks passed
```

The warning is the known MDAnalysis 2.10 deprecation notice. Coverage was reviewed
by publication-facing area rather than gated by an arbitrary global percentage:

| Area | Coverage |
|---|---:|
| CLI | 95% |
| Validate | 88% |
| Screen | 89% |
| Oracle | 84% |
| Input validation | 83% |
| Ligand utilities | 82% |
| CCD cache tool | 77% |

## Boltz-2 acceptance environment

A new isolated conda prefix was cloned from the previously validated Boltz-2
environment and the release commit was reinstalled into it with the
`acceptance`, `analysis`, and `boltz2` extras. `pip check` passed. The environment
used Python 3.12.14, COFOLDER 1.0.0, Boltz 2.2.1, PyTorch 2.14.0+cu130, an NVIDIA
L40S, and seed 20260921. GPU 3 was selected explicitly. The validated shared
Boltz cache was reused; no backend package namespaces were mixed.

The following real CLI routes were executed with one diffusion sample and
`confidence_metrics`. `ACCEPTANCE_ROOT` was
`/tmp/cofolder-issue26-boltz2.peRbCw`:

```text
cofolder validate -s src/cofolder/acceptance/data/system.yaml \
  -o src/cofolder/acceptance/data/options_acceptance.yaml \
  -w $ACCEPTANCE_ROOT/runs/validate-smiles --runner boltz2 --repeats 1 \
  --seed 20260921 --scoring_functions confidence_metrics

cofolder validate -s examples/system_screen.yaml \
  -o src/cofolder/acceptance/data/options_acceptance.yaml \
  --conformers sdf --sdf_file examples/ethanol.sdf \
  -w $ACCEPTANCE_ROOT/runs/validate-sdf --runner boltz2 --seed 20260921 \
  --scoring_functions confidence_metrics

cofolder screen -s examples/system_screen.yaml \
  -o src/cofolder/acceptance/data/options_acceptance.yaml \
  -c $ACCEPTANCE_ROOT/one-ligand.csv --ligand_chain B \
  --smiles_column SMILES --col_id Name --merge_data pIC50 \
  -w $ACCEPTANCE_ROOT/runs/screen-csv --runner boltz2 --seed 20260921 \
  --scoring_functions confidence_metrics

cofolder screen -s $ACCEPTANCE_ROOT/runs/screen-csv/compound_000001/raw/screen_system.yaml \
  -o src/cofolder/acceptance/data/options_acceptance.yaml \
  -c examples/ethanol.sdf --ligand_chain B --id_property ID \
  -w $ACCEPTANCE_ROOT/runs/screen-sdf --runner boltz2 --seed 20260921 \
  --scoring_functions confidence_metrics

cofolder screen -s $ACCEPTANCE_ROOT/runs/screen-csv/compound_000001/raw/screen_system.yaml \
  -o src/cofolder/acceptance/data/options_acceptance.yaml \
  -c examples/ethanol.mol --ligand_chain B \
  -w $ACCEPTANCE_ROOT/runs/screen-mol --runner boltz2 --seed 20260921 \
  --scoring_functions confidence_metrics

cofolder oracle -s $ACCEPTANCE_ROOT/runs/screen-csv/compound_000001/raw/screen_system.yaml \
  -o src/cofolder/acceptance/data/options_acceptance.yaml \
  --input_mol_file examples/ethanol.mol \
  --output_metric system__confidence_score --aggregate first \
  -w $ACCEPTANCE_ROOT/runs/oracle-mol --runner boltz2 --seed 20260921 \
  --scoring_functions confidence_metrics

cofolder-tools populate-ccd-cache --sdf examples/ethanol.sdf --property-id ID \
  --cache-path cache/.boltz --on-conflict overwrite
cofolder validate -s examples/system_custom_ccd.yaml \
  -o src/cofolder/acceptance/data/options_acceptance.yaml \
  -w $ACCEPTANCE_ROOT/runs/validate-ccd --runner boltz2 --seed 20260921 \
  --scoring_functions confidence_metrics
```

The later commands reused the Screen-generated system containing the exact staged
protein MSA from the successful CSV route. This prevented unnecessary repeated
MSA-service calls for the SDF, MOL, and Oracle checks.

## Results and issue closure evidence

| Route | Status | Public metric records | Successful executions | Failures |
|---|---|---:|---:|---:|
| Validate, system-YAML SMILES | success | 55 | 1 | 0 |
| Validate, external SDF conformer | success | 55 | 1 | 0 |
| Screen, CSV | success | 86 | 1 | 0 |
| Screen, SDF | success | 86 | 1 | 0 |
| Screen, MOL | success | 86 | 1 | 0 |
| Oracle, MOL file | success | 2 | 1 | 0 |
| Validate, preprocessed CCD | success | 45 | 1 | 0 |

Every route produced a successful manifest and its documented primary output.
CCD population created a readable `ET5.pkl`, and the subsequent real Validate run
consumed that CCD. Automated regressions additionally cover direct single/multiple
SMILES-to-SDF conversion, invalid-SMILES rejection without partial output, cache
identifier/conflict/failure behavior, Oracle MOL and SDF ingestion, canonical
protein-sequence transfer to the runner request, real mmCIF structural diagnostics,
and successful plus contradictory Screen/Oracle CLI dispatch.

These results satisfy every unchecked item in the consolidated issue #26 release
checklist.
