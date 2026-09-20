# Issue 16 tutorial acceptance

Acceptance was completed on 2026-09-19 against the COFOLDER 1.0 working tree.
This record covers the tutorial-specific release gate; the broader release evidence
remains in the [v1.0 claim matrix](v1.0-claim-matrix.md).

## Environment

- Python 3.12.4 in a fresh virtual environment created with access to the host's
  CUDA packages, followed by `pip install --no-deps -e
  '.[analysis,tutorials,boltz2]'` and a successful `pip check`.
- COFOLDER 1.0.0, Boltz 2.2.1, PyTorch 2.10.0, RDKit 2025.9.4, and marimo
  0.6.26.
- NVIDIA L40S (46,068 MiB), NVIDIA driver 610.57.04.
- Seed 20260919 for runner-backed acceptance.

The first Validate attempt used an already occupied GPU and stopped with CUDA
out-of-memory before producing an accepted result. The recorded run was repeated on
an available GPU and completed successfully. This is an environmental retry, not a
discarded workflow failure.

## Commands

The commands below were run from the repository root. `ACCEPTANCE_ROOT` was
`/tmp/cofolder-issue16-acceptance`; it can be changed to any writable directory.

```bash
export ACCEPTANCE_ROOT=/tmp/cofolder-issue16-acceptance
export CUDA_VISIBLE_DEVICES=3

cofolder validate -s examples/system.yaml -o examples/options.yaml \
  -w "$ACCEPTANCE_ROOT/validate-gpu3" --runner boltz2 --repeats 1 \
  --seed 20260919 --scoring_functions confidence_metrics

cofolder validate -s examples/system_screen.yaml -o examples/options.yaml \
  --conformers sdf --sdf_file examples/ethanol.sdf \
  -w "$ACCEPTANCE_ROOT/validate-sdf" --runner boltz2 --seed 20260919 \
  --scoring_functions confidence_metrics

cofolder screen -s examples/system_screen.yaml -o examples/options.yaml \
  -c examples/ligand_screen.csv --ligand_chain B --smiles_column SMILES \
  --col_id Name --merge_data pIC50 -w "$ACCEPTANCE_ROOT/screen-csv" \
  --runner boltz2 --seed 20260919 --scoring_functions confidence_metrics

cofolder screen -s examples/system_screen.yaml -o examples/options.yaml \
  -c examples/ethanol.sdf --ligand_chain B --id_property ID \
  -w "$ACCEPTANCE_ROOT/screen-sdf" --runner boltz2 --seed 20260919 \
  --scoring_functions confidence_metrics

cofolder screen -s examples/system_screen.yaml -o examples/options.yaml \
  -c examples/ethanol.mol --ligand_chain B \
  -w "$ACCEPTANCE_ROOT/screen-mol" --runner boltz2 --seed 20260919 \
  --scoring_functions confidence_metrics

cofolder oracle -s examples/system.yaml -o examples/options.yaml \
  --input_smiles CCO --output_metric system__confidence_score --aggregate first \
  -w "$ACCEPTANCE_ROOT/oracle" --runner boltz2 --seed 20260919 \
  --scoring_functions confidence_metrics

python -m cofolder.tools.cli setup-boltz2-cache --cache-path cache/.boltz
python -m cofolder.tools.cli populate-ccd-cache --sdf examples/ethanol.sdf \
  --property-id ID --cache-path cache/.boltz --on-conflict overwrite
cofolder validate -s examples/system_custom_ccd.yaml -o examples/options.yaml \
  -w "$ACCEPTANCE_ROOT/validate-ccd" --runner boltz2 --seed 20260919 \
  --scoring_functions confidence_metrics
```

The Bias tutorial's shared fixture builder was executed with its generated command:

```text
cofolder bias --system_path $ACCEPTANCE_ROOT/bias-inputs/bias_system.yaml \
  --wrk_dir $ACCEPTANCE_ROOT/bias \
  --custom_protein_reference_path $ACCEPTANCE_ROOT/bias-inputs/custom_protein.csv \
  --custom_ligand_reference_path $ACCEPTANCE_ROOT/bias-inputs/custom_ligand.csv
```

## Results

Every row below has a `results/manifest.json` with `status: success`, a public
`results/metrics.csv`, and at least one computed metric.

| Route | Public metric rows | Computed rows |
|---|---:|---:|
| Bias | 8 | 4 |
| Validate, system-YAML SMILES | 180 | 74 |
| Validate, external SDF conformer | 180 | 74 |
| Screen, CSV | 800 | 120 |
| Screen, SDF | 200 | 30 |
| Screen, MOL | 200 | 30 |
| Oracle | 2 | 2 |
| Validate, custom preprocessed CCD | 140 | 64 |

Expected structure artifacts were present for each runner-backed workflow. The CCD
population command created `cache/.boltz/mols/ET5.pkl`, which the custom-CCD system
then consumed.

## Automated gates

- `ruff check src tests scripts tutorials examples`: passed.
- `python scripts/run_test_lane.py contracts-tutorial`: 27 passed on Python
  3.12.4. The CI lane is matrixed over Python 3.11 and 3.12.
- Focused ligand, input, packaging, example, tutorial, and Screen regression set:
  117 passed.
- Aggregation, public-contract, Validate, Screen, ligand, and cache regressions:
  102 passed.
- `python scripts/run_test_lane.py all`: 668 passed.
- `mkdocs build --strict`: passed.

The tutorial suite executes the complete offline Bias workflow, validates all
packaged ligand fixtures, preflights tutorial-derived Validate/Screen/Oracle routes
without services, and rejects references to removed or fictional public functions.
