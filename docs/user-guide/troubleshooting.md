# Troubleshooting

Run the relevant command with `--preflight_only` first. It validates inputs, runner
availability, capabilities, metric dependencies, and planned outputs without model
execution, downloads, or external searches.

| Symptom | Meaning and resolution |
| --- | --- |
| Selected runner is unavailable | Install exactly the matching backend extra and run its setup steps from the [installation guide](../getting-started/installation.md). Keep conflicting Boltz package lines in separate environments. |
| Unsupported entity, constraint, or metric group | The selected runner does not advertise that capability. Use the [capability matrix](../tutorials/runners.md#supported-runner-capabilities) or change the input/request; COFOLDER does not simulate unsupported outputs. |
| External MSA request fails | Supply a valid local `.a3m` or paired-MSA CSV on the protein, verify service/network settings, or retry after the service recovers. Screen reuses successfully resolved fixed-protein MSAs. |
| Invalid SMILES or malformed SDF/MOL | Correct the source record. Screen retains malformed structure records as source-mapped failures while continuing valid records. Validate requires an SDF conformer to match the selected system ligand. |
| Missing or conflicting ligand representation | Each system ligand must provide exactly one of `smiles`, `ccd`, or `ccd_codes`; SDF/MOL paths are workflow inputs, not system-ligand keys. |
| Missing custom CCD | Initialize and populate the same cache configured in the options YAML, then verify the requested code and atom names. |
| Some Screen members fail | Inspect `results/failures.csv`, `results/records.jsonl`, and the member work directory. Independent failures do not remove successful members. |
| Metric is empty or unsupported | Inspect the metric `status` and `reason`. Confirm the scoring group was requested, its evidence was supplied, and the runner supports it. Empty is not numerical zero. |
| Optional analysis import fails | Install `cofolder[analysis]`; the error identifies the feature that requested the dependency. |
| CUDA out of memory | Use an idle GPU, lower `runtime.diffusion_samples`, reduce workflow concurrency, or select a smaller backend configuration. Retrying an environmental OOM should be recorded. |
| Cache or model setup fails | Verify cache permissions and free space, then use the installed `cofolder-tools` setup command. Do not mix backend caches or assume a partial download is complete. |

Exit status `2` denotes command or input validation, `1` denotes workflow execution
failure, and `0` denotes success. Tracebacks and backend diagnostics remain in logs;
public failure records contain sanitized, machine-readable stages and error codes.
