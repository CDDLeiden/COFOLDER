# Public output contract

COFOLDER workflows publish schema version `1.0.0` under `<wrk_dir>/results`.
The authoritative output is `records.jsonl`; each line is one independently
validatable success, metric, or failure record. The files `successes.csv`,
`metrics.csv`, and `failures.csv` are flattened views of the same records.
`manifest.json` identifies the invocation, evidence inputs, requested metrics,
artifacts, completion status, and record counts.

The former wide public files (`system_metrics.csv`, `chain_metrics.csv`,
`robustness_metrics.csv`, `screen_results.csv`,
`screen_results_with_scores.csv`, and `oracle_result.csv`) are not produced.
Runner-native normalized CSV files may still appear below `raw/`; they are
private backend/debugging artifacts and are not a supported public interface.

## Identity

Every record contains `workflow`, `run_id`, and `system_id`. Runner-backed
workflows also contain `runner_id`. Screen records carry the exact source
compound identifier. Repeat/model records carry `repeat_id`, `model_id`, and
`sample_id`, while chain records additionally carry `entity_id`, `entity_type`,
and `chain_id`. Pair metrics use `related_entity_id` and `related_chain_id`
instead of encoding chain identifiers in the metric name.

`record_id` is a deterministic UUID derived from the complete identity plus the
metric/statistic or failure stage/code. Duplicate record IDs are rejected.

## Metric states and evidence

Metric records use one of these states:

- `computed`: a finite, catalog-valid value is present;
- `missing`: the metric was applicable and requested, but no value was produced;
- `unsupported`: runner capability or available evidence is incompatible;
- `failed`: computation was attempted and failed;
- `not_requested`: the metric belongs to the workflow profile but was not requested.

Non-computed records always serialize `value` as JSON `null` and include a
machine-readable reason where applicable. Each metric also records its class,
unit, optimization direction, statistic, evidence regime, and evidence
provenance. Reference-coordinate comparisons use `reference_structure`, custom
pocket comparisons use `custom_pocket`, and model/structure/internal-consistency
diagnostics use `reference_free`.

The canonical metric definitions are exposed through
`cofolder.modules.contracts.METRIC_CATALOG`. Unregistered metrics cannot enter
public output.

## Failures

Failures are ordinary public records and preserve all identity known at the
failure boundary. Stages distinguish input validation, preparation, backend
execution, output validation, gathering, analytics, aggregation, and
serialization. Independent repeat/model/compound failures do not discard
successful records. A manifest is `partial` when both successes and failures
exist and `failed` when no usable result exists.

Tracebacks remain in logs. Public failures contain the exception type, stable
error code, actionable message, retryability, and sanitized details.

## Python API

The contract types, validators, catalog, DataFrame adapters, aggregation helper,
and serializer are importable from `cofolder.modules.contracts`. Recipe
constructor signatures and Oracle's scalar return value are unchanged.

