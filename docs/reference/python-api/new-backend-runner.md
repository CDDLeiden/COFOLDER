# Constructing a new backend runner

This reference explains in words what a new prediction backend must provide to work
with COFOLDER's Python recipes. The current contract is defined by
`modules/runners/base.py`, `modules/runners/contracts.py`, and
`modules/runners/validators.py`. Boltz2 is the most complete in-tree example;
OpenFold3 is useful when the backend requires a different input and output format.

## The runner's job

COFOLDER owns workflow validation, repeats, seed planning, partial-failure handling,
cross-backend analytics, and public result serialization. A runner owns five things:

1. declaring what its backend supports;
2. validating and translating backend options;
3. checking that the external backend is installed and ready;
4. preparing and executing one repeat; and
5. normalizing backend output into COFOLDER's private runner bundle.

The central rule is that recipes must not learn the backend's filenames, command
syntax, or native metric schema. That translation ends inside `runner.run()`.

## 1. Add a discoverable module

Create a non-private module under `src/cofolder/modules/runners/`, normally named
`<backend>_runner.py`. Export exactly one ready-to-use object named `RUNNER`:

```python
from cofolder.modules.runners.base import BaseRunner

class ExampleRunner(BaseRunner):
    name = "example"
    # implement the contract below

RUNNER = ExampleRunner()
```

`discover_runners()` scans every non-underscore module in that package, imports it,
and registers any module-level `RUNNER` by `RUNNER.name`. A new runner therefore
appears in `list_runner_names()` and `get_runner("example")` without editing a fixed
registry. Keep module import side effects light: discovery imports runner modules
even when another backend is selected.

## 2. Declare identity, capabilities, and input limits

Set these class attributes:

- `name`: COFOLDER's stable selector, used by Python and CLI callers.
- `backend_name`: scientific/backend identity written to provenance.
- `backend_distribution`: installed Python distribution used for version detection.
- `capabilities`: runner-produced metric groups, currently such values as
  `confidence_metrics`, `affinity_metrics`, and `affinity_metrics_ext`.
- `input_capabilities`: a `RunnerInputCapabilities` value listing faithfully
  supported entity types and constraint types, plus any pocket/force restrictions.
- `options_schema`: the `RunnerOptionsSchema` used to validate this backend's YAML.
- `ligand_preparation_capabilities`: whether native SMILES works and which conformer
  modes (`2D`, `3D`, `sdf`) are actually supported.
- `supports_msa_reuse`: opt in only after implementing the reuse methods correctly.

Do not claim a capability because a similarly named backend offers it. These
declarations control early validation and whether missing result payloads are treated
as unsupported or as backend failures.

## 3. Define the options contract

Create a `RunnerOptionsSchema` describing shared runtime fields and backend-specific
runner fields. `BaseRunner._load_typed_options()` delegates to
`load_runner_options()` and returns immutable, validated `RunnerOptions`.

Implement `load_options(options_path)`. Return `RunnerOptions` directly when that is
enough, or a small runner-owned wrapper when nested backend settings need helper
methods, as OpenFold3 does. Keep original types through command/config construction:
integers such as `0` and `1` must not be treated as booleans. Secrets may reach the
real subprocess input but must be redacted from logs, returned completed-process
objects, and serialized failures.

If the backend has an executable, subcommand, cache, sample-count, or extra-argument
setting, model it explicitly in the typed options/runtime layer rather than searching
arbitrary YAML dictionaries during shared recipe execution.

## 4. Implement availability and setup checks

Implement `check_availability() -> tuple[bool, str | None]`. It must not import a
large backend package merely to discover whether it exists. Prefer installed
distribution metadata and `check_distribution_available()`. Detect incompatible
packages when two distributions provide the same command or module, as the Boltz
adapters do.

Return an actionable installation/setup message on failure. If weights, caches, or
databases are required, distinguish “package missing” from “package installed but
setup incomplete.” `ensure_available()` turns the result into the runtime error used
by recipes.

The default `detect_backend_identity()` records a normalized installed version. Set
`backend_distribution` correctly or override detection when distribution metadata
cannot represent the backend honestly.

## 5. Validate system inputs

Normally inherit `BaseRunner.validate_system()`. It invokes shared source-aware
validation using `input_capabilities`. Override it only for extra native rules that
cannot be expressed by the capability object, and preserve the shared validation
result including its chain-identity mapping.

Validation runs twice around preparation: first with relaxed atom-name checking, then
with full checking. Preparation must not turn a valid system into one the backend
cannot faithfully consume.

## 6. Prepare the system when necessary

Override `prepare_system()` only when the backend needs conversion or enrichment.
It receives the parsed `System`, options object, shared raw directory, conformer
request, optional SDF path, and logger. Return `RunnerPreparationResult` containing:

- the effective system object;
- the effective options object;
- any user-facing warnings; and
- `RunnerRuntime` with authoritative sample count, cache path, and model name.

Do not execute the prediction here. Preparation is the place for ligand/config
translation that must happen before the generated system input is written.

## 7. Describe seeds, models, and optional MSA reuse

The default `resolve_effective_seed()` uses COFOLDER's seed unchanged. Override it
only if the backend requires a transformation, returning explicit
`backend_adjusted` provenance and a reason.

`execution_models(options_obj)` tells Screen how many model/sample terminal slots to
expect before it launches any compound. Override it when the backend has axes that
the default diffusion-sample implementation cannot describe.

For MSA reuse, implement all three behaviors consistently:

- `msa_reuse_settings()` returns non-secret settings that affect cache identity;
- `capture_reusable_msas()` copies valid backend-generated MSAs to the shared cache;
- `inject_reusable_msas()` adds matching cached MSAs to a later system.

Never silently reuse an MSA created under incompatible settings.

## 8. Execute exactly one request

Implement `run(request: RunnerExecutionRequest)`. One call represents one repeat and
one effective seed. Use the request's generated `system_path`, prepared `options_obj`,
`repeat_dir`, backend and seed provenance, chain identities, logger, and timing
context. Write backend-native files beneath `request.repeat_dir`; do not write into
another repeat or directly invent workflow-level public outputs.

Invoke subprocesses with an argument sequence rather than a shell string. Preserve
paths with spaces, check the exit status, retain useful backend output, and redact
credentials from every diagnostic representation.

## 9. Normalize the backend output

Before returning, create this private layout:

```text
<request.repeat_dir>/normalized/
  manifest.json
  system_metrics.csv
  chain_metrics.csv
  structures/
    ... .cif, .mmcif, or .pdb
```

At minimum, `system_metrics.csv` contains `cif_file`, `model_name`, `repeat`, and
`diffusion_sample`. `chain_metrics.csv` additionally contains `conf_chain_id`.
Use the request's chain-identity mapping rather than inferring biological chain IDs
from DataFrame row order. Apply `attach_runner_provenance()` and
`attach_sample_provenance()` so every row/sample agrees with the request.

Copy or convert predicted structures into `normalized/structures/`. Every referenced
structure must exist, and `sample_records` must agree with the files and declared
sample count.

Write `manifest.json` with runner, backend, seed, paths, capabilities, repeat,
runtime, samples, warnings, and companion-artifact information. The manifest is
debug/provenance material; recipes consume the typed result rather than reparsing
loose runtime context.

## 10. Declare metric outcomes explicitly

For every relevant metric group, return a `RunnerMetricOutcome` with one of:

- `computed`: valid required columns/artifacts exist;
- `unsupported`: the backend cannot provide the group;
- `missing`: it claims support but produced no payload;
- `failed`: production was attempted and failed; or
- `not_requested`: supported but not requested for this execution.

List required normalized columns and artifacts in the outcome. Shared validation no
longer guesses outcome state from CSV shape. A required malformed/missing/failed
group invalidates that repeat; unsupported optional groups can continue with explicit
unavailable values.

## 11. Return the typed result

Return `RunnerExecutionResult` with all normalized paths, capabilities, warnings,
`RunnerRuntime`, sample records, metric outcomes, companion artifacts, the request's
chain identities, backend identity, and seed provenance. Built-in runners also use
`build_runner_public_records()` to produce typed per-repeat records.

The recipe immediately passes this result to `validate_runner_bundle()`. That
validator checks layout, manifest, CSV schemas, provenance, identities, structures,
sample cardinality, artifacts, records, capabilities, and metric outcomes before any
shared analytics run.

## 12. Package and test the backend separately

Add the backend dependency to its own optional installation extra. If the backend
conflicts with another implementation, document and test isolated environments.
Package any small adapter resources, but do not place model weights or caches in the
wheel.

Tests should cover:

- discovery and `get_runner()`;
- dependency missing, installed, wrong-version, conflicting, and incomplete-setup
  cases;
- valid and invalid options, including zero, one, floats, booleans, absent fields,
  extra arguments, paths with spaces, and secret redaction;
- every supported and rejected entity/constraint type;
- preparation and conformer capabilities;
- deterministic seed/model/sample planning;
- backend command or native-config construction;
- normalization from representative backend fixtures;
- required CSV columns, structures, sample records, identities, provenance,
  companion artifacts, and every metric-outcome state;
- partial-repeat behavior through Validate;
- Screen cardinality and Oracle metric use; and
- a real clean-environment acceptance run for Validate, Screen, and Oracle before
  claiming release support.

## Component replacement boundary

A runner can be replaced without changing recipes when both adapters accept the
shared validated system contract and return equivalent normalized semantics. Native
backend commands and files may be completely different. What must remain stable is
the runner name/capability truth, typed runtime, seed and backend provenance, chain
identity, planned sample cardinality, explicit metric states, normalized bundle, and
typed failures.
