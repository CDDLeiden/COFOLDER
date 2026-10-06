# COFOLDER Code Audit — 2026-09-07

Scope: PI, IN, and WF requirements from `v1.0-claim-matrix.md`.

Verification: `pytest -q` passed (423 tests); all CLI help commands succeeded.

| Done | ID | Status | Evidence | Notes/Gaps |
|---|---|---|---|---|
| ✓ | PI-01 | [Implemented] | `src/cofolder/cli.py:671`, `README.md:5` | All four commands expose help; invalid subcommand exits 2. |
| ✓ | PI-02 | [Implemented] | `pyproject.toml:68`, `src/cofolder/__main__.py:6` | Both paths call the same CLI `main()`. |
| ✓ | PI-03 | [Implemented] | `src/cofolder/recipes/validate.py:54`, `src/cofolder/modules/runners/contracts.py:57` | Recipes, runner contracts, inputs, analytics, and aggregation are callable without CLI use. |
| ✓ | PI-04 | [Implemented] | `src/cofolder/recipes/bias.py:392`, `docs/user-guide/bias.md:1` | Standalone bias has no runner/options requirement and writes CSV artifacts. |
| ✓ | PI-05 | [Implemented] | `src/cofolder/recipes/validate.py:179`, `src/cofolder/modules/analytics/structure.py:91`, `src/cofolder/modules/analytics/reproduction.py:897` | Optional runner metric groups are explicitly skipped/empty where unsupported. |
| ~ | PI-06 | [Partially Implemented] | `src/cofolder/recipes/screen.py:296` | CSV members run in isolated directories, but output is one consolidated row per compound rather than one normalized row per compound/repeat/model; no SDF-library route. |
| ✓ | PI-07 | [Implemented] | `src/cofolder/recipes/oracle.py:95`, `src/cofolder/recipes/oracle.py:331` | Replaces the SMILES field only, reduces repeats, returns a scalar, and writes `oracle_result.csv`. |
| ~ | PI-08 | [Partially Implemented] | `src/cofolder/modules/runners/__init__.py:22`, `src/cofolder/modules/runners/boltz_runner.py:319`, `src/cofolder/modules/runners/openfold3_runner.py:358` | Four runners are selectable and capabilities are recorded; normalized manifests do not record detected backend version. |
| ✓ | PI-09 | [Implemented] | `src/cofolder/modules/runners/contracts.py:33`, `src/cofolder/recipes/validate.py:487` | Explicit computed/unsupported/missing/failed metric outcomes preserve independent diagnostics. |
| ~ | PI-10 | [Partially Implemented] | `src/cofolder/modules/runners/contracts.py:57`, `src/cofolder/modules/utils/gather.py:108` | Repeat/model/chain data exists, but there is no enforced public schema carrying workflow, system, compound, runner, and entity identifiers in every output. |
| ~ | IN-01 | [Partially Implemented] | `src/cofolder/cli.py:313`, `src/cofolder/modules/input/validation.py:340` | CLI validates paths and runner input before execution, but options YAML has no strict schema and YAML loader logs parse errors instead of consistently raising actionable exceptions. |
| ~ | IN-02 | [Partially Implemented] | `src/cofolder/modules/input/validation.py:82`, `src/cofolder/modules/utils/gather.py:108` | Entity/chain identity is preserved, but validation does not require at least one protein and one ligand. |
| ~ | IN-03 | [Partially Implemented] | `src/cofolder/modules/input/validation.py:119`, `src/cofolder/modules/runners/msa.py:1` | Chain uniqueness and nonempty sequences are checked, but sequence alphabet/content and inaccessible MSA paths are not comprehensively prevalidated for all workflows. |
| ~ | IN-04 | [Partially Implemented] | `src/cofolder/modules/entities/ligand.py:119`, `src/cofolder/recipes/screen.py:769` | SMILES conversion exists for Boltz conformer preparation; invalid SMILES is only explicitly screened in Screen mappings and conversion is not universal. |
| ~ | IN-05 | [Partially Implemented] | `src/cofolder/modules/entities/ligand.py:390`, `src/cofolder/modules/entities/ligand.py:493` | Utilities handle SDF/molblocks, but Screen accepts CSV only and malformed SDF records are filtered rather than represented as source-record failures. |
| ✓ | IN-06 | [Implemented] | `src/cofolder/modules/entities/ligand.py:512`, `src/cofolder/modules/input/validation.py:206` | CCD inputs are accepted and atom references are validated; supported simple inputs can be converted to CCD. |
| ✓ | IN-07 | [Implemented] | `src/cofolder/modules/input/system.py:16`, `src/cofolder/modules/runners/openfold3_runner.py:243` | DNA/RNA are declared and capability-validated before backend launch. |
| ✓ | IN-08 | [Implemented] | `src/cofolder/modules/input/validation.py:340`, `src/cofolder/recipes/screen.py:355` | Screen deep-copies the system, mutates configured fields, and reruns atom/constraint validation. |
| ~ | IN-09 | [Partially Implemented] | `src/cofolder/modules/analytics/reproduction.py:897`, `tests/modules/analytics/test_reproduction.py:177` | PDB/CIF reference structural metrics are supported; no direct reference-complex-to-IFP filtering workflow was found. |
| ✓ | IN-10 | [Implemented] | `src/cofolder/cli.py:176`, `src/cofolder/recipes/oracle.py:280` | Custom pocket references drive coverage, filtering, and Oracle gate metrics without a reference structure. |
| ~ | IN-11 | [Partially Implemented] | `src/cofolder/recipes/screen.py:296` | Tabular CSV SMILES works and row-indexed directories are stable; SDF/molblock libraries and explicit duplicate-ID policy are absent. |
| ✓ | IN-12 | [Implemented] | `src/cofolder/recipes/bias.py:210`, `src/cofolder/modules/analytics/bias.py:459` | Training inputs and provenance fields are validated and retained. |
| ✓ | WF-01 | [Implemented] | `src/cofolder/recipes/validate.py:265`, `src/cofolder/modules/runners/contracts.py:165` | Repeat-specific requests/directories and sample records are retained. |
| ~ | WF-02 | [Partially Implemented] | `src/cofolder/modules/utils/helpers.py:124`, `src/cofolder/recipes/validate.py:198` | Deterministic per-repeat seed generation exists, but effective seeds are logged rather than recorded in normalized output/manifest. |
| ~ | WF-03 | [Partially Implemented] | `src/cofolder/recipes/validate.py:179`, `src/cofolder/modules/analytics/reproduction.py:897` | Inputs enable reference, pocket, and reference-free diagnostics, but no explicit evidence-regime model or output classification exists. |
| ~ | WF-04 | [Partially Implemented] | `src/cofolder/recipes/screen.py:355` | Deep-copy preservation works, but user-provided paths can mutate arbitrary system fields, not only the intended ligand. |
| ✓ | WF-05 | [Implemented] | `src/cofolder/recipes/screen.py:317`, `src/cofolder/recipes/screen.py:429` | Indexed per-member directories prevent collisions; failures are recorded and processing continues. |
| ✓ | WF-06 | [Implemented] | `src/cofolder/recipes/screen.py:153`, `src/cofolder/recipes/validate.py:291` | Shared runner-scoped MSA cache is injected/captured across members and repeats. |
| ✓ | WF-07 | [Implemented] | `src/cofolder/modules/runners/boltz_runner.py:464`, `src/cofolder/modules/runners/openfold3_runner.py:731` | Capability-controlled confidence and affinity groups are normalized separately. |
| ✓ | WF-08 | [Implemented] | `src/cofolder/modules/analytics/reproduction.py:897`, `tests/modules/analytics/test_reproduction.py:484` | RMSD, pose overlap, pocket coverage, and alignment fallback behavior are implemented/tested. |
| ~ | WF-09 | [Partially Implemented] | `src/cofolder/modules/analytics/structure.py:224`, `src/cofolder/modules/analytics/ifp_clustering.py:1` | IFP, custom-pocket filtering, and reference-free clustering exist; direct filtering from a reference-complex interaction pattern is not implemented. |
| ✓ | WF-10 | [Implemented] | `src/cofolder/modules/analytics/structure.py:91`, `src/cofolder/recipes/oracle.py:267` | SASA and normalized heavy-atom SASA are available to workflows and Oracle. |
| ✓ | WF-11 | [Implemented] | `src/cofolder/modules/analytics/bias.py:638`, `docs/user-guide/bias.md:93` | Protein and ECFP ligand proximity retain provenance and are documented as contextual diagnostics. |
| ✓ | WF-12 | [Implemented] | `src/cofolder/modules/utils/gather.py:314`, `tests/modules/utils/test_gather.py:57` | Numeric variation and pairwise structural RMSD are aggregated across runs. |
| ✓ | WF-13 | [Implemented] | `src/cofolder/recipes/screen.py:31`, `src/cofolder/recipes/screen.py:698` | Screen supplies stable requested-metric schema plus statuses/errors and empty unavailable values. |
| ✓ | WF-14 | [Implemented] | `src/cofolder/recipes/oracle.py:98`, `src/cofolder/recipes/oracle.py:484` | Exactly `first`, `mean`, `max`, `min`, and `median` are supported with finite-value validation. |
| ✓ | WF-15 | [Implemented] | `src/cofolder/recipes/oracle.py:85`, `src/cofolder/recipes/oracle.py:386` | Python API supports custom callables and weighted multi-metric scores. |
| ✓ | WF-16 | [Implemented] | `src/cofolder/recipes/oracle.py:34`, `src/cofolder/recipes/oracle.py:408` | Gates use the same collected/aggregated proposed-ligand run metrics. |
| ✓ | WF-17 | [Implemented] | `src/cofolder/recipes/oracle.py:60`, `src/cofolder/recipes/oracle.py:360` | Downweight, fixed-penalty, non-binder actions and audit columns are written. |
| ✓ | WF-18 | [Implemented] | `src/cofolder/recipes/oracle.py:408`, `docs/user-guide/oracle.md:151` | Thresholds are user-configured; documentation rejects universal criteria. |
| ~ | WF-19 | [Partially Implemented] | `src/cofolder/recipes/screen.py:31`, `src/cofolder/modules/utils/helpers.py:259` | Some schemas and affinity conversions are explicit, but no complete public metric catalog defines all units, directionality, classes, and missing-value conventions. |
| ~ | WF-20 | [Partially Implemented] | `src/cofolder/recipes/screen.py:429`, `src/cofolder/modules/runners/validators.py:1` | Screen preserves per-compound failures and runner-boundary validation is specific; Validate/post-hoc failures are not universally wrapped with full stage/system/repeat/model context. |
