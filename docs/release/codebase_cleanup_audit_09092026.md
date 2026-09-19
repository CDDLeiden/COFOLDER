# COFOLDER initial-release cleanup and usability audit

Audit date: 2026-09-09. Audited commit: `e9642f25eb358d44e2a97edcb4203d96d485c14b` (`feature/release_features`).

## Release assessment

The codebase is not yet ready to finalize for publication. The existing architecture has useful boundaries and extensive contract tests, but passing tests conceal a reproducible Boltz argument-serialization defect. Package metadata and distribution contents also need release work. The broken optional UI is not part of the supported initial-release surface under GitHub issue 24 and the user's 2026-09-10 scope decision; align or exclude that development material instead of repairing it for release. Address correctness, security and release-surface work before undertaking broad module moves.

The findings section records the original audit; implementation progress is tracked later in this file. The [v1.0 claim matrix](v1.0-claim-matrix.md) is mandatory. Backward compatibility is unnecessary for this initial release, but **every feature or interface removal needs explicit permission**, including capabilities absent from the matrix. The user approved R01-R10 explicitly on 2026-09-09 and approved withdrawal of the UI from the supported/shipped release surface (R11) on 2026-09-10. Those approvals do not waive the recorded prerequisites or capability safeguards. Issue 24's historical references to v0.1.0 do not change this audit's `1.0.0` release target because the user instructed that all remaining proposed changes be implemented.

Priority definitions: **P1** = resolve before publication; **P2** = recommended cleanup or usability improvement before freezing the code; **P3** = optional simplification. Findings distinguish reproduced defects, static evidence, and design recommendations. No recommendation authorizes weakening a claim or deleting its test coverage.

### Verification performed

| Check | Observed result | Interpretation / limit |
|---|---|---|
| Repository inventory | 66 Python files / 29,313 lines under `src/cofolder`; 47 Python files under `tests`; inspected CLI/UI, recipes, input, entities, runners, contracts, analytics, utilities, packaging, scripts, examples, tutorials, acceptance resources, and release notes | Screening and targeted tracing, not a proof that every function is correct or unused |
| `pytest -q` | **520 passed in 18.28 s**, Python 3.12.4 | Existing configured environment; includes mocked workflow/runner tests. Does not establish actual backend inference or clean dependency installation |
| Console/module smoke checks | Both entry points returned 0 with identical help for root and all four commands; unknown command returned 2 with identical stderr; version outputs matched | PI-01/02 smoke evidence; version is still `0.1.0` |
| Focused Boltz serialization probe | Schema-valid `diffusion_samples: 1`, `devices: 1`, `recycling_steps: 0` produced bare sample/device flags and omitted recycling steps | Reproduced without backend execution; see C01 |
| UI-generated options probes | Current loader rejected old `wrapper/options` YAML and rejected UI Validate-style `runner.cache` | Reproduced without starting Streamlit/inference; launch defects additionally traced statically |
| `ruff check --no-cache --output-format concise src tests scripts tutorials examples` | 301 findings under ambient settings | No repository-owned Ruff configuration found; this count is environment-dependent, not an agreed project standard |
| `ruff check --isolated --no-cache --select E4,E7,E9,F --output-format concise src tests scripts tutorials examples` | **15 findings** | Explicit reproducible baseline; no automatic fixes applied |
| Offline artifact build | Wheel and sdist built successfully with installed `setuptools.build_meta.build_sdist` / `build_wheel`, in a temporary copy of tracked files | `build` frontend unavailable; no isolated build-dependency resolution |
| Wheel install | Installed with `pip --no-index --no-deps` into a temporary `venv --system-site-packages`; confirmed import from its site-packages; console/module help succeeded from outside source tree | Uses host dependencies; **not** PK-01/05 clean-environment acceptance |
| Artifact inspection | Wheel: 84 files, including all 12 acceptance YAML/CSV resources. Sdist: 109 entries, only four top-level test modules | See C05; archive contents explicitly inspected |
| Example validation | `system.yaml`, `system_screen.yaml`, `system_nucleic_acid.yaml` and `options.yaml` passed relevant Boltz2 input/options validation | Covalent example required CCD atom evidence absent from the configured cache; prerequisite gap, not established malformed input |

Commands and probes used temporary directories; no backends, models, or dependencies were downloaded. No GPU inference, external MSA service, live bias-data retrieval, interactive browser session, clean dependency-resolution matrix, or documentation build was performed. Documentation work remains deferred; documentation-related claims remain release gates.

The pre-existing untracked `code_audit_09092026.md` was preserved. It reports 42/42 PI/IN/WF claims met at this commit. This audit independently reproduced its test count but found a gap in runner execution evidence (C01); that earlier conclusion is not sufficient to certify actual execution. Its scope also excluded packaging and UI. Neither earlier audits nor claim-matrix checkboxes were edited.

## Capability and claim guardrails

| Capability to preserve | Claim mapping | Code / existing safeguards | Audit implication |
|---|---|---|---|
| Four commands and matching console/module behavior | PI-01/02, DOC-01/03 | `cli.py`, `__main__.py`, `tests/test_cli.py` | Preserve command functionality while simplifying presentation |
| CLI-independent workflow and component APIs | PI-03, DOC-04 | `recipes/`, `modules/`, recipe and contract tests | Module moves must not strand callers; removal of existing import paths requires permission |
| Standalone bias; Validate diagnostics; complete Screen; scalar Oracle | PI-04–07 | Four recipe modules and `tests/recipes/` | No workflow removal; retain custom Oracle scoring |
| Boltz1, Boltz2, Boltz Community, OpenFold3 and capability/provenance handling | PI-08–10, PK-03/04, DOC-08/09 | Runner registry, runner/public contracts, runner tests | Boltz1 is a required backend, not removable legacy code |
| Separate system/options YAML; sequences, entity identity, DNA/RNA and constraints | IN-01–03/07/08, DOC-05/06 | Input config/validation/system, runner validators | Keep two-file input support, early capability validation, MSA and constraint integrity |
| SMILES, MOL/molblock, SDF and reusable CCD chemistry | IN-04–06/11, DOC-07 | Input ligand/library and entity ligand utilities/tests | Keep conversions, stereochemistry, atom names, source identity and duplicate-ID policies |
| Reference complexes, custom pockets, training references | IN-09/10/12, DOC-10/11 | Reproduction, reference IFP, bias and input tests | Preserve reference-free, pocket-only and reference-complex modes |
| Repeats, seeds, isolated library execution, MSA reuse | WF-01–06, DOC-12 | Recipes, runner contracts/MSA, seed planning tests | C01 directly undermines execution confidence; isolation/seed tests remain essential |
| Confidence, binding, affinity, structural metrics, IFP, SASA, bias and robustness | WF-07–13/18/19, DOC-14 | Analytics, metric catalog, aggregation, workflow tests | Keep all metric classes, provenance, directionality, configurable thresholds and missing-value distinctions |
| Five Oracle reducers, custom scalar, gates and penalties | WF-14–17, DOC-13 | Oracle recipe and its tests | No reducer/gate removal or universal scientific threshold defaults |
| Typed failures with execution identity and stage | WF-20, DOC-18 | Public models/adapters/serialization and workflow tests | Simplified errors must retain machine-readable diagnostic detail |
| Installation, metadata, artifacts, resources and reproducibility | PK-01–08 | `pyproject.toml`, packaging and acceptance tests | C04–C07/C12/C15; complete fresh-environment acceptance later |
| Examples, executable tutorials and complete documentation | DOC-01–20, especially DOC-02/15–20 | Examples/tutorials, MkDocs and prior documentation audit | Code-side resource/test fixes included; full documentation compliance explicitly deferred, not waived |

### Additional capabilities outside explicit matrix requirements

Retain these unless separately approved for removal: automatic in-package runner discovery; old Python aliases and utility wrappers; the misspelled bias CLI alias; archived shell/config entry paths; setup/download/vendor helper scripts and temporary-conda wrapper; bias landscape plots, mixed/same-type pair datasets and enrichment/checkpoint outputs; affinity censoring/correlation/plot utilities; cache-population helpers; and legacy vector-only fingerprint clustering support. Some contribute to broader claims even if their exact interface is not named. “No internal callers found” does not establish absence of external users or permission to delete. The Streamlit interface is governed separately by approved R11: it may remain in the repository only as clearly warned, unsupported development material and must not be installed or advertised as part of the initial release.

## Findings and suggested fixes

### C01 — P1 · Correct Boltz numeric option serialization

**Reproduced defect.** `src/cofolder/modules/input/command.py:147` compares values using membership in tuples containing `False` and `True`. Python treats `0 == False` and `1 == True`. `BoltzRunner.load_options` converts validated typed settings to this old representation (`src/cofolder/modules/runners/boltz_runner.py:240`), and `run` uses it at line 295. The schema permits these values (`src/cofolder/modules/input/config.py:350`). This affects all three Boltz-family adapters.

Reproduction: a valid options file with `runtime.diffusion_samples: 1`, `runner.devices: 1`, and `runner.recycling_steps: 0` becomes:

```text
boltz predict system.yaml --out_dir out --seed 0 --diffusion_samples --devices --model boltz2
```

**Fix:** dispatch on actual value types; serialize integers/floats with explicit values, booleans as flags according to backend semantics, and omit only intended absent values. Make the typed options object the internal source of truth rather than round-tripping through list-of-dictionaries configuration. Preserve any currently exposed legacy `Command` API until R05 is approved.

**Acceptance:** assert complete argv for 0, 1, 2, floating-point values, true, false and absent fields through `BoltzRunner.load_options` and each Boltz-family run adapter. Test default sample count 1 explicitly and verify numeric zero is transmitted. Existing `tests/modules/input/test_command.py:75` only checks selected command tokens and misses this failure. Later backend acceptance must verify execution and resulting sample cardinality.

**Claims:** IN-01, PI-08, WF-01/02/07. No feature removal required.

### C02 — P1 · Remove the unsupported UI from the initial-release surface

**Reproduced schema failures and static launch defects.** `src/cofolder/ui/app.py:463` reads/saves old `wrapper/options` YAML. `run_screen` (606) and `run_oracle` (683) also write it. `run_validate` (561) writes the new outer shape but puts obsolete and runtime keys from `DEFAULT_OPTIONS` (33) into `runner`. Validate/Oracle argv use removed `-b` instead of `-o`; Oracle supplies neither required query input nor `--output_metric` (`src/cofolder/cli.py:647`). All run functions overwrite the selected options file and ignore `save_config` failure. They use literal `python`, not the current interpreter.

**Scope decision.** GitHub issue 24 records that the manuscript and Supporting Information require the current Bias, Validate, Screen and Oracle CLI/Python workflows but do not claim a graphical interface. The issue also records that the UI's Predict/Evaluate model and command construction are obsolete. The user confirmed on 2026-09-10 that the UI will not be shipped and that the remaining audit work should be implemented.

**Fix:** remove UI launch instructions from supported release documentation and direct users to `bias`, `validate`, `screen`, `oracle`, and the Python API. Retain `src/cofolder/ui`, `run_ui.sh`, `.streamlit`, and their shallow tests only as clearly labelled unsupported development material, or move them outside the installable package; do not add an installed UI launcher. Explicitly exclude the UI package and UI-only dependencies from the wheel, supported extras, release test lanes and release claims. Record future UI work in a separate post-release issue based on current recipe/runner contracts.

**Acceptance:** supported documentation contains no unqualified UI launch path; any retained source UI and launcher display an explicit unsupported/development-only warning; built wheel metadata, import inventory, entry points and supported extras contain no UI surface or Streamlit dependency; core release tests do not import Streamlit; all four current CLI commands and Python recipes remain documented and tested. Do not describe issue 24 as completed feature work; it is closed or tracked as not planned for the initial release.

**Claims:** UI itself is additional functionality and is not required by the manuscript, Supporting Information or claim matrix. Underlying PI-03/04–08, IN-01 and WF-20 remain applicable through the CLI and Python API. Supported-surface withdrawal is explicitly authorized by R11; deletion of unrelated workflow capabilities is not.

### C03 — P1 · Redact credentials from command logs and recorded failures

**Static defect.** The options schema accepts `msa_server_password` and `api_key_value` (`src/cofolder/modules/input/config.py:380`), and command assembly forwards these as argv. `src/cofolder/modules/runners/boltz_runner.py:130` logs the entire command at INFO; its failure path includes the original command in `CalledProcessError` (163). Those messages can reach persisted logs and failure records. No real credentials were used or exposed during this audit.

**Fix:** maintain the real argv only for execution and a separately redacted display form for logs/failure details. Apply the same policy to OpenFold3 command reporting (`openfold3_runner.py:240`) and nested option representations. Preserve actionable non-secret context and record credential presence without values.

**Acceptance:** dummy password/token sent unchanged to a mocked process, absent from captured logs and serialized failures on success and failure. **Claims:** PK-08, WF-20. No feature removal.

### C04 — P1 · Finish version and licence metadata

**Confirmed release mismatch.** `pyproject.toml:7`, `src/cofolder/__init__.py:1` and CLI version identify `0.1.0`; `tests/test_packaging.py:47` explicitly enforces that version. PK-02 requires `1.0.0`. The wheel includes MIT `LICENSE`, but has neither `License` nor `License-Expression` metadata because the project omits a licence declaration.

**Fix:** select one authoritative package-version source and align metadata/runtime/CLI at release finalization; update the test to validate consistency and intended release version. Declare MIT in package metadata using a form supported by the chosen build-tool floor. Keep public output schema version independent of package-version implementation. Release tag and documentation alignment are later explicit gates; do not create a tag during cleanup.

**Acceptance:** built wheel metadata, installed runtime, both CLI version outputs and final release identifier agree at `1.0.0`; installed metadata identifies MIT and includes the licence text. **Claims:** PK-02/06, DOC-19.

### C05 — P1 · Make release artifact contents deliberate and usable

**Confirmed artifact gaps.** `pyproject.toml:72` discovers only `cofolder*`; line 76 includes acceptance YAML/CSV data. Wheel and sdist omit `THIRD_PARTY_SOFTWARE.md`. The wheel has no root examples, tutorials, setup scripts or launcher. The sdist contains only `tests/test_cli.py`, `test_packaging.py`, `test_tutorials.py`, and `test_ui.py`, without `conftest.py`, nested tests, referenced tutorials, `.streamlit` config or launcher. Consequently the shipped tests are not a self-contained suite. `run_ui.sh:26` requires a source-relative app path. This does not mean every repository file belongs in a wheel; it means ownership and distribution intent need explicit decisions.

**Fix:** include licence/attribution in both artifacts; include the full supported test suite and required repository resources in the sdist. Keep notebooks/source documentation in the source release. Place the small runnable example set in package resources, expose a Python resource locator, and provide an additive command to copy examples to a user workspace. Provide installed entry points for backend setup and user-facing data preparation, backed by packaged implementations; retain source wrappers unless removal is approved. Do not provide an installed UI entry point or package the unsupported `cofolder.ui` implementation in the wheel. Include the PDB/CIF references needed by examples in the distributable example set, not model weights or caches.

**Acceptance:** build wheel from an unpacked sdist; assert an explicit artifact inventory and UI exclusion; install outside the checkout, locate/copy fixtures, and run lightweight examples plus setup/data-tool `--help` or dry-run paths. Execute the supported sdist tests with the declared test extras. No development-tree fallback via `sys.path` may be required. **Claims:** PK-01/05–08, DOC-15–17/20.

### C06 — P1 · Define and verify the supported dependency environment

**Static inconsistencies / validation gap.** README supports Python 3.11/3.12 but `pyproject.toml:14` admits every Python >=3.11; Boltz1/2 extras silently omit their backend dependency on >=3.13. Boltz2 has no major-version constraint despite its runtime check requiring major 2 (`boltz2_runner.py:42`). Boltz Community references a moving Git branch (`pyproject.toml:39`). The base install includes optional scientific stacks such as ProLIF/pdb2pqr and plotting packages; the code already lazily loads some corresponding functionality (`structure.py:32`, `reference_ifp.py:320`).

**Fix:** keep the advertised support target at 3.11/3.12 until broader support is demonstrated. Prefer explicit explanatory refusal for unsupported backend/Python combinations rather than an apparently successful backend-extra install with no backend. Bound the Boltz2 dependency to the supported major line; pin the Community source to a tested immutable revision or tested released distribution. Derive tested dependency floors from APIs actually imported, not arbitrary blanket pins. Prepare separate validated environments for the four backends; their shared `boltz` namespace conflicts are intentional and must remain detectable.

Separate core, optional interaction/plotting analysis, testing and development tooling dependency groups. Do not publish a supported UI dependency group for the initial release. Preserve all analytical capabilities and provide installation guidance when an optional group is absent. **Changing base-install availability or narrowing advertised Python support is an interface change: R08 requires permission before applying those portions.** Bounds within the already-supported backend major and immutable provenance can be fixed independently.

**Acceptance:** fresh 3.11/3.12 base installs run backend-independent paths; a dedicated full-test environment includes analysis/tutorial dependencies but not UI dependencies; each backend environment passes availability and smoke checks; omitted analysis extras fail with explicit installation guidance. Test unsupported Python behavior. Index-specific publication eligibility of the direct Git requirement was not verified and remains a separate release check. **Claims:** PK-01/03/04/08, PI-08/09.

### C07 — P2 · Put reproducible quality checks in the repository

**Confirmed tooling/test gap.** No tracked CI workflow or lint/test configuration was found. Ruff/Black live under the docs extra (`pyproject.toml:56`); the test extra does not include UI dependencies, although `tests/test_ui.py:18` imports them unconditionally. Baseline Ruff reports 15 issues, including a duplicate `logging` import (`input/system.py:12`), unused imports, unused locals in plots/structure, and compound one-line statements. `tests/tmp_test.py` is a tracked zero-byte placeholder, not a missing feature.

**Fix:** add project-owned Ruff configuration starting with the explicit baseline used above and expand rules deliberately; put development tools in a development extra. Define supported core, contracts, tutorial, packaging and isolated backend lanes. UI tests are development-only under R11 and must not be imported or counted by a supported release lane. Add CI checks for lint, unit/contracts, artifact inventory and installed CLI smoke tests. Add meaningful tests for C01 and the C02 release-surface exclusions. Remove the empty placeholder and unused imports/locals after checking side effects; do not apply unsafe bulk fixes or delete supported tests to obtain a clean run.

**Acceptance:** a fresh documented test environment can collect and run the supported suite without Streamlit; baseline lint passes; tests exercise argument values, schema round-trips, partial failures, wheel resources and UI exclusion. Any retained UI tests are clearly outside the release gate. **Claims:** PK-01/03/05, DOC-15/20; cross-cutting support for PI/IN/WF.

### C08 — P2 · Give legacy interfaces a finite, permission-gated removal list

**Confirmed compatibility surface.** `legacy/README.md:3` labels historical environments/configs/wrappers. `runners/base.py:26` defines aliases re-exported by `runners/__init__.py:14`; `utils/helpers.py:40` declares legacy `set_dir`, and line 309 onwards contains forwarding wrappers. `cli.py:130` retains `--asess_bias`. `scripts/build_bias_training_data.py:2` is explicitly a compatibility wrapper. `screen.py:1345` accepts vector-only IFPs by manufacturing `_legacy` receptor/residue identities.

**Fix:** use R01–R06 below to approve exact removals. Update every internal caller and tutorial/example invocation to the retained canonical implementation before removing an approved entry point. Move tests for underlying capabilities to their canonical module rather than discarding them with wrappers.

**Acceptance:** repository searches show no dangling references to each approved removal; supported API and workflow tests still pass. Distance-contact clustering must continue with real identities, not disappear with a legacy fallback. **Claims:** PI-03, WF-09; interface spellings/archive paths are mostly additional capabilities. This finding itself authorizes no removal.

### C09 — P2 · Keep typed configuration through runner execution

**Structural recommendation grounded in C01.** `BoltzRunner.load_options` (240) validates `RunnerOptions` then rebuilds old `Command.options`; `run` (289) copies that untyped structure. `Command` still lives under generic input despite encoding Boltz executable semantics. `modules/entities/ligand.py:122` consumes the command-shaped options during conformer preparation, mixing chemistry and backend cache management.

**Fix:** separate pure validated input, ligand preparation and backend argv construction. Let each runner translate typed settings at its execution boundary; carry cache/sample identity through existing `RunnerRuntime` contracts. Isolate CCD storage from chemistry conversion. Retain all SMILES/SDF/MOL/CCD and conformer routes. Keep old public facades until separately approved under R05/R09.

**Acceptance:** the same typed request supplies capability validation, preparation, argv and provenance; no lossy conversion through legacy YAML occurs internally. Existing chemical-identity/constraint and seed/model tests plus C01 numeric cases pass. **Claims:** PI-03/08, IN-01/04–08, WF-01/02/04.

### C10 — P2 · Split workflows and bias analytics by responsibility

**Structural recommendation.** `screen.py` is 1,947 lines: `_run_impl` (512), `_write_public_results` (999), clustering (1331), filtering (1568), and metric selection (1813) combine distinct concerns. `validate.py:337` holds orchestration and diagnostic collection. `analytics/bias.py` is 3,834 lines and includes reference loading (425), network retrieval (652), similarity (884), pair enrichment (1798), dataset materialization (3021), plots (3086) and scoring (3575). Size is supporting evidence; the mixed responsibilities are the reason to split.

**Fix:** extract workflow planning/execution, diagnostic collection, screen postprocessing, and public result assembly into cohesive internal modules. Split bias reference loading, similarity computation, remote enrichment, pair datasets and artifacts; keep `apply_bias_metrics` as a stable orchestrator. Preserve the network-enrichment, mixed/same-type datasets, plotting and checkpoints even where not individually named in the matrix. Keep one canonical metric catalog and common failure/identity construction.

**Acceptance:** fixture-based before/after comparisons preserve execution counts, identities, seeds, metric states/values, provenance, gates, cluster/filter decisions and all current artifact families. Independent analytics imports must not execute a workflow or network request. **Claims:** PI-03–07/09/10, IN-12, WF-03/05/07–20. No feature removal.

### C11 — P2 · Make the common CLI journey easier without reducing controls

**Usability recommendation.** Help is 93 lines for Validate, 138 for Screen and 110 for Oracle. `cli.py:75` puts backend, bias-building, evidence and execution controls into one common argument section. `--ligand_training_data_path` changes from an input to an output depending on build mode (194). Screen always requires `--ligand_chain` (526), while Oracle can infer a unique target (674). This creates unnecessary decisions for single-ligand cases.

**Current journey:** assemble two YAMLs, inspect a long help page, choose chains/columns and several reference paths, then discover incompatibilities during setup. **Proposed journey:** copy schema-valid packaged examples, select a backend, provide the required system/options, and see grouped help plus a preflight summary. Infer the ligand only when unique and validate before execution; ambiguity still requires explicit selection.

**Fix:** group help into required inputs, execution, optional evidence/scoring, bias preparation, and advanced controls. Preserve every accepted option and default unless separately approved. Add a preflight-only path using the same validation as execution, reporting selected runner/capabilities, planned repeats/samples, external MSA needs and result locations without invoking services. Add explicit bias-reference input/output fields while retaining the old dual-purpose option until a separately enumerated removal is approved; reject contradictory settings clearly. Keep separate system/options YAML support (IN-01).

**Acceptance:** users can reach a valid single-ligand dry run from bundled examples with no backend/service execution; ambiguous selectors and missing prerequisites fail early; CLI/module parity and all advanced scientific controls remain. **Claims:** PI-01/02/08, IN-01/10–12, WF-18/20. No existing option removal is proposed here beyond R04.

### C12 — P2 · Consolidate executable and cache discovery

**Static duplicate implementation.** MMseqs lookup appears in `scripts/fetch_bias_training_data.py:28`, `analytics/bias_training.py:22`, and `analytics/build_bias_training_data.py:547`. Repository-relative parent calculations differ; the packaged builder uses `parents[1]`, which points under analytics' package ancestry rather than repository root, whereas the runner helper uses `parents[4]`. Existing home/PATH routes can mask this discrepancy.

**Fix:** one shared resolver with explicit argument first, then `COFOLDER_MMSEQS_BIN`, user vendor directory, supported source vendor locations when applicable, then PATH. Report rejected non-executable candidates and the chosen path. Keep all currently supported lookup routes; move maintainer wrappers to orchestration only. Ensure the same resolver works in a wheel installation without assuming a checkout.

**Acceptance:** parameterized tests cover precedence, paths with spaces, missing/non-executable files, source/vendor and installed layouts; setup helpers and analytics choose the same executable. **Claims:** PK-08, IN-12, WF-11. If implementation proposes dropping a lookup route, add it to the removal register first.

### C13 — P2 · Simplify progress and results discovery

**Usability recommendation.** `utils/log.py:15` includes filename, line number and function in every console message and uses stdout (81). CLI completion messages (`cli.py:503`, 697, 713) do not consistently identify the canonical manifest and result tables. Public serialization already provides a coherent bundle (`contracts/serialization.py:196`), so this is a presentation issue rather than justification to remove files.

**Current journey:** read detailed implementation logs and browse output folders to find the relevant tables or Oracle value. **Proposed journey:** concise stage progress followed by success/failure counts, resolved output directory, primary table, manifest, and the Oracle scalar when applicable. Retain detailed logs for debugging.

**Fix:** use a concise console formatter and detailed file/debug formatter; expose a common completion summary for the four supported commands and Python recipes. Keep every current output, identity and machine-readable failure state. Do not silently change stdout/stderr contracts; offer explicit machine-readable/quiet modes additively before considering any stream change.

**Acceptance:** success, partial Screen failure, total failure and unavailable Oracle metrics each produce a useful summary pointing to existing files; test paths with spaces. **Claims:** PI-07/09/10, WF-05/13/17/20, DOC-09/18.

### C14 — P2 · Remove test-framework knowledge from production logging after approval

**Confirmed coupling.** `utils/log.py:91` traverses private `logging._handlerList` and reattaches `_pytest.logging` handlers by module name. `tests/modules/utils/test_log.py:63` explicitly locks in recovery after tests clear root handlers. This is compatibility for a test setup, not scientific functionality.

**Fix:** preserve ordinary application handlers and configure only COFOLDER-owned handlers; keep pytest capture setup/teardown in test fixtures using supported fixture behavior. R07 covers removal of the currently tested special recovery behavior.

**Acceptance:** repeated configuration does not duplicate logs; file/console capture works; unrelated application handlers survive; production source no longer references pytest/private handler registries after approval. **Claims:** no specific scientific claim; supports WF-20 diagnostics. Do not drop logging coverage.

### C15 — P2 · Correct the cache-initialization helper without losing CCD utilities

**Static unreachable branch.** `modules/entities/ligand.py:746` creates `<cache>/mols` before checking whether `<cache>` exists. Its subsequent download branch cannot run after successful directory creation. `input/command.py:14` implements the intended cache download via a minimal prediction in a fixed `./tmp` directory. These utilities are callable Python functionality, despite limited repository callers.

**Fix:** validate source input before creating cache artifacts, check the actual required cache components, and isolate explicit cache setup from reusable CCD conversion. Use a unique temporary workspace for any setup prediction. Keep conversion reusable without automatically invoking inference; because this changes the documented automatic setup behavior, that change requires R10 approval. Independently fix validation order and workspace collisions while preserving the intended capability.

**Acceptance:** mocked empty/partial/complete caches, invalid SDF and concurrent setup requests behave predictably; conversion preserves names/connectivity/stereochemistry; no accidental model execution occurs during tests. **Claims:** IN-06, PI-03, PK-08. The R10 change is approved but remains gated on its explicit setup replacement and retained conversion contract.

## Approved feature/interface removals — all prerequisite-gated

These are finite proposals, not execution instructions. Even replacement by an equivalent interface counts as a removal. No matrix-required capability is proposed for deletion.

| ID | Exact proposed removal | Reason / retained alternative | Claim impact and required safeguards | Status |
|---|---|---|---|---|
| R01 | Tracked `legacy/` tree: README, environment/requirements, old YAMLs, `run.sh`, `screening.sh` | Remove archived first-run paths; retain current CLI, examples and Git history | Historical invocation/config reproduction disappears; not itself required. Check CCD/input/backend capabilities remain and record later doc-link repairs | **Approved; migration recorded in S0** |
| R02 | `RunnerPreparation`, `RunnerRequest`, `RunnerResult` aliases from runner base and package exports | Keep `RunnerPreparationResult`, `RunnerExecutionRequest`, `RunnerExecutionResult` | Old Python imports disappear; PI-03 preserved through canonical contracts; update alias tests/callers without losing execution coverage | **Approved; migration recorded in S0** |
| R03 | `helpers.set_dir` and forwarding exports `read_yaml`, `read_csv`, `read_sdf`, `delete_last_line`, `parse_censored_affinity`, `remove_censored_affinity`, `strip_censoring_signs`, `prepare_affinity_dataframe`, `convert_boltz_affinity_to_ic50`, `calculate_affinity_correlations`, `plot_affinity_correlation` | Keep directory creation and implementations in read/write/dataset/stats/plots modules | Removes import paths, not I/O, affinity conversion, censoring, correlation or plotting capabilities; PI-03/WF-19 retained; move tests to canonical APIs | **Approved; migration recorded in S0** |
| R04 | `--asess_bias` and generated boolean-negation alias | Keep correctly spelled `--assess_bias` / `--no-assess_bias` | Only misspelled invocation support disappears; preserve bias flags/functionality and CLI/module parity | **Approved; ready** |
| R05 | Old public `Command(options={"options": [...]})` / old-YAML loading interface and recursive mutation/query methods, after internal migration | Use typed `RunnerOptions` plus backend command builder | Non-CLI old Python configuration support disappears; PI-03 and IN-01 retain current preparation/execution APIs. Enumerate any additional caller-facing replacement signatures in the implementation plan | **Approved; gated by typed migration** |
| R06 | Screen's vector-only clustering fallback that fabricates `_legacy` identities (`screen.py:1359`) | Canonical identity-bearing interaction fingerprints, including distance taxonomy | Old vector-only clustering inputs disappear; retain distance clustering and pocket bitstring/residue-list input functionality. Prove current runs generate canonical fingerprints before removal (WF-09) | **Approved; gated by identity proof** |
| R07 | Automatic pytest-handler recovery after arbitrary root-handler clearing | Test fixtures manage capture; application logging retained | Removes explicitly tested test-harness behavior, no scientific feature; keep normal logging/capture tests | **Approved; gated by capture proof** |
| R08 | Always-installed optional interaction/plotting stacks; any withdrawal of metadata-advertised Python >=3.13 support | Explicit analysis extras; tested 3.11/3.12 backend environments | Changes installation availability, not analytical capabilities. PK-03/04 must pass; apply only the exact S0 dependency/Python proposal | **Approved; exact proposal recorded in S0** |
| R09 | Exact old import/module paths listed in the S0 path-migration table, and the source builder compatibility wrapper | Keep canonical APIs and installed tool entry points | PI-03 requires usable Python APIs; each path waits for its replacement and migrated callers | **Approved; exact list recorded in S0** |
| R10 | Implicit cache-download/prediction behavior of the CCD cache utility | Explicit setup followed by backend-independent conversion | Setup convenience disappears; CCD conversion (IN-06) and optional setup remain available. Current branch is unreachable, but documented intended behavior must still be disclosed | **Approved; gated by replacement contract** |
| R11 | Streamlit UI as an installed, documented or supported initial-release interface, including any planned installed UI launcher and supported UI dependency/test lane | Retain the four current CLI commands and Python recipes. UI source may remain only as explicitly unsupported development material and must be excluded from the wheel | Removes no manuscript/claim-matrix capability. Documentation, artifacts, dependencies and tests must consistently describe the same release surface | **Approved; ready for scope alignment** |

Not recommended for removal: any backend (including Boltz1); standalone bias; extra bias datasets/plots/enrichment; affinity analyses; SDF/MOL/CCD routes; IFP taxonomies, filters or clustering; setup/download/vendor tooling; maintainer temporary-conda wrapper; custom scores; five Oracle reducers; gates or penalties; current public output files. These are retained in the cleanup design. R11 withdraws only the UI's supported and installed release status; it does not authorize removal of any underlying recipe or runner capability.

Do not automatically delete dependency adapters because they contain the word “legacy.” `structure.py:37` supports two pdb2pqr programmatic entry points; preserve both until a tested supported-version policy proves one unnecessary. Likewise, old public-file cleanup in `contracts/serialization.py:290`, `screen.py:1319` and `recipes/bias.py:524` can protect against stale mixed outputs on reruns; its removal needs an explicit output-lifecycle decision, not a keyword search.

## Target structure and dependency rules

Prefer responsibility extraction over a wholesale cosmetic rename of `modules/`. The first stage can preserve public imports while moving private implementation into the following layout. Add module-level responsibility docstrings during later implementation; full API/user documentation remains separate.

```text
src/cofolder/
  cli.py                       # parsing, request assembly, presentation
  ui/                          # unsupported development material; excluded from wheel
  recipes/                     # public Validate/Screen/Oracle/Bias orchestration
    _execution.py              # shared planning and lifecycle operations
    _diagnostics.py            # applicable diagnostics and failure collection
    _screen_postprocess.py     # ranking/filtering/clustering orchestration
    _results.py                # workflow-to-public-record assembly
  modules/
    contracts/                 # identities, metric definitions, validation, serialization
    input/                     # YAML schemas, system/ligand/library validation
    entities/                  # chemistry representations and conversions
    runners/                   # registry, typed backend adapters, MSA reuse
      _boltz_command.py        # typed Boltz argv serialization
    analytics/
      bias/                    # reference loading, similarity, enrichment, datasets, artifacts
      ...                      # structure, reproduction, IFP, affinity, aggregation
    utils/                     # small I/O, timing, logging, executable discovery
  tools/                       # packaged setup/data-fetch/build implementations
  resources/examples/          # distributable runnable inputs and small references
  acceptance/                  # keep backend acceptance apps/resources
tests/                         # mirrors responsibilities; optional lanes explicit
  scripts/                       # thin source launchers; remove only the R09-listed builder after migration
examples/                      # source examples, synchronized with packaged fixtures
tutorials/                     # executable teaching material
docs/                          # documentation and audit reports
```

| Current responsibility | Proposed home / public-surface treatment |
|---|---|
| `input/command.py` backend serialization | `runners/_boltz_command.py`; retain `Command` until the approved R05 readiness gate passes |
| `entities/ligand.py` chemistry versus backend cache operations | Keep chemistry in entities; extract cache/setup orchestration into runners/tools; apply only the exact approved R09/R10 migrations |
| Large Screen/Validate orchestration | Internal recipe execution/diagnostic/postprocessing/results modules; retain public recipe classes and signatures |
| `analytics/bias.py` many responsibilities | `analytics/bias/` package with explicit submodules and re-exported public orchestrator; preserve private helpers used by current tests until tests cover the new boundaries |
| `analytics/build_bias_training_data.py`, fetch/setup scripts | Packaged tools implementations; migrate and remove only the R09-listed module/source builder paths |
| Three MMseqs resolvers | One utility resolver shared by tools and analytics |
| Analytical aggregation in `utils/gather.py` | Analytics aggregation implementation; retain import facade until an exact path-removal proposal is approved |
| Streamlit UI | Outside the initial-release architecture under issue 24/R11; retain only as warned source development material, exclude it from installed artifacts and do not refactor it during this plan |

Release dependency direction: CLI → recipes → input/runners/analytics → foundational contracts and small I/O utilities. The unsupported development UI may call the Python API but is not part of this release dependency graph. Foundation model/metric definitions must not import workflows or backends. DataFrame adapters/serialization may depend on pandas, but model-only imports should not eagerly load every adapter; the current contracts package eagerly imports adapters (`contracts/__init__.py:1`). Analytics should not route substantive work through a catch-all helpers module. Tools reuse library functions; library code must not rely on a source script path. Backend modules may expose capabilities without importing the external backend package. Preserve runner discovery behavior; do not replace it with a fixed list and silently remove extension support.

## Implementation plan and persistent progress record

This section is the source of truth for implementing this audit. It is deliberately
kept in this file so another session can resume without reconstructing intent from
Git history or chat. Work is organized into dependency-ordered phases; the existing
`S<n>` identifiers are retained as stable phase/task IDs so the historical work log
does not need to be rewritten. Update the control block, phase ledger, task ledger
and work log as work proceeds; do not maintain a separate unchecked checklist.

### Plan control block

```yaml
plan_version: 3
overall_status: COMPLETE
active_phase: S8
active_task: S8.5
last_completed_task: S8.5
next_action: "S8 handoff complete; final documentation alignment, release tagging and publication remain separate authorized work."
blocked_on: []
implementation_base_commit: e9642f25eb358d44e2a97edcb4203d96d485c14b
last_updated: 2026-09-19
updated_by: Codex
```

Allowed `overall_status` values are `NOT_STARTED`, `IN_PROGRESS`, `BLOCKED`,
`VERIFYING`, and `COMPLETE`. Allowed task/phase values are `NOT_STARTED`,
`IN_PROGRESS`, `BLOCKED`, `DONE`, `SKIPPED`, and `NOT_APPLICABLE`. There must be at
most one `IN_PROGRESS` task. `SKIPPED` requires a reason in the work log;
`NOT_APPLICABLE` requires either a rejected removal proposal or superseding evidence.

### Resume and update protocol

At the start of every implementation session:

1. Read the control block, the active task row, its acceptance criteria, the removal
   decision ledger, and the latest work-log entry. Do not simply start at the first
   unchecked item.
2. Run `git status --short`, record the current commit, and distinguish prior task
   changes from unrelated user changes. Never discard or overwrite unrelated work.
3. If a task is already `IN_PROGRESS`, inspect its listed files and diffs, then
   continue from its `resume_from` note. Otherwise set the next task and its phase to
   `IN_PROGRESS`, set `overall_status: IN_PROGRESS`, and update `next_action` before
   editing implementation files.
4. Work on one task at a time. Add or tighten a failing regression test before the
   implementation when practical. Do not combine permission-gated removal with a
   behavior-preserving task.
5. Before stopping, update the task's status, evidence, changed paths and
   `resume_from`. If its acceptance checks pass, mark it `DONE`, set
   `last_completed_task`, and point `active_task`/`next_action` at the next eligible
   task. If interrupted, leave it `IN_PROGRESS` and make `resume_from` exact enough
   to continue. If blocked, identify the missing decision, dependency or resource.
6. A phase becomes `DONE` only after all of its required tasks are `DONE`,
   `NOT_APPLICABLE`, or explicitly `SKIPPED`, and its exit check is recorded. Do not
   infer completion from a commit message or a green broad test run.

For each completed task, append one work-log row containing the date, task ID, commit
or `working-tree`, paths changed, commands run with results, decisions made, and the
next task. Keep verbose command output outside this document; record the exact command
and a concise result here. Commits are recommended at phase boundaries, but this plan
does not authorize commits, tags, publication, deletion of user work, network/service
use, or feature/interface removal beyond approved R01-R11 and their gates.

### Phase ledger

| Phase | Scope | Depends on | Status | Exit evidence |
|---|---|---|---|---|
| S0 | Baseline, ownership and release-scope freeze | None | **DONE** | 2026-09-09 baseline/map/R01-R10 records plus 2026-09-10 issue-24/R11 scope decision and Python API flow references |
| S1 | Command correctness and secret-safe failures (C01, C03, narrow C09 seam) | S0 | **DONE** | 51 focused tests and 530-test full suite pass; all seven console/module parity cases match |
| S2 | Release-surface alignment and API orientation (C02, R11) | S0 | **DONE** | Supported docs point to four CLI/Python workflows; clean 84-file wheel excludes UI/Streamlit; 525 supported tests plus 6 development-only UI checks pass |
| S3 | Version, licence and deliberate artifacts (C04, C05) | S1, S2 | **DONE** | Clean 228-file sdist produced the 104-file wheel; installed-outside-checkout smoke and all 551 supported sdist tests passed without Streamlit |
| S4 | Reproducible quality and supported environments (C06, C07) | S3 | **DONE** | Repository-owned lint/test/package lanes pass in documented environments |
| S5 | Additive CLI/Python usability and discovery improvements (C11–C13) | S3, S4 | **DONE** | 71 focused usability tests and 641-test full lane pass; eight console/module journeys match, all four preflights are no-service/non-mutating, and artifact/docs/lint gates pass |
| S6 | Behavior-preserving responsibility extraction (C09, C10) | S1–S5 | **DONE** | Pre-S6/current semantic digests match for all four workflows; import, focused, full, lint, docs and artifact gates pass |
| S7 | Permission-gated cleanup and cache/logging changes (C08, C14, C15; R01–R11) | S0 decision gate, relevant earlier phases | **DONE** | Per-ID reference/import checks, 227 focused and 657 full supported tests, lint, strict docs and clean artifact verification pass at `9d40404` |
| S8 | Release-candidate verification and claim reconciliation | S1–S7 | **DONE** | Fresh source/artifact/Python 3.11/3.12/four-backend gates pass; claims and removals are reconciled below, with tag/publication and final documentation alignment handed off explicitly |

S0.4 documentation work was completed while S1 code was active; only one
implementation task may be active at a time. S5 and non-conflicting S7 preparation
may be developed before S6. S8 cannot be marked complete while a
required task is merely `BLOCKED`. Documentation completion, release tagging and
publication remain separate work and are not authorized by this plan.

### S0 baseline evidence

Captured on 2026-09-09 before implementation changes:

- Repository: `feature/release_features` at
  `e9642f25eb358d44e2a97edcb4203d96d485c14b`; the only worktree entries were the
  pre-existing untracked `docs/release/code_audit_09092026.md` and this audit. The
  former had SHA-256
  `48bdb9ce09541a86ab5b6636e467ff6b6a66869a97b445c8a17bcabfdd9c82ab`
  and is outside the implementation scope.
- Environment: CPython 3.12.4 at `/home/remco/apps/miniconda3/bin/python` on
  Linux 5.14.0-687.44.1.el9_8.x86_64 with glibc 2.34.
- `pytest -q`: **520 passed in 18.72 s**, matching the audited count.
- `ruff check --isolated --no-cache --select E4,E7,E9,F --output-format concise
  src tests scripts tutorials examples`: **15 findings**, matching the audited
  baseline (command exit 1 because findings exist); no fixes were applied.
- Console `/home/remco/apps/miniconda3/bin/cofolder` and
  `python -m cofolder` had byte-identical stdout/stderr and equal exit codes for root
  help, help for `validate`, `screen`, `oracle`, and `bias`, version, and an unknown
  command. Help/version cases returned 0; the unknown command returned 2. Version
  remained `cofolder 0.1.0`. Help lengths remained root 19, Validate 93, Screen 138,
  Oracle 110 and Bias 61 lines.

No baseline drift was found. These checks used the existing environment and did not
perform backend inference, dependency installation, network/service access or source
changes; they do not replace the later clean-environment and backend release gates.

### S0 finding ownership and capability map

The refreshed baseline is the audited commit, so no C01-C15 finding is stale. Status
below distinguishes evidence already reproduced/confirmed by the audit from design
recommendations and checks still requiring later environments.

| Finding | Current classification | Owning code | Existing/required test ownership | Claims and retained safeguards |
|---|---|---|---|---|
| C01 | **Reproduced defect; current** | `input/config.py`, `input/command.py`, all Boltz runner adapters | `test_command.py`, `test_boltz_runner.py`, Boltz1/Community adapter tests; add complete typed argv cases | IN-01, PI-08, WF-01/02/07; keep all three Boltz-family backends, samples and numeric controls |
| C02 | **Defects confirmed; superseded as implementation work by issue 24/R11 release-scope decision** | supported docs, packaging configuration, `src/cofolder/ui`, `run_ui.sh`, `.streamlit`, `tests/test_ui.py` | Test supported docs/artifact/dependency exclusion; retain recipe/CLI coverage; any source UI test is development-only | UI is not a release claim; retain Bias/Validate/Screen/Oracle CLI and Python APIs |
| C03 | **Statically confirmed security defect; current** | `input/config.py`, `boltz_runner.py`, `openfold3_runner.py`, failure serialization | Boltz/OpenFold3 runner tests and public-contract failure tests; add dummy-secret success/failure cases | PK-08, WF-20; preserve real credential delivery and actionable non-secret diagnostics |
| C04 | **Confirmed release mismatch; current** | `pyproject.toml`, package `__init__.py`, CLI version | `test_packaging.py` plus installed-wheel metadata/version checks | PK-02/06, DOC-19; keep schema version independent and ship MIT text/metadata |
| C05 | **Confirmed artifact gaps; current, installed behavior awaits S3** | packaging metadata, examples/tutorials/scripts and acceptance resources | `test_packaging.py`, `test_tutorials.py`, acceptance resource tests; add sdist-to-wheel inventory and UI exclusion | PK-01/05-08, DOC-15-17/20; retain runnable examples, PDB/CIF inputs, setup/data tools, supported tests and attribution |
| C06 | **Static inconsistencies; current, clean/backend verification pending** | `pyproject.toml`, backend availability/version checks, lazy analytics imports | packaging/registry/runner tests plus four isolated backend acceptance apps and fresh 3.11/3.12 lanes | PK-01/03/04/08, PI-08/09; retain four backends and every analysis capability via explicit extras |
| C07 | **Confirmed tooling gap; current** | packaging extras, repository test/lint/CI configuration | Supported suite without Streamlit, artifact/installed CLI checks; UI tests remain development-only if retained | PK-01/03/05, DOC-15/20 and cross-cutting claims; do not delete supported coverage to obtain green checks |
| C08 | **Confirmed compatibility surface; current; removals approved but gated** | `legacy/`, runner aliases, `utils/helpers.py`, CLI misspelling, builder wrapper, Screen fallback | registry, helpers, CLI, packaging and Screen clustering tests; move capability tests to canonical APIs | PI-03, WF-09; preserve underlying APIs/workflows and identity-bearing distance clustering |
| C09 | **Structural recommendation supported by C01; current** | typed input schemas, `input/command.py`, Boltz runner boundary, ligand preparation/cache code | input, ligand, runner contract/adapter and seed/model tests plus C01 regressions | PI-03/08, IN-01/04-08, WF-01/02/04; retain SMILES/MOL/SDF/CCD/conformer routes and runtime provenance |
| C10 | **Structural recommendation; current** | Validate/Screen recipes, `analytics/bias.py`, contracts/metrics/serialization | all recipe, bias, IFP/filter/cluster and public-contract tests; add before/after fixtures | PI-03-07/09/10, IN-12, WF-03/05/07-20; preserve scientific values/states, identities, decisions, network behavior and artifact families |
| C11 | **Usability recommendation; current** | `cli.py`, shared validation, recipes and packaged examples | `test_cli.py` plus Validate/Screen/Oracle preflight and ambiguity tests | PI-01/02/08, IN-01/10-12, WF-18/20; retain all advanced flags/defaults and two-file YAML support except approved R04 |
| C12 | **Statically confirmed duplication; current** | fetch script, `analytics/bias_training.py`, `analytics/build_bias_training_data.py`, future shared resolver | bias-training/build/bias analytics tests; add precedence, executable and installed-layout cases | PK-08, IN-12, WF-11; retain every existing MMseqs lookup route and report resolution |
| C13 | **Usability recommendation; current** | `utils/log.py`, CLI completion paths, contract serialization and recipe results | log, CLI, public-contract and workflow success/partial/total-failure tests | PI-07/09/10, WF-05/13/17/20, DOC-09/18; preserve outputs, identities, machine-readable failures and stream contracts |
| C14 | **Confirmed production/test coupling; current; R07 approved but gated** | `utils/log.py`, pytest logging fixtures | `test_log.py` and supported capture/handler-preservation cases | WF-20 support; retain normal console/file capture and unrelated application handlers |
| C15 | **Statically confirmed unreachable branch; current; R10 approved but gated** | `entities/ligand.py`, `input/command.py`, future packaged setup tool | `test_ligand.py`, `test_command.py`; add empty/partial/complete cache, invalid input and concurrency cases | IN-06, PI-03, PK-08; retain CCD conversion and explicit optional setup without test inference |

Additional protected surfaces remain runner auto-discovery; standalone Bias; all
Validate/Screen/Oracle modes; setup/download/vendor and temporary-conda tooling; bias
enrichment, plots, mixed/same-type datasets and checkpoints; affinity censoring,
conversion, correlation and plotting; all ligand formats and IFP taxonomies; five
Oracle reducers, custom scoring, gates and penalties; cache population; and every
current public output. Approval of R01-R11 is limited to the exact interfaces listed
in the removal register and is not permission to remove these underlying capabilities.
The UI is intentionally absent from this protected release list under R11.

### S0 removal prerequisite record

The user's 2026-09-09 correction, **"accept R01-R10 (all)"**, supersedes the
earlier R01-R02-only instruction and explicitly approves every ID in the finite
removal register then present. On 2026-09-10 the user additionally confirmed that the
UI will not be shipped and supplied issue 24's initial-release scope assessment;
this is recorded as R11. Approval settles intent but does not waive capability
safeguards, migrations, tests or the readiness gates below. S0 executes no removal.

#### R01-R04 inventories

- **R01:** no production code, test fixture or required runtime resource resolves a
  file under `legacy/`. Three current documentation surfaces do rely on the directory
  as an archival link/description: `README.md`, `docs/index.md`, and
  `docs/getting-started/repository-tour.md`. Remove or replace those references in
  the same S7 change. The `/legacy/cache` strings in `test_validate.py` are synthetic
  cache paths, not repository dependencies.
- **R02:** the aliases are defined only in `runners/base.py`, imported/re-exported in
  `runners/__init__.py`, and asserted in `test_registry.py`. No production caller,
  example or tutorial authors against an alias. `docs/tutorials/runners.md` mentions
  them only as legacy compatibility names. Replace
  `RunnerPreparation` -> `RunnerPreparationResult`, `RunnerRequest` ->
  `RunnerExecutionRequest`, and `RunnerResult` -> `RunnerExecutionResult`; retain
  the canonical contract tests.
- **R03:** repository callers of every listed helper compatibility wrapper are
  confined to `tests/modules/utils/test_helpers.py`; production code imports only
  the unrelated canonical `resolve_seed_plan` from `helpers`. Move wrapper capability
  coverage to the canonical APIs before removal: `set_dir` -> `create_dir`;
  `read_yaml`/`read_csv`/`read_sdf` -> `utils.read`; `delete_last_line` ->
  `utils.write`; censoring helpers -> `analytics.dataset`; affinity conversion and
  correlations -> `analytics.stats`; affinity plotting -> `analytics.plots`.
- **R04:** all user-facing repository invocations found in tutorials, guides and
  acceptance code use `--assess_bias`; only `cli.py` registers `--asess_bias` (and
  `argparse.BooleanOptionalAction` consequently exposes `--no-asess_bias`). Remove
  both misspelled forms together and retain `--assess_bias`/`--no-assess_bias` plus
  console/module parity coverage.

#### R08 exact environment proposal

- Change `requires-python` from `>=3.11` to **`>=3.11,<3.13`**, matching the stated
  and tested Python 3.11/3.12 support policy.
- Keep `numpy`, `pandas`, `pyyaml`, `biopython`, `rdkit`, `gemmi`, and `packaging` in
  base dependencies.
- Move exactly `scipy`, `scikit-learn`, `prolif`, `pdb2pqr`, `matplotlib`, and
  `seaborn` from base dependencies into a new **`analysis`** extra. This extra owns
  IFP clustering, affinity statistics, ProLIF/pdb2pqr interaction processing and
  plotting. Tutorials, acceptance, test and development groups remain separate; the
  full supported test lane installs `analysis` but no UI extra.
- Before the move, make imports at each optional feature boundary lazy enough for a
  base install to import COFOLDER and run backend-independent paths. Missing analysis
  packages must raise feature-specific installation guidance, never silently skip a
  requested metric. S4.4 applies this approved proposal and S4.5 verifies it.

#### R09 exact path migration proposal

Only the following old paths are approved for withdrawal; new private recipe/bias
submodules that do not displace an import are not additional removals.

| Old path/interface to remove | Canonical replacement and treatment |
|---|---|
| `cofolder.modules.input.command.Command` and `download_cache` | Typed `RunnerOptions` loaded through `cofolder.modules.input.config`; internal Boltz argv building in `cofolder.modules.runners._boltz_command`; explicit cache setup in `cofolder.tools.cache`. R05/R10 gates apply. |
| Backend-facing `cofolder.modules.entities.ligand.handle_conformers` | Internal `cofolder.modules.runners._ligand_preparation.prepare_ligand_conformers`; public chemistry conversion functions remain in `entities.ligand`. |
| `cofolder.modules.entities.ligand.cache_mols_from_sdf` | `cofolder.tools.cache.populate_ccd_cache_from_sdf` plus the installed cache/data-preparation entry point; `mol_to_ccd` and pure conversion routes remain public in `entities.ligand`. R10 gate applies. |
| `cofolder.modules.utils.gather` and its nine public functions | `cofolder.modules.analytics.aggregation` with the same behaviors for gathering, merging, chain identity, variance, bitstring similarity, robustness and metric classification; migrate recipe callers and tests, then remove the old module facade. |
| `cofolder.modules.analytics.build_bias_training_data` | `cofolder.tools.build_bias_training_data`; update the subprocess module target, packaging test and builder tests. |
| `scripts/build_bias_training_data.py` | Installed data-builder entry point backed by `cofolder.tools.build_bias_training_data`; update installation docs before deleting the source wrapper. |
| Monolithic `cofolder.modules.analytics.bias` implementation file and patched private-helper locations | Convert the same import path into a package; re-export `apply_bias_metrics` at `cofolder.modules.analytics.bias` while moving private reference/similarity/enrichment/dataset/artifact helpers to package submodules. No compatibility promise remains for underscored helper paths. |

Recipe extraction under C10 retains the public `Validate`, `Screen`, `Oracle`, and
`Bias` classes and their signatures. Runner discovery remains dynamic. Source setup,
download/vendor and temporary-conda wrappers other than the explicitly listed bias
builder are not removed by R09.

#### R11 exact UI release-surface proposal

- Keep Bias, Validate, Screen and Oracle as the four supported CLI commands and
  Python recipes. The new files under `docs/reference/python-api/` describe their
  current Python-only execution paths and runner extension boundary.
- Remove unqualified UI launch instructions and claims from supported release docs.
  Direct users to the four current commands and Python APIs.
- Do not add the UI launcher proposed by the original C05. Exclude `cofolder.ui` and
  Streamlit-only dependencies from the wheel, supported extras, installed entry
  points and release test lanes.
- The existing UI implementation, `.streamlit` data, `run_ui.sh` and shallow UI tests
  may remain in the source repository only as explicitly warned, unsupported
  development material. They are not repaired or refactored in this release plan.
- Do not describe issue 24 as successfully implemented. Its initial-release outcome
  is `not planned`; future UI work requires a new issue based on current recipe and
  runner contracts.
- This decision withdraws only the UI surface. It does not authorize removal or
  weakening of any CLI command, Python recipe, runner, metric, input route or public
  output contract.

#### Later readiness gates

- **R05:** finish S1's typed Boltz command migration, remove all internal `Command`
  callers, and document the retained `load_runner_options`/typed runner request
  signatures before withdrawing the legacy object.
- **R06:** during S6, prove normal Screen paths produce identity-bearing
  `InteractionFingerprint` values for every supported taxonomy, including distance,
  before deleting the vector-only `_legacy` fallback.
- **R07:** first replace the production handler-recovery expectation with supported
  pytest fixture/capture behavior and prove repeated application logging preserves
  unrelated handlers.
- **R10:** document and test `cofolder.tools.cache` explicit setup plus retained
  `mol_to_ccd`/conversion signatures before removing implicit prediction/download.

### Removal decision ledger

No `Rxx` item may be implemented unless its decision is `APPROVED` and its readiness
gate is satisfied. `APPROVED; GATED` means intent is settled but the named later
evidence is still required. Execute removals independently in S7 or in the explicitly
linked gated task, never as incidental cleanup.

| ID | Decision | Readiness / prerequisite evidence | Decision note |
|---|---|---|---|
| R01 | **APPROVED; APPLIED S7.2** | Documentation dependencies migrated and tracked archive deleted | User: "accept R01-R10 (all)", 2026-09-09 |
| R02 | **APPROVED; APPLIED S7.2** | Canonical runner contracts and discovery tests retained | User: "accept R01-R10 (all)", 2026-09-09 |
| R03 | **APPROVED; APPLIED S7.2** | Capability tests migrated to canonical modules before wrappers were removed | User: "accept R01-R10 (all)", 2026-09-09 |
| R04 | **APPROVED; APPLIED S7.2** | Correct positive/negative spelling and CLI parity retained | User: "accept R01-R10 (all)", 2026-09-09 |
| R05 | **APPROVED; APPLIED S7.2** | Typed options/runner request path documented and regression-tested | User: "accept R01-R10 (all)", 2026-09-09 |
| R06 | **APPROVED; APPLIED S7.2** | Identity-bearing distance/ProLIF clustering retained; vector-only rows are not evaluable | User: "accept R01-R10 (all)", 2026-09-09 |
| R07 | **APPROVED; APPLIED S7.2** | Supported capture, repeated setup and unrelated-handler preservation tests pass | User: "accept R01-R10 (all)", 2026-09-09 |
| R08 | **APPROVED; APPLIED S4.4** | Exact Python bound and six-dependency `analysis` move verified again in S7.3 | User: "accept R01-R10 (all)", 2026-09-09 |
| R09 | **APPROVED; APPLIED S7.2** | Exact finite paths migrated to canonical implementations and callers | User: "accept R01-R10 (all)", 2026-09-09 |
| R10 | **APPROVED; APPLIED S7.2** | Explicit setup/population APIs and commands replace implicit prediction | User: "accept R01-R10 (all)", 2026-09-09 |
| R11 | **APPROVED; APPLIED S2/S3** | UI exclusion from supported docs, dependencies, entry points, tests and wheel verified again in S7.3 | User: "The UI will not be shipped", 2026-09-10 |

### S7.3 cleanup verification evidence

Verification at clean commit `9d40404` found no live dangling reference in source,
scripts, supported documentation, tutorials or examples. References in the dated
release audits remain historical descriptions of their audited revisions. The only
current search matches were deliberate negative assertions and packaging exclusions.

| Removal | Negative and retained-capability evidence |
|---|---|
| R01 | The `legacy/` tree is absent; supported documentation and artifact tests do not depend on it. |
| R02 | The three aliases are absent from runner exports; canonical preparation/request/result contracts import and pass registry/contract tests. |
| R03 | All approved helper forwarders are absent; canonical read/write, dataset, stats, plots and aggregation coverage passes. |
| R04 | Neither misspelled flag is registered or documented; the correctly spelled positive/negative CLI paths pass. |
| R05 | The old command module is absent and non-importable; typed configuration, `_boltz_command` and all numeric serialization cases pass. Commit `3009c87` established the typed boundary before `9d40404` removed the facade. |
| R06 | Production code no longer fabricates `_legacy` identities; vector-only inputs stay unclustered and distance/ProLIF identity-bearing paths pass. Commit `ba45810` established and proved the postprocessing boundary before removal. |
| R07 | Production logging contains no pytest/private-handler lookup; capture, unrelated-handler preservation and repeated setup pass. Commit `d7a6735` established normal console/file/capture behavior before the special recovery path was removed. |
| R08 | Metadata is exactly `>=3.11,<3.13`; the six approved packages are absent from base and present in `analysis`, with packaging and optional-dependency tests passing. |
| R09 | The finite old modules/script are absent and non-importable; canonical aggregation, bias package, builder, ligand preparation and installed tools import and pass. |
| R10 | Old implicit cache/conformer APIs are absent; `mol_to_ccd`, explicit cache setup/population, invalid/partial/concurrent cases and installed help pass without inference. Commit `6ad502f` hardened cache behavior first; the explicit API, docs/tests and withdrawal then landed atomically in `9d40404`. |
| R11 | Streamlit/UI remain only warned source-development material; supported extras, entry points, test lane and built wheel exclude them while all four workflows remain covered. |

### Task ledger

#### S5.1 reopened user-story verification and cache design

S5.1 was reopened on 2026-09-14 because its earlier unit, documentation and artifact
checks did not establish the complete training-data journey requested by the user.
The existing S5.1 working-tree implementation was retained and extended. The
reopened work is complete with the evidence in the ledger and work log below; S5.2
is now the active task.

Current capability assessment:

| Scenario | How it works now | Current limit / S5.1 decision |
|---|---|---|
| One protein plus one ligand | Prepare one protein and one ligand source bundle, then run Bias/Validate/Screen/Oracle with `--assess_bias` and a selected `--bias_release_cutoff`. Both query chains receive their own reference rows and mixed-pair artifacts. | Supported in principle and by smaller component tests; add one end-to-end, no-network fixture that asserts chain identities, strict cutoff, metrics and artifact names. |
| Two non-identical proteins plus two non-identical ligands | One source-bundle pair is searched for every distinct selected protein sequence and canonical ligand. Outputs preserve `query_chain_id`; downstream bias analytics produces pair-specific mixed, protein-protein and ligand-ligand datasets. | The present database test proves two distinct protein searches but only one distinct ligand shared by two chains. Add an explicit 2P/2L fixture with two non-identical ligand chemotypes and prove there is no cross-chain contamination or accidental query collapse. |
| PDB comparison before 2023-06-01 | Reuse one immutable source snapshot and supply the exact date through `--bias_release_cutoff`. Both protein metadata and ligand occurrences use the strict rule `release_date < 2023-06-01`. | Retain a separate result directory and record the exact comparison cutoff and source snapshot. The date is not presented as a backend training cutoff. |
| Whole prepared PDB snapshot | Use a cutoff later than every release date represented by the prepared snapshot. | There is no explicit `whole`/unbounded selector, and a guessed future date is an opaque workaround. Add a documented `whole snapshot` mode or a manifest-derived inclusive snapshot boundary, with manifest provenance and tests, while retaining exact-date mode. |
| PDB slice supplemented with a custom crystal-structure set | Database-backed public sources can be combined with `--custom_protein_reference_path` and `--custom_ligand_reference_path`. Custom rows are appended as `source=custom` and are not removed by the public PDB cutoff, which matches supplementation semantics. A shared `pdb_id`, `source_structure_path` or `source_reference_path` can associate protein and ligand sides of a complex. | Partly possible: the workflow accepts derived protein CSV (`sequence` required) plus ligand CSV/SDF (`smiles` required), but it cannot ingest a directory/list of raw PDB/mmCIF protein-ligand complexes and build both sides automatically. Add a preparation command for a custom complex manifest/directory, stable complex identity/provenance, duplicate policy, extraction diagnostics and a generated supplement bundle or paired reference files. Document that custom structures are supplements, not members of the date-filtered public slice, unless a future explicit custom cutoff policy is selected. |

Required comparison layout and provenance:

```text
<assessment_root>/
├── pre-2023-06-01/
├── whole-snapshot/
├── pre-2023-06-01-plus-custom/
└── custom-only/
```

Each scenario must have its own result directory and record: exact cutoff policy,
source-bundle fingerprints/snapshot dates, custom-set fingerprint when present,
selected query chains and query hashes. A small documented matrix runner/config
should run these four scenarios without requiring four hand-written commands. It must
call the same public recipes rather than introduce a separate scientific path.

Current cache behavior is narrower than the intended behavior:

- A source bundle itself is reusable and read-only across systems, Screen compounds,
  Oracle candidates and cutoffs.
- During one reference-materialization call, identical protein sequences share one
  MMseqs invocation and identical canonical ligand molecules share one fingerprint
  calculation. Distinct 2P/2L queries are intentionally computed separately.
- `<wrk_dir>/results/bias_train/reference_manifest.json` reuses the derived protein
  and ligand CSVs only when the complete request matches: both bundle fingerprints,
  all protein/ligand query hashes and chain IDs, cutoff, ligand threshold and selected
  chains. This safely reuses an identical rerun in the same output directory.
- Screen places every compound under a separate execution directory and invokes
  Validate there. Its invariant protein query is therefore searched again for every
  compound. Oracle writes under `oracle_run`; an identical candidate rerun can reuse
  its references, but changing the ligand changes the complete manifest and causes
  the unchanged protein side to be searched again. Independent systems/work
  directories do not share derived searches. Files under `_protein_search` are
  working files, not a validated persistent cross-run cache.
- Consequently the requested protein-similarity reuse in Screen/Oracle and between
  systems is **not yet satisfied**. Tests must count MMseqs calls; the presence of
  similarly named files is not sufficient evidence of reuse.

Streamline this with a content-addressed, concurrency-safe bias query cache separate
from result directories. Protein cache keys must include the protein source-bundle
fingerprint, normalized sequence hash, MMseqs version/search parameters and cache
schema; store raw hits before release-date filtering so one search can serve the
pre-2023-06-01 and whole-snapshot comparisons. Ligand cache keys must include the ligand-bundle
fingerprint, canonical molecule/fingerprint parameters and cache schema; likewise
apply cutoff and reporting threshold after loading reusable raw similarities. Cache
writes must be atomic and locked for concurrent Screen/Oracle workers, validate
checksums/schema before reuse, preserve `query_chain_id` only when materializing a
run, and expose hit/miss/rebuild provenance in `reference_manifest.json`. Provide a
documented shared-cache location/override and an explicit safe invalidation command;
never use only a chain ID, output directory or cutoff as a cache identity.

S5.1 acceptance for these stories:

1. Run mocked/no-network 1P/1L and genuinely non-identical 2P/2L fixtures through
   standalone Bias and the bias path used by Validate; assert all expected pair
   artifacts, chain-qualified metrics and exact source/query provenance.
2. Exercise strict pre-2023-06-01 and whole-snapshot scenarios against one fixed
   bundle; prove strict boundary behavior, monotonic eligible reference sets,
   separate result retention and raw-search reuse across the two policies.
3. Exercise pre-2023-06-01 plus a two-complex custom supplement and custom-only;
   prove public/custom attribution, paired complex identity and unchanged public
   cutoff semantics. Include raw PDB/mmCIF custom-complex preparation coverage.
4. In Screen, assess at least two different ligands against an invariant protein and
   prove one shared protein search rather than one per compound. In Oracle, prove
   repeated and different candidates reuse the invariant protein search while
   ligand searches are reused only for matching canonical molecules.
5. Run a second system sharing one protein and one ligand with the first; prove cache
   hits for shared query components and misses for non-identical components. Also
   test stale bundle fingerprints, parameter/schema changes, corrupt entries and two
   concurrent writers.
6. Document the scenario matrix, exact cutoff versus backend-training-cutoff
   distinction, custom input schema/preparation, output comparison, shared-cache
   identity/lifetime/invalidation, and Screen/Oracle behavior. Preflight must report
   scenario sources and cache readiness without searching, inference or mutation.

Evidence fields start as `—`. Replace them with concise test results, artifact paths,
or a work-log reference as tasks finish.

| Task | Deliverable | Acceptance / required evidence | Status | Evidence / resume_from |
|---|---|---|---|---|
| S0.1 | Fresh baseline and worktree reconciliation | Record commit, Python/platform, tracked/untracked changes, `pytest -q`, focused Ruff baseline, CLI/module smoke results; identify pre-existing changes without modifying them | **DONE** | Baseline section and 2026-09-09 work-log row: commit/worktree owned, 520 tests passed, 15 expected Ruff findings, all seven console/module smoke cases matched. |
| S0.2 | Scope and capability map | Map C01–C15 to owning files/tests and confirm matrix-required plus additional retained capabilities; note findings that have become stale | **DONE** | S0 map records owner/test/claim/safeguard for C01-C15; audited commit and baselines match, so none is stale. |
| S0.3 | Permission freeze | Present R01–R10 individually; fill prerequisites for R08/R09; record each as approved/rejected/pending. Pending decisions do not block S1–S6 behavior-preserving work | **DONE** | R01-R10 explicitly approved; R01-R04 inventories, exact R08 environment proposal, exact finite R09 path table and R05-R07/R10 later gates recorded. No removal executed. |
| S0.4 | UI scope and Python API orientation | Record issue 24/R11 exactly; create one implementation-flow reference per recipe plus a new-runner construction reference before changing later phases | **DONE** | UI excluded from supported/shipped initial release by user decision; five references plus index added under `docs/reference/python-api/`; no UI code repaired or removed. |
| S1.1 | Typed Boltz argv regression tests | Cover integer `0/1/2`, floats, booleans and absent values through `load_options` and Boltz1/Boltz2/Community adapters; assert full argv and default sample cardinality | **DONE** | `test_boltz_command.py` covers all three adapters, typed 0/1/2/float/true/false/absent behavior and default sample count; included in 2026-09-10 focused 51-pass run. |
| S1.2 | Correct typed command construction | Introduce the narrow `_boltz_command` boundary or equivalent; preserve the public `Command` facade; remove internal truth-value ambiguity and lossy typed→legacy→argv conversion | **DONE** | `_boltz_command.py` builds argv from `RunnerOptions`; Boltz runners retain typed options through preparation/execution; legacy `Command` remains and its numeric regression passes. |
| S1.3 | Secret-safe execution reporting | Use real argv only for subprocess execution and a recursively redacted display form for Boltz/OpenFold3 logs and typed failures; test success/failure without leaking dummy secrets | **DONE** | `_command_reporting.py` plus Boltz/OpenFold3 integrations preserve real subprocess argv and redact logs/results/failures; success/failure tests included in focused 51-pass run. |
| S1.4 | Phase verification | Run focused command/runner/contract tests, CLI smoke tests and `pytest -q`; record exact counts and remaining known failures | **DONE** | Focused runner/input tests: 51 passed; `pytest -q`: 530 passed in 26.17s; console/module root/four-command help and version matched at exit 0, unknown command matched at exit 2. |
| S2.1 | Supported documentation scope | Remove unqualified UI launch/availability claims; mark any retained references unsupported and development-only; direct users to Bias/Validate/Screen/Oracle CLI and Python APIs | **DONE** | README UI installation/launch section removed; repository tour marks retained files unsupported; README, docs home and MkDocs navigation link all five Python API flow references; supported-doc search passes. |
| S2.2 | Development-material boundary | Add explicit warnings to retained UI source/launcher material and keep it out of the supported architecture; do not repair obsolete workflow or command behavior | **DONE** | UI package/app, visible header, launcher, Streamlit config and renamed checker warn that behavior is obsolete/unsupported; no workflow command behavior repaired. |
| S2.3 | Installed/test surface boundary | Exclude `cofolder.ui`, UI entry points and Streamlit-only dependencies from the wheel and supported extras/test collection while preserving supported CLI/recipe coverage | **DONE** | `ui` extra removed; setuptools excludes `cofolder.ui`; checker renamed outside pytest's pattern; packaging regressions pass and default collection contains 525 supported tests with no UI nodes. |
| S2.4 | UI-scope verification | Search supported docs for launch claims; inspect package metadata/wheel/import inventory; prove supported tests collect without Streamlit and all four CLI/Python workflows remain represented | **DONE** | Clean 84-file wheel has no UI/Streamlit and only the `cofolder` entry point; 525 supported tests and 6 optional UI checks pass; all seven console/module parity cases and link/search checks pass. |
| S3.1 | Single version source and MIT metadata | Runtime, wheel metadata and both CLI entry points resolve to `1.0.0`; output schema version remains independent; wheel metadata and licence file identify MIT | **DONE** | `cofolder.__version__` is the setuptools dynamic source; focused 38-pass and full 526-pass suites; clean 84-file wheel and outside-checkout install report `1.0.0`, `License-Expression: MIT`, packaged `LICENSE`, no UI files and 12 acceptance resources. |
| S3.2 | Explicit artifact inventory | Define wheel/sdist allowlists; include licence/attribution in both, complete supported tests and required resources in sdist, exclude UI from the wheel, and package runnable examples plus required small PDB/CIF references without caches/weights | **DONE** | Explicit manifest/package-data policies and 10 synchronized examples added; 10 packaging regressions pass. Clean setuptools 84 build produced an exact 217-file sdist with all 47 supported test files and an exact 97-file wheel with 10 examples, 12 acceptance resources and both legal notices; UI/generated/model material absent. |
| S3.3 | Installed access paths | Add resource locator/copy command plus installed backend-setup/data-preparation entry points; retain source wrappers until their approved R09 migration and S7 removal step is ready; add no UI entry point | **DONE** | `cofolder-tools` exposes copy-examples plus all four setup/data helpers; resource APIs and retained wrappers are tested; 60 focused and 551 full tests pass. |
| S3.4 | Artifact verification | Build sdist, build wheel from unpacked sdist, inspect both inventories, install outside checkout, run import/version/help/resource-copy/dry-run checks, then run sdist tests with declared extras | **DONE** | Clean `HEAD` `6627de3` built a 228-file sdist and its 104-file wheel; external install metadata/resources/help passed; 25 focused tool/setup tests and all 551 supported sdist tests passed without Streamlit. |
| S4.1 | Repository-owned lint/dev configuration | Add a minimal explicit Ruff baseline (`E4,E7,E9,F` initially), place development tools in a development extra, fix the 15 audited findings deliberately, and remove `tests/tmp_test.py` only as non-interface dead test data | **DONE** | Ruff passes with the project-owned `E4/E7/E9/F` configuration; 208 focused and 552 full supported tests pass; optional UI checks remain 6/6 development-only; metadata and diff checks pass. |
| S4.2 | Test lanes and CI | Define core, contracts/tutorial, artifact and acceptance lanes with correct extras and no Streamlit requirement; add CI for supported Python versions; keep backend jobs isolated because their `boltz` namespaces conflict | **DONE** | Four exhaustive lanes pass independently (482 core, 22 contracts/tutorial, 14 artifact, 42 lightweight acceptance; 560 total); Ruff and diff checks pass; sdist-to-wheel installed smoke passes; routine CI defines six Python 3.11/3.12 test jobs without Streamlit or backend extras. |
| S4.3 | Safe dependency bounds/provenance | Bound Boltz2 to supported major 2 and replace moving Community provenance with a tested immutable source/release; add explanatory failures for unsupported backend/Python combinations without applying R08 | **DONE** | Boltz2 is bounded to `>=2.0.0,<3`, Community is pinned to PyPI `2.10.12`, unsupported Python/backend versions fail explicitly; 42 focused and all 568 supported tests plus artifact verification pass without applying R08. |
| S4.4 | R08-dependent environment changes | Apply the approved exact `>=3.11,<3.13` bound and six-package base→`analysis` move after adding lazy feature boundaries and installation guidance | **DONE** | Exact Python/base/analysis metadata applied; fresh-process blocked-import regression covers analysis and Boltz packages; 115 focused and 580 full supported tests, four lanes, artifact verification, Ruff, diff and seven CLI/module parity cases pass. |
| S4.5 | Environment verification | Fresh 3.11/3.12 base installs, full supported test environment without Streamlit, optional-analysis guidance, and four separate backend availability/smoke environments pass or have resource-specific blockers recorded | **DONE** | Fresh wheel SHA-256 `967cb43f45133d868580b588d0399516e4b7e862a164256e5e67b1e263e847cd`; Python 3.11.16/3.12.14 base smokes and both 580-test full matrices passed without Streamlit; isolated real inference passed for Boltz 1.0.0, Boltz 2.2.1, Community 2.10.12 and OpenFold3 0.5.0. Evidence: `/tmp/cofolder-s4.5.uUjsHm`. |
| S5.1 | Easier CLI journey and bias training-data scenarios | Group help without removing flags; add shared preflight-only validation, unique-ligand inference, database-backed bias sources, contradiction checks, provisioning documentation and bundled-example path. Verify 1P/1L and non-identical 2P/2L through the pre-2023-06-01, whole-snapshot, pre-2023-06-01-plus-custom and custom-only scenarios; implement component-level shared cache reuse for Bias/Validate/Screen/Oracle and cross-system runs as specified above | **DONE** | 609-test full lane passed. Database-backed Validate 1P/1L and standalone 2P/2L, the exact four-scenario runner, raw PDB/mmCIF custom preparation, shared component cache/corruption/concurrency/invalidation, and offline MMseqs backfill/plot boundaries are covered. Strict docs, Ruff, diff and sdist-to-wheel artifact verification pass; real local MMseqs version `01683a607f83878e95436632d73e1d7d9ae30955` returned `mmseqs_pident` plotting value `0.000 < 0.25`. |
| S5.2 | Shared executable discovery | Implement one tested MMseqs resolver with documented precedence, executable rejection reasons, space-safe paths, source and installed layouts; preserve all lookup routes | **DONE** | Shared structured resolver and three retained forwarding seams pass 42 focused and 621 full supported tests; Ruff, strict docs, diff and installed-artifact verification pass. |
| S5.3 | Progress and result summaries | Add concise console versus detailed debug/file formatting and shared completion summaries for success, partial/total failure and unavailable Oracle metrics; preserve output files and stream contracts | **DONE** | Default console output is concise while debug/file output remains detailed; an outermost-only reporter reads canonical public records and covers all four workflows, partial Screen, total failure and available/unavailable Oracle values. Focused 85 and full 630-test lanes, Ruff, strict docs, diff, artifact verification and seven console/module parity cases pass. |
| S5.4 | Usability verification | Test console/module parity, no-service preflight, ambiguity failures, resolver precedence, paths with spaces and canonical result/manifest pointers | **DONE** | Added durable real-entry-point and four-workflow preflight regressions; 71 focused and 641 full tests pass with resolver, ambiguity, space-containing path and completion-pointer coverage. Ruff, diff, strict docs and sdist-to-wheel installed verification pass. |
| S6.1 | Workflow shared internals | Extract planning/lifecycle, diagnostics and public-result assembly in small diffs; preserve public recipe classes/signatures and failure/identity semantics | **DONE** | Private execution, diagnostic and result seams plus import/semantic characterization tests; 154 focused and 646 full supported tests pass with lint, docs and artifact gates. |
| S6.2 | Screen postprocessing | Extract ranking/filtering/clustering while preserving every taxonomy, decision and output family; remove the approved vector fallback only after the R06 identity proof | **DONE** | Frozen private postprocess configuration and delegated metric/filter/cluster implementation; 66 focused and 654 full tests pass. PDB/mmCIF, distance/ProLIF and public-record/CSV identity proofs satisfy R06 readiness without executing the removal. |
| S6.3 | Bias analytics package | Split reference loading, similarity, enrichment, pair datasets and artifacts behind a stable `apply_bias_metrics`; imports must not trigger workflows/network | **DONE** | Stable facade plus seven acyclic responsibility modules; 88 focused and 655 full supported tests pass with import, lint, docs and artifact gates. |
| S6.4 | Chemistry/cache and aggregation boundaries | Separate chemistry from backend cache orchestration and move aggregation to analytics; follow the approved R05/R09/R10 path table only after each readiness gate | **DONE** | Typed runner preparation and canonical analytics aggregation boundaries pass 59 focused and 662 full supported tests; compatibility facades and all gated cache/command interfaces remain. |
| S6.5 | Structural parity verification | Compare fixture outputs before/after for counts, identities, seeds, values/states, provenance, gates, filters/clusters and artifact families; run import-boundary checks and full suite | **DONE** | Identical normalized digests from `6d5c48b` and `6175f41`: Validate `c05c6ee2`, Screen `087af4ef`, Oracle `ca7787d2`, Bias `6245cbe7`; 3 parity/import tests, 273 focused tests and all 665 supported tests pass with lint, docs and artifact gates. |
| S7.1 | Safe cache corrections independent of R10 | Validate input before writes, inspect required cache components, use unique temporary workspaces, and test empty/partial/complete/concurrent/error cases without changing the documented automatic behavior | **DONE** | Source validation precedes writes; cache-component inspection, serialized setup, unique temporary workspaces and destination forwarding pass 50 focused and 682 full supported tests plus lint, docs and artifact gates. |
| S7.2 | Execute approved removals | Implement only `APPROVED` R01–R11 items, one ID per reviewable change; update canonical callers/tests/docs and retain the underlying capability safeguards | **DONE** | R01–R07/R09/R10 applied; R08/R11 reconciled as already applied. Focused 139, full 656 and final tool/packaging 35 tests pass; Ruff, strict docs, artifact verification and diff checks pass. |
| S7.3 | Cleanup verification | Search for dangling references per executed ID; prove gated interfaces remained until their prerequisites passed; run relevant capability tests and full suite | **DONE** | Per-ID searches/import probes found no live dangling references; prerequisite history and atomic replacement evidence recorded above. Focused 227 and full 657 tests pass with only the two known warnings; Ruff, strict MkDocs, artifact verification and diff checks pass. |
| S8.1 | Static and unit release gate | Repository lint, supported core/contract/tutorial tests and CLI/module parity pass without Streamlit in the declared test environment | **DONE** | Independent initial lanes passed 577 core, 22 contracts/tutorial, 42 lightweight acceptance, 5 explicit redaction and 8 real console/module parity tests; after the S8.3 correction, Ruff and the complete 658-test supported suite passed again |
| S8.2 | Distribution release gate | Clean sdist→wheel build, explicit inventory, installed metadata/licence, external-checkout smoke and sdist test execution pass | **DONE** | Verifier passed again after the correction; retained handoff build `/tmp/cofolder-s8-handoff.PZc2KT` contains a 312-entry sdist and its exclusively derived 125-entry wheel with SHA-256 `55674a8e...a0b9b` and `bb399581...45e9` |
| S8.3 | Environment/backend release gate | Supported Python/base/analysis matrix and actual Boltz1/Boltz2/Community/OpenFold3 acceptance run in isolated environments; record unavailable hardware/service/model resources without claiming a pass | **DONE** | Python 3.11.15/3.12.13 base/full installs pass `pip check`, installed smokes, no-Streamlit checks and both 658-test full lanes; isolated actual inference passes one sample/zero failures for all four backends with 45/55/57/43 metrics |
| S8.3-FIX | Fresh-environment verification correction | Scope any defect exposed by S8 separately, add a regression and rerun every affected source/artifact/environment/backend gate | **DONE** | Fixed pandas view mutation and duplicate-index alignment in protein pairwise similarity; isolated structural-parity fixtures from ambient backend availability/version. New regression and corrected semantic fixture pass under pandas 2/3, followed by complete source, artifact, environment and backend reruns |
| S8.4 | Claim and removal reconciliation | Re-evaluate every affected PI/IN/WF/PK claim, all additional retained capabilities and R01–R11; link evidence and list any residual risk | **DONE** | Every PI/IN/WF/PK ID is accounted for in `v1.0-claim-matrix.md`; retained additional capabilities remain represented; R01–R11 remain approved/applied with no new removal; PK-02's tag/documentation portion is an explicit handoff rather than a false pass |
| S8.5 | Handoff | Set `overall_status` accurately; summarize shipped changes, deferred documentation/tag/publication work and exact blockers. Mark `COMPLETE` only if all required gates pass | **DONE** | All required S8 verification-matrix rows pass; residual warnings/limits and the precise documentation/tag/publication handoff are recorded below. No tag or publication was performed |

### Phase exit checks

- **S0:** baseline evidence is current, unrelated work is identified, permissions
  include R01-R11, and the recipe/runner Python flow references capture the current
  architecture before refactoring.
- **S1:** numeric values retain type and value end-to-end, credentials never appear in
  display/log/failure representations, and existing public command APIs still work.
- **S2:** supported documentation, dependencies, installed artifacts and release
  tests contain no UI claim or requirement; any retained UI source is unmistakably
  unsupported development material, while all four current CLI/Python workflows
  remain supported.
- **S3:** artifacts are intentional, self-describing and usable without repository
  paths; version/licence evidence comes from the installed artifact.
- **S4:** reproducible project-owned checks exist, supported environments are stated
  honestly, and no R08 interface change slipped into dependency cleanup.
- **S5:** a basic CLI/Python user can preflight bundled inputs and find results, while advanced
  controls, stream behavior and executable lookup routes remain available. The
  1P/1L and non-identical 2P/2L bias matrices work for exact dated and whole-snapshot
  PDB sources with/without custom crystal-complex supplements, and content-addressed
  component searches are demonstrably reused across cutoffs, Screen/Oracle
  candidates and systems without losing chain identity or provenance.
- **S6:** module boundaries improve without changing public workflows, scientific
  semantics, outputs or network behavior; compatibility facades remain unless
  specifically approved.
- **S7:** every deletion/interface withdrawal cites an approved R ID and evidence;
  unapproved items remain implemented and covered.
- **S8:** acceptance evidence distinguishes passing mocked/unit checks from actual
  backend inference and clean dependency installation. No unavailable check is
  silently treated as passed.

### Final verification matrix

Fill this table during S8. A row is `PASS` only when its command/evidence was produced
from the release candidate, not inherited from the audit baseline.

| Gate | Required check | Status | Evidence |
|---|---|---|---|
| Source | Repository-owned Ruff baseline | **PASS** | `ruff check src tests scripts tutorials examples`: clean on the corrected S8 working-tree candidate derived from `fcae176` |
| Source | Full supported core/contracts/tutorial test lanes without Streamlit | **PASS** | Initial independent lanes passed 577 core, 22 contracts/tutorial and 42 lightweight acceptance tests; final `python scripts/run_test_lane.py all` passed all 658 after the S8.3 correction |
| CLI | Console/module help, version, invalid command and preflight parity | **PASS** | Eight real subprocess cases passed: root help/version, four command helps, invalid command and copied-example preflight |
| Security | Success/failure secret-redaction tests for all command reporters | **PASS** | `tests/modules/runners/test_command_reporting.py`: 5 passed; Boltz/OpenFold3 real argv remains usable while success/failure diagnostics are redacted |
| Artifact | Sdist inventory and wheel built from unpacked sdist | **PASS** | Final verifier passed; retained 312-entry sdist SHA-256 `55674a8ee6e1fd6085cc1cc1d8699764eae629694bbb87f15ebde72a0d9a0b9b` exclusively produced the 125-entry wheel SHA-256 `bb399581c7b762d0b5ebcab7421e3dcb15efaddd28db394c804298327f4b45e9` |
| Artifact | Installed-outside-checkout resources, entry points, metadata and licence | **PASS** | Clean verifier passed installed import/version, MIT metadata/two legal files, `cofolder`/`cofolder-tools` entry points, tool help, 12 acceptance resources and 12 copied example resources |
| Scope | UI absent from wheel, supported extras/entry points/docs and release test requirements | **PASS** | Explicit inventory reports zero `cofolder.ui` files; verifier confirms no UI/Streamlit dependency, entry point or supported test requirement |
| Environment | Fresh Python 3.11 base/full lanes | **PASS** | CPython 3.11.15 base/full wheel installs pass `pip check`, import/version, installed entry points/resources/options/input smokes and the 658-test full lane; Streamlit is absent and omitted analysis imports give declared guidance |
| Environment | Fresh Python 3.12 base/full lanes | **PASS** | CPython 3.12.13 base/full wheel installs pass the same checks and independent 658-test full lane; package imports resolve to the isolated site-packages and Streamlit is absent |
| Backend | Boltz1 actual acceptance | **PASS** | Isolated Boltz 1.0.0 actual GPU run: manifest `success`, one success, zero failures, 45 metric records, detected runner/backend provenance |
| Backend | Boltz2 actual acceptance | **PASS** | Isolated Boltz 2.2.1 actual GPU run: manifest `success`, one success, zero failures, 55 metric records, detected runner/backend provenance |
| Backend | Boltz Community actual acceptance | **PASS** | Isolated Boltz Community 2.10.12 actual GPU run: manifest `success`, one success, zero failures, 57 metric records, detected runner/backend provenance |
| Backend | OpenFold3 actual acceptance | **PASS** | Isolated OpenFold3 0.4.1 actual GPU run: manifest `success`, one success, zero failures, 43 metric records, detected runner/backend provenance |
| Claims | PI/IN/WF/PK and additional-capability reconciliation | **PASS** | Every ID is explicitly reconciled in `v1.0-claim-matrix.md`; runtime/code gates pass, while PK-02 accurately hands off its forbidden tag and deferred final-documentation portions. Retained additional capabilities are accounted for below |
| Removals | Approved/rejected/pending R01–R11 reconciled | **PASS** | R01–R11 are all approved and applied with their recorded replacements/prerequisites intact; fresh full, artifact, environment and real-backend checks expose no additional or unapproved removal |

### S8 claim, capability and handoff reconciliation

The claim-by-claim S8 status is recorded in the mandatory
[`v1.0-claim-matrix.md`](v1.0-claim-matrix.md). PI-01–PI-10,
IN-01–IN-12, WF-01–WF-20 and the code/runtime/artifact portions of PK-01–PK-08
have fresh passing evidence. PK-02 is intentionally split: installed metadata,
runtime and both CLI entry points report `1.0.0`, but a release tag was not created
and final documentation-wide release alignment remains handoff work.

The additional capability guardrails also remain accounted for:

- automatic runner discovery and selected-runner provenance pass in contracts and
  in four isolated actual backend runs;
- supported setup, download, vendor and temporary-environment paths remain exposed
  through packaged `cofolder-tools`, source wrappers and the installed help smokes;
- Bias landscape plots, mixed/same-type pair data, enrichment and checkpoint output
  remain covered by fresh recipe/structural-parity tests;
- affinity censoring, conversion, correlation and plotting remain in their canonical
  dataset/stats/plots modules with fresh full-suite coverage;
- explicit cache setup/population, shared Bias query-cache inspection/invalidation,
  MSA reuse and executable discovery remain covered and packaged;
- the aliases, misspelled flag, archived entry paths and vector-only fingerprint
  fallback named by the old guardrail text were removed only under approved
  R01–R06/R09. The UI boundary remains the separately approved R11 unsupported
  source-development surface, absent from installed artifacts.

Residual limits and handoff items are not release-gate passes by implication:

- the complete source suite still reports the known MDAnalysis deprecation and
  pandas concatenation future warnings in the ambient environment; fresh full
  environments report the MDAnalysis warning, and backend logs contain upstream
  Lightning/Biotite multiprocessing or deprecation notices;
- actual inference is one small reference-free sample per backend using retained
  local model caches. OpenFold3 intentionally used its query-only dummy MSA with
  `--use-msa-server=False`; no fresh model download, live MSA service, credentialed
  service, exhaustive diagnostic combination or remote CI run is claimed;
- the OpenFold3 optional dependency remains unbounded in metadata; the detected
  and passing S8 backend version is 0.4.1, so later dependency upgrades require the
  same isolated backend gate;
- broad documentation completion, final version/repository/citation/archive
  alignment, release tagging and artifact publication remain outside this plan and
  require separate authorization. These are the precise resume point for release
  finalization; S8 itself has no unavailable required check or blocker.

### Post-audit reconciliation — issue 24 UI retirement

On 2026-09-19, the user explicitly approved complete removal of the obsolete
Streamlit prototype for the v1.0 release. This supersedes only R11's earlier option
to retain warned source-development material; the historical S2-S8 entries below
remain unchanged as evidence of the state that was verified at the time. The release
continues to expose only the supported `bias`, `validate`, `screen`, and `oracle`
CLI and Python workflows. Issue 24 is to be closed as not planned after this removal
is merged and verified, without creating a replacement UI issue.

### Work log

Append entries; do not rewrite history except to correct a factual error and note the
correction. The initial row records plan creation, not implementation progress.

| Date | Task | Commit/state | Paths changed | Verification/result | Decisions and resume point |
|---|---|---|---|---|---|
| 2026-09-09 | PLAN | `e9642f2`; audit files untracked | `docs/release/codebase_cleanup_audit_09092026.md` | Converted the audit's outline into a staged, resumable plan; no implementation or verification run | Start S0.1 by refreshing the baseline and identifying ownership of all worktree changes |
| 2026-09-09 | S0.1 | `e9642f2`; working tree with two pre-existing untracked audits | `docs/release/codebase_cleanup_audit_09092026.md` | `pytest -q`: 520 passed; isolated Ruff: 15 findings; root/four-command help, version and invalid-command console/module outputs and exits all matched | Baseline matches the audit; other audit preserved by hash. Start S0.2 capability/owner mapping |
| 2026-09-09 | S0.2 | `e9642f2`; working tree | `docs/release/codebase_cleanup_audit_09092026.md` | Inspected C01-C15 owners/callers/tests and claim guardrails; no finding is stale at the unchanged audited commit | Recorded all matrix and additional retained capabilities. Start S0.3 approval/prerequisite record |
| 2026-09-09 | S0.3 | `e9642f2`; working tree | `docs/release/codebase_cleanup_audit_09092026.md` | Recorded explicit approval for R01-R10; enumerated R01-R04, fixed exact R08 Python/dependency policy, fixed finite R09 path migrations, and linked remaining technical gates; no removal executed | S0 exit satisfied. Advance to S1.1 typed Boltz argv regression tests |
| 2026-09-10 | S0.4 | `e9642f2`; working tree | `docs/reference/python-api/*.md`, this audit | Traced all four recipe paths and runner contract; `git diff --check` passed | Issue 24 and user decision authorize R11 supported/shipped UI withdrawal. Python API flow references precede plan adjustment |
| 2026-09-10 | S1.1 | `e9642f2`; working tree | `tests/modules/runners/test_boltz_command.py`, runner/input tests | Focused input/runner selection: 51 passed, including all three Boltz adapters and typed numeric/boolean/absent cases | Regression coverage complete; record existing typed implementation as S1.2 |
| 2026-09-10 | S1.2 | `e9642f2`; working tree | `_boltz_command.py`, `command.py`, Boltz runner modules and tests | Focused input/runner selection: 51 passed | Typed `RunnerOptions` now reaches Boltz argv construction without legacy round-trip; legacy facade retained. Record redaction work as S1.3 |
| 2026-09-10 | S1.3 | `e9642f2`; working tree | `_command_reporting.py`, Boltz/OpenFold3 runner modules and tests | Focused input/runner selection: 51 passed; dummy secrets absent from success/failure logs and serialized failure evidence | Redaction implementation complete at focused scope. Resume S1.4 with CLI parity and full suite |
| 2026-09-10 | PLAN-2 | `e9642f2`; working tree | `docs/release/codebase_cleanup_audit_09092026.md` | Reconciled issue 24, R11, phase dependencies, task ledger and final gates; no release implementation added by this plan edit | S1 remains active at S1.4; S2 replaces UI repair with supported-surface exclusion; all non-UI findings remain scheduled |
| 2026-09-10 | S1.4 | `e9642f2`; working tree | S1 implementation/tests and this audit | `pytest -q`: 530 passed in 26.17s; console/module root and four-command help/version outputs and exits matched; unknown command matched at exit 2 | S1 exit satisfied. Advance to S2.1 supported documentation scope |
| 2026-09-10 | S2.1 | `4d1fbd6`; working tree | `README.md`, documentation home/tour/navigation | Link/search checks pass for all five Python API flows; no supported UI installation or launch recommendation remains | Supported surface directs users to Bias/Validate/Screen/Oracle CLI and Python APIs. Advance to S2.2 |
| 2026-09-10 | S2.2 | `4d1fbd6`; working tree | UI package/app, launcher, Streamlit config and development checker | `python tests/ui_development_check.py`: 6/6 passed with unsupported-development result | Retained UI material is visibly unsupported and obsolete; its workflow behavior was not repaired. Advance to S2.3 |
| 2026-09-10 | S2.3 | `4d1fbd6`; working tree | `pyproject.toml`, `tests/test_packaging.py`, renamed UI checker | Packaging tests: 5 passed; default collection: 525 tests and no UI nodes | Removed supported `ui` extra; excluded `cofolder.ui` from discovery; preserved the checker outside the release lane. Advance to S2.4 |
| 2026-09-10 | S2.4 | `4d1fbd6`; working tree | S2 changes and this audit | Clean wheel: 84 files, no UI/Streamlit, only `cofolder` entry point; `pytest -q`: 525 passed; seven CLI/module cases matched; link/search and `git diff --check` passed | S2 exit satisfied. Advance to S3.1 version and MIT metadata |
| 2026-09-10 | S3.1 | `31b3cf4`; working tree | `pyproject.toml`, `src/cofolder/__init__.py`, `tests/test_packaging.py`, `tests/test_cli.py`, this audit | `pytest -q tests/test_packaging.py tests/test_cli.py`: 38 passed; `pytest -q`: 526 passed in 25.86s; `git diff --check` passed; locally cached setuptools 83/packaging 25 built an 84-file wheel whose metadata/version/licence, UI exclusion and 12 acceptance resources passed inspection; no-network outside-checkout install and both CLI entry points reported `1.0.0` | `cofolder.__version__` is the sole version source via dynamic setuptools metadata; SPDX MIT metadata requires setuptools `>=77.0.3`; `include-package-data = false` preserves the approved UI exclusion under the newer backend. Output schema version remains independent. Advance to S3.2 explicit artifact inventory |
| 2026-09-10 | S3.2 | `e173b35`; working tree | `MANIFEST.in`, `pyproject.toml`, `src/cofolder/resources/`, `tests/test_packaging.py`, this audit | `pytest -q tests/test_packaging.py`: 10 passed; `git diff --check` passed; isolated setuptools 84 build produced exact allowlisted 217-file sdist/97-file wheel; sdist has 47 supported test files, wheel has 10 synchronized examples, 12 acceptance resources and 2 legal notices; importlib resource inventory passed | Excluded unsupported UI/legacy and ignored generated docs plus caches/weights; retained all supported source docs/tutorials/scripts/tests. Advance to S3.3 installed access paths |
| 2026-09-10 | S3.3 | `89d7e2a`; working tree | `pyproject.toml`, `src/cofolder/tools/`, example resource API, retained scripts, setup/data docs and focused tests | Focused tools/resources/packaging/OpenFold3/builder selection: 60 passed; `pytest -q`: 551 passed in 29.31s; targeted isolated Ruff and `git diff --check` passed; existing console/module help stayed identical and all six tools help paths returned 0 without setup/network work | Added only `cofolder-tools` with five subcommands; all four source helpers remain as wrappers, the old analytics builder remains behind the new tool facade, and no UI entry point was added. Advance to S3.4 clean artifact verification |
| 2026-09-10 | S3.4 | `6627de3`; working tree | This audit; temporary artifacts under `/tmp/cofolder-s3.4.8FTbTH` | Python 3.12.12/setuptools 80.10.2 built `cofolder-1.0.0.tar.gz` from a clean `git archive`, then `cofolder-1.0.0-py3-none-any.whl` only from the unpacked sdist; inventory assertions passed at 228/104 files, 51 supported test files, 10 examples and 12 acceptance resources; outside-checkout import/version/MIT/legal files/two entry points, CLI/module parity, all six tool help paths and example copy passed; focused mocked tool/setup tests: 25 passed; supported sdist suite: 551 passed in 29.63s with Streamlit absent | The first temporary dependency projection preferred the intentionally incompatible Boltz1 namespace and failed six collection imports; rebuilding the temporary environment with the established full-test dependency set first resolved that environment-only conflict. No network, backend, GPU or service work ran. S3 exit satisfied; advance to S4.1. |
| 2026-09-10 | S4.1 | `eb83b54`; working tree | `pyproject.toml`, `MANIFEST.in`, development docs, 11 lint-touched source/test files, removed `tests/tmp_test.py`, this audit | `ruff check src tests scripts tutorials examples`: clean; focused packaging/analytics/runner/input/entity/utils suite: 208 passed; `python tests/ui_development_check.py`: 6/6 passed; `pytest -q`: 552 passed in 29.80s; TOML metadata assertions and `git diff --check` passed | Added only the `development` extra while retaining Black/Ruff in `docs`; resolved all 15 findings manually without suppressions or broad formatting; removed only the approved zero-byte placeholder. S4 remains active; advance to S4.2 test lanes and CI. |
| 2026-09-10 | S4.2 | `4c8537e`; working tree | `.github/workflows/quality.yml`, `scripts/run_test_lane.py`, `scripts/verify_release_artifacts.py`, test/dependency/manifest/contributor updates, this audit | Four independent lane commands: 482 core, 22 contracts/tutorial, 14 artifact and 42 lightweight acceptance passed (560 total); `python scripts/verify_release_artifacts.py` built an sdist then its wheel and passed inventory, UI-exclusion, outside-checkout install, metadata/licence, entry-point/help and example-copy checks; Ruff, workflow regression/YAML parsing and `git diff --check` passed | The `test` extra is pytest-only, Marimo remains in `tutorials`/`acceptance`, and `development` owns `build`; routine GitHub Actions covers both supported Python versions without Streamlit, backend extras, downloads, GPU/service work or inference. Remote CI and actual backend acceptance are not claimed by this local evidence. Advance to S4.3 dependency bounds/provenance. |
| 2026-09-10 | S4.3 | `3b971ae`; working tree | `pyproject.toml`, Boltz-family runner availability checks, focused runner/packaging tests, installation guide, this audit | Focused packaging/Boltz runner tests: 42 passed; four lanes: 489 core, 22 contracts/tutorial, 15 artifact, 42 acceptance (568 total); full supported lane: 568 passed; `python scripts/verify_release_artifacts.py`, Ruff and `git diff --check` passed | Bounded Boltz2 to major 2, replaced moving Community Git provenance with exact PyPI 2.10.12, and added pre-discovery Python 3.13 diagnostics for Boltz1/2 plus exact Community validation. R08 remained unapplied; actual backend installs/inference remain S4.5. Advance to S4.4. |
| 2026-09-11 | S4.4 | `4f71d2b`; working tree | Optional-dependency boundaries, R08 metadata, supported CI/acceptance install hints, installation/tutorial guidance, packaging and analysis regressions, this audit | Focused dependency/analytics/packaging/acceptance selection: 115 passed; lanes: 500 core, 22 contracts/tutorial, 16 artifact, 42 acceptance; full supported lane: 580 passed; `python scripts/verify_release_artifacts.py`, Ruff, `git diff --check` and seven console/module parity cases passed | Applied only approved R08: Python `>=3.11,<3.13`, exact seven-package base and six-package `analysis` extra. Optional feature use now raises installation guidance; base-safe imports also required deferring the existing Boltz-only CCD parser import. No capability removed and no clean environment/backend result claimed. Advance to S4.5. |
| 2026-09-11 | S4.5 | `0aa5aa0`; working tree after verification fixes | pandas censor-sign preservation, artifact verifier isolation, Oracle optional-backend test isolation, backend metric catalog/normalizers and regressions, this audit | Final wheel SHA-256 `967cb43f45133d868580b588d0399516e4b7e862a164256e5e67b1e263e847cd`; fresh Python 3.11.16 and 3.12.14 base installs passed imports/version/two entry points/example copy/config and input validation with all analysis deps absent and actionable `cofolder[analysis]` guidance. Both full installs proved Streamlit absent and independently passed 500 core + 22 contracts/tutorial + 16 artifact + 42 acceptance tests, Ruff and sdist-to-wheel installed verification. Separate Python 3.12 environments passed `pip check`, discovery, availability, options and one-sample/repeat real GPU validation: Boltz 1.0.0 (45 metric records), Boltz 2.2.1 (55), Community 2.10.12 (57), OpenFold3 0.5.0 (43); all manifests report one success/zero failures. Complete inventories, commands, GPU data, logs and outputs: `/tmp/cofolder-s4.5.uUjsHm`. | Fresh execution exposed and resolved four release-boundary defects rather than weakening checks: pandas 3 nullable censor signs, Oracle unit coupling to Boltz2, artifact verifier same-version install reuse, and real backend confidence fields/scales (including OpenFold3 0-100 pLDDT). OpenFold setup downloaded its required default checkpoint and CCD database successfully; shared `/home/remco/.boltz` and `/home/remco/.openfold3` caches retained. S4 exit satisfied; advance to S5.1. |
| 2026-09-11 | S5.1 | `cbf429e`; working tree | CLI and four recipes; bias database/materialization and query-scoping analytics; shared preflight and ligand resolution; fetch tool; Bias/Validate/Screen/Oracle/install/configuration/quickstart/Python documentation; focused regressions; this audit | `python scripts/run_test_lane.py all`: 587 passed in 34.65s (one third-party deprecation warning); focused fetch/bundle regressions: 7 passed; Ruff and `git diff --check` passed; `python -m mkdocs build --strict` passed; `python scripts/verify_release_artifacts.py` passed sdist-to-wheel inventory, external install, entry-point/help and example-copy checks. No live RCSB/MMseqs fetch or inference was used by tests. | New source bundles are the default while explicit legacy CSV/build inputs retain legacy mode; mixing routes is rejected. Preparation supports already-built expanded sources and has bounded retries, batched persistent RCSB checkpoints, checksummed manifests, validation-only mode and atomic publication. Workflow preflight is offline/non-mutating; Screen and Oracle share exact-one-ligand inference. S5 remains active; advance to S5.2 shared executable discovery. |
| 2026-09-14 | S5.1-REOPEN | `cbf429e`; existing S5.1 working tree preserved | This audit only | Reconciled the requested 1P/1L, non-identical 2P/2L, dated/whole PDB, custom-complex and Screen/Oracle reuse stories against the current implementation. No tests were rerun and no implementation was changed during this plan correction. | The earlier S5.1 checks remain evidence but were insufficient for completion. Source bundles and exact identical-output reruns are reusable; component-level derived searches are not shared across Screen rows, changed Oracle candidates or systems. Set S5.1 back to IN_PROGRESS and S4.5 as last completed; complete the six acceptance groups before S5.2. |
| 2026-09-14 | S5.1-REOPEN-COMPLETE | `cbf429e`; working tree | Bias database/materialization/cache, all four recipes and CLI/preflight, custom-complex/matrix/cache tools, packaged examples, user/Python docs and regressions | `python scripts/run_test_lane.py all`: 609 passed in 24.37s (one third-party deprecation warning); Ruff and `git diff --check` passed; strict MkDocs and `python scripts/verify_release_artifacts.py` passed. Installed real MMseqs `01683a607f83878e95436632d73e1d7d9ae30955` ran the `createdb/search/convertalis` backfill adapter offline and produced `sequence_similarity=0.000`, plot value `0.000 < 0.25`, method `mmseqs_pident`. | Implemented only the requested pre-2023-06-01, whole, pre-2023-06-01-plus-custom and custom-only matrix. Protein and ligand thresholds are independently configurable on the 0–1 scale and applied after raw-cache reuse. No live RCSB call, download, backend inference, commit or tag was performed. Advance to S5.2. |
| 2026-09-14 | S5.2 | `c16a06d`; working tree | Shared executable resolver, MMseqs fetch/builder/analysis forwarding seams and required-failure diagnostics, installation/bias-data guidance, focused regressions, this audit | Focused resolver/tool/analytics selection: 42 passed; `python scripts/run_test_lane.py all`: 621 passed in 24.26s (one third-party deprecation warning); Ruff and `git diff --check` passed; `python -m mkdocs build --strict` passed in 2.63s; `python scripts/verify_release_artifacts.py` passed sdist-to-wheel inventory, installed import/version/tool help and example-copy checks. | Explicit argument now precedes environment, managed user vendor, checkout-only vendor and PATH candidates; only executable regular files qualify, structured rejection/selection diagnostics are retained, and space-containing paths remain single argv values. Existing private resolver seams, installer behavior, S5.1 cache identities, network behavior and subprocess commands remain intact. Advance to S5.3. |
| 2026-09-14 | S5.3 | `a36b447`; working tree | Root logging formatters, private recipe completion reporter, four recipe/CLI completion paths, focused logging/completion regressions, this audit | Focused logging/completion/CLI/Screen/Bias selection: 85 passed; `python scripts/run_test_lane.py all`: 630 passed in 25.68s (one third-party deprecation warning); Ruff, `git diff --check`, strict MkDocs, sdist-to-wheel artifact verification and seven console/module parity cases passed. | Normal console logs retain stdout with concise `LEVEL | message` formatting; debug console and files retain detailed source context. Only the outermost recipe reports canonical status/counts/result/manifest paths, including Oracle scalar or unavailable state; nested Validate runs are suppressed. No public signature, flag, exit code, schema or output file changed. Advance to S5.4. |
| 2026-09-14 | S5.4 | `d7a6735`; working tree | CLI and completion-summary usability regressions; this audit | Focused CLI/completion/resolver/example selection: 71 passed; `python scripts/run_test_lane.py all`: 641 passed in 40.46s (one third-party deprecation warning); Ruff and `git diff --check` passed; `python -m mkdocs build --strict` passed in 2.60s; `python scripts/verify_release_artifacts.py` passed the clean sdist-to-wheel inventory and installed smoke. | Real console/module root help, version, four command helps, invalid command and copied-example preflight have identical output/status. All four CLI preflights prohibit run/subprocess/network/materialization calls and leave output/cache paths absent; Screen/Oracle infer a unique ligand and reject ambiguity before writes. Canonical completion pointers and space-containing paths are covered. No inference, service, download, interface change, removal, refactor, commit, tag or publication occurred. S5 exit satisfied; advance to S6.1. |
| 2026-09-15 | S6.1 | `6d5c48b`; working tree | Private recipe execution/diagnostic/result helpers, Bias/Validate/Screen/Oracle delegations, characterization tests and this audit | New helper tests: 5 passed; focused recipe/contract selection: 154 passed; `python scripts/run_test_lane.py all`: 646 passed in 38.34s with one third-party deprecation warning; Ruff, `git diff --check`, strict MkDocs (2.47s) and `python scripts/verify_release_artifacts.py` passed, including sdist-to-wheel build and installed smoke. | Preserved public signatures, terminal bundle schemas/paths, failure identity/stage/code, seed/backend provenance, Screen's historically unvalidated runner seed adjustment and Validate's strict validation. Imports perform no subprocess/network/write work. No postprocessing move, removal, inference, service, tag or publication occurred. Advance to S6.2 identity proof and Screen postprocessing extraction. |
| 2026-09-15 | S6.2 | `800ddfa`; working tree | Private Screen postprocessor, Screen delegation wrappers, structured identity/import regressions and this audit | Focused Screen/IFP/shared-internal selection: 66 passed; `python scripts/run_test_lane.py all`: 654 passed in 34.53s with one third-party deprecation warning; Ruff, `git diff --check`, strict MkDocs (4.09s) and `python scripts/verify_release_artifacts.py` passed through clean sdist-to-wheel build and installed smoke. | Frozen resolved settings now drive metric shaping, filtering and clustering without changing public constructors, records, columns, decisions or artifacts. PDB/mmCIF structures through public records and CSV fallback produce real distance/ProLIF ligand, chain and residue identities. R06 is ready but its `_legacy` fallback remains for isolated S7.2 removal. Advance to S6.3 Bias analytics extraction. |
| 2026-09-15 | S6.3 | `ba45810`; working tree | `modules/analytics/bias/` package and compatibility facade, bias analytics/recipe tests, this audit | Focused bias/database/Bias/Validate selection: 88 passed; `python scripts/run_test_lane.py all`: 655 passed in 40.10s with one third-party deprecation warning; Ruff, `git diff --check`, strict MkDocs (2.46s) and `python scripts/verify_release_artifacts.py` passed through clean sdist-to-wheel build and installed smoke. | Split constants, references, similarity, datasets, cached enrichment, artifacts and orchestration into an acyclic package while preserving all 104 former helper names through the facade, the `apply_bias_metrics` signature, logger identity, output/artifact behavior and passive imports. Updated test patch targets to their owning modules only. No network, inference, removal, tag or publication occurred. Advance to S6.4 chemistry/cache and aggregation boundaries. |
| 2026-09-15 | S6.4 | `bde6c73`; working tree | Runner ligand-preparation boundary, analytics aggregation module, retained entity/utility facades, Validate import and characterization tests; this audit | Focused runner/entity/aggregation/Validate selection: 59 passed; `python scripts/run_test_lane.py all`: 662 passed in 34.82s with one third-party deprecation warning; Ruff and `git diff --check` passed; strict MkDocs passed in 2.24s; `python scripts/verify_release_artifacts.py` passed clean sdist-to-wheel inventory and installed smoke. | Boltz preparation now receives the typed runtime cache path and owns conformer/CCD orchestration; all nine aggregation functions live under analytics and Validate uses that boundary. Legacy `handle_conformers` and `utils.gather` remain compatible, including shared alignment patching. `Command`, cache download/population and every R05/R09/R10 removal remain untouched. Imports perform no subprocess/network/write work. No inference, service, tag or publication occurred. Advance to S6.5 structural parity verification. |
| 2026-09-15 | S6.5 | `6175f41`; working tree | Structural parity collector/fixture, sdist fixture inventory, import/dependency-boundary regressions and this audit | The same public-workflow fixture passed against temporary `6d5c48b` and current `6175f41` archives with exact normalized digests: Validate `c05c6ee2`, Screen `087af4ef`, Oracle `ca7787d2`, Bias `6245cbe7`. New parity/import tests: 3 passed; focused recipe/contract/analytics/preparation/aggregation selection: 273 passed in 59.13s; `python scripts/run_test_lane.py all`: 665 passed in 45.77s. Ruff, `git diff --check`, strict MkDocs (2.51s) and `python scripts/verify_release_artifacts.py` passed clean sdist-to-wheel inventory and installed smoke. | Only temporary roots and record timestamps are normalized; CSV/JSON order, values/states, identities, seeds/provenance, decisions and text artifacts remain exact, while PNG/PDF content is represented by required filenames and Bias checkpoint callbacks are compared explicitly. The known MDAnalysis deprecation and pandas concatenation future warning remain; the pandas warning occurs in both revisions. No semantic drift, production change, removal, inference, network/service use, commit, tag or publication occurred. S6 exit satisfied; advance to S7.1. |
| 2026-09-15 | S7.1 | `c3cc48a`; working tree | Legacy Boltz cache setup/population helpers, focused input/entity tests and this audit | Focused input/entity/ligand-preparation selection: 50 passed; `python scripts/run_test_lane.py all`: 682 passed in 45.66s with the two known MDAnalysis/pandas warnings. Ruff and `git diff --check` passed; strict MkDocs passed in 3.99s; `python scripts/verify_release_artifacts.py` passed clean sdist-to-wheel inventory and installed smoke. Existing `/home/remco/.boltz` satisfied the read-only component inspection. | Input and identifiers are validated before writes; incomplete Boltz2 caches use a cache-local lock, unique cleaned setup workspaces and post-setup verification; requested CCD destinations are honored. Tests mocked setup, so no inference, network, download, interface removal, commit, tag or publication occurred. R10 remains gated; advance to S7.2 with a fresh per-ID readiness/caller check. |

| 2026-09-15 | S7.2 | `6ad502f`; working tree | Approved compatibility removals, explicit Boltz2 cache tools, canonical aggregation/bias-builder/bias-submodule migrations, tests/docs and this audit | Focused removal/migration selection: 139 passed; `python scripts/run_test_lane.py all`: 656 passed with the two known MDAnalysis/pandas warnings; final tool/packaging selection: 35 passed; Ruff, strict MkDocs, `python scripts/verify_release_artifacts.py`, three installed-tool help smokes and `git diff --check` passed. | Applied only R01–R07/R09/R10; R08/R11 were already applied. Vector-only Screen fixtures now remain explicitly unclustered while identity-bearing distance/ProLIF coverage remains. Cache population performs no implicit setup. No network, inference, commit, tag or publication occurred. Advance to S7.3. |
| 2026-09-17 | S7.3 | `9d40404`; working tree | This audit; generated `site/` verification output trashed after the check | Per-ID path/reference searches and import probes passed; focused cleanup/capability selection: 227 passed in 18.87s with the known MDAnalysis warning; `python scripts/run_test_lane.py all`: 657 passed in 39.13s with the known MDAnalysis and pandas warnings; Ruff, strict MkDocs (2.12s), `python scripts/verify_release_artifacts.py` and `git diff --check` passed. | No live dangling references or lost retained capabilities found. Dated audits retain historical citations; current matches are negative assertions/exclusions. R05/R06/R07 prerequisites predate removal; R10 cache hardening predates its atomic explicit-API replacement. No network, inference, commit, tag or publication occurred. S7 exit satisfied; advance to S8.1. |
| 2026-09-17 | S8.1 | `fcae176`; audit working tree | This audit only | Ruff clean; independent lanes passed: core 577 with the two known warnings, contracts/tutorial 22, lightweight acceptance 42; explicit command-reporting security tests 5 passed; eight real console/module help/version/error/preflight parity cases passed; `git diff --check` passed. | `fcae176` differs from product-code candidate `9d40404` only by the S7 audit record. No product behavior, network, inference, tag or publication action occurred. Advance to S8.2 clean distribution verification. |
| 2026-09-17 | S8.2 | `fcae176`; audit working tree | This audit; temporary artifacts `/tmp/cofolder-s8.2.NmZmzR` | `python scripts/verify_release_artifacts.py` passed; retained isolated builds produced a 312-entry sdist (`28906ba...4f0f9`) and, only from its unpacked source, a 125-entry wheel (`5a2ee5...08497`). Installed wheel metadata/licence, both entry points, tool help, example copy, 12 acceptance and 12 example resources passed; UI inventory was zero. A wheel-installed environment ran all 657 supported tests from the unpacked sdist with the two known warnings. | The first ambient sdist run resolved the checkout's editable install, so it was not used as isolation evidence; the repeated run proved imports came from the S8 wheel environment. No product change, inference, tag or publication occurred. Advance to S8.3 fresh environments and real backends. |
| 2026-09-19 | S8.3-FIX | `fcae176`; corrected working-tree candidate | Protein pairwise-similarity implementation/regression, structural-parity isolation/fixture, this audit | Fresh pandas 2 and pandas 3 focused runs pass the new duplicate-index/non-identical-query regression and three structural-parity tests. Final `python scripts/run_test_lane.py all`: 658 passed; Ruff, explicit 5-test redaction lane, 8-test CLI parity selection and `git diff --check` passed. `python scripts/verify_release_artifacts.py` passed after the correction. | Fresh environments exposed a real pandas view-mutation leak and two fixture dependencies on ambient Boltz2 availability/version. Copied positional assignment now prevents source mutation and duplicate-index alignment errors; parity mocks only environmental availability/provenance, and the corrected Bias digest changes from `6245cbe7...` to `3aed2e50...`. All affected gates were rerun rather than waiving the failures. |
| 2026-09-19 | S8.3 | corrected S8 working tree | Temporary environments, backend logs/results under `/tmp/cofolder-s8.3`; no tracked environment output | CPython 3.11.15 and 3.12.13 base/full installs pass `pip check`, package-origin/version/entry-point/resource/config/input/no-Streamlit checks; both full environments pass Ruff and all 658 tests. Isolated real GPU runs pass Boltz1 1.0.0 (45 metrics), Boltz2 2.2.1 (55), Community 2.10.12 (57), OpenFold3 0.4.1 (43); every manifest reports `success`, one success and zero failures. | Separate prefixes prevent conflicting `boltz` namespaces. Existing validated `/home/remco/.boltz` and `/home/remco/.openfold3` caches were reused; OpenFold3 used no live MSA server. No missing hardware/model/service result was counted as a pass. Advance to S8.4. |
| 2026-09-19 | S8.4 | corrected S8 working tree | `docs/release/v1.0-claim-matrix.md`, this audit | PI-01–PI-10, IN-01–IN-12, WF-01–WF-20 and PK-01–PK-08 were each reconciled against the fresh 658-test, artifact, environment and real-backend evidence. R01–R11 and all retained additional capability groups were rechecked. `python -m mkdocs build --strict` passed with the known informational Material/MkDocs 2 notice and unnavlisted-page report. | PK-02's installed/runtime/artifact version passes at 1.0.0; tag and final documentation alignment remain explicit handoff items. No removal decision changed and no new removal was made. Advance to S8.5. |
| 2026-09-19 | S8.5 | corrected S8 working tree | This audit; retained handoff artifacts `/tmp/cofolder-s8-handoff.PZc2KT` | Final `python scripts/verify_release_artifacts.py` passed. A retained 312-entry sdist (`55674a8e...a0b9b`) exclusively produced the 125-entry wheel (`bb399581...45e9`); installed metadata/licence, resources, UI exclusion and entry-point smokes pass. Final Ruff, complete 658-test lane, strict docs and diff checks pass as recorded in S8. | Every required S8 matrix row is PASS, so S8 and `overall_status` are COMPLETE. Known warnings, one-sample/cache/service limits, unbounded OpenFold3 provenance and the exact documentation/tag/publication resume point are recorded above. No commit, tag or publication was performed. |
| 2026-09-19 | ISSUE-24 | post-S8 working tree | Removed `.streamlit/`, `run_ui.sh`, `src/cofolder/ui/`, `tests/ui_development_check.py`; updated repository tour, artifact guard, packaging regression and this audit | Packaging: 16 passed; complete supported lane: 658 passed with the two previously recorded warnings; Ruff, strict MkDocs, clean sdist-to-wheel installed verification and `git diff --check` passed. Built artifacts contain no UI package, Streamlit dependency, UI extra or UI entry point. | The user approved full prototype retirement, superseding only R11's permission to retain warned development material. The four supported CLI/Python workflows are unchanged. Close issue 24 as `not planned` only after this change is merged; do not create a replacement UI issue. |

Completion means the P1 defects are resolved with meaningful regressions, the
intended runtime/optional environments and artifact contents are tested, every
retained capability is accounted for, and every implemented removal has specific
approval. A clean lint run, fewer files, or a green unit suite alone is not sufficient.
