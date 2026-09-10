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
plan_version: 2
overall_status: IN_PROGRESS
active_phase: S4
active_task: S4.4
last_completed_task: S4.3
next_action: "Start S4.4 by adding lazy analysis feature boundaries and installation guidance, then apply only the approved R08 Python bound and six-package base-to-analysis move."
blocked_on: []
implementation_base_commit: e9642f25eb358d44e2a97edcb4203d96d485c14b
last_updated: 2026-09-10
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
| S4 | Reproducible quality and supported environments (C06, C07) | S3 | **IN_PROGRESS** | Repository-owned lint/test/package lanes pass in documented environments |
| S5 | Additive CLI/Python usability and discovery improvements (C11–C13) | S3, S4 | **NOT_STARTED** | CLI/Python journeys and resolver/result-summary tests pass |
| S6 | Behavior-preserving responsibility extraction (C09, C10) | S1–S5 | **NOT_STARTED** | Before/after contract fixtures and full suite show semantic parity |
| S7 | Permission-gated cleanup and cache/logging changes (C08, C14, C15; R01–R11) | S0 decision gate, relevant earlier phases | **NOT_STARTED** | Only approved removals applied; each has specific regression evidence |
| S8 | Release-candidate verification and claim reconciliation | S1–S7 | **NOT_STARTED** | Required matrix below is complete; unresolved external checks are explicit |

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
| R01 | **APPROVED; READY WITH MIGRATION** | Three documentation dependencies enumerated above; repair them with deletion | User: "accept R01-R10 (all)", 2026-09-09 |
| R02 | **APPROVED; READY WITH MIGRATION** | Definitions, export, alias test, documentation mention and canonical types enumerated | User: "accept R01-R10 (all)", 2026-09-09 |
| R03 | **APPROVED; READY WITH MIGRATION** | Every wrapper, repository caller and canonical implementation enumerated | User: "accept R01-R10 (all)", 2026-09-09 |
| R04 | **APPROVED; READY** | Correct spelling confirmed across user invocations; both generated misspelled forms identified | User: "accept R01-R10 (all)", 2026-09-09 |
| R05 | **APPROVED; GATED** | Await S1 typed migration and published replacement signatures | User: "accept R01-R10 (all)", 2026-09-09 |
| R06 | **APPROVED; GATED** | Await S6 proof of identity-bearing fingerprints on all current workflow paths | User: "accept R01-R10 (all)", 2026-09-09 |
| R07 | **APPROVED; GATED** | Await supported pytest capture and handler-preservation evidence | User: "accept R01-R10 (all)", 2026-09-09 |
| R08 | **APPROVED; READY FOR S4** | Exact Python bound and six-dependency `analysis` move recorded above | User: "accept R01-R10 (all)", 2026-09-09 |
| R09 | **APPROVED; GATED BY RELATED PHASES** | Exact finite old-to-new path table recorded; each row waits for its canonical implementation and migrated callers | User: "accept R01-R10 (all)", 2026-09-09 |
| R10 | **APPROVED; GATED** | Await explicit setup documentation/tests and retained CCD conversion contract | User: "accept R01-R10 (all)", 2026-09-09 |
| R11 | **APPROVED; READY FOR S2/S3** | Issue 24 scope supplied; exact documentation, wheel, dependency and release-test exclusions recorded above | User: "The UI will not be shipped", 2026-09-10 |

### Task ledger

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
| S4.4 | R08-dependent environment changes | Apply the approved exact `>=3.11,<3.13` bound and six-package base→`analysis` move after adding lazy feature boundaries and installation guidance | **NOT_STARTED** | — |
| S4.5 | Environment verification | Fresh 3.11/3.12 base installs, full supported test environment without Streamlit, optional-analysis guidance, and four separate backend availability/smoke environments pass or have resource-specific blockers recorded | **NOT_STARTED** | — |
| S5.1 | Easier CLI journey | Group help without removing flags; add shared preflight-only validation, unique-ligand inference, explicit bias input/output fields, contradiction checks and bundled-example path | **NOT_STARTED** | — |
| S5.2 | Shared executable discovery | Implement one tested MMseqs resolver with documented precedence, executable rejection reasons, space-safe paths, source and installed layouts; preserve all lookup routes | **NOT_STARTED** | — |
| S5.3 | Progress and result summaries | Add concise console versus detailed debug/file formatting and shared completion summaries for success, partial/total failure and unavailable Oracle metrics; preserve output files and stream contracts | **NOT_STARTED** | — |
| S5.4 | Usability verification | Test console/module parity, no-service preflight, ambiguity failures, resolver precedence, paths with spaces and canonical result/manifest pointers | **NOT_STARTED** | — |
| S6.1 | Workflow shared internals | Extract planning/lifecycle, diagnostics and public-result assembly in small diffs; preserve public recipe classes/signatures and failure/identity semantics | **NOT_STARTED** | — |
| S6.2 | Screen postprocessing | Extract ranking/filtering/clustering while preserving every taxonomy, decision and output family; remove the approved vector fallback only after the R06 identity proof | **NOT_STARTED** | — |
| S6.3 | Bias analytics package | Split reference loading, similarity, enrichment, pair datasets and artifacts behind a stable `apply_bias_metrics`; imports must not trigger workflows/network | **NOT_STARTED** | — |
| S6.4 | Chemistry/cache and aggregation boundaries | Separate chemistry from backend cache orchestration and move aggregation to analytics; follow the approved R05/R09/R10 path table only after each readiness gate | **NOT_STARTED** | — |
| S6.5 | Structural parity verification | Compare fixture outputs before/after for counts, identities, seeds, values/states, provenance, gates, filters/clusters and artifact families; run import-boundary checks and full suite | **NOT_STARTED** | — |
| S7.1 | Safe cache corrections independent of R10 | Validate input before writes, inspect required cache components, use unique temporary workspaces, and test empty/partial/complete/concurrent/error cases without changing the documented automatic behavior | **NOT_STARTED** | — |
| S7.2 | Execute approved removals | Implement only `APPROVED` R01–R11 items, one ID per reviewable change; update canonical callers/tests/docs and retain the underlying capability safeguards | **NOT_STARTED** | — |
| S7.3 | Cleanup verification | Search for dangling references per executed ID; prove gated interfaces remained until their prerequisites passed; run relevant capability tests and full suite | **NOT_STARTED** | — |
| S8.1 | Static and unit release gate | Repository lint, supported core/contract/tutorial tests and CLI/module parity pass without Streamlit in the declared test environment | **NOT_STARTED** | — |
| S8.2 | Distribution release gate | Clean sdist→wheel build, explicit inventory, installed metadata/licence, external-checkout smoke and sdist test execution pass | **NOT_STARTED** | — |
| S8.3 | Environment/backend release gate | Supported Python/base/analysis matrix and actual Boltz1/Boltz2/Community/OpenFold3 acceptance run in isolated environments; record unavailable hardware/service/model resources without claiming a pass | **NOT_STARTED** | — |
| S8.4 | Claim and removal reconciliation | Re-evaluate every affected PI/IN/WF/PK claim, all additional retained capabilities and R01–R11; link evidence and list any residual risk | **NOT_STARTED** | — |
| S8.5 | Handoff | Set `overall_status` accurately; summarize shipped changes, deferred documentation/tag/publication work and exact blockers. Mark `COMPLETE` only if all required gates pass | **NOT_STARTED** | — |

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
  controls, stream behavior and executable lookup routes remain available.
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
| Source | Repository-owned Ruff baseline | **NOT_RUN** | — |
| Source | Full supported core/contracts/tutorial test lanes without Streamlit | **NOT_RUN** | — |
| CLI | Console/module help, version, invalid command and preflight parity | **NOT_RUN** | — |
| Security | Success/failure secret-redaction tests for all command reporters | **NOT_RUN** | — |
| Artifact | Sdist inventory and wheel built from unpacked sdist | **NOT_RUN** | — |
| Artifact | Installed-outside-checkout resources, entry points, metadata and licence | **NOT_RUN** | — |
| Scope | UI absent from wheel, supported extras/entry points/docs and release test requirements | **NOT_RUN** | — |
| Environment | Fresh Python 3.11 base/full lanes | **NOT_RUN** | — |
| Environment | Fresh Python 3.12 base/full lanes | **NOT_RUN** | — |
| Backend | Boltz1 actual acceptance | **NOT_RUN** | — |
| Backend | Boltz2 actual acceptance | **NOT_RUN** | — |
| Backend | Boltz Community actual acceptance | **NOT_RUN** | — |
| Backend | OpenFold3 actual acceptance | **NOT_RUN** | — |
| Claims | PI/IN/WF/PK and additional-capability reconciliation | **NOT_RUN** | — |
| Removals | Approved/rejected/pending R01–R11 reconciled | **NOT_RUN** | — |

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

Completion means the P1 defects are resolved with meaningful regressions, the
intended runtime/optional environments and artifact contents are tested, every
retained capability is accounted for, and every implemented removal has specific
approval. A clean lint run, fewer files, or a green unit suite alone is not sufficient.
