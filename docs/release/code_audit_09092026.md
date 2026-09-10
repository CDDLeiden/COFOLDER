# COFOLDER Code Audit — 2026-09-09

Scope: PI, IN, and WF requirements from `v1.0-claim-matrix.md`.

Verification performed on `feature/release_features` at `e9642f2`:

- `pytest -q` — **520 passed**.
- `cofolder -h`, `cofolder {validate,screen,oracle,bias} -h`, and their
  `python -m cofolder` counterparts all succeeded.
- Both entry points rejected an unknown command with exit status 2 and emitted
  identical top-level help.

The 2026-09-07 findings have been re-evaluated against the current public
contracts and regression tests. All requirements in this scope are met.

| ID | Assessment | Evidence | Audit conclusion |
|---|---|---|---|
| PI-01 | Met | `src/cofolder/cli.py`; `tests/test_cli.py` | The four documented commands have help and invalid invocation fails non-zero. |
| PI-02 | Met | `pyproject.toml`; `src/cofolder/__main__.py` | Console script and module invocation call the same `main()`. |
| PI-03 | Met | `src/cofolder/recipes/`; `src/cofolder/modules/contracts/__init__.py` | Recipe, input, runner, analytics, contract, and aggregation APIs are importable independently of the CLI. |
| PI-04 | Met | `src/cofolder/recipes/bias.py`; `tests/recipes/test_bias.py` | Bias writes machine-readable outputs and does not construct or run a co-folding backend. |
| PI-05 | Met | `src/cofolder/recipes/validate.py`; `src/cofolder/modules/contracts/adapters.py` | Validate publishes applicable diagnostics and represents unavailable capability/evidence groups explicitly without suppressing independent diagnostics. |
| PI-06 | Met | `src/cofolder/recipes/screen.py`; `tests/recipes/test_screen.py` | Screen emits validated execution-level public records for every compound/repeat/model and supports CSV, SDF, and MOL sources. |
| PI-07 | Met | `src/cofolder/recipes/oracle.py`; `tests/recipes/test_oracle.py` | Oracle replaces only the selected ligand, reduces repeat results, applies gates, and publishes/returns its scalar. |
| PI-08 | Met | `src/cofolder/modules/runners/__init__.py`; `src/cofolder/modules/contracts/models.py`; `tests/modules/runners/test_contracts.py` | All four runners are selectable; public manifests and identities retain detected/unavailable/unparseable backend version and capabilities. |
| PI-09 | Met | `src/cofolder/modules/contracts/metrics.py`; `docs/user-guide/public-output-contract.md` | Metric records distinguish computed, missing, unsupported, failed, and not-requested states with reasons. |
| PI-10 | Met | `src/cofolder/modules/contracts/models.py`; `tests/modules/contracts/test_public_contracts.py` | Validated public identities carry workflow, system, compound, repeat, model/sample, entity/chain, and runner fields. |
| IN-01 | Met | `src/cofolder/modules/input/config.py`; `tests/modules/input/test_config.py`; `tests/recipes/test_validate.py` | System and options YAML are schema-validated before runner execution; malformed/missing input yields actionable typed failures. |
| IN-02 | Met | `src/cofolder/modules/input/validation.py`; `src/cofolder/modules/contracts/models.py` | Runner-backed workflows require protein and ligand entities and preserve entity/chain identity in public records. |
| IN-03 | Met | `src/cofolder/modules/input/validation.py`; `tests/modules/input/test_validation.py` | Protein alphabet, position, collisions, empty sequences, MSA format/content, and MSA-query matching are checked before execution. |
| IN-04 | Met | `src/cofolder/modules/input/ligand.py`; `tests/modules/entities/test_ligand.py` | SMILES is RDKit-validated/normalized and conversion failures are explicit. |
| IN-05 | Met | `src/cofolder/modules/input/compound_library.py`; `tests/modules/input/test_compound_library.py` | MOL/SDF and CSV converge on the ligand contract; each malformed structural source record becomes a source-mapped failure. |
| IN-06 | Met | `src/cofolder/modules/entities/ligand.py`; `src/cofolder/modules/input/validation.py` | CCD input and atom constraints are retained and validated; supported simpler ligand input can be converted. |
| IN-07 | Met | `src/cofolder/modules/input/system.py`; `src/cofolder/modules/runners/validators.py`; `tests/modules/runners/test_openfold3_runner.py` | DNA/RNA are modeled as supported entity types and runner capabilities are checked before launch. |
| IN-08 | Met | `src/cofolder/modules/input/ligand.py`; `tests/recipes/test_screen.py`; `tests/modules/input/test_system.py` | Replacement asserts that only selected ligand chemistry changes; unrelated constraints round-trip and invalidated constraints fail precisely. |
| IN-09 | Met | `src/cofolder/modules/analytics/reproduction.py`; `src/cofolder/modules/analytics/reference_ifp.py`; `tests/modules/analytics/test_reference_ifp.py` | PDB/CIF reference structures support structural metrics and reference-complex IFP extraction/filtering. |
| IN-10 | Met | `src/cofolder/recipes/{validate,screen,oracle}.py`; `tests/recipes/test_oracle.py` | Custom pockets operate without a reference complex across applicable diagnostics, filters, and gates. |
| IN-11 | Met | `src/cofolder/modules/input/compound_library.py`; `docs/user-guide/screen.md`; `tests/modules/input/test_compound_library.py` | CSV, SDF, and MOL inputs use stable source/execution identities with documented deterministic duplicate-ID policy. |
| IN-12 | Met | `src/cofolder/recipes/bias.py`; `src/cofolder/modules/analytics/bias.py`; `tests/recipes/test_bias.py` | Training references are validated and their provenance is retained. |
| WF-01 | Met | `src/cofolder/modules/runners/contracts.py`; `src/cofolder/recipes/validate.py` | Execution plans and public output identities retain repeat and model/sample dimensions. |
| WF-02 | Met | `src/cofolder/modules/utils/helpers.py`; `src/cofolder/modules/contracts/models.py`; `tests/modules/utils/test_helpers.py` | Seed plans deterministically map distinct generated/user seeds to repeats and publish the effective seed. |
| WF-03 | Met | `src/cofolder/modules/contracts/models.py`; `src/cofolder/modules/contracts/metrics.py`; `docs/user-guide/public-output-contract.md` | Metrics are explicitly classified as reference-structure, custom-pocket, or reference-free evidence. |
| WF-04 | Met | `src/cofolder/modules/input/ligand.py`; `tests/recipes/test_screen.py` | Screen changes only the selected ligand chemistry in a deep-copied fixed system and protects unrelated fields/constraints. |
| WF-05 | Met | `src/cofolder/recipes/screen.py`; `tests/recipes/test_screen.py` | Source-indexed isolated directories prevent collisions; failed members become records while unrelated members continue. |
| WF-06 | Met | `src/cofolder/modules/runners/msa.py`; `tests/recipes/test_screen.py` | Precomputed MSAs are preserved and generated MSAs are cached/reused by eligible screen members/repeats. |
| WF-07 | Met | `src/cofolder/modules/runners/{boltz_runner,openfold3_runner}.py`; `src/cofolder/modules/contracts/metrics.py` | Capability-controlled confidence, binding, and affinity metrics are collected in distinct metric classes. |
| WF-08 | Met | `src/cofolder/modules/analytics/reproduction.py`; `tests/modules/analytics/test_reproduction.py` | Structural metrics retain alignment/mapping behavior and explicit unavailable results. |
| WF-09 | Met | `src/cofolder/modules/analytics/{reference_ifp,ifp_filtering,ifp_clustering}.py`; `tests/modules/analytics/test_reference_ifp.py` | Interpretable IFPs support custom/reference filtering and reference-free clustering. |
| WF-10 | Met | `src/cofolder/modules/analytics/structure.py`; `src/cofolder/recipes/oracle.py` | Ligand SASA and heavy-atom-normalized SASA are workflow metrics and available to Oracle gates. |
| WF-11 | Met | `src/cofolder/modules/analytics/bias.py`; `docs/user-guide/bias.md` | Protein and ECFP4 proximity carry provenance and remain contextual diagnostics. |
| WF-12 | Met | `src/cofolder/modules/utils/gather.py`; `tests/modules/utils/test_gather.py` | Aggregation computes metric variation and pairwise structural agreement with correct directionality. |
| WF-13 | Met | `src/cofolder/recipes/screen.py`; `src/cofolder/modules/contracts/models.py` | Requested workflow metrics have stable records and states that distinguish zero from absent/unsupported/failed/not-requested values. |
| WF-14 | Met | `src/cofolder/recipes/oracle.py`; `tests/recipes/test_oracle.py` | Oracle supports exactly `first`, `mean`, `max`, `min`, and `median`, validating unavailable/non-numeric inputs. |
| WF-15 | Met | `src/cofolder/recipes/oracle.py`; `tests/recipes/test_oracle.py` | Python callers may supply custom multi-component scalar functions/weighted composites. |
| WF-16 | Met | `src/cofolder/recipes/oracle.py`; `src/cofolder/modules/analytics/reference_ifp.py` | IFP, SASA, custom-pocket, and reference-derived gates use the same proposed-ligand execution context. |
| WF-17 | Met | `src/cofolder/recipes/oracle.py`; `tests/recipes/test_oracle.py` | Down-weight, fixed-penalty, and non-binder policies write base/reduced values, gate evidence, and applied action. |
| WF-18 | Met | `src/cofolder/recipes/oracle.py`; `docs/user-guide/oracle.md` | Thresholds are caller-configured; no metric is universalized as a pass/fail rule. |
| WF-19 | Met | `src/cofolder/modules/contracts/metrics.py`; `docs/user-guide/public-output-contract.md`; `tests/modules/contracts/test_public_contracts.py` | `METRIC_CATALOG` enforces stable names, definitions, class, direction, units/scales, evidence compatibility, and null-with-reason missing conventions. |
| WF-20 | Met | `src/cofolder/modules/contracts/models.py`; `src/cofolder/modules/contracts/validation.py`; `tests/modules/contracts/test_public_contracts.py`; `tests/recipes/test_validate.py` | Typed failure records retain known identity context, explicit lifecycle stage, stable error code, message, and details without discarding partial successes. |

Result: **42/42 PI, IN, and WF requirements met; 0 partial or unmet findings.**
