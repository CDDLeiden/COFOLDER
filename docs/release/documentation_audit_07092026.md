# COFOLDER Documentation Audit — 2026-09-07

Scope: PK and DOC requirements from `v1.0-claim-matrix.md`.

Verification: `mkdocs build --strict` completed successfully; no wheel/sdist release test was found.

| Done | ID | Status | Evidence | Notes/Gaps |
|---|---|---|---|---|
| ~ | PK-01 | [Partially Implemented] | `pyproject.toml:68`, `docs/getting-started/installation.md:11` | Installation documentation and entry points exist, but no clean-environment install acceptance test was found. |
| — | PK-02 | [Missing] |  | Package/runtime version is `0.1.0`; no v1.0.0 tag or documentation version was found. |
| ✓ | PK-03 | [Implemented] | `pyproject.toml:15`, `pyproject.toml:30`, `docs/getting-started/installation.md:33` | Core and runner extras are declared and installation guidance is provided. |
| ✓ | PK-04 | [Implemented] | `src/cofolder/__init__.py:1`, `src/cofolder/modules/runners/base.py:22` | Core import does not require optional runner distributions; availability is checked on use. |
| — | PK-05 | [Missing] |  | No wheel/sdist build, clean-install, or artifact smoke-test automation was found. |
| ~ | PK-06 | [Partially Implemented] | `LICENSE`, `README.md:143` | MIT license file and README statement exist, but package metadata lacks a license field. |
| — | PK-07 | [Missing] |  | Only acceptance data is declared as package data; examples/tutorial inputs are not declared for distribution and no artifact test exists. |
| ~ | PK-08 | [Partially Implemented] | `docs/getting-started/installation.md:46`, `src/cofolder/modules/input/command.py:82` | Runner setup and external MSA behavior are documented, but no comprehensive release audit covers undeclared paths/resources/credentials. |
| ✓ | DOC-01 | [Implemented] | `README.md:1`, `README.md:5` | README distinguishes standalone bias from runner-backed workflows. |
| ✓ | DOC-02 | [Implemented] | `docs/getting-started/installation.md:3`, `README.md:29` | Covers supported Python, extras, runner setup, and verification. |
| ~ | DOC-03 | [Partially Implemented] | `src/cofolder/cli.py:671`, `docs/user-guide` | Per-workflow docs exist, but no complete CLI reference documents every default/exit behavior or states module-entry equivalence. |
| ~ | DOC-04 | [Partially Implemented] | `docs/api`, `src/cofolder/modules/runners/__init__.py:1` | Recipe APIs are generated, but runner selection/execution, analytics, and aggregation lack a clearly defined supported public API surface. |
| ~ | DOC-05 | [Partially Implemented] | `docs/getting-started/configuration.md:1` | Components and constraints are documented, but MSA/reference schema and invalid-rule examples are incomplete. |
| ~ | DOC-06 | [Partially Implemented] | `docs/getting-started/configuration.md:118`, `src/cofolder/cli.py:82` | Basic options are documented, but the complete runner/options, diagnostics, Oracle-gate, reducer, and penalty contract is fragmented/incomplete. |
| ~ | DOC-07 | [Partially Implemented] | `docs/tutorials/ligands.md:9`, `src/cofolder/modules/entities/ligand.py:390` | Conversion/SDF/CCD material exists, but current documented system routes overstate PDB/CIF support and do not fully define multi-record policy. |
| ~ | DOC-08 | [Partially Implemented] | `docs/getting-started/configuration.md:65`, `docs/tutorials/runners.md:309` | Matrix covers entities/constraints, but not all MSA, confidence, binding, affinity, and limitation fields for each runner. |
| ~ | DOC-09 | [Partially Implemented] | `docs/user-guide/screen.md:193`, `docs/tutorials/runners.md:202` | Partial normalized and Screen schemas exist; no complete field-level public output schema. |
| ✓ | DOC-10 | [Implemented] | `docs/user-guide/bias.md:1` | Covers references, provenance, standalone execution, interpretation, and non-correctness meaning. |
| ~ | DOC-11 | [Partially Implemented] | `docs/user-guide/validate.md:1`, `src/cofolder/modules/analytics/reproduction.py:897` | Inputs are described, but the three evidence regimes and applicability matrix are not explicitly documented. |
| ~ | DOC-12 | [Partially Implemented] | `docs/user-guide/screen.md:1`, `src/cofolder/recipes/screen.py:296` | Strong coverage of CSV behavior, MSA, failures, and outputs; docs claim SDF support that the workflow does not implement. |
| ✓ | DOC-13 | [Implemented] | `docs/user-guide/oracle.md:1`, `src/cofolder/recipes/oracle.py:60` | Documents reducers, custom scoring, gates, penalties, and audit output. |
| ~ | DOC-14 | [Partially Implemented] | `README.md:11`, `docs/user-guide/oracle.md:151` | Correct non-universal framing exists, but no consolidated documentation explains all listed metric categories and scales. |
| ~ | DOC-15 | [Partially Implemented] | `examples`, `tests/test_tutorials.py:17` | Examples exist, but tests only import tutorial modules; they do not schema-validate or execute all example configurations. |
| ~ | DOC-16 | [Partially Implemented] | `tutorials`, `tests/test_tutorials.py:31` | Bias/Validate/Screen/Oracle/API tutorial files exist, but no end-to-end execution test validates them. |
| ~ | DOC-17 | [Partially Implemented] | `tutorials/README.md:1`, `tests/test_tutorials.py:31` | Marimo modules are import-tested, not sequentially executed from a clean kernel. |
| ~ | DOC-18 | [Partially Implemented] | `docs/getting-started/installation.md:223`, `docs/tutorials/screening.md:206` | Troubleshooting is distributed across pages; it does not comprehensively cover all required failure categories. |
| — | DOC-19 | [Missing] |  | Documentation and metadata identify 0.1.0, with no v1.0.0 release identity, citation, or software archive reference. |
| ~ | DOC-20 | [Partially Implemented] | `mkdocs.yml:1`, `tests/test_tutorials.py:31` | `mkdocs build --strict` passed during audit, but automated checks do not verify documented commands/API signatures against the installed package. |
