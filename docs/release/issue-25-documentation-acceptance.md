# Issue 25 documentation acceptance

Issue: [#25 — COFOLDER Documentation Plan](https://github.com/CDDLeiden/COFOLDER/issues/25)

This record covers the publication-facing documentation gate for COFOLDER 1.0.0.
The final candidate was revalidated on 2026-10-06 at
`b8906443ad3851f9359c420d64465ccb9e2821a6`. The acceptance-record update that
follows that commit changes documentation evidence only.

## Checklist reconciliation

| Area | Status | Evidence |
| --- | --- | --- |
| MkDocs, Material, search, themes, navigation, highlighting, and mkdocstrings | Complete | `mkdocs.yml`; strict build |
| Getting started, workflow guides, tutorials, API reference, contributing, and changelog | Complete | MkDocs navigation and public pages |
| Current workflows | Complete | Public documentation covers `bias`, `validate`, `screen`, and `oracle`; former Predict/Evaluate wording is superseded |
| Installation and configuration | Complete | Clean-install guidance, strict system/options separation, ligand XOR rule, MSA and constraint documentation |
| CLI reference | Complete | Consolidated command reference, parser help, exit statuses, output locations, and console/module parity tests |
| Screen input scope | Superseded | The current implementation and docs support CSV, SDF, and MOL libraries rather than the older CSV-only proposal |
| Ligand handling | Complete | SMILES, Validate SDF conformers, Screen SDF/MOL libraries, standard/custom CCD, identifiers, stereochemistry, and conversion utilities |
| Tutorials and advanced examples | Complete | Issue [#16 acceptance](issue-16-tutorial-acceptance.md), packaged examples, preflight tests, and fictional-API regression scan |
| Metrics and outputs | Complete | Public output contract and catalog-derived metric reference document state, scope, units, direction, evidence, and unavailable-value semantics |
| Runner capabilities | Complete | Boltz1, Boltz2, Boltz Community, and OpenFold3 matrix plus backend acceptance documentation |
| Troubleshooting | Complete | Consolidated failure and recovery guide |
| Release identity | Complete with explicit follow-up | Version 1.0.0, repository, MIT licence, changelog, and `CITATION.cff`; final publication/archive identifiers are deferred to [#40](https://github.com/CDDLeiden/COFOLDER/issues/40) |
| `setuptools` minimum | No change required | `setuptools >= 77.0.3` is a compatible minimum, not a pin |
| GitHub Pages | Deferred | A hosted deployment is not required by the manuscript or the v1.0.0 release gate |

## Automated release gates

The documentation job in `.github/workflows/quality.yml` installs `.[docs]` and
runs `python -m mkdocs build --strict`. The contracts/tutorial lane additionally:

- validates every packaged `system*.yaml` with the packaged options file;
- accepts each single ligand representation and rejects missing/conflicting forms;
- executes the backend-free Bias example;
- preflights Validate, Screen CSV/SDF/MOL, Oracle, and custom-CCD commands;
- checks packaged/canonical example parity;
- checks the metric table against `METRIC_CATALOG`;
- rejects legacy names and removed or fictional tutorial APIs;
- checks installed console/module help and exit behavior through the CLI suite.

## Candidate verification

The complete gate was rerun from the repository root on 2026-10-06 with Python
3.12.4 after the final mapped-screening, ProLIF/clustering, and structure-gated
Oracle documentation changes.

| Command | Result |
| --- | --- |
| `ruff check src tests scripts tutorials examples` | PASS |
| `python scripts/run_test_lane.py all` | PASS — 744 tests, 2 known third-party/future warnings |
| `python -m mkdocs build --strict` | PASS — 57 HTML output files |
| `python scripts/verify_release_artifacts.py` | PASS — sdist-to-wheel build, inventory, outside-checkout install, metadata, entry points, tool help, and example copy |
| `git diff --check` | PASS — checked after this evidence update |

Routine CI matrices the documentation/example contracts over Python 3.11 and 3.12.
The existing issue #16 record supplies real tutorial workflow evidence. This change
alters input validation and documentation, not backend execution or dependencies, so
it does not require repeating GPU inference.

## Closure conditions

Issue #25 can close when this evidence update is published. The final
publication/archive metadata remains tracked by
[#40](https://github.com/CDDLeiden/COFOLDER/issues/40) and does not block this
documentation gate. Tagging, artifact publication, and GitHub Pages deployment are
separate release actions.
