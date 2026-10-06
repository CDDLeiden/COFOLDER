# COFOLDER v1.0 Documentation Claim Audit — 2026-09-21

Scope: documentation claims DOC-01 through DOC-20 in the
[v1.0 claim matrix](v1.0-claim-matrix.md). This audit does not re-audit source-code,
backend, packaging, or release-tag claims; those remain for a separate commit window.

## Result

Nineteen documentation claims are verified. DOC-19 remains open because a final
publication/software-archive identifier has not yet been issued. The public pages
correctly identify version `1.0.0`, the repository, MIT licence, and citation metadata,
but they deliberately do not invent an archive identifier.

| Claims | Status | Primary evidence |
| --- | --- | --- |
| DOC-01 | Verified | [README](https://github.com/CDDLeiden/COFOLDER/blob/main/README.md), [workflow overview](../user-guide/overview.md) |
| DOC-02 | Verified | [Installation](../getting-started/installation.md), [backend acceptance](../tutorials/backend-acceptance.md) |
| DOC-03 | Verified | [CLI reference](../reference/cli.md) |
| DOC-04 | Verified | [Python API execution flows](../reference/python-api/index.md), generated API reference |
| DOC-05–DOC-06 | Verified | [Configuration](../getting-started/configuration.md), [CLI reference](../reference/cli.md), [Oracle](../user-guide/oracle.md) |
| DOC-07 | Verified | [Ligand handling](../tutorials/ligands.md) |
| DOC-08 | Verified | [Runner capability matrix](../tutorials/runners.md#supported-runner-capabilities) |
| DOC-09 | Verified | [Public output contract](../user-guide/public-output-contract.md), [metric reference](../user-guide/metric-reference.md) |
| DOC-10 | Verified | [Bias](../user-guide/bias.md), [bias training data](../user-guide/bias-training-data.md) |
| DOC-11 | Verified | [Validate evidence regimes](../user-guide/validate.md#evidence-regimes) |
| DOC-12 | Verified | [Screen](../user-guide/screen.md) |
| DOC-13 | Verified | [Oracle](../user-guide/oracle.md) |
| DOC-14 | Verified | [Metric reference](../user-guide/metric-reference.md), [public output contract](../user-guide/public-output-contract.md) |
| DOC-15–DOC-17 | Verified | [Issue 16 tutorial acceptance](issue-16-tutorial-acceptance.md) |
| DOC-18 | Verified | [Troubleshooting](../user-guide/troubleshooting.md) |
| DOC-19 | Open | [Home](../index.md), [README](https://github.com/CDDLeiden/COFOLDER/blob/main/README.md), [`CITATION.cff`](https://github.com/CDDLeiden/COFOLDER/blob/main/CITATION.cff); archive identifier pending |
| DOC-20 | Verified | Strict MkDocs build completed on 2026-09-21; [Issue 25 documentation acceptance](issue-25-documentation-acceptance.md) records the automated documentation-contract gates |

## Verification

`python -m mkdocs build --strict` completed successfully on 2026-09-21. MkDocs
reported its normal informational list of documentation files outside the configured
navigation; it did not report broken links or build errors.

## DOC-19 release handoff

Before DOC-19 can be checked, publish the approved release/archive record and add its
stable identifier to the citation and public release documentation. Historical release
notes may retain earlier-version entries; they are not current-version claims.
