# Changelog

All notable changes to COFOLDER will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Python Oracle composite objectives, structured custom scoring callbacks, and
  auditable structure gates with down-weight, fixed-penalty, or non-binder outcomes.
- Qualified Oracle metric selectors plus documented ligand-bias and custom-pocket
  coverage examples.
- Screen-level reuse of fixed-protein MSAs for Boltz-family runners, plus a stable
  manuscript-facing consolidated output schema and a DataFrame return value from
  `Screen.run()`.
- Opt-in deterministic average-linkage clustering of binary distance IFPs, with
  stable row annotations and cluster medoid/consensus summaries.
- Runner-specific DNA, RNA, and constraint input contracts with pre-execution
  validation, OpenFold3 nucleic-acid/pocket translation, and manual backend
  acceptance fixtures.
- Comprehensive documentation with MkDocs Material
- API reference documentation
- User guides and tutorials
- Contributing guidelines

### Changed
- Package name from `boltz-eval` to `cofolder`
- Repository structure overhaul
- Improved logging system

### Fixed
- Import paths updated for new package structure
- Screening and validation now warn and continue with empty affinity columns when
  default affinity outputs are unavailable but affinity was not explicitly requested
  or activated in the system YAML.
- Protein sequence similarities in generated combined bias tables now consistently use
  MMseqs `pident`, retain below-threshold hits for PDB lookup, and record method
  provenance instead of substituting PairwiseAligner scores. Downstream bias outputs
  also keep PairwiseAligner scores in a separate `sequence_similarity_pairwise` column.

## [0.1.0] - 2026-08-06

### Added
- Initial development version
- `validate` command for single system co-folding and validation
- `screen` command for virtual screening
- `oracle` command for oracle function usage
- Support for SMILES, SDF, PDB, and CIF input formats
- 2D and 3D conformer generation
- RMSD calculation against reference structures
- Interaction fingerprint analysis
- CSV and SDF input/output handling
- Comprehensive CLI interface
- Example configurations and tutorials

[Unreleased]: https://github.com/CDDLeiden/COFOLDER/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/CDDLeiden/COFOLDER/releases/tag/v0.1.0
