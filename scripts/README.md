# Scripts

This folder contains helper utilities that support optional setup and data-preparation tasks around COFOLDER.

These scripts are not the primary public interface. The primary public interface is the `cofolder` CLI.

## Installed User-Facing Helpers

After installing COFOLDER, run `cofolder-tools --help`. It provides:

- `copy-examples`: copies the packaged runnable examples to a workspace
- `setup-openfold3`: prepares the OpenFold3 cache and checkpoints
- `install-mmseqs`: installs a verified vendored `mmseqs2` binary
- `fetch-bias-training-data`: downloads CCD data and optional MMseqs2 databases
- `build-bias-training-data`: builds protein and ligand training-reference CSVs

The similarly named scripts in this directory are retained source-tree wrappers.
They delegate to the packaged implementations and remain available for existing
checkout-based workflows.

## Maintainer-Oriented Helper

- `fetch_bias_training_data_tmp_mmseqs_env.sh`: convenience wrapper that creates a temporary conda environment, installs `mmseqs2`, runs `fetch_bias_training_data.py`, and removes the temporary environment afterward

This wrapper can be useful for maintainers, but it is intentionally more operational and less reproducible than installing the required tooling in your target environment up front.

## Recommended User Path

Most users should:

1. install COFOLDER
2. use the `cofolder` CLI directly
3. use `cofolder-tools` when they need examples, optional backend setup, or bias-data preparation
