# Scripts

This folder contains helper utilities that support optional setup and data-preparation tasks around COFOLDER.

These scripts are not the primary public interface. The primary public interface is the `cofolder` CLI.

## Supported User-Facing Helpers

- `setup_openfold3.sh`: prepares the OpenFold3 cache and checkpoints after `python -m pip install -e ".[openfold3]"`
- `install_mmseqs_vendor.sh`: installs a vendored `mmseqs2` binary when you do not want to rely on conda-forge
- `fetch_bias_training_data.py`: downloads CCD data and, optionally, MMseqs2 databases used for bias-data preparation
- `build_bias_training_data.py`: builds public protein and ligand training-reference CSVs for the bias workflow

## Maintainer-Oriented Helper

- `fetch_bias_training_data_tmp_mmseqs_env.sh`: convenience wrapper that creates a temporary conda environment, installs `mmseqs2`, runs `fetch_bias_training_data.py`, and removes the temporary environment afterward

This wrapper can be useful for maintainers, but it is intentionally more operational and less reproducible than installing the required tooling in your target environment up front.

## Recommended User Path

Most users should:

1. install COFOLDER
2. use the `cofolder` CLI directly
3. return to this folder only when they need optional OpenFold3 setup or bias-data preparation
