#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  fetch_bias_training_data_tmp_mmseqs_env.sh <output_root> [fetch_bias_training_data.py args...]

Examples:
  fetch_bias_training_data_tmp_mmseqs_env.sh /path/to/training_data --overwrite
  fetch_bias_training_data_tmp_mmseqs_env.sh /path/to/training_data --overwrite --mmseqs_db_name PDB

What it does:
  1) Creates a temporary conda env
  2) Installs mmseqs2
  3) Runs scripts/fetch_bias_training_data.py
  4) Returns to the previous conda env
  5) Deletes the temporary conda env
USAGE
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" || $# -lt 1 ]]; then
  usage
  exit 0
fi

OUTPUT_ROOT="$1"
shift || true
FETCH_ARGS=("$@")

if ! command -v conda >/dev/null 2>&1; then
  echo "[error] conda not found in PATH" >&2
  exit 1
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FETCH_SCRIPT="${REPO_ROOT}/scripts/fetch_bias_training_data.py"

if [[ ! -f "${FETCH_SCRIPT}" ]]; then
  echo "[error] fetch script not found: ${FETCH_SCRIPT}" >&2
  exit 1
fi

ORIG_ENV="${CONDA_DEFAULT_ENV:-}"
CONDA_BASE="$(conda info --base)"
# shellcheck source=/dev/null
source "${CONDA_BASE}/etc/profile.d/conda.sh"

TMP_ENV="cofolder-mmseqs-tmp-$(date +%s)-${RANDOM}"

cleanup() {
  local exit_code=$?
  set +e

  if [[ "${CONDA_DEFAULT_ENV:-}" == "${TMP_ENV}" ]]; then
    conda deactivate >/dev/null 2>&1 || true
  fi

  conda env remove -y -n "${TMP_ENV}" >/dev/null 2>&1 || true

  if [[ -n "${ORIG_ENV}" && "${CONDA_DEFAULT_ENV:-}" != "${ORIG_ENV}" ]]; then
    conda activate "${ORIG_ENV}" >/dev/null 2>&1 || true
  fi

  if [[ ${exit_code} -eq 0 ]]; then
    echo "[done] cleaned temporary conda env: ${TMP_ENV}"
    if [[ -n "${ORIG_ENV}" ]]; then
      echo "[done] restored conda env: ${ORIG_ENV}"
    fi
  fi

  exit ${exit_code}
}
trap cleanup EXIT INT TERM

echo "[info] creating temporary conda env: ${TMP_ENV}"
conda create -y -n "${TMP_ENV}" -c conda-forge python=3.12

echo "[info] installing mmseqs2 in ${TMP_ENV}"
conda install -y -n "${TMP_ENV}" -c conda-forge -c bioconda mmseqs2

echo "[info] activating temporary env: ${TMP_ENV}"
conda activate "${TMP_ENV}"

echo "[info] running fetch script"
python "${FETCH_SCRIPT}" --output_root "${OUTPUT_ROOT}" "${FETCH_ARGS[@]}"

echo "[done] fetch complete"
