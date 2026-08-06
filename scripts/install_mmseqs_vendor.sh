#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Install a vendored MMseqs2 binary for COFOLDER (no conda-forge required).

Usage:
  install_mmseqs_vendor.sh [--url <tarball_url> --sha256 <sha256>] [--install-dir <dir>]

Example:
  install_mmseqs_vendor.sh

  Custom source:
  install_mmseqs_vendor.sh \
    --url https://example.org/mmseqs-linux-avx2.tar.gz \
    --sha256 0123456789abcdef... \
    --install-dir ~/.cofolder/vendor/mmseqs

After install, either:
  1) export COFOLDER_MMSEQS_BIN=<install-dir>/bin/mmseqs
  2) or add <install-dir>/bin to PATH
USAGE
}

URL=""
SHA256_SUM=""
INSTALL_DIR="${HOME}/.cofolder/vendor/mmseqs"
DEFAULT_URL="https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz"
DEFAULT_SHA256="ea05b7c706a166874f56ff0cbdf14791b0de1b7ae3f44e59d2a78383ff8f8d62"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --url)
      URL="${2:-}"
      shift 2
      ;;
    --sha256)
      SHA256_SUM="${2:-}"
      shift 2
      ;;
    --install-dir)
      INSTALL_DIR="${2:-}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "[error] unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

if [[ -z "${URL}" && -z "${SHA256_SUM}" ]]; then
  URL="${DEFAULT_URL}"
  SHA256_SUM="${DEFAULT_SHA256}"
elif [[ -z "${URL}" || -z "${SHA256_SUM}" ]]; then
  echo "[error] --url and --sha256 must be provided together" >&2
  exit 1
fi

echo "[info] source URL: ${URL}"
echo "[info] expected SHA256: ${SHA256_SUM}"

if command -v curl >/dev/null 2>&1; then
  DL_CMD=(curl -L --fail --silent --show-error)
elif command -v wget >/dev/null 2>&1; then
  DL_CMD=(wget -q -O -)
else
  echo "[error] neither curl nor wget found" >&2
  exit 1
fi

TMP_DIR="$(mktemp -d)"
cleanup() {
  rm -rf "${TMP_DIR}"
}
trap cleanup EXIT INT TERM

ARCHIVE="${TMP_DIR}/mmseqs.tar.gz"
EXTRACT_DIR="${TMP_DIR}/extract"
mkdir -p "${EXTRACT_DIR}"

echo "[info] downloading MMseqs2 archive"
if [[ "${DL_CMD[0]}" == "curl" ]]; then
  "${DL_CMD[@]}" "${URL}" -o "${ARCHIVE}"
else
  "${DL_CMD[@]}" "${URL}" > "${ARCHIVE}"
fi

echo "[info] verifying SHA256"
ACTUAL_SHA="$(sha256sum "${ARCHIVE}" | awk '{print $1}')"
if [[ "${ACTUAL_SHA}" != "${SHA256_SUM}" ]]; then
  echo "[error] sha256 mismatch" >&2
  echo "expected: ${SHA256_SUM}" >&2
  echo "actual:   ${ACTUAL_SHA}" >&2
  exit 1
fi

echo "[info] extracting archive"
tar -xzf "${ARCHIVE}" -C "${EXTRACT_DIR}"

MMSEQS_SRC="$(find "${EXTRACT_DIR}" -type f -name mmseqs | head -n 1 || true)"
if [[ -z "${MMSEQS_SRC}" ]]; then
  echo "[error] mmseqs executable not found in archive" >&2
  exit 1
fi

mkdir -p "${INSTALL_DIR}/bin"
install -m 0755 "${MMSEQS_SRC}" "${INSTALL_DIR}/bin/mmseqs"

echo "[done] installed mmseqs to ${INSTALL_DIR}/bin/mmseqs"
echo "[next] export COFOLDER_MMSEQS_BIN='${INSTALL_DIR}/bin/mmseqs'"
echo "[next] or add '${INSTALL_DIR}/bin' to PATH"
