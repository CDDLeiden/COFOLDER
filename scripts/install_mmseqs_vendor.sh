#!/usr/bin/env bash
set -euo pipefail

# Compatibility wrapper; the installed interface is `cofolder-tools install-mmseqs`.
exec "${PYTHON:-python}" -m cofolder.tools install-mmseqs "$@"
