#!/usr/bin/env bash
set -euo pipefail

# Compatibility wrapper; the installed interface is `cofolder-tools setup-openfold3`.
exec "${PYTHON:-python}" -m cofolder.tools setup-openfold3 -- "$@"
