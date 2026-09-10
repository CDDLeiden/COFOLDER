#!/usr/bin/env bash

set -euo pipefail

cd "$(dirname "$0")" || exit 1

cat >&2 <<'EOF'
WARNING: run_ui.sh starts an unsupported, obsolete development prototype.
It is not part of the COFOLDER release, and its configuration and command
construction may not match the current recipe contracts. Use the bias, validate,
screen, or oracle CLI commands or the Python API for supported workflows.
EOF

if ! python -c "import streamlit, streamlit_option_menu" >/dev/null 2>&1; then
    cat >&2 <<'EOF'
COFOLDER UI dependencies are not installed in this environment.

There is no supported COFOLDER UI extra. For source-tree experimentation only,
install the development dependencies directly:
  python -m pip install "streamlit>=1.28.0" "streamlit-option-menu>=0.3.6"

Then rerun:
  ./run_ui.sh
EOF
    exit 1
fi

echo "Starting unsupported COFOLDER development UI..."
echo "Access the UI at: http://localhost:8501"
echo ""

exec python -m streamlit run src/cofolder/ui/app.py
