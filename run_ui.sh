#!/usr/bin/env bash

set -euo pipefail

cd "$(dirname "$0")" || exit 1

if ! python -c "import streamlit, streamlit_option_menu" >/dev/null 2>&1; then
    cat >&2 <<'EOF'
COFOLDER UI dependencies are not installed in this environment.

Install them explicitly with:
  python -m pip install -e ".[ui]"

Then rerun:
  ./run_ui.sh
EOF
    exit 1
fi

echo "Starting COFOLDER UI..."
echo "Access the UI at: http://localhost:8501"
echo ""

exec python -m streamlit run src/cofolder/ui/app.py
