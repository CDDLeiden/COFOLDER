#!/bin/bash
# Start the Boltz-Lab Streamlit UI

cd "$(dirname "$0")" || exit

# Check if streamlit is installed
if ! command -v streamlit &> /dev/null; then
    echo "Streamlit is not installed. Installing dependencies..."
    pip install streamlit streamlit-option-menu pyyaml
fi

# Run the Streamlit app
echo "Starting Boltz-Lab UI..."
echo "Access the UI at: http://localhost:8501"
echo ""

streamlit run src/boltz_eval/ui/app.py
