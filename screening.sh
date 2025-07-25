#!/usr/bin/env bash

set -euo pipefail

# Defaults
DEFAULT_YAML="examples/vgfr2_options.yaml"
DEFAULT_CONDA_ENV="boltz"

YAML_FILE="${1:-$DEFAULT_YAML}"
CONDA_ENV="${2:-$DEFAULT_CONDA_ENV}"

if [[ ! -f "$YAML_FILE" ]]; then
    echo "YAML file not found: $YAML_FILE"
    echo "Usage: $0 [yaml_file] [conda_env_name]"
    echo "Defaults: yaml_file=$DEFAULT_YAML, conda_env_name=$DEFAULT_CONDA_ENV"
    exit 1
fi

# Extract run_dir from YAML (robustly finds - run_dir: under wrapper, allowing for indentation)
RUN_DIR=$(awk '/wrapper:/,0' "$YAML_FILE" | grep -m1 -E '^[[:space:]]*- run_dir:' | sed 's/^[[:space:]]*- run_dir:[[:space:]]*//;s/["\'']//g')
RUN_DIR=$(echo "$RUN_DIR" | xargs)

# Fallback if not found
if [[ -z "$RUN_DIR" ]]; then
    echo "Could not extract run_dir from $YAML_FILE"
    exit 2
fi

mkdir -p "$RUN_DIR"
LOGFILE="${RUN_DIR}/screening_run_$(date +%Y%m%d_%H%M%S).log"

{
    echo "===== Screening Run Log ====="
    echo "Start time: $(date)"
    echo "YAML file: $YAML_FILE"
    echo "Run dir: $RUN_DIR"
    echo "Conda env: $CONDA_ENV"
    echo "Hostname: $(hostname)"
    echo "nvidia-smi:"
    nvidia-smi || echo "nvidia-smi not found"
    echo "which conda: $(which conda || echo 'not found')"
    echo "which python (before env): $(which python || echo 'not found')"
    echo "-----------------------------"
    echo "Sourcing ~/.bashrc and activating conda env..."
} | tee -a "$LOGFILE"

START_TIME=$(date +%s)

{
    source ~/.bashrc
    conda activate "$CONDA_ENV"
    echo "which python (after env): $(which python || echo 'not found')"
    echo "-----------------------------"
    echo "Running screening.py..."
    python screening.py "$YAML_FILE"
} 2>&1 | tee -a "$LOGFILE"

EXIT_CODE=${PIPESTATUS[0]}
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))

{
    echo "-----------------------------"
    echo "End time: $(date)"
    echo "Total runtime (seconds): $RUNTIME"
    if [[ $EXIT_CODE -ne 0 ]]; then
        echo "ERROR: screening.py exited with code $EXIT_CODE"
    else
        echo "screening.py completed successfully."
    fi
    echo "===== End of Log ====="
} | tee -a "$LOGFILE"

exit $EXIT_CODE 