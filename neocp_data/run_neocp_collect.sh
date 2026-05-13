#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-/opt/anaconda3/bin/python}"

cd "$SCRIPT_DIR"
mkdir -p logs

"$PYTHON_BIN" "$SCRIPT_DIR/collect_neocp.py" "$@" >> "$SCRIPT_DIR/logs/neocp_collect.log" 2>&1
