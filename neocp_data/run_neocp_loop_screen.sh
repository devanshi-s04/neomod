#!/usr/bin/env bash
# Run the NEOCP collector in a screen session on Arnor.
#
# Start:   screen -S neocp
#          cd ~/digest3/neocp_data && ./run_neocp_loop_screen.sh
#          Ctrl+A, D  to detach
#
# Reattach: screen -r neocp
# Stop:     screen -r neocp, then Ctrl+C

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INTERVAL_SEC=${NEOCP_INTERVAL_SEC:-14400}   # default: every 4 hours

mkdir -p "$SCRIPT_DIR/logs"
LOG="$SCRIPT_DIR/logs/neocp_collect.log"

echo "=== NEOCP loop started $(date -u '+%Y-%m-%dT%H:%M:%SZ') interval=${INTERVAL_SEC}s ===" | tee -a "$LOG"

while true; do
    echo "--- $(date -u '+%Y-%m-%dT%H:%M:%SZ') collecting ---" | tee -a "$LOG"
    python3 "$SCRIPT_DIR/collect_neocp.py" 2>&1 | tee -a "$LOG" || true
    echo "--- $(date -u '+%Y-%m-%dT%H:%M:%SZ') sleeping ${INTERVAL_SEC}s ---" | tee -a "$LOG"
    sleep "$INTERVAL_SEC"
done
