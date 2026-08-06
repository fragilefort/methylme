#!/usr/bin/env bash
set -euo pipefail

# Drop-in threaded runner for the existing unified analysis script.
# It creates a temporary edited copy that adds THREADS and -p "$THREADS"
# to every computeMatrix command, then runs that copy.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_SCRIPT="$SCRIPT_DIR/run_unified_deeptools_option1.sh"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
export THREADS

[[ -s "$SOURCE_SCRIPT" ]] || {
  echo "Cannot find source script: $SOURCE_SCRIPT" >&2
  exit 1
}

TEMP_SCRIPT="$(mktemp --tmpdir run_unified_deeptools_option1_threaded_XXXXXX.sh)"
trap 'rm -f "$TEMP_SCRIPT"' EXIT

awk '
  /^OUTROOT=/ {
    print
    print "THREADS=\"${SLURM_CPUS_PER_TASK:-4}\""
    next
  }
  /^computeMatrix / {
    sub(/^computeMatrix /, "computeMatrix -p \\"$THREADS\\" ")
    print
    next
  }
  { print }
' "$SOURCE_SCRIPT" > "$TEMP_SCRIPT"

chmod +x "$TEMP_SCRIPT"
echo "Running unified deepTools workflow with THREADS=$THREADS"
bash "$TEMP_SCRIPT"
