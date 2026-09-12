#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Initialize a SCENIC+ analysis directory.
#
# Copies the config.yaml template from $SCENICPLUS_PATH/config/ into the
# current (or specified) directory, so the user has a per-analysis config
# they can edit without touching the central pipeline.
#
# Usage:
#   scenicplus_init.sh                # initialize current directory
#   scenicplus_init.sh /path/to/dir   # initialize a specific directory
#   scenicplus_init.sh -f             # overwrite an existing config.yaml
#
# Then edit config/config.yaml and run the workflow with scenicplus.run.sh --
# see quickstart.md. The bash launchers this script used to name are obsolete.
# -----------------------------------------------------------------------------
set -euo pipefail

FORCE=0
TARGET_DIR=""
for arg in "$@"; do
    case "$arg" in
        -f|--force) FORCE=1 ;;
        -h|--help)
            # Derived from the comment block, not a hardcoded line range. The
            # runner shipped `sed -n '2,26p'` and an edit to its header
            # truncated the help mid-sentence with nothing to notice it; this
            # file had the same `sed -n '2,15p'` and had already grown past it.
            awk 'NR>1 && /^#/ {print} NR>1 && !/^#/ {exit}' "$0"
            exit 0
            ;;
        *) TARGET_DIR="$arg" ;;
    esac
done

if [[ -z "${SCENICPLUS_PATH:-}" ]]; then
    echo "[scenicplus_init] ERROR: SCENICPLUS_PATH is not set." >&2
    echo "  Point it at the central pipeline directory containing config/, workflow/, scripts/." >&2
    exit 1
fi
if [[ ! -d "$SCENICPLUS_PATH" ]]; then
    echo "[scenicplus_init] ERROR: SCENICPLUS_PATH=$SCENICPLUS_PATH does not exist." >&2
    exit 1
fi

TEMPLATE="$SCENICPLUS_PATH/config/config.yaml"
if [[ ! -f "$TEMPLATE" ]]; then
    echo "[scenicplus_init] ERROR: template not found at $TEMPLATE" >&2
    exit 1
fi

TARGET_DIR="${TARGET_DIR:-$PWD}"
mkdir -p "$TARGET_DIR/config"
DEST="$TARGET_DIR/config/config.yaml"

if [[ -f "$DEST" && "$FORCE" -ne 1 ]]; then
    echo "[scenicplus_init] $DEST already exists. Pass -f to overwrite." >&2
    exit 1
fi

cp "$TEMPLATE" "$DEST"
echo "[scenicplus_init] Wrote $DEST"
echo "[scenicplus_init] Edit this file, then run:"
echo "    scenicplus.run.sh -n            # dry run: show the plan first"
echo "    scenicplus.run.sh -j 8          # run here, on 8 cores"
echo "    scenicplus.run.sh --lsf -j 20   # one LSF job per rule, 20 at once"
echo "[scenicplus_init] A run ends with report.html; open that first."
echo "[scenicplus_init] Set input.reduction before the first run -- changing it"
echo "                  later re-runs the whole pipeline (RUNBOOK section 5)."
# This used to print scenicplus_run_workstation.sh / scenicplus_run_lsf.sh, so
# the FIRST thing a new user was told to do was start the obsolete driver --
# one job for twenty steps, sized for the heaviest, and no report or provenance
# bundle at the end. Those launchers still work and are kept as a fallback;
# they are not what anyone should be pointed at.
