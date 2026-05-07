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
# -----------------------------------------------------------------------------
set -euo pipefail

FORCE=0
TARGET_DIR=""
for arg in "$@"; do
    case "$arg" in
        -f|--force) FORCE=1 ;;
        -h|--help)
            sed -n '2,15p' "$0"
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
echo "    scenicplus_run_workstation.sh   # local workstation"
echo "    scenicplus_run_lsf.sh           # LSF cluster"
