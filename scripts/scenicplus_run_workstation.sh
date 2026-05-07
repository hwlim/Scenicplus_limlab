#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Run the SCENIC+ pipeline on a regular workstation.
#
# Run this from inside an analysis directory that was initialized with
# `scenicplus_init.sh`. Thin wrapper around scenicplus_run_pipeline.sh.
#
# Usage:
#   scenicplus_run_workstation.sh                # run everything stale
#   scenicplus_run_workstation.sh --dry-run
#   scenicplus_run_workstation.sh --from 4
#   scenicplus_run_workstation.sh --only 4
#   scenicplus_run_workstation.sh --force
# -----------------------------------------------------------------------------
set -euo pipefail

if [[ -z "${SCENICPLUS_PATH:-}" ]]; then
    echo "[scenicplus_run_workstation] ERROR: SCENICPLUS_PATH is not set." >&2
    exit 1
fi

exec "$SCENICPLUS_PATH/scripts/scenicplus_run_pipeline.sh" "$@"
