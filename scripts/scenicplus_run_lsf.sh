#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Submit the SCENIC+ pipeline to an LSF cluster as a single bsub'd driver job.
#
# The driver runs scenicplus_run_pipeline.sh, which walks the 9 steps sequentially. The
# heavy step (07 run_scenicplus) internally invokes its own snakemake with
# {LSF_CORES} cores, so the bsub'd job needs that many slots.
#
# Run from inside an analysis directory initialized with scenicplus_init.sh.
#
# Customize via env vars:
#   LSF_QUEUE  LSF_PROJECT  LSF_CORES  LSF_MEM_MB  LSF_WALLTIME  LOG_DIR
#   SCENICPLUS_ENV         (conda env to activate inside the bsub'd shell)
#   SCENICPLUS_SKIP_CHECK  (set to 1 to skip the package preflight)
#
# Extra args are forwarded to scenicplus_run_pipeline.sh (e.g. --from 4, --force).
# -----------------------------------------------------------------------------
set -euo pipefail

if [[ -z "${SCENICPLUS_PATH:-}" ]]; then
    echo "[scenicplus_run_lsf] ERROR: SCENICPLUS_PATH is not set." >&2
    exit 1
fi

CONFIG="$PWD/config/config.yaml"
DRIVER="$SCENICPLUS_PATH/scripts/scenicplus_run_pipeline.sh"

if [[ ! -f "$CONFIG" ]]; then
    echo "[scenicplus_run_lsf] ERROR: $CONFIG not found." >&2
    echo "  Run scenicplus_init.sh in this directory first." >&2
    exit 1
fi

LSF_QUEUE="${LSF_QUEUE:-normal}"
LSF_PROJECT="${LSF_PROJECT:-scenicplus}"
LSF_CORES="${LSF_CORES:-16}"
LSF_MEM_MB="${LSF_MEM_MB:-64000}"
LSF_WALLTIME="${LSF_WALLTIME:-72:00}"
LOG_DIR="${LOG_DIR:-$PWD/logs/lsf}"
mkdir -p "$LOG_DIR"

if [[ -n "${SCENICPLUS_ENV:-}" ]]; then
    CONDA_PRELUDE="source \$(conda info --base)/etc/profile.d/conda.sh && conda activate ${SCENICPLUS_ENV} && "
else
    CONDA_PRELUDE=""
fi

# Preflight on the submitting host (skip if cluster compute env differs).
if [[ "${SCENICPLUS_SKIP_CHECK:-0}" != "1" ]]; then
    "$SCENICPLUS_PATH/scripts/scenicplus_check.sh"
fi

# Forward CLI flags to the driver inside the bsub.
DRIVER_ARGS=""
for arg in "$@"; do
    DRIVER_ARGS+=" $(printf %q "$arg")"
done

echo "[scenicplus_run_lsf] queue=${LSF_QUEUE} cores=${LSF_CORES} mem=${LSF_MEM_MB}MB walltime=${LSF_WALLTIME}"

bsub \
    -q "${LSF_QUEUE}" \
    -P "${LSF_PROJECT}" \
    -J "scenicplus_driver" \
    -n "${LSF_CORES}" \
    -W "${LSF_WALLTIME}" \
    -R "rusage[mem=${LSF_MEM_MB}] span[hosts=1]" \
    -o "${LOG_DIR}/driver_%J.out" \
    -e "${LOG_DIR}/driver_%J.err" \
    /bin/bash -c "${CONDA_PRELUDE}${DRIVER}${DRIVER_ARGS}"
