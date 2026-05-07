#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# CCHMC LSF launcher for the SCENIC+ pipeline.
#
# Same as scenicplus_run_lsf.sh, but the bsub'd shell sets up the environment
# the CCHMC way: the site's `scenicplus` conda env does NOT bundle R, so R is
# brought in via the module system *after* the conda env is activated.
#
# Inside the bsub'd shell:
#   module purge
#   module load anaconda3
#   conda activate ${SCENICPLUS_ENV:-scenicplus}
#   module load ${R_MODULE:-R/4.4.0-R0}
#
# Run from inside an analysis directory initialized with scenicplus_init.sh.
# Extra args are forwarded to scenicplus_run_pipeline.sh (e.g. --from 4).
#
# Customize via env vars:
#   LSF_QUEUE  LSF_PROJECT  LSF_CORES  LSF_MEM_MB  LSF_WALLTIME  LOG_DIR
#   SCENICPLUS_ENV         (conda env name; default: scenicplus)
#   R_MODULE               (Lmod module spec; default: R/4.4.0-R0)
# -----------------------------------------------------------------------------
set -euo pipefail

if [[ -z "${SCENICPLUS_PATH:-}" ]]; then
    echo "[scenicplus_run_lsf_cchmc] ERROR: SCENICPLUS_PATH is not set." >&2
    exit 1
fi

CONFIG="$PWD/config/config.yaml"
DRIVER="$SCENICPLUS_PATH/scripts/scenicplus_run_pipeline.sh"

if [[ ! -f "$CONFIG" ]]; then
    echo "[scenicplus_run_lsf_cchmc] ERROR: $CONFIG not found." >&2
    echo "  Run scenicplus_init.sh in this directory first." >&2
    exit 1
fi

LSF_QUEUE="${LSF_QUEUE:-normal}"
LSF_PROJECT="${LSF_PROJECT:-scenicplus}"
LSF_CORES="${LSF_CORES:-16}"
LSF_MEM_MB="${LSF_MEM_MB:-64000}"
LSF_WALLTIME="${LSF_WALLTIME:-72:00}"
LOG_DIR="${LOG_DIR:-$PWD/logs/lsf}"
SCENICPLUS_ENV="${SCENICPLUS_ENV:-scenicplus}"
R_MODULE="${R_MODULE:-R/4.4.0-R0}"
mkdir -p "$LOG_DIR"

# The submitting host probably doesn't have R loaded, so the standard
# preflight (which does `Rscript -e ...`) would falsely fail. Skip by default;
# the driver still runs scenicplus_check.sh inside the bsub'd shell after the
# modules are loaded (unless SCENICPLUS_SKIP_CHECK=1).
PRELUDE="module purge && module load anaconda3 && conda activate ${SCENICPLUS_ENV} && module load ${R_MODULE}"

# Forward CLI flags to the driver inside the bsub.
DRIVER_ARGS=""
for arg in "$@"; do
    DRIVER_ARGS+=" $(printf %q "$arg")"
done

echo "[scenicplus_run_lsf_cchmc] queue=${LSF_QUEUE} cores=${LSF_CORES} mem=${LSF_MEM_MB}MB walltime=${LSF_WALLTIME}"
echo "[scenicplus_run_lsf_cchmc] env=${SCENICPLUS_ENV}  R_module=${R_MODULE}"

# `bash -lc` so the login shell defines the `module` function.
bsub \
    -q "${LSF_QUEUE}" \
    -P "${LSF_PROJECT}" \
    -J "scenicplus_driver" \
    -n "${LSF_CORES}" \
    -W "${LSF_WALLTIME}" \
    -R "rusage[mem=${LSF_MEM_MB}] span[hosts=1]" \
    -o "${LOG_DIR}/driver_%J.out" \
    -e "${LOG_DIR}/driver_%J.err" \
    /bin/bash -lc "${PRELUDE} && ${DRIVER}${DRIVER_ARGS}"
