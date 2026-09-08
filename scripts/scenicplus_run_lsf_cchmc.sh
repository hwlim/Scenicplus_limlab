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

# -M vs rusage[mem]: `rusage` RESERVES memory for scheduling; it is not a
# ceiling. Without -M the ENFORCED ceiling is the queue's MEMLIMIT default,
# which nobody here chose -- so a job can reserve 128 GB, use 8, and still be
# killed with TERM_MEMLIMIT. Both flags carry the same number on purpose.
#
# Both are read in LSF_UNIT_FOR_LIMITS, whose documented default when unset is
# KB, not MB. LSF_MEM_MB is megabytes by its name, so print the site's unit
# rather than let a 1024x error look like a memory bug in the pipeline.
lsf_unit_note() {
    local u
    u="$(grep -hs '^[[:space:]]*LSF_UNIT_FOR_LIMITS' \
         "${LSF_ENVDIR:-/etc/lsf}"/lsf.conf 2>/dev/null | tail -1 | cut -d= -f2 | tr -d '[:space:]')"
    if [[ -z "$u" ]]; then
        echo "[lsf] WARNING: LSF_UNIT_FOR_LIMITS not found in ${LSF_ENVDIR:-/etc/lsf}/lsf.conf." >&2
        echo "[lsf]          LSF defaults it to KB, in which case -M ${LSF_MEM_MB} means" >&2
        echo "[lsf]          $((LSF_MEM_MB / 1024)) MB, not ${LSF_MEM_MB} MB. Check with your site." >&2
    elif [[ "${u^^}" != "MB" ]]; then
        echo "[lsf] WARNING: LSF_UNIT_FOR_LIMITS=$u, but LSF_MEM_MB is in MB." >&2
        echo "[lsf]          -M ${LSF_MEM_MB} will be read as ${LSF_MEM_MB} $u." >&2
    else
        echo "[lsf] LSF_UNIT_FOR_LIMITS=$u -> -M ${LSF_MEM_MB} = ${LSF_MEM_MB} MB"
    fi
}

echo "[scenicplus_run_lsf_cchmc] queue=${LSF_QUEUE} cores=${LSF_CORES} mem=${LSF_MEM_MB}MB walltime=${LSF_WALLTIME}"
echo "[scenicplus_run_lsf_cchmc] env=${SCENICPLUS_ENV}  R_module=${R_MODULE}"
lsf_unit_note

# `bash -lc` so the login shell defines the `module` function.
bsub \
    -q "${LSF_QUEUE}" \
    -P "${LSF_PROJECT}" \
    -J "scenicplus_driver" \
    -n "${LSF_CORES}" \
    -W "${LSF_WALLTIME}" \
    -M "${LSF_MEM_MB}" \
    -R "rusage[mem=${LSF_MEM_MB}] span[hosts=1]" \
    -o "${LOG_DIR}/driver_%J.out" \
    -e "${LOG_DIR}/driver_%J.err" \
    /bin/bash -lc "${PRELUDE} && ${DRIVER}${DRIVER_ARGS}"
