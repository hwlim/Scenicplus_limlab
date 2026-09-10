#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Runner for the SCENIC+ Snakemake workflow.
#
# Run from inside an analysis directory created by scenicplus_init.sh, the same
# place the bash driver runs from and reading the same config/config.yaml.
#
#   scenicplus.run.sh -n                # dry run: print the plan, do nothing
#   scenicplus.run.sh -j 8              # run locally on 8 cores
#   scenicplus.run.sh --lsf             # submit each rule via profiles/lsf
#   scenicplus.run.sh --lsf -j 20       # ... at most 20 cluster jobs at once
#   scenicplus.run.sh -p                # also print each shell command
#   scenicplus.run.sh -f 7              # force from step 7 (rule R07_*) onward
#   scenicplus.run.sh -- --forceall     # anything after -- goes to snakemake
#
# WHY THE PREFLIGHT IS HERE AND NOT A RULE. scenicplus_check.sh has no outputs,
# so as a rule it would either run on every invocation or need a sentinel that
# lies about when it last passed. More importantly, half of what it checks --
# which libstdc++ actually loaded, what pip's config says -- is a property of
# the environment the work will run in, and it is worth knowing BEFORE anything
# is queued. Under Snakemake that matters more than it did under the bash
# driver: a bad environment would otherwise be discovered once per rule, in
# twenty separate jobs, each after its own queue wait.
#
# STATUS: increment I0. `--lsf` uses profiles/lsf via the cluster-generic
# executor, which is proven on CCHMC but has no step rules to schedule yet. I5
# is where the per-rule resource numbers arrive.
# -----------------------------------------------------------------------------
set -euo pipefail

if [[ -z "${SCENICPLUS_PATH:-}" ]]; then
    echo "[scenicplus.run] ERROR: SCENICPLUS_PATH is not set." >&2
    exit 1
fi

SNAKEFILE="$SCENICPLUS_PATH/Snakefile"
CONFIG="$PWD/config/config.yaml"

if [[ ! -f "$CONFIG" ]]; then
    echo "[scenicplus.run] ERROR: $CONFIG not found." >&2
    echo "  Run scenicplus_init.sh in this directory first." >&2
    exit 1
fi
if [[ ! -f "$SNAKEFILE" ]]; then
    echo "[scenicplus.run] ERROR: $SNAKEFILE not found." >&2
    echo "  SCENICPLUS_PATH should point at the repository root." >&2
    exit 1
fi

DRY=0; JOBS=1; PRINT=0; FROM=""; LSF=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        -n|--dry-run) DRY=1; shift;;
        -j)           JOBS="$2"; shift 2;;
        -p)           PRINT=1; shift;;
        -f|--from)    FROM="$2"; shift 2;;
        -l|--local)   LSF=0; shift;;           # already the default; explicit is fine
        --lsf)        LSF=1; shift;;
        --)           shift; break;;
        -h|--help)    sed -n '2,26p' "$0"; exit 0;;
        *) echo "[scenicplus.run] unknown arg: $1" >&2; exit 2;;
    esac
done

PROFILE="$SCENICPLUS_PATH/profiles/lsf"
if [[ "$LSF" -eq 1 ]]; then
    # snakemake 8 moved cluster submission into executor PLUGINS and ships
    # none, so this can be missing even in a healthy environment. Ask snakemake
    # rather than checking for the package: a plugin can be installed and still
    # fail to register.
    if ! snakemake --executor cluster-generic --help >/dev/null 2>&1; then
        echo "[scenicplus.run] ERROR: the cluster-generic executor is not installed." >&2
        echo "  Add it to this environment:" >&2
        echo "    SCP_FROM=3 ./install_cchmc.sh <prefix>" >&2
        echo "  or directly:" >&2
        echo "    PIP_USER=0 pip install snakemake-executor-plugin-cluster-generic==1.0.8" >&2
        echo "  1.0.9 wants an interface version scenicplus pins away from; see" >&2
        echo "  profiles/lsf/config.yaml." >&2
        exit 2
    fi
    [[ -f "$PROFILE/config.yaml" ]] || {
        echo "[scenicplus.run] ERROR: $PROFILE/config.yaml not found." >&2; exit 2; }
fi

# --- Preflight ---------------------------------------------------------------
if [[ "${SCENICPLUS_SKIP_CHECK:-0}" != "1" ]]; then
    "$SCENICPLUS_PATH/scripts/scenicplus_check.sh"
fi

# NOT passing --configfile on purpose. The Snakefile already declares
# `configfile: "config/config.yaml"`, and a --configfile on the command line
# does not REPLACE that, it EXTENDS it: snakemake merges the two, so a key
# missing from the second one silently keeps the first one's value. Passing the
# same path twice is harmless but teaches the wrong habit, and passing a
# different one is a trap. The workspace's config is the config.
ARGS=(--snakefile "$SNAKEFILE")
if [[ "$LSF" -eq 1 ]]; then
    # --profile carries the executor, the bsub template, the status command and
    # the default resources. `jobs` comes from the profile too, so -j here means
    # concurrent CLUSTER JOBS rather than local cores.
    ARGS+=(--profile "$PROFILE")
    [[ "$JOBS" -ne 1 ]] && ARGS+=(--jobs "$JOBS")
else
    ARGS+=(--cores "$JOBS")
fi
[[ "$DRY"   -eq 1 ]] && ARGS+=(--dry-run)
[[ "$PRINT" -eq 1 ]] && ARGS+=(--printshellcmds)

# --- Step number -> rule name -------------------------------------------------
# Rules are named R01_* .. R20_*, so the step number survives as the interface
# while Snakemake gets real identifiers. It has to be a lookup rather than a
# guess, because the suffix carries the stage's name: `--forcerun R07` is not a
# rule. `--list` is the authority, so a renamed rule surfaces here as "no rule
# matches" rather than as a flag that silently forces nothing.
if [[ -n "$FROM" ]]; then
    if [[ "$FROM" =~ ^[0-9]+$ ]]; then
        want="$(printf 'R%02d_' "$FROM")"
        # `--list` prints "R01_seurat_export (the rule's docstring...)" for any
        # rule that has one, so take the FIRST FIELD. Without that the whole
        # line, docstring and all, was handed to --forcerun as a rule name.
        rule="$(snakemake --snakefile "$SNAKEFILE" --list 2>/dev/null \
                | grep -m1 "^${want}" | awk '{print $1}' || true)"
        if [[ -z "$rule" ]]; then
            echo "[scenicplus.run] ERROR: no rule named ${want}* in this workflow." >&2
            echo "  Available: $(snakemake --snakefile "$SNAKEFILE" \
                                  --list 2>/dev/null | tr '\n' ' ')" >&2
            exit 2
        fi
    else
        rule="$FROM"
    fi
    echo "[scenicplus.run] forcing from ${rule} onward"
    ARGS+=(--forcerun "$rule")
fi

echo "[scenicplus.run] snakemake ${ARGS[*]} $*"
exec snakemake "${ARGS[@]}" "$@"
