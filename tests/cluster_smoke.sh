#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Prove the cluster-generic executor and profiles/lsf work, before any real step
# depends on them.
#
#   tests/cluster_smoke.sh              # local: no LSF, no cluster needed
#   tests/cluster_smoke.sh --lsf        # the real thing, ON THE CLUSTER
#
# Needs SCENICPLUS_PATH set and snakemake on PATH. Writes into a throwaway
# directory and removes it, except on --lsf, where the workspace is kept so the
# bsub logs can be read afterwards.
#
# THE LOCAL MODE IS NOT A SUBSTITUTE. It swaps bsub for a shell that runs the
# job script directly, so it exercises the executor, the profile's parsing, the
# resource plumbing and the failure path -- but not LSF. Queue names, the
# memory unit, slot allocation and job-state polling are only tested by --lsf.
# Run local first anyway: it catches the mistakes that would otherwise waste a
# queue wait to discover.
# -----------------------------------------------------------------------------
set -uo pipefail

: "${SCENICPLUS_PATH:?set SCENICPLUS_PATH to the repository root}"
command -v snakemake >/dev/null || { echo "snakemake not on PATH"; exit 1; }

MODE=local
[[ "${1:-}" == "--lsf" ]] && MODE=lsf

SMK="$SCENICPLUS_PATH/tests/cluster_smoke/Snakefile"
FAIL=0
say() { printf '  %-4s %s\n' "$1" "$2"; [[ "$1" == FAIL ]] && FAIL=1; return 0; }

if ! snakemake --executor cluster-generic --help >/dev/null 2>&1; then
    echo "the cluster-generic executor is not installed in this environment."
    echo "  pip install snakemake-executor-plugin-cluster-generic==1.0.8"
    echo "  (1.0.9 requires an interface version scenicplus pins away from; see"
    echo "   profiles/lsf/config.yaml)"
    exit 1
fi

if [[ "$MODE" == lsf ]]; then
    WORK="$PWD/cluster_smoke.$$"
    command -v bsub >/dev/null || { echo "--lsf but no bsub on PATH"; exit 1; }
    EXEC=(--profile "$SCENICPLUS_PATH/profiles/lsf")
    # Long enough to cover a queue wait, short enough that a hang is a result
    # rather than an afternoon.
    LIMIT=1800
else
    WORK="$(mktemp -d)"
    export XDG_CACHE_HOME="$WORK/cache"
    # `bash` as the submit command: the plugin appends the job script path, so
    # this runs the job right here and prints nothing. With no status command,
    # snakemake falls back to watching for the outputs, which is precisely the
    # arrangement that hangs on a dead job -- and that is what makes the
    # deliberate failure below worth timing out.
    EXEC=(--executor cluster-generic --cluster-generic-submit-cmd bash --jobs 4)
    LIMIT=180
fi

mkdir -p "$WORK"
cd "$WORK" || exit 1
echo "mode:      $MODE"
echo "workspace: $WORK"
echo

# --- 0. the profile itself ---------------------------------------------------
# A dry run READS the profile and turns every entry into a command-line flag,
# so this catches a typo'd key or a malformed submit template without a cluster
# and without submitting anything. Worth doing first: the alternative is
# discovering it after a queue wait.
timeout 120 snakemake --snakefile "$SMK" --profile "$SCENICPLUS_PATH/profiles/lsf" -n \
    >"$WORK/profile.log" 2>&1 \
    && say ok "profiles/lsf parses and plans" \
    || { say FAIL "profiles/lsf does not parse"; tail -n 5 "$WORK/profile.log"; }

# The same check pointed at a deliberately broken copy, so a pass above means
# something. A profile key snakemake does not recognise must be refused, not
# ignored -- an ignored `latency-wait` on a shared filesystem is a run that
# fails on missing outputs that were merely late.
mkdir -p "$WORK/badprofile"
sed 's/^latency-wait:/latency_wait:/' "$SCENICPLUS_PATH/profiles/lsf/config.yaml" \
    > "$WORK/badprofile/config.yaml"
if timeout 120 snakemake --snakefile "$SMK" --profile "$WORK/badprofile" -n >/dev/null 2>&1; then
    say FAIL "a typo'd profile key was ACCEPTED -- the check above proves nothing"
else
    say ok "a typo'd profile key is refused"
fi
echo

# --- 1. the happy path -------------------------------------------------------
timeout "$LIMIT" snakemake --snakefile "$SMK" "${EXEC[@]}" >"$WORK/run1.log" 2>&1
rc=$?
case "$rc" in
    0)   say ok   "three jobs submitted, ran, and joined" ;;
    124) say FAIL "TIMED OUT after ${LIMIT}s -- submission or polling is stuck"; sed -n '1,25p' "$WORK/run1.log" ;;
    *)   say FAIL "exit $rc"; sed -n '1,40p' "$WORK/run1.log" ;;
esac

if [[ -f smoke/summary.txt ]]; then
    say ok "summary.txt exists"
    # The point of b_bigger: LSF must give it the four slots the rule asked
    # for. `unset` here means the executor never told LSF about `threads`,
    # which is the 16x oversubscription bug from the sibling repo, in advance.
    got="$(grep -h '^slots_allocated=' smoke/b.txt 2>/dev/null | cut -d= -f2)"
    if [[ "$MODE" == lsf ]]; then
        [[ "$got" == "4" ]] && say ok "b_bigger was allocated 4 slots" \
                            || say FAIL "b_bigger asked for 4 slots, LSF gave '${got:-nothing}'"
    else
        say ok "slots_allocated=${got:-unset} (local mode: LSF is not involved)"
    fi
else
    say FAIL "summary.txt was not produced"
fi

# --- 2. the case that matters: a job that dies must be REPORTED, not awaited --
echo
timeout "$LIMIT" snakemake --snakefile "$SMK" "${EXEC[@]}" smoke/fail.txt \
    >"$WORK/run2.log" 2>&1
rc=$?
case "$rc" in
    124) say FAIL "a failing job HUNG the workflow -- this is the failure the status command exists to prevent"
         tail -n 15 "$WORK/run2.log" ;;
    0)   say FAIL "a job that exits 1 was reported as success" ;;
    *)   say ok   "a failing job was detected and reported (exit $rc)" ;;
esac

echo
if [[ "$MODE" == lsf ]]; then
    echo "workspace kept for inspection: $WORK"
    echo "  bsub logs:  $WORK/logs/lsf/"
    echo "  what ran:   cat $WORK/smoke/*.txt"
else
    rm -rf "$WORK"
fi

[[ "$FAIL" -ne 0 ]] && { echo "FAILED"; exit 1; }
echo "all checks passed"
