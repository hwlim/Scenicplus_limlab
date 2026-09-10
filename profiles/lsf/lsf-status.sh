#!/usr/bin/env bash
# =============================================================================
# LSF job-status probe for snakemake's cluster-generic executor.
#
#   lsf-status.sh <what bsub printed>      ->  success | failed | running
#
# WHY THIS EXISTS
#
# Without a status command snakemake cannot ask LSF anything: it submits, then
# waits for the rule's output files to appear. A job killed by LSF -- walltime,
# memory, a dead node -- produces no outputs and no notification, so the
# orchestrator waits FOREVER on a job that is already gone. Observed in
# scRNA_LimLab_Snake, the sibling pipeline this file is adapted from: a markers
# job hit TERM_RUNLIMIT at that profile's `-W 120` and the workflow hung with
# nothing else runnable, because every remaining job depended on its output.
# Nothing here has been through an LSF failure yet -- tests/cluster_smoke.sh
# exists to put it through one on purpose.
#
# It also makes `restart-times` and `keep-going` meaningful. Both need a
# DETECTED failure; without this they are inert.
#
# ARGUMENT SHAPE. The plugin hands back whatever the submit command printed on
# stdout, and bsub prints "Job <12345> is submitted to queue <normal>." rather
# than a bare id -- so pull the first run of digits instead of assuming a clean
# job id. That keeps the submit command untouched, which matters because
# submission already works and is the riskier half to change.
#
# Must print EXACTLY one of the three words and exit 0. Anything else is taken
# as a protocol error by snakemake.
# =============================================================================
set -uo pipefail

raw="${1:-}"
jid="$(printf '%s' "$raw" | grep -oE '[0-9]+' | head -1)"

# No id to query. Report failed rather than running: a job we cannot name will
# never be observed to finish, and reporting `running` would reproduce exactly
# the hang this script exists to remove.
if [ -z "$jid" ]; then
    echo failed
    exit 0
fi

# `-o stat -noheader` needs LSF 10+; fall back to parsing plain bjobs (STAT is
# field 3) so this works on older sites too.
stat="$(bjobs -o stat -noheader "$jid" 2>/dev/null | tr -d '[:space:]')"
if [ -z "$stat" ]; then
    stat="$(bjobs -a -w "$jid" 2>/dev/null | awk 'NR==2 {print $3}' | tr -d '[:space:]')"
fi

# Finished jobs leave bjobs after the cluster's CLEAN_PERIOD; bhist still knows
# them. Without this fallback a job that completed while the poll interval
# elapsed would be reported failed.
if [ -z "$stat" ]; then
    hist="$(bhist -n 0 -l "$jid" 2>/dev/null)"
    case "$hist" in
        *"Done successfully"*) echo success; exit 0 ;;
        *Exited*|*TERM_*)      echo failed;  exit 0 ;;
    esac
    # Neither tool knows this job. It is not going to produce output.
    echo failed
    exit 0
fi

case "$stat" in
    DONE)                                  echo success ;;
    EXIT|ZOMBI)                            echo failed  ;;
    # UNKWN is a transient loss of contact with the execution host, not a
    # failure -- LSF recovers it. Treating it as failed would abort real work.
    PEND|RUN|PROV|WAIT|PSUSP|USUSP|SSUSP|UNKWN) echo running ;;
    *)                                     echo running ;;
esac
exit 0
