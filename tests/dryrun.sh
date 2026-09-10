#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# I0 gate for the Snakemake workflow: does it parse, and does it refuse what it
# is supposed to refuse.
#
#   tests/dryrun.sh
#
# Needs snakemake on PATH and SCENICPLUS_PATH set to the repository root. Builds
# a throwaway workspace from the shipped config template, so it touches nothing
# real and needs no data: at I0 there are no step rules, so nothing reads the
# .rds the config points at.
#
# Every case here is a negative one except the first. That is deliberate: a
# workflow that parses proves only that it is not broken, while a workflow that
# REFUSES a typo'd key is the thing being built. If one of these stops failing,
# the contract has a hole.
# -----------------------------------------------------------------------------
set -uo pipefail

: "${SCENICPLUS_PATH:?set SCENICPLUS_PATH to the repository root}"
command -v snakemake >/dev/null || { echo "snakemake not on PATH"; exit 1; }

# Snakemake writes a source cache under $XDG_CACHE_HOME, and on a box where
# $HOME is read-only that is the first thing to fail, with a traceback that
# looks like a workflow error rather than a permissions one.
WORK="$(mktemp -d)"
export XDG_CACHE_HOME="$WORK/cache"
trap 'rm -rf "$WORK"' EXIT

mkdir -p "$WORK/ws/config"
cp "$SCENICPLUS_PATH/config/config.yaml" "$WORK/ws/config/config.yaml"
cd "$WORK/ws"

# R01 declares the .rds as an INPUT, so the DAG cannot be built while
# input.seurat_rds is the template's placeholder path -- which is correct
# behaviour, and the reason this stand-in exists. Its contents never matter:
# every check here is a dry run or a refusal, so no rule ever opens it.
touch "$WORK/ws/stand-in.rds"
python3 - <<PY
import yaml, pathlib
c = yaml.safe_load(open("config/config.yaml"))
c["input"]["seurat_rds"] = "$WORK/ws/stand-in.rds"
pathlib.Path("config/config.yaml").write_text(yaml.safe_dump(c))
PY

FAIL=0
say() { printf '  %-4s %s\n' "$1" "$2"; [[ "$1" == FAIL ]] && FAIL=1; return 0; }

run_sm() { snakemake --snakefile "$SCENICPLUS_PATH/Snakefile" -n --cores 1 >"$WORK/out" 2>&1; }

edit_cfg() {   # edit_cfg <python expression on `c`>
    python3 - "$1" <<'PY'
import sys, yaml, pathlib
c = yaml.safe_load(open("config/config.yaml"))
exec(sys.argv[1])
pathlib.Path("config/config.yaml").write_text(yaml.safe_dump(c))
PY
}

cp config/config.yaml "$WORK/pristine.yaml"
restore() { cp "$WORK/pristine.yaml" config/config.yaml; }

echo "workflow: $SCENICPLUS_PATH/Snakefile"
echo

# 1. the positive case
run_sm && say ok "the shipped template parses and dry-runs" \
       || { say FAIL "the shipped template parses and dry-runs"; sed -n '1,12p' "$WORK/out"; }

# 2. a typo'd key must stop the DAG BUILD, before any job is scheduled
edit_cfg 'c["cistopic"]["n_topic"] = [10]'
if run_sm; then say FAIL "a typo'd key (n_topic) is refused"
elif grep -q "n_topic" "$WORK/out"; then say ok "a typo'd key (n_topic) is refused, and named"
else say FAIL "a typo'd key is refused, but the message does not name it"; fi
restore

# 3. a value of the wrong type, which used to surface inside python after LSF
#    had already scheduled the job
edit_cfg 'c["resources"]["n_cpu"] = "x"'
run_sm && say FAIL "a non-numeric n_cpu is refused" || say ok "a non-numeric n_cpu is refused"
restore

# 4. a species with no reference table. The schema allows it, because the bash
#    driver passes species straight through; the workflow must not.
edit_cfg 'c["input"]["species"] = "dmelanogaster"'
if run_sm; then say FAIL "a species with no reference row is refused"
elif grep -q "common.smk" "$WORK/out"; then say ok "a species with no reference row is refused, pointing at the table"
else say FAIL "refused, but the message does not say where the table is"; fi
restore

# 5. another pipeline's config must not be runnable here by accident
edit_cfg 'c["Pipeline"] = "scRNA_LimLab_Snake"'
run_sm && say FAIL "another pipeline's config is refused" || say ok "another pipeline's config is refused"
restore

# 6. the runner's own refusals
RUN="$SCENICPLUS_PATH/scripts/scenicplus.run.sh"

# `--lsf` is environment-dependent, so assert the contract rather than one
# outcome: with the executor plugin present it must build a profile invocation,
# and without it it must refuse WITH the install line rather than fall back to
# running everything on the submit host. Checking only one of those would pass
# vacuously on whichever machine happens to run this.
if snakemake --executor cluster-generic --help >/dev/null 2>&1; then
    SCENICPLUS_SKIP_CHECK=1 "$RUN" --lsf -n >/dev/null 2>&1
    [[ $? -eq 0 ]] && say ok "--lsf plans through profiles/lsf (executor present)" \
                   || say FAIL "--lsf failed even though the executor is installed"
else
    out="$(SCENICPLUS_SKIP_CHECK=1 "$RUN" --lsf 2>&1)"; rc=$?
    if [[ $rc -eq 2 ]] && grep -q "cluster-generic" <<<"$out"; then
        say ok "--lsf refuses and names the missing executor (not installed here)"
    else
        say FAIL "--lsf did not refuse cleanly without the executor (exit $rc)"
    fi
fi

SCENICPLUS_SKIP_CHECK=1 "$RUN" -f 7 -n >/dev/null 2>&1
[[ $? -eq 2 ]] && say ok "-f on a step with no rule refuses instead of forcing nothing" \
               || say FAIL "-f on a step with no rule refuses instead of forcing nothing"

# The other direction, which only became testable once rules existed: a step
# number must RESOLVE to its rule. Checking only the refusal would leave the
# lookup itself unexercised, and a `-f` that silently forces nothing is exactly
# what the runner's rule lookup exists to prevent.
out="$(SCENICPLUS_SKIP_CHECK=1 "$RUN" -f 1 -n 2>&1)"
if grep -q "forcing from R01_seurat_export onward" <<<"$out"; then
    say ok "-f 1 resolves to R01_seurat_export"
else
    say FAIL "-f 1 did not resolve to a rule name"
fi

SCENICPLUS_SKIP_CHECK=1 "$RUN" -n >/dev/null 2>&1
[[ $? -eq 0 ]] && say ok "the runner's dry run succeeds" \
               || say FAIL "the runner's dry run succeeds"

echo
if [[ "$FAIL" -ne 0 ]]; then echo "FAILED"; exit 1; fi
echo "all checks passed"
