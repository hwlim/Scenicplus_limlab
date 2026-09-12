#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Gate for the provenance bundle (I7).
#
#   tests/provenance.sh
#
# Needs snakemake and python3. Runs a REAL miniature workflow -- twice, once
# succeeding and once failing -- rather than calling the script by hand.
#
# WHY A REAL RUN. The bundle is assembled by snakemake's `onstart` / `onsuccess`
# / `onerror` handlers, and those are the part most likely to be wrong: they do
# not fire on a dry run, an exception inside one aborts the workflow, and the
# `log` variable they read is not documented anywhere this repo controls.
# Calling `scenicplus_provenance.py` directly would exercise the easy half and
# skip the half that has to work on the cluster.
#
# THE FAILING RUN IS THE POINT. A successful run is already described by
# report.html; the bundle earns its keep on the run that stopped at rule 9 of 21
# and has nothing else. So the plan's stated gate for I7 is "a bundle from a
# FAILED run contains the failing log", and that is checked here against a job
# that genuinely dies.
# -----------------------------------------------------------------------------
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$HERE")"
PROV="$ROOT/scripts/scenicplus_provenance.py"
[[ -f "$PROV" ]] || { echo "missing $PROV"; exit 1; }
command -v snakemake >/dev/null || { echo "snakemake not on PATH"; exit 1; }

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
export XDG_CACHE_HOME="$WORK/cache"
FAIL=0
say() { printf '  %-4s %s\n' "$1" "$2"; [[ "$1" == FAIL ]] && FAIL=1; return 0; }

# A miniature workflow with the SAME handler wiring as the real Snakefile: the
# helper, `--mode start` at onstart, `--mode finish` at both outcomes. Copied
# rather than imported because the real Snakefile needs a full config, a genome
# pair and six stand-in databases -- and this is a test of the HANDLERS, not of
# the DAG.
mk_ws() {                       # mk_ws <dir> <shell-body-for-the-rule>
    local d="$1" body="$2"
    mkdir -p "$d/config" "$d/logs/lsf"
    printf 'Pipeline: "ScenicPlus"\ninput: {species: hsapiens}\n' \
        > "$d/config/config.yaml"
    # An LSF epilogue, so lsf_jobs.tsv has something real to parse.
    cat > "$d/logs/lsf/R09_cistarget.1.out" <<'LSF'
Job was executed on host(s) <16*bmi-200m5-04>, in queue <normal>, as user <x>
    Run time :                                   682 sec.
    Max Memory :                                 152566 MB
    Max Processes :                              28
    Max Threads :                                1151
LSF
    mkdir -p "$d/QC"; printf '{"assembly_detected": "hg38"}\n' > "$d/QC/assembly.json"
    # A report to copy. Without one the bundle's report.html path is untested,
    # which is exactly how it shipped unchecked the first time.
    printf '<html><body>THE-REPORT-BODY</body></html>\n' > "$d/report.html"
    # `( ... )` and not `{ ... }`: snakemake formats the shell string, so braces
    # are placeholders and a rule body wrapped in them dies with
    # "NameError: the name ' echo hello; touch out' is unknown". Parentheses are
    # inert to the formatter and give the same grouping. Learned here.
    cat > "$d/Snakefile" <<SMK
import os, subprocess
PROV = "$PROV"
CONFIG_FILE = os.path.abspath("config/config.yaml")
PIPELINE = "$ROOT"
def _provenance(mode, **kw):
    cmd = ["python", PROV, "--mode", mode, "--workspace", ".",
           "--pipeline", PIPELINE]
    for k, v in kw.items():
        if v:
            cmd += ["--" + k.replace("_", "-"), str(v)]
    try:
        subprocess.run(cmd, check=True)
    except Exception as e:
        print("provenance", mode, "failed", e)
onstart:
    _provenance("start")
onsuccess:
    _provenance("finish", config=CONFIG_FILE, status="success", snakemake_log=log)
onerror:
    _provenance("finish", config=CONFIG_FILE, status="error", snakemake_log=log)
# The real Snakefile sets this too, and without it `false | tee` exits 0 and the
# deliberately-failing case below would silently pass.
shell.prefix("set -o pipefail; ")
rule all:
    output: "out.txt"
    log: "logs/R99_thing.log"
    shell: "( $body ) 2>&1 | tee {log}"
SMK
}

bundle_of() { ls -d "$1"/provenance/*_"$2" 2>/dev/null | tail -1; }

# --- 1. a SUCCESSFUL run -----------------------------------------------------
WS="$WORK/ok"
mk_ws "$WS" 'echo hello; touch out.txt'
( cd "$WS" && snakemake -c1 >"$WORK/ok.out" 2>&1 )
rc=$?
[[ $rc -eq 0 ]] && say "ok" "a successful run completes with the handlers wired" \
   || { say FAIL "the successful run failed (rc=$rc)"; tail -12 "$WORK/ok.out"; }

B="$(bundle_of "$WS" success)"
[[ -n "$B" && -d "$B" ]] && say "ok" "it leaves a *_success bundle" \
                         || say FAIL "no *_success bundle"

if [[ -n "$B" ]]; then
    for f in manifest.txt config.used.yaml lsf_jobs.tsv assembly.json; do
        [[ -f "$B/$f" ]] && say "ok" "bundle carries $f" \
                         || say FAIL "bundle is missing $f"
    done
    [[ -f "$B/logs/R99_thing.log" ]] && say "ok" "bundle carries this run's rule log" \
                                     || say FAIL "bundle is missing the rule log"
    grep -q "bmi-200m5-04" "$B/lsf_jobs.tsv" \
        && say "ok" "lsf_jobs.tsv is parsed, not just created" \
        || say FAIL "lsf_jobs.tsv has no accounting in it"
    grep -q "^status: *success" "$B/manifest.txt" \
        && say "ok" "the manifest records the status" \
        || say FAIL "the manifest does not record the status"
    # Content, not existence: a zero-byte copy passes `-f` and is useless. The
    # bundle should be self-contained -- evidence and readable summary together.
    grep -q "THE-REPORT-BODY" "$B/report.html" 2>/dev/null \
        && say "ok" "the bundle carries report.html, with its content" \
        || say FAIL "report.html is missing from the bundle or empty"
    grep -q "^report.html: *included" "$B/manifest.txt" \
        && say "ok" "...and the manifest records that it is included" \
        || say FAIL "the manifest does not record report.html"
fi

# --- 2. THE COMMIT IS CAPTURED AT ONSTART, NOT AT FINISH ---------------------
# The whole reason `--mode start` exists. Recorded early and reported verbatim,
# so a checkout during the run cannot make the bundle name a commit that
# produced nothing. Simulated by rewriting the recorded meta between the two
# phases and confirming the bundle reports what was RECORDED.
WS2="$WORK/commit"
mk_ws "$WS2" 'echo hi; touch out.txt'
python3 "$PROV" --mode start --workspace "$WS2" --pipeline "$ROOT" >/dev/null 2>&1
python3 - "$WS2" <<'PY'
import json, sys
p = sys.argv[1] + "/logs/.run_meta.json"
m = json.load(open(p))
m["commit"] = "c0ffee1234567890"          # what HEAD was AT ONSTART
m["branch"] = "the-branch-that-ran"
json.dump(m, open(p, "w"))
PY
python3 "$PROV" --mode finish --workspace "$WS2" --pipeline "$ROOT" \
    --config "$WS2/config/config.yaml" --status success >/dev/null 2>&1
B2="$(bundle_of "$WS2" success)"
if grep -q "c0ffee1234567890" "$B2/manifest.txt" 2>/dev/null; then
    say "ok" "the bundle reports the commit recorded AT ONSTART"
else
    say FAIL "the bundle re-read HEAD at finish time instead of using the record"
fi

# The other half: with NO record, it must say UNKNOWN rather than substitute a
# fresh rev-parse. A plausible wrong commit is the failure being designed out.
WS3="$WORK/nostart"
mk_ws "$WS3" 'true'
python3 "$PROV" --mode finish --workspace "$WS3" --pipeline "$ROOT" \
    --config "$WS3/config/config.yaml" --status success >/dev/null 2>&1
B3="$(bundle_of "$WS3" success)"
if grep -q "^commit: *unknown" "$B3/manifest.txt" 2>/dev/null \
   && grep -q "did not go through" "$B3/manifest.txt"; then
    say "ok" "with no onstart record the commit is UNKNOWN, and it says why"
else
    say FAIL "a missing onstart record was filled in rather than reported"
fi

# --- 3. THE FAILING RUN, which is the plan's stated gate ---------------------
WS4="$WORK/bad"
mk_ws "$WS4" 'echo ValueError: nope; false'
# A real failed run has NO report: snakemake does not build a target whose
# inputs failed. The fixture seeds one for the success case, so remove it here
# BEFORE the run -- otherwise the failure path is tested against a workspace
# that cannot occur.
rm -f "$WS4/report.html"
( cd "$WS4" && snakemake -c1 >"$WORK/bad.out" 2>&1 )
rc=$?
[[ $rc -ne 0 ]] && say "ok" "the failing run really fails (rc=$rc)" \
                || say FAIL "the deliberately-failing run exited 0"

# rc!=0 alone does not prove the RULE failed -- a Snakefile that dies at parse
# also exits non-zero, and then no handler ever runs. That is exactly what a
# double quote in the rule body did here on the first attempt, and it read as
# "onerror produced no bundle". Require evidence the rule actually EXECUTED.
grep -q "ValueError" "$WS4/logs/R99_thing.log" 2>/dev/null \
    && say "ok" "...by running the rule, not by failing to parse" \
    || { say FAIL "the run died before the rule executed"; tail -6 "$WORK/bad.out"; }

B4="$(bundle_of "$WS4" error)"
[[ -n "$B4" && -d "$B4" ]] && say "ok" "a FAILED run still leaves a bundle" \
                           || say FAIL "no *_error bundle from the failed run"

if [[ -n "$B4" ]]; then
    if [[ -f "$B4/logs/R99_thing.log" ]] && grep -q "ValueError" "$B4/logs/R99_thing.log"; then
        say "ok" "the bundle contains the FAILING log, with the error in it"
    else
        say FAIL "the failing log is absent or empty -- I7's stated gate"
    fi
    grep -q "R99_thing.log" "$B4/manifest.txt" \
        && say "ok" "the manifest NAMES the log with the error signature" \
        || say FAIL "the manifest does not name the failing log"
    [[ -f "$B4/snakemake.log" ]] && say "ok" "the snakemake log is bundled too" \
                                 || say FAIL "the snakemake log is missing"
    [[ ! -f "$B4/report.html" ]] \
        && say "ok" "a failed run's bundle has no report.html" \
        || say FAIL "a report.html appeared in a failed run's bundle"
    # Absence has to be STATED. A bundle that merely lacks the file leaves a
    # reader guessing between "never built" and "lost on the way here".
    grep -q "^report.html: *NOT INCLUDED" "$B4/manifest.txt" \
        && say "ok" "...and the manifest says WHY it is absent" \
        || say FAIL "the manifest does not explain the missing report"
fi

# --- 4. logs are SCOPED to this run ------------------------------------------
# bsub -o APPENDS and {jobid} restarts at 0, so a stale log next to a fresh one
# is the normal case. Bundling an older run's log as this run's misattributes a
# failure, which is worse than omitting it.
WS5="$WORK/scope"
mk_ws "$WS5" 'echo new; touch out.txt'
printf 'FROM-AN-OLDER-RUN\n' > "$WS5/logs/R00_stale.log"
touch -t 202001010000 "$WS5/logs/R00_stale.log"
( cd "$WS5" && snakemake -c1 >/dev/null 2>&1 )
B5="$(bundle_of "$WS5" success)"
if [[ -n "$B5" ]]; then
    [[ ! -f "$B5/logs/R00_stale.log" ]] \
        && say "ok" "a log older than the run marker is NOT bundled" \
        || say FAIL "a stale log was bundled as this run's"
    [[ -f "$B5/logs/R99_thing.log" ]] \
        && say "ok" "...while this run's own log still is" \
        || say FAIL "scoping dropped this run's log too"
fi

# --- 5. the cap SKIPS rather than truncates ----------------------------------
# A half log is a trap: the interesting part of a traceback is at the END.
WS6="$WORK/cap"
mk_ws "$WS6" 'true'
python3 "$PROV" --mode start --workspace "$WS6" --pipeline "$ROOT" >/dev/null 2>&1
head -c 200000 /dev/zero | tr '\0' 'x' > "$WS6/logs/R98_big.log"
python3 "$PROV" --mode finish --workspace "$WS6" --pipeline "$ROOT" \
    --config "$WS6/config/config.yaml" --status success --cap-kb 100 >/dev/null 2>&1
B6="$(bundle_of "$WS6" success)"
if [[ -n "$B6" ]]; then
    if [[ ! -d "$B6/logs" ]] && grep -q "exceed the" "$B6/manifest.txt"; then
        say "ok" "over the cap NO logs are copied, and the manifest says so"
    else
        say FAIL "the cap truncated or silently copied anyway"
    fi
fi

# --- 6. a broken provenance script must NOT fail the run ---------------------
# The handler wrapper exists for this: provenance is bookkeeping, and a finished
# run reported as failed because its bundle broke would be a worse bug than no
# bundle at all.
WS7="$WORK/broken"
mk_ws "$WS7" 'echo hi; touch out.txt'
sed -i 's|^PROV = .*|PROV = "/nonexistent/provenance.py"|' "$WS7/Snakefile"
( cd "$WS7" && snakemake -c1 >"$WORK/broken.out" 2>&1 )
rc=$?
[[ $rc -eq 0 && -f "$WS7/out.txt" ]] \
    && say "ok" "a broken provenance script does not fail the run" \
    || { say FAIL "a broken provenance script took the run down (rc=$rc)"; }

# --- 7. THE REAL SNAKEFILE, not a copy of its shape --------------------------
# Everything above uses a miniature workflow that MIRRORS the wiring. That
# proves the script and the handler contract, and proves nothing about whether
# the real Snakefile calls them -- a typo there would leave every check above
# green. `dryrun.sh` cannot cover it either: handlers do not fire on a dry run.
#
# So: the real Snakefile, a real (non-dry) run, in a stand-in workspace where
# R01 is guaranteed to fail on an empty .rds. Failing is the point -- it is the
# path that has to work, and it is cheap because it fails immediately.
if [[ -n "${SCENICPLUS_PATH:-}" && -f "$ROOT/Snakefile" ]]; then
    WS8="$WORK/real"
    mkdir -p "$WS8/config"
    cp "$ROOT/config/config.yaml" "$WS8/config/" 2>/dev/null
    for f in a.rds ga.tsv cs.tsv c.f d.f m.tbl; do : > "$WS8/$f"; done
    python3 - "$WS8" <<'PY'
import yaml, pathlib, sys
w = sys.argv[1]
c = yaml.safe_load(open(w + "/config/config.yaml")); i = c["input"]
for k, v in dict(seurat_rds="a.rds", genome_annotation="ga.tsv",
                 chromsizes="cs.tsv", ctx_db="c.f", dem_db="d.f",
                 motif_annotations="m.tbl").items():
    i[k] = w + "/" + v
pathlib.Path(w + "/config/config.yaml").write_text(yaml.safe_dump(c))
PY
    ( cd "$WS8" && timeout 300 snakemake --snakefile "$ROOT/Snakefile" -c1 \
        --until R07_genome_annot >"$WORK/real.out" 2>&1 )
    B8="$(bundle_of "$WS8" error)"
    if [[ -n "$B8" && -f "$B8/manifest.txt" ]]; then
        say "ok" "the REAL Snakefile's handlers produce a bundle on failure"
        grep -q "^commit: *[0-9a-f]\{7,\}" "$B8/manifest.txt" \
            && say "ok" "...with a real commit captured at its onstart" \
            || say FAIL "the real run recorded no commit"
        [[ -s "$B8/lsf_jobs.tsv" ]] \
            && say "ok" "...and lsf_jobs.tsv, even with no LSF accounting to find" \
            || say FAIL "lsf_jobs.tsv missing from the real run's bundle"
    else
        say FAIL "the REAL Snakefile produced no bundle"
        tail -8 "$WORK/real.out"
    fi
else
    say "ok" "(real-Snakefile check skipped: SCENICPLUS_PATH unset)"
fi

echo
if [[ $FAIL -eq 0 ]]; then echo "all checks passed"; else echo "FAILED"; fi
exit $FAIL
