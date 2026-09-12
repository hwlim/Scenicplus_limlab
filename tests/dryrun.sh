#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Gate for the Snakemake workflow (I0-I4): does it parse, and does it refuse what it
# is supposed to refuse.
#
#   tests/dryrun.sh
#
# Needs snakemake on PATH and SCENICPLUS_PATH set to the repository root. Builds
# a throwaway workspace from the shipped config template, so it touches nothing
# real and needs no data: every check is a dry run or a refusal, so no rule ever
# opens the files the config points at.
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

# The GRN rules require the genome pair, the cisTarget databases and the motif
# table as real INPUTS, so the DAG cannot be built while the template holds
# placeholder paths. Stand-ins for all of them: nothing here runs a rule, so no
# file is ever opened. Requiring them is the point -- a config that names a
# missing database should fail at DAG build, not after a queue wait.
for f in stand-in.rds genome_annotation.tsv chromsizes.tsv ctx.feather dem.feather motifs.tbl; do
    touch "$WORK/ws/$f"
done
python3 - <<PY
import yaml, pathlib
c = yaml.safe_load(open("config/config.yaml"))
i = c["input"]
i["seurat_rds"]        = "$WORK/ws/stand-in.rds"
i["genome_annotation"] = "$WORK/ws/genome_annotation.tsv"
i["chromsizes"]        = "$WORK/ws/chromsizes.tsv"
i["ctx_db"]            = "$WORK/ws/ctx.feather"
i["dem_db"]            = "$WORK/ws/dem.feather"
i["motif_annotations"] = "$WORK/ws/motifs.tbl"
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

# 6. half a genome pair, which is worse than none: the download that would fill
#    the gap cannot produce chromsizes for anyone, so a config with only the
#    annotation set looks configured and is not.
edit_cfg 'c["input"]["chromsizes"] = ""'
if run_sm; then say FAIL "only ONE of the genome pair is refused"
elif grep -q "BOTH or NEITHER" "$WORK/out"; then
     say ok "only ONE of the genome pair is refused, saying both or neither"
else say FAIL "refused, but not with the both-or-neither reason"; fi
restore

# 7. the DAG's SHAPE, not just that it builds.
#
# The whole argument for I3 is that declaring each stage's real inputs recovers
# the parallelism flattening gave up. That is a claim about edges, so check
# edges. Properties rather than a snapshot of the whole graph: a snapshot breaks
# on every legitimate change and teaches people to re-bless it.
snakemake --snakefile "$SCENICPLUS_PATH/Snakefile" --dag >"$WORK/dag.dot" 2>/dev/null
python3 - "$WORK/dag.dot" <<'PY' >"$WORK/dag.txt" 2>&1
import re, sys
t = open(sys.argv[1]).read()
label = dict(re.findall(r'(\d+)\[label = "([^"\\]+)', t))
deps = {}
for a, b in re.findall(r'(\d+) -> (\d+)', t):          # a -> b : b needs a
    deps.setdefault(label.get(b, b), set()).add(label.get(a, a))
rules = {r for r in label.values() if r.startswith("R")}
checks = [
    # 20 steps + R21_report, which is not a step: the bash driver has no
    # equivalent and `rule all` targets it (I6).
    ("all 20 step rules plus R21_report are in the graph",
     len(rules) == 21 and "R21_report" in rules, sorted(rules)),
    # The non-obvious edge: tf_to_gene reads tf_names.txt, which prepare_menr
    # writes. The driver satisfied this by being sequential, not by declaring it.
    ("R12 needs R11, for tf_names.txt",
     "R11_prepare_menr" in deps.get("R12_tf_to_gene", set()), None),
    ("R08 needs R07, so the genome pair is checked first",
     "R07_genome_annot" in deps.get("R08_search_space", set()), None),
    # The edge this suite MISSED. dem takes --genome_annotation under
    # `dem_balance_number_of_promoters` (shipped default: true) and the rule did
    # not declare it, so a changed genome pair re-ran R07, R08 and R13 while
    # R10's promoter balancing kept using the old annotation -- no error, other
    # motifs. Nothing above could see it: every check here was about the edges
    # that WERE declared. 7c asserts the other half, that the edge disappears
    # when the branch is off.
    ("R10 needs R07, because balanced dem reads the annotation",
     "R07_genome_annot" in deps.get("R10_dem", set()),
     sorted(deps.get("R10_dem", set()))),
    ("R13 needs R08, for the search space",
     "R08_search_space" in deps.get("R13_region_to_gene", set()), None),
    # The parallelism the increment exists for.
    ("R09 and R10 do not depend on each other",
     "R10_dem" not in deps.get("R09_cistarget", set())
     and "R09_cistarget" not in deps.get("R10_dem", set()), None),
    ("R14 and R15 do not depend on each other",
     "R15_egrn_extended" not in deps.get("R14_egrn_direct", set())
     and "R14_egrn_direct" not in deps.get("R15_egrn_extended", set()), None),
    ("R16 and R17 do not depend on each other",
     "R17_aucell_extended" not in deps.get("R16_aucell_direct", set())
     and "R16_aucell_direct" not in deps.get("R17_aucell_extended", set()), None),
]
bad = [(n, d) for n, ok, d in checks if not ok]
for n, d in bad:
    print(f"{n} :: FAILED" + (f" :: {d}" if d else ""))
sys.exit(1 if bad else 0)
PY
if [[ $? -eq 0 ]]; then
    say ok "the DAG has the right shape: 21 rules, and the four pairs stay parallel"
else
    say FAIL "the DAG's shape is wrong"; sed -n '1,5p' "$WORK/dag.txt"
fi

# 7c. ...and the R10->R07 edge is CONDITIONAL, because the argument is.
#
# Without this, declaring the annotation unconditionally would pass 7 just as
# well, and the false setting would carry a rerun trigger for a file dem never
# opens -- re-running a 12.9 GB-database step whenever the genome pair changed
# for someone who had switched the balancing off. Asserting an edge's ABSENCE
# is also what proves 7's version of it is measuring something: both cases
# build the same graph and differ only here.
edit_cfg 'c["scenicplus"]["dem_balance_number_of_promoters"] = False'
snakemake --snakefile "$SCENICPLUS_PATH/Snakefile" --dag >"$WORK/dag_nb.dot" 2>/dev/null
python3 - "$WORK/dag_nb.dot" <<'PY'
import re, sys
t = open(sys.argv[1]).read()
label = dict(re.findall(r'(\d+)\[label = "([^"\\]+)', t))
deps = {}
for a, b in re.findall(r'(\d+) -> (\d+)', t):
    deps.setdefault(label.get(b, b), set()).add(label.get(a, a))
dem = deps.get("R10_dem", set())
# The graph must still be a graph -- an empty parse would satisfy the assertion
# below for the wrong reason, which is how a check stops being one.
assert "R10_dem" in label.values(), "R10_dem absent from the DAG; parse failed"
assert dem, f"R10_dem has no dependencies at all: {sorted(label.values())[:5]}"
sys.exit(1 if "R07_genome_annot" in dem else 0)
PY
if [[ $? -eq 0 ]]; then
    say ok "...and with balancing OFF the R10->R07 edge is gone"
else
    say FAIL "R10 still depends on R07 with dem_balance_number_of_promoters false"
fi
restore

# 7d. THE DOCUMENTED `-f N` EXAMPLES, against the real DAG.
#
# `-f N` becomes `--forcerun R<NN>_*`, so snakemake re-runs that rule and its
# DEPENDENTS -- not every step numbered N or higher, because this workflow
# forks. Four documents said "from step N onward", copying the obsolete
# driver's semantics, where `--from N` really does compare numbers. Following
# that, `-f 9` looks like it rebuilds dem and region_to_gene and does not, and
# the run still ends green.
#
# The fix was prose, so the prose is what can rot. These examples are now
# derived from the resolved graph and compared against what the documents
# claim, which also means a DAG change surfaces here rather than in a reader's
# wrong expectation. It already would have: adding R10's genome-annotation edge
# (7 above) changed R07's skip set from {R09, R10} to {R09}.
snakemake --snakefile "$SCENICPLUS_PATH/Snakefile" --dag >"$WORK/dag2.dot" 2>/dev/null
python3 - "$WORK/dag2.dot" "$SCENICPLUS_PATH" <<'PY' >"$WORK/fdoc.txt" 2>&1
import os, re, sys
dot, root = sys.argv[1], sys.argv[2]
t = open(dot).read()
label = dict(re.findall(r'(\d+)\[label = "([^"\\]+)', t))
kids = {}
for a, b in re.findall(r'(\d+) -> (\d+)', t):
    kids.setdefault(label.get(a, a), set()).add(label.get(b, b))
rules = sorted(r for r in set(label.values()) if re.match(r"R\d\d_", r))
assert len(rules) == 21, f"expected 21 rules, parsed {len(rules)}"


def skipped(n):
    """Numbers >= n that --forcerun R<n> does NOT reach."""
    start = next(r for r in rules if r.startswith(f"R{n:02d}_"))
    seen, stack = {start}, [start]
    while stack:
        for k in kids.get(stack.pop(), ()):
            if k not in seen:
                seen.add(k); stack.append(k)
    return sorted(int(r[1:3]) for r in rules
                  if int(r[1:3]) >= n and r not in seen)


# What the DAG says, for the three fork points the documents use as examples.
actual = {n: skipped(n) for n in (7, 9, 14)}
bad = []
for n, want in actual.items():
    if not want:
        bad.append(f"-f {n} skips nothing, so it is a poor example to document")

# What the documents say. Each must name exactly the rules the DAG skips.
docs = ["RUNBOOK.md", "quickstart.md", "scripts/scenicplus.run.sh",
        "scripts/scenicplus_run_pipeline.sh", "README.md", "CLAUDE.md"]
for d in docs:
    text = open(os.path.join(root, d)).read()
    # The retired wording. It is what every one of these files used to say.
    for phrase in ("from step 7 onward", "from step 9 onward", "force from step"):
        if phrase in text:
            bad.append(f"{d}: still says {phrase!r}, which is the driver's "
                       f"semantics, not --forcerun's")
    # The generic form too, but only on a line that is about the RUNNER. The
    # obsolete driver's own header says "force re-run from step N onward" and
    # is CORRECT to: it walks a line and compares numbers. Flagging that would
    # be a check punishing the one place the phrase is true.
    for i, ln in enumerate(text.splitlines(), 1):
        if "from step N onward" in ln and ("-f " in ln or "scenicplus.run.sh" in ln):
            bad.append(f"{d}:{i}: describes `-f N` as 'from step N onward'")
    # Where a file gives the skip list for a fork point, it must be right.
    for n, want in actual.items():
        m = re.search(rf"-f {n}\D{{0,40}}?skips? (R\d\d(?:[,\s]+(?:and )?R\d\d)*)",
                      text)
        if not m:
            continue
        got = sorted(int(x) for x in re.findall(r"R(\d\d)", m.group(1)))
        if got != want:
            bad.append(f"{d}: says -f {n} skips {got}, DAG says {want}")

print("\n".join(bad) if bad else
      "  ".join(f"-f {n} skips {v}" for n, v in actual.items()))
sys.exit(1 if bad else 0)
PY
if [[ $? -eq 0 ]]; then
    say ok "the documented -f examples match the DAG ($(cat "$WORK/fdoc.txt"))"
else
    say FAIL "the -f documentation disagrees with the DAG"; sed -n '1,6p' "$WORK/fdoc.txt"
fi

# 7e. scenicplus_init.sh must point a new user at the RUNNER.
#
# It printed "Edit this file, then run: scenicplus_run_workstation.sh /
# scenicplus_run_lsf.sh" -- so the very first instruction a new user received
# was to start the obsolete driver, months after it was retired. Checked by
# RUNNING it, not by grepping the source: what matters is what reaches the
# terminal.
_INIT_WS="$WORK/initws"; mkdir -p "$_INIT_WS"
( cd "$_INIT_WS" && bash "$SCENICPLUS_PATH/scripts/scenicplus_init.sh" ) \
    >"$WORK/init.out" 2>&1
if grep -q "scenicplus.run.sh" "$WORK/init.out"; then
    say ok "scenicplus_init.sh points a new user at scenicplus.run.sh"
else
    say FAIL "scenicplus_init.sh does not name the runner"; cat "$WORK/init.out"
fi
if grep -qE 'scenicplus_run_(workstation|lsf|pipeline)\.sh' "$WORK/init.out"; then
    say FAIL "scenicplus_init.sh still tells a new user to run the obsolete driver"
else
    say ok "...and does not name the obsolete launchers"
fi

# 7f. the ORCHESTRATOR is not itself an LSF job, and the docs must say so.
#
# `--lsf` submits one bsub per RULE. Snakemake stays in the calling shell and is
# the only thing polling LSF and scheduling the next rule, so a dropped ssh
# session halts the run halfway -- the queued jobs finish, nothing starts what
# follows, and NOTHING FAILS. Every document said "--lsf submits each rule as
# its own job" and stopped there, which reads as though the whole thing were
# submitted; the retired `scenicplus_run_lsf.sh` really did bsub the driver, so
# the conversion tables were promising a property the replacement had dropped.
_BSUB_TMPL="$SCENICPLUS_PATH/scripts/scenicplus.bsub.sh"
if bash -n "$_BSUB_TMPL" 2>"$WORK/bsub.err"; then
    say ok "scenicplus.bsub.sh parses"
else
    say FAIL "scenicplus.bsub.sh has a syntax error"; sed -n '1,3p' "$WORK/bsub.err"
fi
python3 - "$_BSUB_TMPL" "$SCENICPLUS_PATH" <<'PY' >"$WORK/bsub.txt" 2>&1
import os, re, sys
tmpl, root = sys.argv[1], sys.argv[2]
t = open(tmpl).read()
bad = []
# The directives LSF needs. -W above all: an orchestrator that hits its own
# runlimit produces exactly the silent halt this file exists to prevent.
for flag in ("-n", "-W", "-M", "-J", "-oo", "-eo"):
    if not re.search(rf"^#BSUB {re.escape(flag)}\b", t, re.M):
        bad.append(f"the template has no `#BSUB {flag}` line")
# It must drive the RUNNER. Pointing a new wrapper at the retired driver would
# reintroduce the one-oversized-job problem the workflow exists to fix.
if "scenicplus.run.sh" not in t:
    bad.append("the template does not invoke scenicplus.run.sh")
if re.search(r"scenicplus_run_(pipeline|workstation|lsf)\.sh", t):
    bad.append("the template names an obsolete launcher")

# Every document that tells someone to use --lsf must also name the wrapper.
# Coarse on purpose: what must not happen is the warning being dropped while
# the --lsf recipe stays.
for d in ("RUNBOOK.md", "quickstart.md", "README.md",
          "scripts/scenicplus.run.sh", "CLAUDE.md"):
    text = open(os.path.join(root, d)).read()
    if "--lsf" in text and "scenicplus.bsub.sh" not in text:
        bad.append(f"{d}: documents --lsf without pointing at scenicplus.bsub.sh")
print("\n".join(bad) if bad else "template and 5 documents agree")
sys.exit(1 if bad else 0)
PY
if [[ $? -eq 0 ]]; then
    say ok "the bsub template is complete, and every --lsf doc points at it"
else
    say FAIL "the orchestrator-lifetime documentation has a hole"
    sed -n '1,5p' "$WORK/bsub.txt"
fi

# 7b. the per-rule resources SNAKEMAKE RESOLVES, against the table they come
#     from (I5).
#
# tests/test_resources.py reads the same table and the .smk files as text, which
# proves the declarations are right and the tiers clear the measurement. It
# cannot prove snakemake agrees: a rule whose `resources:` never took effect --
# shadowed by a later key, lost to an indentation slip inside `if
# GENOME_SUPPLIED:` -- still reads correctly in the file and still falls through
# to the profile's 8000 MB on the cluster. This is the check that opens the
# resolved job and looks.
snakemake --snakefile "$SCENICPLUS_PATH/Snakefile" -n --cores 1 --verbose \
    >"$WORK/verbose.txt" 2>&1
python3 - "$WORK/verbose.txt" "$SCENICPLUS_PATH/rules/common.smk" <<'PY' >"$WORK/res.txt" 2>&1
import re, sys
txt = open(sys.argv[1]).read()
# Local dry runs emit `localrule`, cluster plans `rule`; accept both.
#
# SPLIT INTO BLOCKS FIRST, rather than one regex spanning from a rule header to
# the next `resources:` line. A rule that declares NO resources has no such line
# in its block, and a spanning pattern then reaches into the NEXT job and
# reports that job's numbers under this rule's name -- a true failure with a
# false explanation. Measured while mutation-testing this very check.
heads = [(m.start(), m.group(1))
         for m in re.finditer(r"^(?:local)?rule (\w+):", txt, re.M)]
got, nores = {}, []
for i, (pos, name) in enumerate(heads):
    block = txt[pos:heads[i + 1][0] if i + 1 < len(heads) else len(txt)]
    m = re.search(r"^    resources:[^\n]*?mem_mb=(\d+)[^\n]*?runtime=(\d+)", block, re.M)
    if m:
        got[name] = (int(m.group(1)), int(m.group(2)))
    elif name != "all":
        nores.append(name)
src, ns = open(sys.argv[2]).read(), {}
for b in ("MEM_TIERS", "TIME_TIERS", "RULE_TIERS"):
    blk = re.findall(rf"^{b} = \{{.*?^\}}", src, re.S | re.M)
    if len(blk) != 1:
        sys.exit(f"common.smk: expected one {b} block, found {len(blk)}")
    exec(blk[0], ns)
steps = {r: v for r, v in got.items() if r != "all"}
bad = []
for r in nores:
    bad.append(f"{r}: resolved with NO mem_mb/runtime -- it would take the "
               f"profile's default-resources on the cluster")
if len(steps) + len(nores) != 21:
    bad.append(f"resolved {len(steps) + len(nores)} step jobs, expected 21: "
               f"{sorted(set(steps) | set(nores))}")
for r, (mb, mn) in sorted(steps.items()):
    want = ns["RULE_TIERS"].get(r)
    if not want:
        bad.append(f"{r}: resolved but absent from RULE_TIERS")
        continue
    wm, wt = ns["MEM_TIERS"][want[0]], ns["TIME_TIERS"][want[1]]
    if (mb, mn) != (wm, wt):
        bad.append(f"{r}: snakemake resolved {mb} MB / {mn} min, "
                   f"table says {wm} / {wt}")
for b in bad:
    print(b)
sys.exit(1 if bad else 0)
PY
if [[ $? -eq 0 ]]; then
    say ok "every rule's resolved mem/runtime matches rules/common.smk"
else
    say FAIL "a rule's resolved resources disagree with the tier table"
    sed -n '1,5p' "$WORK/res.txt"
fi

# 8. the reproducibility pinning actually reaches a job's environment.
#
# Behavioural, not a grep: this includes the REAL common.smk, calls the REAL
# shell_prefix(), and reads the environment from inside a running rule. A rule
# body that does not inherit the pins would produce a run that LOOKS
# reproducible and is not, which is worse than not pinning at all.
mkdir -p "$WORK/pin"
cat > "$WORK/pin/cfg.yaml" <<'YML'
resources:
  n_cpu: 7
YML
cat > "$WORK/pin/Snakefile" <<SMK
configfile: "cfg.yaml"
include: "$SCENICPLUS_PATH/rules/common.smk"
shell.prefix(shell_prefix())
rule all:
    output: "env.txt"
    shell: "env | grep -E 'PYTHONHASHSEED|NUM_THREADS' | sort > {output}; "
           "python -c \\"print(hash('ATF5'))\\" >> {output}; "
           "python -c \\"print(hash('ATF5'))\\" >> {output}"
SMK
( cd "$WORK/pin" && snakemake -c1 >/dev/null 2>&1 )
if [[ -s "$WORK/pin/env.txt" ]]; then
    miss=""
    for v in PYTHONHASHSEED OMP_NUM_THREADS OPENBLAS_NUM_THREADS MKL_NUM_THREADS NUMEXPR_NUM_THREADS; do
        grep -q "^$v=" "$WORK/pin/env.txt" || miss="$miss $v"
    done
    # n_cpu is 7 in that config, so the thread vars must say 7 -- proving the
    # value is derived rather than hardcoded.
    grep -q "^OMP_NUM_THREADS=7$" "$WORK/pin/env.txt" || miss="$miss OMP!=n_cpu"
    # Two separate interpreters must agree on a string's hash.
    h="$(grep -c . <(sort -u <(tail -2 "$WORK/pin/env.txt")))"
    [[ "$h" == "1" ]] || miss="$miss hash-not-stable"
    [[ -z "$miss" ]] && say ok "every rule inherits the hash seed and BLAS thread pinning" \
                     || say FAIL "pinning incomplete:$miss"
else
    say FAIL "the pinning probe did not run"
fi

# 9. the runner's own refusals
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

# A `-f` that silently forces NOTHING is what the lookup exists to prevent, so
# this needs a number with genuinely no rule behind it.
#
# DERIVED, NOT HARDCODED. This check used to pass a literal 21, with a comment
# claiming 21 "stays out of range however far the increments get". I6 added
# R21_report and falsified it, turning a real check into a failure about its own
# assumption. Ask --list what exists and take the first gap after it, so the
# number cannot go stale again.
_MAXN="$(snakemake --snakefile "$SCENICPLUS_PATH/Snakefile" --list 2>/dev/null \
         | sed -n 's/^R\([0-9][0-9]\)_.*/\1/p' | sort -n | tail -1)"
_OOR=$(( 10#${_MAXN:-0} + 1 ))
SCENICPLUS_SKIP_CHECK=1 "$RUN" -f "$_OOR" -n >/dev/null 2>&1
[[ $? -eq 2 ]] && say ok "-f $_OOR, which has no rule at all, refuses instead of forcing nothing" \
               || say FAIL "-f $_OOR did not refuse, though no R$(printf '%02d' $_OOR)_ rule exists"

# The report rule is numbered like the steps, so the number is the interface for
# it too -- and "redraw the report" is the single most likely thing anyone wants
# to force. Checked explicitly, because it is the boundary the case above moved.
out="$(SCENICPLUS_SKIP_CHECK=1 "$RUN" -f 21 -n 2>&1)"
grep -q "forcing R21_report and everything downstream" <<<"$out" \
  && say ok "-f 21 resolves to R21_report, so the report can be redrawn by number" \
  || { say FAIL "-f 21 did not resolve to R21_report"; sed -n '1,4p' <<<"$out"; }

# The other direction, which only became testable once rules existed: a step
# number must RESOLVE to its rule. Checking only the refusal would leave the
# lookup itself unexercised, and a `-f` that silently forces nothing is exactly
# what the runner's rule lookup exists to prevent.
out="$(SCENICPLUS_SKIP_CHECK=1 "$RUN" -f 1 -n 2>&1)"
if grep -q "forcing R01_seurat_export and everything downstream" <<<"$out"; then
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
