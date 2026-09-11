#!/usr/bin/env python3
"""Gate for the per-rule resource tiers (I5).

    python tests/test_resources.py

Needs neither snakemake nor a cluster: it reads `rules/common.smk` and the
`.smk` files as TEXT, and `tests/measured_resources.tsv` as data.

WHY THIS EXISTS. A resource declaration is the one kind of code whose defect
never shows up locally and never shows up in a dry run. It surfaces as a job
that dies at TERM_MEMLIMIT an hour into a cluster run, or -- worse, and this is
the failure this pipeline actually had -- as a job that silently oversubscribes
its node and completes, so nothing reports anything and the number stays wrong.

Four checks, each aimed at a different way the tiers can rot:

1. EVERY RULE HAS A TIER. A rule that declares no `resources:` inherits the
   profile's `default-resources` of 8000 MB, which is under the measured peak of
   twelve of the eighteen measured steps. Falling through is silent.

2. EVERY RULE ASKS FOR ITS OWN. `mem()` takes a rule NAME, so a rule copied from
   its neighbour keeps the neighbour's string and is sized for the wrong work.
   That is invisible to every other check here, because both names are valid.

3. EVERY TIER CLEARS THE HEADROOM POLICY against the measurement. This is the
   one that fails when someone edits a tier down, or adds a rule and guesses.

4. THE MEASUREMENT COVERS THE RULES. A rule with no row in the TSV is reported
   as unmeasured rather than skipped, because a check that quietly skips what it
   cannot verify reports an all-clear it did not earn.
"""
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COMMON = os.path.join(ROOT, "rules", "common.smk")
MEASURED = os.path.join(ROOT, "tests", "measured_resources.tsv")

fails = []


def ok(label, cond, detail=""):
    print(f"{'ok  ' if cond else 'FAIL'}  {label}")
    if not cond:
        if detail:
            print(f"        {detail}")
        fails.append(label)


# --- load the tables out of common.smk --------------------------------------
# By exec'ing the literal blocks rather than importing: common.smk is a
# snakemake file and cannot be imported as a module. Each block is anchored on
# its own name at column 0 through the closing brace, and the match count is
# asserted -- a silently-renamed table would otherwise test nothing.
src = open(COMMON).read()
TABLES = {}
for name in ("MEM_TIERS", "TIME_TIERS", "MEM_HEADROOM", "RULE_TIERS"):
    hits = re.findall(rf"^{name} = \{{.*?^\}}", src, re.S | re.M)
    if len(hits) != 1:
        sys.exit(f"rules/common.smk: expected exactly one {name} block, found {len(hits)}")
    exec(hits[0], TABLES)
MEM_TIERS, TIME_TIERS = TABLES["MEM_TIERS"], TABLES["TIME_TIERS"]
HEADROOM, RULE_TIERS = TABLES["MEM_HEADROOM"], TABLES["RULE_TIERS"]

# --- load the measurement ----------------------------------------------------
measured = {}
with open(MEASURED) as fh:
    for line in fh:
        if line.startswith("#") or line.startswith("rule\t") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        measured[f[0]] = {"runtime_s": int(f[3]), "max_mem_mb": int(f[4])}
ok(f"measured_resources.tsv parsed ({len(measured)} rules)", len(measured) > 0)

# --- find the rules, and what each one asks for ------------------------------
# A rule block runs from `rule NAME:` to the next one. The leading `[ \t]*` is
# load-bearing: R07 is nested inside `if GENOME_SUPPLIED:`, so a column-0 anchor
# silently omits it -- which is exactly what the first version of this file did,
# and the count check below is what caught it.
rules = {}
for smk in sorted(os.listdir(os.path.join(ROOT, "rules"))):
    if not smk.endswith(".smk") or smk == "common.smk":
        continue
    text = open(os.path.join(ROOT, "rules", smk)).read()
    starts = [(m.start(), m.group(1))
              for m in re.finditer(r"^[ \t]*rule (\w+):", text, re.M)]
    for i, (pos, name) in enumerate(starts):
        end = starts[i + 1][0] if i + 1 < len(starts) else len(text)
        rules[name] = {"file": smk, "body": text[pos:end]}
ok(f"found the step rules ({len(rules)})", len(rules) == 20,
   f"got {len(rules)}: {', '.join(sorted(rules))}")

# --- 1. every rule declares memory AND runtime -------------------------------
missing = [n for n, r in rules.items()
           if "mem_mb=" not in r["body"] or "runtime=" not in r["body"]]
ok("every rule declares mem_mb and runtime", not missing,
   f"would inherit the profile's 8000 MB default: {', '.join(sorted(missing))}")

# --- 2. every rule asks for ITS OWN tier -------------------------------------
wrong = []
for name, r in rules.items():
    for fn in ("mem", "rt"):
        for arg in re.findall(rf"\b{fn}\(\s*[\"']([^\"']+)[\"']", r["body"]):
            if arg != name:
                wrong.append(f"{r['file']}: rule {name} calls {fn}({arg!r})")
ok("every rule asks for its own name, not a neighbour's", not wrong,
   "; ".join(wrong))

# --- 3. every rule in the .smk files has a row in RULE_TIERS -----------------
untabled = sorted(set(rules) - set(RULE_TIERS))
orphaned = sorted(set(RULE_TIERS) - set(rules))
ok("RULE_TIERS covers exactly the rules that exist",
   not untabled and not orphaned,
   f"missing from RULE_TIERS: {untabled or 'none'}; "
   f"in RULE_TIERS but no such rule: {orphaned or 'none'}")

# --- 4. the headroom policy, per rule ----------------------------------------
# Runtime is held to 3x for every rule; memory to the factor its growth class
# names. Both are checked against the SAME table the tiers were derived from,
# so editing a tier down without editing the measurement turns this red.
TIME_HEADROOM = 3.0
violations, unmeasured = [], []
for name, (mt, tt, growth) in sorted(RULE_TIERS.items()):
    if mt not in MEM_TIERS:
        violations.append(f"{name}: unknown memory tier {mt!r}")
        continue
    if tt not in TIME_TIERS:
        violations.append(f"{name}: unknown time tier {tt!r}")
        continue
    if growth not in HEADROOM:
        violations.append(f"{name}: unknown growth class {growth!r}")
        continue
    if name not in measured:
        unmeasured.append(name)
        continue
    m = measured[name]
    need_mem = m["max_mem_mb"] * HEADROOM[growth]
    need_min = (m["runtime_s"] / 60.0) * TIME_HEADROOM
    if MEM_TIERS[mt] < need_mem:
        violations.append(
            f"{name}: tier {mt} = {MEM_TIERS[mt]} MB, but the measured peak "
            f"{m['max_mem_mb']} MB at {HEADROOM[growth]}x needs {need_mem:.0f} MB")
    if TIME_TIERS[tt] < need_min:
        violations.append(
            f"{name}: tier {tt} = {TIME_TIERS[tt]} min, but the measured "
            f"{m['runtime_s']} s at {TIME_HEADROOM}x needs {need_min:.1f} min")

ok("every measured rule's tier clears the headroom policy", not violations,
   "\n        ".join(violations))

# Reported, never silently skipped. These are the rules whose numbers are a
# judgement call, and saying so is the whole point.
print(f"\n{len(unmeasured)} rule(s) carry UNMEASURED tiers: {', '.join(unmeasured) or 'none'}")
print("  They are sized by analogy. Replace with figures from the first run "
      "that includes them.")

# --- the table, for the reader -----------------------------------------------
print(f"\n{'rule':<22}{'peak MB':>9}{'tier':>7}{'MB':>8}{'head':>7}   "
      f"{'sec':>6}{'tier':>8}{'min':>6}")
for name, (mt, tt, growth) in sorted(RULE_TIERS.items()):
    if name in measured:
        m = measured[name]
        print(f"{name:<22}{m['max_mem_mb']:>9}{mt:>7}{MEM_TIERS[mt]:>8}"
              f"{MEM_TIERS[mt] / m['max_mem_mb']:>6.1f}x   "
              f"{m['runtime_s']:>6}{tt:>8}{TIME_TIERS[tt]:>6}")
    else:
        print(f"{name:<22}{'--':>9}{mt:>7}{MEM_TIERS[mt]:>8}{'n/a':>7}   "
              f"{'--':>6}{tt:>8}{TIME_TIERS[tt]:>6}")

total = sum(MEM_TIERS[RULE_TIERS[n][0]] for n in RULE_TIERS)
print(f"\nSum of every reservation: {total} MB. Not what the workflow holds at "
      f"once -- that is set by\nthe DAG and by `jobs:`; the single bsub'd driver "
      f"held {20 * 128000} MB for the whole run.")

print()
if fails:
    print(f"FAILED: {len(fails)}")
    for f in fails:
        print(f"  {f}")
    sys.exit(1)
print("all resource checks passed")
