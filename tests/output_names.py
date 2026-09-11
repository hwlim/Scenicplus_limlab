#!/usr/bin/env python
"""Do the rules declare the names the scripts actually write?

A rule that declares an output its script never writes fails with
MissingOutputException -- at the END of the run, after every expensive stage has
already succeeded. A rule that omits one a script always writes loses the check
that caught the plotnine regression, where step 20 crashed before producing the
two heatmap-dotplots.

Neither is visible to a dry run: `-n` never executes a script, so it cannot know
what the script would have written. Nothing else compares these two files.

This compares two INDEPENDENT sources -- the declarations in rules/report.smk and
the write calls in the scripts -- rather than either against itself. Textual on
both sides on purpose: importing the Snakefile needs a config and a workspace,
and the failure being guarded against is a mismatched string.

    python tests/output_names.py
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FAIL = []


def say(ok, msg, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {msg}" + (f"\n          {detail}" if detail and not ok else ""))
    if not ok:
        FAIL.append(msg)


# --- what the figure script writes ------------------------------------------
# save(fig, out_dir, "01_umap_celltype") and the f-string forms. A stem carrying
# a brace is a PATTERN: its literal prefix is all that can be checked.
viz = (ROOT / "scripts" / "scenicplus_08_visualize.py").read_text()
written_stems = set()
for m in re.finditer(r'save\(\s*fig,\s*out_dir,\s*f?"([^"]+)"', viz):
    written_stems.add(m.group(1))
# also the `X if cond else Y` form, whose second branch the regex above misses
for m in re.finditer(r'else\s+f?"([0-9]{2}_[^"]+)"', viz):
    written_stems.add(m.group(1))

literal = {s for s in written_stems if "{" not in s}
patterns = {s.split("{")[0] for s in written_stems if "{" in s}
print(f"figure script writes: {len(literal)} literal, {len(patterns)} patterned")

# --- what the rule declares -------------------------------------------------
smk = (ROOT / "rules" / "report.smk").read_text()
declared = set()
for m in re.finditer(r'_fig\(\s*f?"([^"]+)"', smk):
    declared.add(m.group(1))
print(f"rule declares: {len(declared)}\n")


def matches(stem):
    """Is this declared stem something the script can write?"""
    # `f"08_eGRN_network_top{_NET_TOP}"` in the rule becomes a concrete name; on
    # the script side it is patterned, so compare on the literal prefix.
    base = stem.split("{")[0]
    if stem in literal:
        return True
    return any(base.startswith(p) or p.startswith(base.rstrip("0123456789"))
               for p in patterns)


for stem in sorted(declared):
    say(matches(stem), f"declared {stem!r} is a name the script writes",
        f"script writes: {sorted(literal)} + patterns {sorted(patterns)}")

# --- the other direction: a literal the script always writes, undeclared ----
# Only a WARNING, because two are deliberately undeclared and the file says why:
# 02_* has data-dependent names, 03_rss_per_celltype sits in a try/except that
# prints and continues, so declaring it would turn a warning into a failed run.
EXPECTED_UNDECLARED = {"03_rss_per_celltype"}
undeclared = literal - {s.split("{")[0] for s in declared} - declared
surprise = undeclared - EXPECTED_UNDECLARED
say(not surprise,
    "every literal figure the script writes is declared, or known-exempt",
    f"undeclared and unexplained: {sorted(surprise)}")

# --- the TSV side -----------------------------------------------------------
post = (ROOT / "scripts" / "scenicplus_07_postprocess_tsv.py").read_text()
tsv_written = set(re.findall(r'out_dir / "([^"]+)"', post))
tsv_declared = set(re.findall(r'os\.path\.join\(TSV, "([^"]+)"\)', smk))
print()
for name in sorted(tsv_declared):
    say(name in tsv_written, f"declared TSV {name!r} is written by the script",
        f"script writes: {sorted(tsv_written)}")

# The four unconditional ones must all be declared; the rest are gated on the
# AUC modalities being non-empty or the cell-type column being found.
UNCONDITIONAL_TSV = {"eRegulons_direct.tsv", "eRegulons_extended.tsv",
                     "eRegulons_combined.tsv", "TF_summary.tsv"}
missing = UNCONDITIONAL_TSV - tsv_declared
say(not missing, "every unconditional TSV is declared", f"missing: {sorted(missing)}")

print()
if FAIL:
    print(f"FAILED: {len(FAIL)}")
    sys.exit(1)
print("all checks passed")
