#!/usr/bin/env python
"""Every key the report reads from `assembly.json` must be one R07 writes.

    python3 tests/record_keys.py

Text on both sides, no data and no imports of either file.

WHY THIS EXISTS, and it is a defect that shipped. The report's Genome section
read `assembly_detected` and `assembly_configured`. `scenicplus_genome_
prepare.py` has never written either -- it writes `assembly`, `chr1_bp` and
`chr1_matches`. On a real report both rows rendered `?`, and the "Assembly
mismatch" panel beneath them could never have fired: it compared two keys that
are always absent.

NOTHING CAUGHT IT, because the FIXTURE INVENTED THE SAME KEYS. Consumer and test
were written from the same imagination, so they agreed with each other and with
nothing the producer emits. A fixture derived from the code it checks is not a
check; this file compares the consumer against the PRODUCER instead, which is
the only pair that can disagree usefully.

The same shape appeared twice more this week -- the t-SNE reading an unprefixed
MuData column, and `output_names.py` seeing only literal figure stems -- so the
general lesson is worth stating: WHEN A CONSUMER AND ITS FIXTURE ARE WRITTEN
TOGETHER, THEY TEST NEITHER. Anchor one of them on the producer.
"""
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PRODUCER = os.path.join(ROOT, "scripts", "scenicplus_genome_prepare.py")
CONSUMER = os.path.join(ROOT, "scripts", "scenicplus_09_report.py")

fails = []


def ok(label, cond, detail=""):
    print(f"{'ok  ' if cond else 'FAIL'}  {label}")
    if not cond:
        if detail:
            print(f"        {detail}")
        fails.append(label)


prod = open(PRODUCER).read()
cons = open(CONSUMER).read()

# The producer's record literal. Anchored on `record = {` through its closing
# brace at the same indent, and the match count is asserted -- a renamed
# variable must fail loudly here rather than silently check nothing.
blocks = re.findall(r"^    record = \{.*?^    \}", prod, re.S | re.M)
ok("found R07's record literal", len(blocks) == 1, f"{len(blocks)} matches")
written = set(re.findall(r'^\s+"([a-z0-9_]+)":', blocks[0], re.M)) if blocks else set()
# `source` nests one level; its inner keys are read through a loop, not by name.
ok(f"it writes a plausible number of keys ({len(written)})", len(written) >= 8,
   sorted(written))

# What the report asks for. `rec` is the parsed record in sec_assembly.
read = set(re.findall(r'rec\.get\("([a-z0-9_]+)"', cons))
ok(f"the report reads {len(read)} key(s) from it", bool(read))

missing = sorted(read - written)
ok("every key the report reads is one R07 writes", not missing,
   f"report reads {missing}, which R07 never writes -- these render as '?' and "
   f"any comparison between them is unreachable. R07 writes: {sorted(written)}")

# --- and the VALUES, not only the names --------------------------------------
# The key-name check above shipped, and the same defect came straight back one
# level down: `report_render.sh` got the NAMES from this producer and then wrote
# `"chr1_matches": "hg38"`, a value the producer cannot emit. KNOWN_CHR1 maps a
# length to every accepted spelling of one assembly, joined by "/", so the real
# value is "hg38/GRCh38". The report compared it to the config with `==`, the
# comparison was false for every recognised assembly, and the Genome section
# rendered no panel at all -- while this fixture reported it CONFIRMED.
#
# So: any chr1_matches literal in the test suite must be one the producer's own
# table can produce. Text on both sides, same as everything else here.
known = re.search(r"^KNOWN_CHR1 = \{(.*?)\}", prod, re.S | re.M)
ok("found R07's KNOWN_CHR1 table", known is not None)
spellings = set(re.findall(r'"([^"]+)"', known.group(1))) if known else set()
ok(f"it maps lengths to {len(spellings)} assembly name(s)", len(spellings) >= 3,
   sorted(spellings))

RENDER = os.path.join(ROOT, "tests", "report_render.sh")
# Comment lines dropped first. The prose above that fixture QUOTES the bad value
# while explaining it, and scanning them made this check report the comment as
# the defect -- correct-looking output for the wrong reason, which is the thing
# this file exists to complain about.
render_code = "\n".join(ln for ln in open(RENDER).read().splitlines()
                        if not ln.lstrip().startswith("#"))
fixtures = re.findall(r'"chr1_matches":\s*"([^"]*)"', render_code)
ok(f"report_render.sh pins {len(fixtures)} chr1_matches value(s)", bool(fixtures))
bogus = sorted(set(fixtures) - spellings)
ok("every chr1_matches fixture value is one KNOWN_CHR1 emits", not bogus,
   f"invented {bogus}; the producer can only write {sorted(spellings)}. A "
   f"fixture value that cannot occur makes the assertion beside it meaningless")

# The consumer must not compare that value with `==` -- "hg38" is never equal to
# "hg38/GRCh38". This is a text check because the shape, not the result, is the
# defect: an equality here is wrong no matter what it is fed.
bad_eq = re.search(r'str\(configured\)\s*==\s*str\(detected\)', cons)
ok("the report does not compare configured to detected with ==",
   bad_eq is None,
   "chr1_matches carries \"/\"-separated SYNONYMS, so equality is false for "
   "every recognised assembly; use assembly_agrees()")
ok("...it uses an alias-aware predicate instead",
   "def assembly_agrees(" in cons and "assembly_agrees(configured, detected)" in cons)

# The other direction is informational: a key written and not shown is a choice,
# not a defect, but worth seeing so the choice stays deliberate.
unshown = sorted(written - read)
print(f"\nwritten but not shown in the report ({len(unshown)}): "
      f"{', '.join(unshown) or 'none'}")
print("  Not a failure -- a record may carry more than a page should. Listed so "
      "dropping something useful stays a decision rather than an oversight.")

print()
if fails:
    print(f"FAILED: {len(fails)}")
    sys.exit(1)
print("all record-key checks passed")
