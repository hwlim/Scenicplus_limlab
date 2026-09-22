#!/usr/bin/env python
"""Gate for write_bed()'s minimum-size guard (issue #1).

Run:  python tests/test_region_sets.py

WHAT THIS EXISTS FOR. One degenerate region set aborts the whole
motif-enrichment stage. Measured on a real run: Otsu binarization produced a
29-region topic out of 40, none of whose regions overlapped the cisTarget
database; pycistarget raised ValueError, joblib's loky backend tore down the
worker pool, and BOTH R09_cistarget and R10_dem died -- discarding 86 healthy
sets after many hours of upstream stages.

The guard is a MINIMUM REGION COUNT, and the cases below are written from the
numbers that run produced rather than from round figures: 29 (the set that
crashed), 49 (a set that mapped 34 database regions and survived), and 4409
(the smallest healthy topic). A threshold of 500 sits between them with room on
both sides.

Two things the cases deliberately pin:

  IT IS A PROXY. What fails downstream is "zero DATABASE regions mapped"; a
  count only correlates with that. The 49-region set DID survive, and is
  dropped here anyway because 34 database regions is useless for motif
  enrichment -- so this is a behaviour change beyond fixing the crash, not just
  a crash fix, and the case says so.

  NO SILENT SKIPS. write_bed had two quiet return paths before this: an
  all-unparseable region list, and an empty one. A thin region set must not be
  indistinguishable from a healthy one in the log, so every refusal prints.
"""
import importlib.util
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location(
    "region_sets", ROOT / "scripts" / "scenicplus_05_region_sets.py")
rs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rs)

failed = []
passed = 0


def chk(label, cond):
    global passed
    if cond:
        passed += 1
        print(f"ok    {label}")
    else:
        failed.append(label)
        print(f"FAIL  {label}")


def regions(n):
    return [f"chr1:{i * 1000}-{i * 1000 + 500}" for i in range(n)]


tmp = Path(tempfile.mkdtemp())

# --- the guard, at the sizes the failing run actually produced ---------------
CASES = [
    ("a healthy topic is written (4409 = the smallest real one)", regions(4409), 500, True),
    ("the 29-region topic that aborted the stage is refused",     regions(29),   500, False),
    ("the 49-region topic is refused too, though it SURVIVED --"
     " 34 database regions is not a usable enrichment",           regions(49),   500, False),
    ("exactly at the threshold is kept, not dropped",             regions(500),  500, True),
    ("one below the threshold is dropped",                        regions(499),  500, False),
    ("min_regions = 0 disables the guard entirely",               regions(29),   0,   True),
]
for i, (label, regs, min_regions, want) in enumerate(CASES):
    path = tmp / f"case{i}.bed"
    got = rs.write_bed(regs, path, min_regions)
    chk(label, got == want and path.exists() == want)

# --- the return value is what the caller counts skips with -------------------
chk("write_bed returns True when it writes",
    rs.write_bed(regions(600), tmp / "ret_true.bed", 500) is True)
chk("write_bed returns False when it skips",
    rs.write_bed(regions(5), tmp / "ret_false.bed", 500) is False)

# --- the two paths that used to be silent ------------------------------------
p = tmp / "all_malformed.bed"
chk("an all-unparseable region list is refused and no file appears",
    rs.write_bed(["not-a-region", "also bad"], p, 0) is False and not p.exists())
p = tmp / "empty.bed"
chk("an empty region list is refused", rs.write_bed([], p, 0) is False and not p.exists())
p = tmp / "mixed.bed"
chk("a partly-parseable list still writes the parseable rows",
    rs.write_bed(["chr1:0-100", "garbage", "chr2:0-100"], p, 0) is True
    and len(p.read_text().strip().splitlines()) == 2)

# --- NEGATIVE CONTROL --------------------------------------------------------
# Without the guard every case above would pass by writing everything, so the
# suite would be green against the code it exists to police. Call the same
# function with the guard disabled and the crashing set must come back written.
p = tmp / "control.bed"
chk("negative control: with min_regions = 0 the 29-region set IS written,"
    " so the cases above are testing the guard and not the parser",
    rs.write_bed(regions(29), p, 0) is True and p.exists())

print()
if failed:
    print(f"FAILED {len(failed)} of {passed + len(failed)}:")
    for f in failed:
        print(f"  {f}")
    sys.exit(1)
print(f"all {passed} region-set guard checks passed")
