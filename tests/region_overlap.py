#!/usr/bin/env python
"""Gate for the eRegulon region-overlap heatmap (step 08).

    python tests/region_overlap.py

Needs pandas, seaborn and scipy -- the scenicplus env -- but no MuData, no
cluster and no GRN. The function takes the eRegulon metadata frame directly, so
it can be driven with a constructed one.

WHY A CONSTRUCTED FRAME IS ENOUGH HERE, unlike the model-selection gate. That
one had to drive somebody else's plotting code, so a stub would have tested the
stub. This figure is built in-repo from a groupby and a pairwise loop, so what
needs checking is the ARITHMETIC and the degenerate cases -- and for those, a
frame with known-by-hand overlaps is stronger evidence than real data, because
the right answer is known independently.

THE JACCARD VALUES ARE CHECKED AGAINST HAND-COMPUTED ONES. A heatmap of the
wrong matrix is still a plausible heatmap: every cell in [0,1], a clean diagonal,
believable blocks. Nothing about looking at it would reveal a transposed index or
an intersection counted against the wrong denominator. So the TSV written beside
the figure is read back and compared to numbers worked out here.
"""
import importlib.util
import os
import sys

os.environ.setdefault("MPLCONFIGDIR", os.environ.get("TMPDIR", "/tmp") + "/mpl")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STEP = os.path.join(ROOT, "scripts", "scenicplus_08_visualize.py")

fails = []


def ok(label, cond, detail=""):
    print(f"{'ok  ' if cond else 'FAIL'}  {label}")
    if not cond:
        if detail:
            print(f"        {detail}")
        fails.append(label)


import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

spec = importlib.util.spec_from_file_location("_step08", STEP)
step = importlib.util.module_from_spec(spec)
spec.loader.exec_module(step)
ok("step 08 imports", hasattr(step, "region_overlap"))

import tempfile
WORK = Path(tempfile.mkdtemp())


def meta(rows):
    """rows: {eRegulon_name: [regions]} -> the triplet frame R19 reads."""
    out = []
    for name, regions in rows.items():
        for r in regions:
            out.append({"eRegulon_name": name, "TF": name.split("_")[0],
                        "Region": r, "Gene": "G1"})
    return pd.DataFrame(out)


# Overlaps worked out by hand:
#   A={r1..r4}  B={r3..r6}  -> |A&B|=2 |A|B|=6  -> 1/3
#   A vs C={r1..r4}         -> identical        -> 1
#   A vs D={r9,r10}         -> disjoint         -> 0
FRAME = meta({
    "TFA_direct": ["r1", "r2", "r3", "r4"],
    "TFB_direct": ["r3", "r4", "r5", "r6"],
    "TFC_direct": ["r1", "r2", "r3", "r4"],
    "TFD_direct": ["r9", "r10"],
})

out = WORK / "plots"
plt.close("all")
step.region_overlap(FRAME, out, "direct", 40)

png = out / "09_region_overlap_direct.png"
pdf = out / "09_region_overlap_direct.pdf"
tsv = out / "09_region_overlap_direct.tsv"
ok("the figure is written in both formats", png.exists() and pdf.exists())
ok("...and the data sheet beside it", tsv.exists())
ok("...in the SAME directory, not R19's tsv/",
   not (out.parent / "tsv").exists(),
   "step 08 must not write into the rule that owns tsv/")

if tsv.exists():
    m = pd.read_csv(tsv, sep="\t", index_col=0)
    ok("the matrix is square and named by eRegulon",
       m.shape == (4, 4) and "TFA_direct" in m.index, f"{m.shape}, {list(m.index)}")
    ok("the diagonal is 1", all(abs(m.loc[k, k] - 1.0) < 1e-9 for k in m.index))
    ok("it is symmetric",
       all(abs(m.loc[a, b] - m.loc[b, a]) < 1e-9 for a in m.index for b in m.index))
    # The three values that matter, each independently known.
    ok("half-overlapping sets give Jaccard 1/3, hand-computed",
       abs(m.loc["TFA_direct", "TFB_direct"] - 1 / 3) < 1e-9,
       f'got {m.loc["TFA_direct", "TFB_direct"]}')
    ok("identical sets give 1", abs(m.loc["TFA_direct", "TFC_direct"] - 1.0) < 1e-9,
       f'got {m.loc["TFA_direct", "TFC_direct"]}')
    ok("disjoint sets give 0", abs(m.loc["TFA_direct", "TFD_direct"]) < 1e-9,
       f'got {m.loc["TFA_direct", "TFD_direct"]}')

# --- the top-N cut, and that the figure SAYS what it cut ---------------------
plt.close("all")
big = meta({f"TF{i}_direct": [f"r{i}", f"r{i+1}", f"r{i+2}"] for i in range(30)})
out2 = WORK / "capped"
step.region_overlap(big, out2, "direct", 8)
m2 = pd.read_csv(out2 / "09_region_overlap_direct.tsv", sep="\t", index_col=0)
ok("top_n caps the matrix", m2.shape == (8, 8), f"{m2.shape}")

# --- the degenerate cases must SKIP, not crash and not emit an empty figure --
for label, frame, why in (
        ("an empty frame", pd.DataFrame(), "no metadata at all"),
        ("a frame with no Region column",
         pd.DataFrame({"eRegulon_name": ["a"], "Gene": ["g"]}), "wrong shape"),
        ("a single eRegulon", meta({"TFA_direct": ["r1", "r2"]}),
         "nothing to compare against")):
    d = WORK / label.replace(" ", "_")
    plt.close("all")
    try:
        step.region_overlap(frame, d, "direct", 40)
        crashed = False
    except Exception as e:
        crashed = repr(e)
    ok(f"{label} skips cleanly ({why})", crashed is False, str(crashed))
    ok(f"...and writes NO figure for {label}",
       not (d / "09_region_overlap_direct.png").exists())

# --- extended goes to its own file, not over direct's ------------------------
plt.close("all")
out3 = WORK / "both"
step.region_overlap(FRAME, out3, "direct", 40)
step.region_overlap(FRAME, out3, "extended", 40)
names = sorted(p.name for p in out3.iterdir())
ok("direct and extended write DIFFERENT stems",
   "09_region_overlap_direct.png" in names
   and "10_region_overlap_extended.png" in names, names)

# --- the figure has pixels, not just a filename ------------------------------
import struct


def png_size(p):
    with open(p, "rb") as fh:
        head = fh.read(24)
    return struct.unpack(">II", head[16:24]) if head[:8] == b"\x89PNG\r\n\x1a\n" else None


sz = png_size(png)
ok("the PNG is a real image", sz and sz[0] > 200 and sz[1] > 200, f"IHDR {sz}")
ok("no figure is left open", plt.get_fignums() == [], f"{plt.get_fignums()}")

print()
if fails:
    print(f"FAILED: {len(fails)}")
    sys.exit(1)
print("all region-overlap checks passed")
