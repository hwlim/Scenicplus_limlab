#!/usr/bin/env python
"""Gate for the eRegulon-activity t-SNE (step 08).

    python tests/eregulon_tsne.py

Needs scanpy, anndata and mudata -- the scenicplus env -- but no GRN and no
cluster. A MuData is constructed with the two AUC modalities the function reads.

WHAT NEEDS PROVING HERE IS NOT THE t-SNE. sklearn's embedding is not this
repo's code and testing it would test sklearn. What is this repo's, and what can
silently be wrong, is everything around it:

  * that the SEED is honoured, so two runs of the same data give the same
    picture. This figure is the only stochastic one in the pipeline, and the
    whole reproducibility argument for the track -- PYTHONHASHSEED, BLAS thread
    pinning, drivers agreeing to 1e-9 -- would be undone by one unseeded plot
    that differs every run while every table matches.
  * that a DIFFERENT seed actually changes it, because "reproducible" is
    worthless if the parameter does nothing at all.
  * that perplexity is clamped rather than allowed to raise. sklearn REQUIRES
    perplexity < n_samples and errors otherwise, and `input.celltype_scope` can
    legitimately leave few cells.
  * that the two feature spaces produce DIFFERENT files, not one overwriting
    the other.

THE SEED PAIR IS THE POINT. Either check alone is satisfiable the wrong way: a
function that ignores its input entirely passes "same seed, same output", and
one that ignores the seed passes nothing but would still look fine in a figure.
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
import numpy as np
import pandas as pd
import anndata as ad
import mudata as mu
from pathlib import Path

spec = importlib.util.spec_from_file_location("_step08", STEP)
step = importlib.util.module_from_spec(spec)
spec.loader.exec_module(step)
ok("step 08 imports", hasattr(step, "eregulon_tsne"))

import tempfile
WORK = Path(tempfile.mkdtemp())

# Three cell types with genuinely different regulon activity, so the embedding
# has something to find and a degenerate result would be visible.
rng = np.random.default_rng(0)
N_PER, N_ER = 120, 24
types, blocks = [], []
for t in range(3):
    m = rng.normal(0.1, 0.02, size=(N_PER, N_ER))
    m[:, t * 8:(t + 1) * 8] += 0.6          # this type's own regulons fire
    blocks.append(m)
    types += [f"celltype_{t}"] * N_PER
X = np.vstack(blocks)
cells = [f"c{i}" for i in range(X.shape[0])]
obs = pd.DataFrame({"cell_type": types}, index=cells)


def make_md(with_region=True):
    mods = {"direct_gene_based_AUC": ad.AnnData(
        X.copy(), obs=obs.copy(),
        var=pd.DataFrame(index=[f"g_eReg{i}" for i in range(N_ER)]))}
    if with_region:
        mods["direct_region_based_AUC"] = ad.AnnData(
            X.copy() + rng.normal(0, 0.01, X.shape), obs=obs.copy(),
            var=pd.DataFrame(index=[f"r_eReg{i}" for i in range(N_ER)]))
    m = mu.MuData(mods)
    m.obs = obs.copy()
    return m


MD = make_md()


def png_bytes(p):
    return open(p, "rb").read()


# --- 1. it renders, per feature space, into its own file ---------------------
out = WORK / "a"
plt.close("all")
step.eregulon_tsne(MD, "cell_type", out, "gene", 555)
step.eregulon_tsne(MD, "cell_type", out, "region", 555)
names = sorted(p.name for p in out.iterdir())
ok("gene and region write DIFFERENT stems, in both formats",
   names == ["11_tsne_eRegulon_gene_based.pdf", "11_tsne_eRegulon_gene_based.png",
             "12_tsne_eRegulon_region_based.pdf", "12_tsne_eRegulon_region_based.png"],
   names)
ok("no figure is left open", plt.get_fignums() == [], f"{plt.get_fignums()}")

# --- 2. THE SEED PAIR --------------------------------------------------------
out_same = WORK / "same"
plt.close("all")
step.eregulon_tsne(MD, "cell_type", out_same, "gene", 555)
a1 = png_bytes(out / "11_tsne_eRegulon_gene_based.png")
a2 = png_bytes(out_same / "11_tsne_eRegulon_gene_based.png")
ok("the SAME seed gives a byte-identical figure", a1 == a2,
   f"{len(a1)} vs {len(a2)} bytes")

out_diff = WORK / "diff"
plt.close("all")
step.eregulon_tsne(MD, "cell_type", out_diff, "gene", 999)
a3 = png_bytes(out_diff / "11_tsne_eRegulon_gene_based.png")
ok("a DIFFERENT seed gives a different figure, so the seed is not ignored",
   a1 != a3, "identical output for seeds 555 and 999")

# --- 3. perplexity is clamped, not allowed to raise --------------------------
small = make_md()[:12].copy() if hasattr(make_md(), "__getitem__") else None
try:
    tiny_obs = obs.iloc[:12].copy()
    tiny = mu.MuData({"direct_gene_based_AUC": ad.AnnData(
        X[:12].copy(), obs=tiny_obs,
        var=pd.DataFrame(index=[f"g_eReg{i}" for i in range(N_ER)]))})
    tiny.obs = tiny_obs
    out_t = WORK / "tiny"
    plt.close("all")
    step.eregulon_tsne(tiny, "cell_type", out_t, "gene", 555)
    made = (out_t / "11_tsne_eRegulon_gene_based.png").exists()
    ok("12 cells still renders -- perplexity clamped below n_samples", made)
except Exception as e:
    ok("12 cells still renders -- perplexity clamped below n_samples", False, repr(e))

# --- 4. the degenerate cases SKIP, and write nothing -------------------------
for label, md, why in (
        ("a MuData with no AUC modality",
         mu.MuData({"something_else": ad.AnnData(X.copy(), obs=obs.copy())}),
         "nothing to embed"),
        ("a missing cell-type column", MD, "cannot colour it")):
    d = WORK / label.replace(" ", "_").replace("-", "_")
    plt.close("all")
    col = "cell_type" if "no AUC" in label else "not_a_column"
    try:
        step.eregulon_tsne(md, col, d, "gene", 555)
        crashed = False
    except Exception as e:
        crashed = repr(e)
    ok(f"{label} skips cleanly ({why})", crashed is False, str(crashed))
    ok(f"...and writes nothing for {label}",
       not (d / "11_tsne_eRegulon_gene_based.png").exists())

# --- 5. the caption says it is a DIFFERENT layout ----------------------------
# Every other figure in the report is drawn on the Seurat reduction. A reader
# who has scrolled past six of those will read this as the same coordinates
# unless the figure says otherwise, and `input.reduction` already cost this
# project a run over that confusion.
import re
smk = open(os.path.join(ROOT, "scripts", "scenicplus_09_report.py")).read()
ok("the report's caption for 11_ names the layout as its own",
   re.search(r'\("11_",.*?activity', smk, re.S) is not None,
   "FIGURE_ORDER's 11_ entry must say the layout is eRegulon activity")

print()
if fails:
    print(f"FAILED: {len(fails)}")
    sys.exit(1)
print("all eRegulon-tSNE checks passed")
