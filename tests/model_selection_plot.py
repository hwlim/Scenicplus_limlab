#!/usr/bin/env python
"""Gate for the topic-model selection figures (step 04).

    python tests/model_selection_plot.py

Needs pycisTopic and matplotlib -- i.e. the scenicplus env -- but NO data, no
cluster and no LDA fitting. `evaluate_models` touches exactly two attributes of
a model, `.n_topic` and `.metrics`, so the sweep can be stood in for and the
REAL library function driven end to end.

WHY THAT MATTERS. The thing being tested is not arithmetic, it is an
interaction with somebody else's plotting code: that `plot=True` leaves the
figure open rather than closing it, that the per-metric figures carry titles to
name files from, and that what lands on disk is a readable image. Every one of
those is a property of pycisTopic, not of this repo, and a stub would test the
stub.

THE EMPTY-ARTIFACT TRAP IS THE POINT. `savefig` produces a file whether or not
anything was drawn, exactly as R's `dev.off()` writes a valid ZERO-page PDF --
a trap that bit the sibling repo three separate times. So these checks open
what was written and assert it has pages and pixels, never that the path
exists.
"""
import importlib.util
import os
import re
import sys

os.environ.setdefault("MPLCONFIGDIR", os.environ.get("TMPDIR", "/tmp") + "/mpl")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STEP = os.path.join(ROOT, "scripts", "scenicplus_04_topic_modeling.py")

fails = []


def ok(label, cond, detail=""):
    print(f"{'ok  ' if cond else 'FAIL'}  {label}")
    if not cond:
        if detail:
            print(f"        {detail}")
        fails.append(label)


spec = importlib.util.spec_from_file_location("_step04", STEP)
step = importlib.util.module_from_spec(spec)
spec.loader.exec_module(step)
ok("step 04 imports without running the sweep", hasattr(step, "save_open_figures"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from pycisTopic.lda_models import evaluate_models


class FakeModel:
    """Only what evaluate_models reads: `.n_topic` and `.metrics`.

    Verified by reading the function, not assumed -- it indexes
    `models[i].metrics.loc["Metric", <name>]` for the four metrics and sorts on
    `.n_topic`. If pycisTopic starts reading more, this construction raises
    rather than silently drawing something meaningless.
    """

    def __init__(self, n, arun, cao, mimno, ll):
        self.n_topic = n
        self.metrics = pd.DataFrame(
            {"Arun_2010": [arun], "Cao_Juan_2009": [cao],
             "Mimno_2011": [mimno], "loglikelihood": [ll]},
            index=["Metric"])


# A sweep with a real optimum in the middle, so the chosen-topic marker has
# somewhere meaningful to land.
MODELS = [FakeModel(*row) for row in [
    (5,  0.90, 0.80, -6.0, -900.0),
    (10, 0.55, 0.50, -4.2, -840.0),
    (15, 0.30, 0.35, -3.1, -805.0),
    (20, 0.42, 0.48, -3.8, -815.0),
    (25, 0.61, 0.66, -4.9, -830.0),
]]

import tempfile
WORK = tempfile.mkdtemp()


def png_size(path):
    """(width, height) from the PNG header. A file that exists and decodes to
    0x0 is the failure mode being guarded, so read the IHDR rather than stat."""
    import struct
    with open(path, "rb") as fh:
        head = fh.read(24)
    if head[:8] != b"\x89PNG\r\n\x1a\n":
        return None
    return struct.unpack(">II", head[16:24])


def pdf_pages(path):
    """Page count, without a PDF library: count the /Type /Page objects."""
    data = open(path, "rb").read()
    return len(re.findall(rb"/Type\s*/Page[^s]", data))


# --- 1. the premise: plot=True leaves figures open ---------------------------
plt.close("all")
best = evaluate_models(MODELS, select_model=None, return_model=True, plot=True)
open_now = plt.get_fignums()
ok("evaluate_models(plot=True) leaves the figure OPEN to be saved",
   len(open_now) == 1, f"figures open: {len(open_now)}")
ok("...and still selects a model", getattr(best, "n_topic", None) in
   [m.n_topic for m in MODELS], f"got {getattr(best, 'n_topic', None)}")

# The other half, because the whole design rests on the difference.
plt.close("all")
evaluate_models(MODELS, select_model=None, return_model=True, plot=False)
ok("evaluate_models(plot=False) CLOSES it, which is why plot=True is passed",
   plt.get_fignums() == [], f"figures open: {plt.get_fignums()}")

# --- 2. the combined figure lands, in both formats, with pixels --------------
plt.close("all")
evaluate_models(MODELS, select_model=None, return_model=True, plot=True)
out1 = os.path.join(WORK, "one")
written = step.save_open_figures(out1, "topic_model_selection")
names = sorted(os.path.basename(p) for p in written)
ok("one plot per file, both formats",
   names == ["topic_model_selection.pdf", "topic_model_selection.png"], names)
sz = png_size(os.path.join(out1, "topic_model_selection.png"))
ok("the PNG is a real image with non-zero dimensions", sz and sz[0] > 200 and sz[1] > 100,
   f"IHDR says {sz}")
ok("the PDF has at least one page (not a zero-page file)",
   pdf_pages(os.path.join(out1, "topic_model_selection.pdf")) >= 1)
ok("no figure is left open afterwards", plt.get_fignums() == [])

# --- 3. per-metric figures are named from their TITLES, not their order ------
plt.close("all")
evaluate_models(MODELS, select_model=None, return_model=True, plot=True,
                plot_metrics=True)
out2 = os.path.join(WORK, "metrics")
written = step.save_open_figures(out2, "topic_model_selection")
stems = sorted({os.path.splitext(os.path.basename(p))[0] for p in written})
ok("the combined figure keeps the bare stem",
   "topic_model_selection" in stems, stems)
# The real titles are "Arun_2010 - Minimize" and so on, so the slug carries the
# OPTIMISATION DIRECTION as well as the metric -- `mimno_2011_maximize` tells a
# reader which way is better without opening the figure. Asserted, because it
# is the reason to take names from titles rather than from a hardcoded list.
for metric, direction in (("arun_2010", "minimize"),
                          ("cao_juan_2009", "minimize"),
                          ("mimno_2011", "maximize"),
                          ("loglikelihood", "maximize")):
    hit = [x for x in stems if metric in x]
    ok(f"...and {metric} is named from its title, with its direction",
       len(hit) == 1 and direction in hit[0], f"{hit or 'no file'} in {stems}")
pngs = [str(p) for p in written if str(p).endswith(".png")]
ok(f"every per-metric PNG has pixels too ({len(pngs)} of them)",
   len(pngs) == 5 and all((png_size(p) or (0, 0))[0] > 100 for p in pngs),
   f"sizes: {[png_size(p) for p in pngs]}")

# --- 3b. the library's own multi-page PDF is written too ---------------------
# Kept alongside the per-figure output, never instead of it: a multi-page PDF
# cannot be embedded in the report, so on its own the section reports the
# figures as never drawn -- verified by rendering a workspace holding only that
# file. Here we check the extra exists and really has every panel.
import subprocess, tempfile as _tf
out3 = os.path.join(WORK, "allpages")
os.makedirs(out3, exist_ok=True)
plt.close("all")
from matplotlib.backends.backend_pdf import PdfPages
allp = os.path.join(out3, "topic_model_selection_all_pages.pdf")
evaluate_models(MODELS, select_model=None, return_model=True, plot=True,
                plot_metrics=True, save=allp)
ok("the multi-page PDF is written", os.path.exists(allp))
ok("...and holds every panel, not just the first",
   pdf_pages(allp) == 5, f"pages: {pdf_pages(allp)}")
ok("...while save= does NOT suppress the individual capture",
   len(plt.get_fignums()) == 5, f"figures open: {len(plt.get_fignums())}")
plt.close("all")

# --- 4. the declared name and the written name are the SAME string ----------
# `topic_model_selection` is a literal in two files: the stem this script passes
# to save_open_figures, and the output rules/prepare.smk declares. Drift between
# them is a MissingOutputException at the end of the most expensive rule in the
# pipeline -- which is precisely what tests/output_names.py exists to prevent for
# report.smk, and R04 is outside its scope.
smk = open(os.path.join(ROOT, "rules", "prepare.smk")).read()
step_src = open(STEP).read()
declared = set(re.findall(r'stage_path\("qc", f?"([\w.{}]+)"', smk))
passed = set(re.findall(r'save_open_figures\(\s*[\w.]+,\s*"([\w]+)"', step_src))
ok("the script passes exactly one figure stem", len(passed) == 1, f"{passed}")
stem = next(iter(passed), None)
ok("prepare.smk declares that same stem",
   any(d.startswith(stem or "\0") for d in declared),
   f"script writes {stem!r}; prepare.smk declares {sorted(declared)}")
ok("...and declares BOTH formats, per the output contract",
   any("{ext}" in d or d.endswith(".pdf") for d in declared)
   and any("{ext}" in d or d.endswith(".png") for d in declared),
   f"declared: {sorted(declared)}")

print()
if fails:
    print(f"FAILED: {len(fails)}")
    sys.exit(1)
print("all model-selection plot checks passed")
