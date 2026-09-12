#!/usr/bin/env python
"""Build report.html: one page that says what this run did and what it found.

    python scenicplus_09_report.py --workspace . --config config/config.yaml \
        --out report.html

WHY THIS EXISTS. Until now a finished run left ~20 loose PNGs and 4 TSVs in
`5.analysis/`, and nothing said which of them mattered or whether the run that
produced them was healthy. This pipeline has already shipped a "successful" run
whose RSS figure was 301 megapixels and unopenable, and nobody noticed for a
day, because nothing made looking mandatory. `rule all` now targets this file,
so an ordinary run is not finished until the run is readable.

TWO CHOICES WORTH KNOWING.

**Standard library only.** No jinja2, no markdown, no pandas. jinja2 IS
guaranteed here (snakemake requires `jinja2 <4.0,>=3.0`), and pandas is a
scenicplus dependency -- so this is not about availability. It is about the
failure mode a template introduces: the sibling repo's report drivers pass
params to an .Rmd, and an undeclared param is a HARD RENDER FAILURE that no dry
run catches, discovered at the end of a cluster run after every expensive rule
has already succeeded. A function that returns a string cannot fail that way.

**Self-contained, within a budget.** PNGs are embedded as base64 data URIs so
the file survives being copied off the cluster or emailed, which is what makes
it the thing people actually look at. But the figures here are not small -- one
is 12 megapixels -- so an unbudgeted embed produces a report no browser will
open, which is the same class of failure as the 301-megapixel figure it exists
to surface. Anything over --max-embed-mb is LINKED instead, and the page says
which and why rather than silently dropping it.

A MISSING INPUT IS REPORTED, NEVER SKIPPED. Every section renders even when its
source is absent, and says what is missing. A report that quietly omits the
section whose data failed to appear is worse than no report: it looks complete.
"""
import argparse
import base64
import csv
import datetime as dt
import html
import json
import mimetypes
import os
import re
import subprocess
import sys

# Figures, in the order a reader should meet them: what the cells are, then
# what the eRegulons look like across them, then the TF-level summaries, then
# the network. Stems are matched as PREFIXES, so `08_eGRN_network_top50`
# matches `08_` without this file having to know the config value.
FIGURE_ORDER = [
    ("01_", "Cell types on the embedding",
     "The layout every other figure is read against. If this is not the "
     "reduction you chose, nothing below is either -- check input.reduction."),
    ("02_", "Per-eRegulon activity",
     "One panel per top eRegulon. Data-dependent, so the count varies."),
    ("03_", "Regulon specificity (RSS) per cell type",
     "Computed inside a try/except upstream: a cohort with sparse cell types "
     "can legitimately leave this out."),
    ("04_", "eRegulon heatmap-dotplot, direct",
     "Direct = the TF's own motif is enriched in the region."),
    ("05_", "eRegulon heatmap-dotplot, extended",
     "Extended = the motif is linked to the TF by annotation rather than "
     "found directly. More regulons, weaker evidence per regulon."),
    ("06_", "Targets per TF", ""),
    ("07_", "TF importance distribution", ""),
    ("08_", "eGRN network", "Top TFs by eRegulon count."),
]

# Tables R19 declares. Everything else it writes is gated upstream and may
# legitimately be absent.
TABLES = [
    ("eRegulons_direct.tsv", "Direct eRegulons"),
    ("eRegulons_extended.tsv", "Extended eRegulons"),
    ("eRegulons_combined.tsv", "All eRegulons"),
    ("TF_summary.tsv", "TF summary"),
]

CSS = """
:root{--fg:#1a1a1a;--mut:#666;--line:#ddd;--bg:#fff;--accent:#0b5d9e;
      --warnbg:#fff8e1;--warnln:#e6c200;--okbg:#eef7ee;--okln:#8bc34a}
*{box-sizing:border-box}
body{margin:0;font:15px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI",
     Roboto,Helvetica,Arial,sans-serif;color:var(--fg);background:var(--bg)}
.wrap{max-width:1100px;margin:0 auto;padding:2rem 1.5rem 5rem}
h1{font-size:1.8rem;margin:0 0 .2rem}
h2{font-size:1.25rem;margin:2.5rem 0 .6rem;padding-bottom:.3rem;
   border-bottom:2px solid var(--line)}
h3{font-size:1rem;margin:1.6rem 0 .4rem}
.sub{color:var(--mut);margin:0 0 1.5rem}
.note{color:var(--mut);font-size:.9rem;margin:.2rem 0 .8rem}
table{border-collapse:collapse;width:100%;font-size:.88rem;margin:.5rem 0 1rem}
th,td{border:1px solid var(--line);padding:.35rem .55rem;text-align:left;
      vertical-align:top}
th{background:#f6f6f6;font-weight:600}
td.num,th.num{text-align:right;font-variant-numeric:tabular-nums}
.scroll{overflow-x:auto;max-width:100%}
figure{margin:0 0 2rem}
figure img{max-width:100%;height:auto;border:1px solid var(--line)}
figcaption{color:var(--mut);font-size:.88rem;margin-top:.4rem}
.miss{background:var(--warnbg);border-left:4px solid var(--warnln);
      padding:.6rem .8rem;margin:.5rem 0 1.2rem;font-size:.9rem}
.ok{background:var(--okbg);border-left:4px solid var(--okln);
    padding:.6rem .8rem;margin:.5rem 0 1.2rem;font-size:.9rem}
code,pre{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:.85rem}
pre{background:#f6f6f6;border:1px solid var(--line);padding:.7rem;
    overflow-x:auto;max-height:26rem}
a{color:var(--accent)}
.toc{background:#fafafa;border:1px solid var(--line);padding:.8rem 1.2rem;
     margin:1.5rem 0}
.toc ul{margin:.3rem 0;padding-left:1.2rem}
"""


# --- small helpers -----------------------------------------------------------
def esc(x):
    return html.escape(str(x), quote=True)


def read_tsv(path, limit=None):
    """Header + rows. Returns (header, rows, total_rows) or None if unreadable.

    csv.reader rather than pandas: this must not fail on a ragged line written
    by a step that died halfway, and a report that crashes on its own input is
    worse than one that reports the input as unreadable.
    """
    try:
        with open(path, newline="") as fh:
            r = csv.reader(fh, delimiter="\t")
            header = next(r, None)
            if header is None:
                return None
            rows, total = [], 0
            for row in r:
                total += 1
                if limit is None or len(rows) < limit:
                    rows.append(row)
            return header, rows, total
    except (OSError, csv.Error):
        return None


def table_html(header, rows, total=None, limit_note=True):
    num = [i for i in range(len(header))]
    for row in rows:
        for i, cell in enumerate(row):
            if i < len(num) and num[i] is not None:
                try:
                    float(cell)
                except (TypeError, ValueError):
                    num[i] = None
    out = ['<div class="scroll"><table><thead><tr>']
    for i, h in enumerate(header):
        out.append(f'<th class="{"num" if num[i] is not None else ""}">{esc(h)}</th>')
    out.append("</tr></thead><tbody>")
    for row in rows:
        out.append("<tr>")
        for i in range(len(header)):
            cell = row[i] if i < len(row) else ""
            cls = "num" if i < len(num) and num[i] is not None else ""
            out.append(f'<td class="{cls}">{esc(cell)}</td>')
        out.append("</tr>")
    out.append("</tbody></table></div>")
    if limit_note and total is not None and total > len(rows):
        out.append(f'<p class="note">Showing {len(rows)} of {total:,} rows. '
                   f'The full table is the TSV beside this file.</p>')
    return "\n".join(out)


def missing(what, why=""):
    return (f'<div class="miss"><strong>Not in this run:</strong> {esc(what)}'
            + (f' &mdash; {esc(why)}' if why else "") + "</div>")


def embed_or_link(path, rel, max_bytes):
    """(src, note). Embeds under the cap, links over it, and SAYS WHICH.

    A silently dropped figure and a silently linked one look the same in a page
    that does not mention it, and only one of them still works after the file
    is copied somewhere else.
    """
    size = os.path.getsize(path)
    if size <= max_bytes:
        mime = mimetypes.guess_type(path)[0] or "image/png"
        with open(path, "rb") as fh:
            b64 = base64.b64encode(fh.read()).decode("ascii")
        return f"data:{mime};base64,{b64}", ""
    return rel, (f"Linked, not embedded: {size / 1e6:.1f} MB exceeds the "
                 f"{max_bytes / 1e6:.0f} MB budget. This image will not appear "
                 f"if the HTML is moved away from {esc(os.path.dirname(rel))}/.")


def run_git(repo, *args):
    try:
        out = subprocess.run(("git", "-C", repo) + args, capture_output=True,
                             text=True, timeout=10)
        return out.stdout.strip() if out.returncode == 0 else ""
    except (OSError, subprocess.SubprocessError):
        return ""


# --- LSF accounting ----------------------------------------------------------
# Same fields tests/measured_resources.tsv carries, read from this run's own
# logs so the report says what THIS run cost rather than what a past one did.
_HOST = re.compile(r"host\(s\) <(?:(\d+)\*)?([\w\-]+)>")
# re.M is load-bearing: without it `^` anchors to the START OF THE FILE, not to
# each line, so every field silently reads 0 while the host regex (unanchored)
# keeps working -- a compute table full of zeroes next to correct node names.
# Caught by tests/report_render.sh asserting a humanised wall clock.
_FIELD = re.compile(r"^\s*(Run time|Max Memory|Max Processes|Max Threads)"
                    r"\s*:\s*(\d+)", re.M)


def lsf_accounting(lsf_dir):
    """One row per rule, from the LAST job for that rule.

    `bsub -o` APPENDS and snakemake's {jobid} restarts at 0 each run, so one
    file can hold several runs' epilogues -- take the last block, not the first.
    """
    if not os.path.isdir(lsf_dir):
        return []
    per_rule = {}
    for name in sorted(os.listdir(lsf_dir)):
        if not name.endswith(".out"):
            continue
        rule = name.split(".")[0]
        try:
            text = open(os.path.join(lsf_dir, name), errors="replace").read()
        except OSError:
            continue
        rec = {"rule": rule}
        for m in _HOST.finditer(text):
            rec["slots"], rec["node"] = int(m.group(1) or 1), m.group(2)
        for m in _FIELD.finditer(text):
            rec[m.group(1)] = int(m.group(2))
        if len(rec) > 1:
            per_rule[rule] = rec
    return [per_rule[k] for k in sorted(per_rule)]


def fmt_hms(sec):
    h, rem = divmod(int(sec), 3600)
    m, s = divmod(rem, 60)
    return f"{h}:{m:02d}:{s:02d}" if h else f"{m}:{s:02d}"


# --- sections ----------------------------------------------------------------
def sec_run(cfg, cfg_path, ws, repo):
    inp = cfg.get("input", {}) or {}
    desc = run_git(repo, "describe", "--always", "--dirty", "--tags")
    head = run_git(repo, "rev-parse", "--short", "HEAD")
    branch = run_git(repo, "rev-parse", "--abbrev-ref", "HEAD")
    rows = [
        ("Generated", dt.datetime.now().astimezone().strftime("%Y-%m-%d %H:%M %Z")),
        ("Workspace", os.path.abspath(ws)),
        ("Config", cfg_path),
        ("Pipeline", os.path.abspath(repo)),
        ("Pipeline version", " ".join(x for x in (desc, f"({branch})" if branch else "") if x)
         or head or "not a git checkout"),
        ("Species", inp.get("species", "?")),
        ("Assembly (configured)", inp.get("assembly", "not set")),
        ("Cell-type column", inp.get("celltype_column", "?")),
        ("Reduction", inp.get("reduction") or "none named -- figures use a "
                                              "UMAP of eRegulon activity"),
    ]
    scope = inp.get("celltype_scope") or []
    if scope:
        rows.append(("Cell-type scope", ", ".join(map(str, scope))))
    # The two artifacts should be discoverable from each other. Deliberately
    # NOT naming this run's bundle: it is written by `onsuccess`, after this
    # page, so any bundle visible from here belongs to an EARLIER run. Naming
    # one would be the "an artifact written during a run cannot describe that
    # run" mistake in its other direction.
    prov = os.path.join(ws, "provenance")
    n_prev = len([d for d in os.listdir(prov)
                  if os.path.isdir(os.path.join(prov, d))]) \
        if os.path.isdir(prov) else 0
    rows.append(("Provenance", f"provenance/ \u2014 one bundle per run "
                               f"(commit, config, logs, LSF accounting). "
                               f"This run's is written after this page; "
                               f"{n_prev} earlier one(s) present."))
    body = "".join(f"<tr><th>{esc(k)}</th><td>{esc(v)}</td></tr>" for k, v in rows)
    return f"<h2 id=run>Run</h2><table><tbody>{body}</tbody></table>"


def sec_assembly(path):
    h = "<h2 id=genome>Genome</h2>"
    if not os.path.exists(path):
        return h + missing("QC/assembly.json",
                           "the genome checks did not run or did not finish")
    try:
        rec = json.load(open(path))
    except (OSError, ValueError) as e:
        return h + missing("QC/assembly.json", f"unreadable: {e}")
    unann = rec.get("peak_chromosomes_unannotated") or []
    rows = [
        ("Detected assembly", rec.get("assembly_detected", "?")),
        ("Configured assembly", rec.get("assembly_configured", "?")),
        ("Chromosomes", rec.get("n_chromosomes", "?")),
        ("Annotation rows", f'{rec.get("n_annotation_rows", 0):,}'),
        ("Peak chromosomes", rec.get("peak_chromosomes", "?")),
        ("...also in the annotation", rec.get("peak_chromosomes_annotated", "?")),
    ]
    body = "".join(f"<tr><th>{esc(k)}</th><td>{esc(v)}</td></tr>" for k, v in rows)
    out = h + f"<table><tbody>{body}</tbody></table>"
    det, conf = rec.get("assembly_detected"), rec.get("assembly_configured")
    if det and conf and str(det) != str(conf):
        out += (f'<div class="miss"><strong>Assembly mismatch:</strong> the '
                f'chromsizes say {esc(det)}, the config says {esc(conf)}. '
                f'Coordinates will not line up and nothing downstream will say '
                f'so.</div>')
    if unann:
        out += (f'<p class="note">{len(unann)} peak chromosome(s) absent from '
                f'the annotation: {esc(", ".join(map(str, unann[:12])))}'
                f'{" ..." if len(unann) > 12 else ""}. Peaks there get no '
                f'search space.</p>')
    src = rec.get("source", {})
    if src:
        body = "".join(
            f"<tr><th>{esc(k)}</th><td><code>{esc(v.get('path',''))}</code><br>"
            f"<code>{esc((v.get('sha256') or '')[:16])}…</code></td></tr>"
            for k, v in src.items())
        out += f"<h3>Inputs</h3><table><tbody>{body}</tbody></table>"
    return out


def sec_tables(tsv_dir, head_rows):
    out = ["<h2 id=tables>eRegulon tables</h2>"]
    if not os.path.isdir(tsv_dir):
        return "\n".join(out) + missing(tsv_dir, "R19 did not produce its tables")
    counts = []
    for fname, label in TABLES:
        got = read_tsv(os.path.join(tsv_dir, fname), limit=0)
        counts.append((label, fname, f"{got[2]:,}" if got else "absent"))
    out.append(table_html(["Table", "File", "Rows"],
                          [[a, b, c] for a, b, c in counts], limit_note=False))
    for fname, label in TABLES:
        got = read_tsv(os.path.join(tsv_dir, fname), limit=head_rows)
        out.append(f"<h3>{esc(label)} &mdash; <code>{esc(fname)}</code></h3>")
        if not got:
            out.append(missing(fname))
            continue
        header, rows, total = got
        out.append(table_html(header, rows, total))
    return "\n".join(out)


def sec_figures(plots_dir, ws, max_bytes):
    out = ["<h2 id=figures>Figures</h2>"]
    if not os.path.isdir(plots_dir):
        return "\n".join(out) + missing(plots_dir, "R20 produced no figures")
    pngs = sorted(f for f in os.listdir(plots_dir) if f.endswith(".png"))
    if not pngs:
        return "\n".join(out) + missing("any .png in " + plots_dir)
    used = set()
    for prefix, title, blurb in FIGURE_ORDER:
        group = [f for f in pngs if f.startswith(prefix)]
        out.append(f"<h3>{esc(title)}</h3>")
        if blurb:
            out.append(f'<p class="note">{esc(blurb)}</p>')
        if not group:
            out.append(missing(f"figure {prefix}*",
                               "not produced by this run"))
            continue
        for f in group:
            used.add(f)
            full = os.path.join(plots_dir, f)
            rel = os.path.relpath(full, ws)
            src, note = embed_or_link(full, rel, max_bytes)
            pdf = os.path.splitext(rel)[0] + ".pdf"
            has_pdf = os.path.exists(os.path.join(ws, pdf))
            cap = [f"<code>{esc(rel)}</code>"]
            if has_pdf:
                cap.append(f'<a href="{esc(pdf)}">PDF</a>')
            if note:
                cap.append(note)
            out.append(f'<figure><img src="{src}" alt="{esc(f)}">'
                       f'<figcaption>{" &middot; ".join(cap)}</figcaption>'
                       f"</figure>")
    extra = [f for f in pngs if f not in used]
    if extra:
        out.append("<h3>Other figures</h3>")
        out.append('<p class="note">Produced by this run but not in the '
                   'expected set. Listed rather than dropped.</p>')
        for f in extra:
            full = os.path.join(plots_dir, f)
            rel = os.path.relpath(full, ws)
            src, note = embed_or_link(full, rel, max_bytes)
            out.append(f'<figure><img src="{src}" alt="{esc(f)}">'
                       f'<figcaption><code>{esc(rel)}</code>'
                       f'{" &middot; " + note if note else ""}</figcaption>'
                       f"</figure>")
    return "\n".join(out)


def sec_compute(lsf_dir):
    out = ["<h2 id=compute>Compute</h2>"]
    rows = lsf_accounting(lsf_dir)
    if not rows:
        out.append('<p class="note">No LSF accounting found in '
                   f'<code>{esc(lsf_dir)}</code>. Expected for a local run; '
                   "for a cluster run it means the epilogues were not written "
                   "where the profile says.</p>")
        return "\n".join(out)
    header = ["Rule", "Node", "Slots", "Wall clock", "Peak MB", "Processes", "Threads"]
    body, tot = [], 0
    for r in rows:
        sec = r.get("Run time", 0)
        tot += sec
        body.append([r["rule"], r.get("node", "?"), r.get("slots", "?"),
                     fmt_hms(sec), f'{r.get("Max Memory", 0):,}',
                     r.get("Max Processes", "?"), r.get("Max Threads", "?")])
    out.append(table_html(header, body, limit_note=False))
    peak = max((r.get("Max Memory", 0) for r in rows), default=0)
    nodes = {r.get("node") for r in rows if r.get("node")}
    out.append(f'<p class="note">{len(rows)} rule{"" if len(rows) == 1 else "s"}, '
               f'{fmt_hms(tot)} of summed '
               f'job time across {len(nodes)} node(s); heaviest single job '
               f'{peak:,} MB. Summed time is not wall clock &mdash; independent '
               f'rules run at once. Peak memory is LSF\'s figure for the whole '
               f'process TREE, while <code>-M</code> is enforced per process, '
               f'so a total above the reservation is not a breach.</p>')
    if len(nodes) > 1:
        out.append('<p class="note">More than one node: floating-point '
                   'differences between rules are expected and are '
                   'node-class-determined, not stochastic.</p>')
    return "\n".join(out)


def sec_config(cfg_path):
    out = ["<h2 id=config>Configuration</h2>",
           '<p class="note">The config as it was read. This is the record of '
           "what the run was ASKED to do; it is not proof that every key was "
           "used.</p>"]
    try:
        out.append(f"<pre>{esc(open(cfg_path).read())}</pre>")
    except OSError as e:
        out.append(missing(cfg_path, str(e)))
    return "\n".join(out)


def sec_logs(logs_dir, ws, self_log=None):
    """The run's logs, and which of them are suspiciously empty.

    THE REPORT'S OWN LOG IS ALWAYS EMPTY HERE, and that is not a finding. The
    rule pipes through `tee`, so the file exists from the moment the job starts
    and receives this script's output only as it is printed -- which happens
    AFTER the page has been built and this directory read. Reported as an empty
    log it is a false alarm on every single run, and a check that cries wolf
    every time is one people learn to skip.

    So it is separated rather than suppressed: the row is annotated, the alarm
    excludes it, and the page says why. Suppressing it silently would be the
    other error -- a reader who counts the rules and finds one log unaccounted
    for deserves the explanation in the page, not in the source.

    This is the same shape as the sibling repo's "an artifact written DURING a
    run cannot describe that run completely", where a report's own missing log
    was once reported as evidence the run had been local.
    """
    out = ["<h2 id=logs>Logs</h2>"]
    if not os.path.isdir(logs_dir):
        return "\n".join(out) + missing(logs_dir)
    files = sorted(f for f in os.listdir(logs_dir) if f.endswith(".log"))
    if not files:
        return "\n".join(out) + missing("any .log in " + logs_dir)
    self_name = os.path.basename(self_log) if self_log else None
    rows, empty, self_seen = [], [], False
    for f in files:
        st = os.stat(os.path.join(logs_dir, f))
        is_self = (f == self_name)
        note = ""
        if is_self:
            self_seen = True
            note = "this rule's own log \u2014 written after this page"
        elif st.st_size == 0:
            empty.append(f)
        rows.append([f, f"{st.st_size:,}",
                     dt.datetime.fromtimestamp(st.st_mtime).strftime("%m-%d %H:%M"),
                     note])
    out.append(table_html(["Log", "Bytes", "Modified", "Note"], rows,
                          limit_note=False))
    if self_seen:
        out.append(f'<p class="note"><code>{esc(self_name)}</code> reads as '
                   f'0 bytes above and is NOT counted as an empty log. The rule '
                   f'pipes through <code>tee</code>, so the file exists from the '
                   f'start of the job and receives this script\u2019s output only '
                   f'after the page has been written \u2014 an artifact written '
                   f'during a run cannot describe that run completely. Read it on '
                   f'disk afterwards for the real contents.</p>')
    if empty:
        out.append(f'<div class="miss"><strong>Empty log(s):</strong> '
                   f'{esc(", ".join(empty))}. A rule that wrote nothing at all '
                   f'is worth a look even when the run was green.</div>')
    return "\n".join(out)


def build(ws, cfg_path, repo, head_rows, max_bytes, self_log=None):
    try:
        import yaml
        cfg = yaml.safe_load(open(cfg_path)) or {}
    except Exception:
        cfg = {}
    j = lambda *p: os.path.join(ws, *p)
    parts = [
        sec_run(cfg, cfg_path, ws, repo),
        sec_assembly(j("QC", "assembly.json")),
        sec_tables(j("5.analysis", "tsv"), head_rows),
        sec_figures(j("5.analysis", "plots"), ws, max_bytes),
        sec_compute(j("logs", "lsf")),
        sec_logs(j("logs"), ws, self_log),
        sec_config(cfg_path),
    ]
    toc = """<div class="toc"><strong>On this page</strong><ul>
<li><a href="#run">Run</a></li><li><a href="#genome">Genome</a></li>
<li><a href="#tables">eRegulon tables</a></li><li><a href="#figures">Figures</a></li>
<li><a href="#compute">Compute</a></li><li><a href="#logs">Logs</a></li>
<li><a href="#config">Configuration</a></li></ul></div>"""
    name = os.path.basename(os.path.abspath(ws))
    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>SCENIC+ &mdash; {esc(name)}</title>
<style>{CSS}</style></head><body><div class="wrap">
<h1>SCENIC+ run report</h1>
<p class="sub">{esc(name)} &middot; generated
{esc(dt.datetime.now().strftime("%Y-%m-%d %H:%M"))}</p>
{toc}
{"".join(parts)}
</div></body></html>
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workspace", default=".", help="analysis directory")
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--pipeline", default=os.environ.get("SCENICPLUS_PATH", ""),
                    help="install root, for the git version line")
    ap.add_argument("--head-rows", type=int, default=15,
                    help="rows shown per table [%(default)s]")
    ap.add_argument("--max-embed-mb", type=float, default=4.0,
                    help="embed a PNG up to this size; link it above "
                         "[%(default)s]")
    # Passed by the rule as {log}. Not guessed from the rule name: the report
    # must not depend on a naming convention it cannot see, and being wrong here
    # means either a false empty-log alarm every run or a real one suppressed.
    ap.add_argument("--self-log", default=None, metavar="PATH",
                    help="this rule's own log. It is always empty while the "
                         "page is being built, so it is annotated rather than "
                         "reported as an empty log.")
    a = ap.parse_args()
    html_text = build(a.workspace, a.config, a.pipeline or ".",
                      a.head_rows, a.max_embed_mb * 1e6, a.self_log)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)
    with open(a.out, "w") as fh:
        fh.write(html_text)
    size = os.path.getsize(a.out)
    print(f"[report] wrote {a.out} ({size / 1e6:.2f} MB)")
    # An empty-ish page means every section reported its data missing, which is
    # a green rule that produced nothing worth opening -- exactly the outcome
    # this report exists to make visible.
    if size < 4000:
        print("[report] WARNING: the page is nearly empty. Every section found "
              "its input missing.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
