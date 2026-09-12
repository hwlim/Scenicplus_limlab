#!/usr/bin/env python
"""The provenance bundle: what ran, from what code, on what, and how it ended.

Two modes, and the split between them is the whole point.

    scenicplus_provenance.py --mode start  --workspace . --pipeline $SCENICPLUS_PATH
    scenicplus_provenance.py --mode finish --workspace . --pipeline $SCENICPLUS_PATH \
        --config config/config.yaml --status success

WHY `--mode start` EXISTS AT ALL, and it is the reason this file is not one
function called from `onerror`.

The sibling repo's `provenance.sh` runs `git rev-parse HEAD` inside its
success/error handler -- at the END of the run. If anyone checks out another
branch while the run is in flight, and the pipeline clone is a shared install so
somebody eventually does, the bundle records the commit they switched TO, which
produced none of the outputs. That is worse than having no provenance: the
artifact whose entire job is to say what ran states something false, and states
it confidently.

So the code identity is captured at ONSTART, written to `logs/.run_meta.json`,
and the finish pass reports what was recorded rather than what it can see. If
the marker is missing the bundle says so instead of substituting a fresh
`rev-parse`, because a plausible wrong answer is the failure being avoided.

THE BUNDLE MUST BE MOST USEFUL WHEN THE RUN FAILED. A successful run is already
described by `report.html`. The bundle earns its keep on the run that stopped
at rule 9 of 21, so:

  * logs are copied for a FAILED run exactly as for a successful one, and the
    manifest names the rules whose logs carry a traceback;
  * the cap SKIPS rather than truncates -- a half log is a trap, because the
    interesting part of a python traceback is at the END;
  * logs are SCOPED to this run by `logs/.run_started`, since `bsub -o` appends
    and snakemake's `{jobid}` restarts at 0, so one file can hold several runs'
    epilogues. Bundling an older run's log as this run's would misattribute a
    failure, which is worse than omitting it.
"""
import argparse
import datetime as dt
import importlib.util
import json
import os
import platform
import re
import shutil
import subprocess
import sys

MARKER = ".run_started"
META = ".run_meta.json"


def _load_report_module():
    """Reuse the report's LSF parser instead of copying it.

    Imported BY PATH from the same directory rather than by name: `scripts/` is
    not a package and is not on sys.path when snakemake invokes this. Two copies
    of an accounting parser would drift, and the drift would be silent -- both
    would keep producing plausible numbers.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, "scenicplus_09_report.py")
    spec = importlib.util.spec_from_file_location("_scp_report", path)
    if spec is None or spec.loader is None:
        return None
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)
        return mod
    except Exception:
        return None


def git(repo, *args):
    try:
        r = subprocess.run(("git", "-C", repo) + args, capture_output=True,
                           text=True, timeout=10)
        return r.stdout.strip() if r.returncode == 0 else ""
    except (OSError, subprocess.SubprocessError):
        return ""


def mode_start(ws, pipeline):
    logs = os.path.join(ws, "logs")
    os.makedirs(logs, exist_ok=True)
    # The marker's MTIME is the scope boundary; its contents are incidental.
    with open(os.path.join(logs, MARKER), "w") as fh:
        fh.write(dt.datetime.now().isoformat() + "\n")
    dirty = git(pipeline, "status", "--porcelain")
    meta = {
        "started_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "pipeline_dir": os.path.abspath(pipeline),
        "commit": git(pipeline, "rev-parse", "HEAD") or "unknown",
        "branch": git(pipeline, "rev-parse", "--abbrev-ref", "HEAD") or "unknown",
        "describe": git(pipeline, "describe", "--always", "--dirty", "--tags") or "unknown",
        # Recorded as a COUNT, not the diff: the diff can be large and is not
        # this file's job. What matters downstream is whether the commit above
        # fully describes the code that ran.
        "dirty_files": len([l for l in dirty.splitlines() if l.strip()]),
        "host": platform.node(),
        "python": platform.python_version(),
        "user": os.environ.get("USER", "?"),
        "lsf_job": os.environ.get("LSB_JOBID", ""),
    }
    with open(os.path.join(logs, META), "w") as fh:
        json.dump(meta, fh, indent=2)
    return meta


def scoped_logs(logs_dir, marker_mtime):
    """Files belonging to THIS run: at or after the marker. The marker and the
    meta file are excluded -- they are bookkeeping, not run output.

    The "at or after" test itself comes from the REPORT module, which is the one
    definition of it. Two copies drifted once already in spirit: the report had
    no scoping at all while this file did, so a `--forcerun R20` produced a
    bundle scoped correctly next to a report listing the whole previous run.
    Sharing the predicate makes that divergence impossible rather than unlikely.
    """
    out = []
    if not os.path.isdir(logs_dir):
        return out
    rep = _load_report_module()
    keep = (rep.in_this_run if rep and hasattr(rep, "in_this_run")
            else lambda p, since: since is None
            or os.path.getmtime(p) >= since - 1)
    for root, _, files in os.walk(logs_dir):
        for f in sorted(files):
            if f in (MARKER, META):
                continue
            p = os.path.join(root, f)
            try:
                if keep(p, marker_mtime):
                    out.append(p)
            except OSError:
                continue
    return out


# A log carrying one of these is worth naming in the manifest. Deliberately
# short and boring: the aim is "look here first", not classification.
_TRACE = re.compile(r"^(Traceback \(most recent call last\)|Error in |"
                    r"\w*Error:|Exception:|Killed|Segmentation fault|"
                    r"TERM_(MEMLIMIT|RUNLIMIT|OWNER))", re.M)


def failing_logs(paths):
    hits = []
    for p in paths:
        try:
            with open(p, errors="replace") as fh:
                if _TRACE.search(fh.read()):
                    hits.append(p)
        except OSError:
            continue
    return hits


def mode_finish(ws, pipeline, cfg_path, status, cap_kb, snakemake_log):
    logs_dir = os.path.join(ws, "logs")
    meta_path = os.path.join(logs_dir, META)
    marker = os.path.join(logs_dir, MARKER)

    meta, meta_note = {}, ""
    if os.path.exists(meta_path):
        try:
            meta = json.load(open(meta_path))
        except (OSError, ValueError) as e:
            meta_note = f"logs/{META} is unreadable ({e})"
    else:
        # NOT filled in with a fresh rev-parse. See the module docstring: a
        # plausible wrong commit is the failure this design exists to avoid.
        meta_note = (f"logs/{META} absent -- this run did not go through "
                     f"`--mode start`, so the code identity below is UNKNOWN "
                     f"rather than guessed at finish time.")

    stamp = dt.datetime.now().strftime("%Y%m%dT%H%M%S")
    out = os.path.join(ws, "provenance", f"{stamp}_{status}")
    os.makedirs(out, exist_ok=True)

    mtime = os.path.getmtime(marker) if os.path.exists(marker) else None
    logs = scoped_logs(logs_dir, mtime)
    failed = failing_logs(logs)

    # --- logs, capped by SKIPPING ---------------------------------------------
    total = sum(os.path.getsize(p) for p in logs if os.path.exists(p))
    logs_note = ""
    if total > cap_kb * 1024:
        logs_note = (f"{len(logs)} log(s) totalling {total/1e6:.1f} MB exceed the "
                     f"{cap_kb/1024:.0f} MB cap, so NONE were copied. They are "
                     f"still in {logs_dir}/. Truncating was rejected: the "
                     f"interesting part of a traceback is at the end.")
    else:
        dest = os.path.join(out, "logs")
        for p in logs:
            rel = os.path.relpath(p, logs_dir)
            d = os.path.join(dest, rel)
            os.makedirs(os.path.dirname(d), exist_ok=True)
            try:
                shutil.copy2(p, d)
            except OSError:
                pass
        logs_note = f"{len(logs)} log(s), {total/1e6:.2f} MB"

    # --- the config AS READ ----------------------------------------------------
    if cfg_path and os.path.exists(cfg_path):
        try:
            shutil.copy2(cfg_path, os.path.join(out, "config.used.yaml"))
        except OSError:
            pass

    # --- side artifacts --------------------------------------------------------
    # report.html is copied so the bundle is self-contained: a run's evidence
    # and its readable summary travel together, and the report embeds its own
    # figures, so the pair survives being moved off the cluster.
    #
    # Its presence or ABSENCE is recorded either way. On a failed run snakemake
    # does not build it, and a bundle that simply lacks the file leaves a reader
    # guessing between "not built" and "lost". The manifest says which.
    copied = {}
    for src, name in ((os.path.join(ws, "QC", "assembly.json"), "assembly.json"),
                      (os.path.join(ws, "report.html"), "report.html"),
                      (snakemake_log or "", "snakemake.log")):
        if src and os.path.exists(src):
            try:
                shutil.copy2(src, os.path.join(out, name))
                copied[name] = os.path.getsize(src)
            except OSError:
                pass

    # --- lsf_jobs.tsv ----------------------------------------------------------
    rep = _load_report_module()
    rows = rep.lsf_accounting(os.path.join(logs_dir, "lsf")) if rep else []
    with open(os.path.join(out, "lsf_jobs.tsv"), "w") as fh:
        fh.write("rule\tnode\tslots\truntime_s\tmax_mem_mb\tmax_processes\t"
                 "max_threads\n")
        for r in rows:
            fh.write("\t".join(str(r.get(k, "")) for k in
                               ("rule", "node", "slots", "Run time", "Max Memory",
                                "Max Processes", "Max Threads")) + "\n")

    # --- the manifest ----------------------------------------------------------
    lines = [
        "# SCENIC+ provenance bundle",
        "",
        f"status:        {status}",
        f"finished_at:   {dt.datetime.now().astimezone().isoformat(timespec='seconds')}",
        f"workspace:     {os.path.abspath(ws)}",
        "",
        "# Code identity, CAPTURED AT ONSTART, not now. A checkout during the run",
        "# would make a finish-time rev-parse name a commit that produced nothing.",
        f"started_at:    {meta.get('started_at', 'unknown')}",
        f"pipeline_dir:  {meta.get('pipeline_dir', 'unknown')}",
        f"commit:        {meta.get('commit', 'unknown')}",
        f"branch:        {meta.get('branch', 'unknown')}",
        f"describe:      {meta.get('describe', 'unknown')}",
        f"dirty_files:   {meta.get('dirty_files', 'unknown')}",
        f"host:          {meta.get('host', 'unknown')}",
        f"lsf_job:       {meta.get('lsf_job', '') or '(not under LSF)'}",
        f"user:          {meta.get('user', 'unknown')}",
        f"python:        {meta.get('python', 'unknown')}",
        "",
        f"logs:          {logs_note}",
        f"lsf_jobs:      {len(rows)} job(s) with accounting",
        f"report.html:   " + (
            f"included, {copied['report.html']/1e6:.2f} MB"
            if "report.html" in copied else
            "NOT INCLUDED -- " + ("snakemake does not build a target whose "
                                  "inputs failed, so a partial run has none"
                                  if status != "success" else
                                  "the run succeeded but no report.html was "
                                  "found, which is worth looking into")),
        f"assembly.json: " + ("included" if "assembly.json" in copied
                              else "absent (the genome checks did not run)"),
    ]
    if meta.get("dirty_files"):
        lines += ["", "# WARNING: the install had uncommitted changes when this run",
                  "# started, so `commit` does not fully describe the code that ran."]
    if meta_note:
        lines += ["", f"# NOTE: {meta_note}"]
    if failed:
        lines += ["", "# Logs carrying an error signature. Read these first:"]
        lines += [f"#   {os.path.relpath(p, ws)}" for p in failed]
    elif status != "success":
        lines += ["", "# The run failed but NO log carries an error signature.",
                  "# Check the snakemake log: the failure may be in scheduling",
                  "# rather than in a rule (a missing input, or a killed job whose",
                  "# epilogue landed elsewhere)."]
    with open(os.path.join(out, "manifest.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return out, failed


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", choices=("start", "finish"), required=True)
    ap.add_argument("--workspace", default=".")
    ap.add_argument("--pipeline", default=os.environ.get("SCENICPLUS_PATH", "."))
    ap.add_argument("--config", default=None)
    ap.add_argument("--status", default="success")
    ap.add_argument("--snakemake-log", default=None)
    ap.add_argument("--cap-kb", type=int, default=51200,
                    help="skip the logs entirely above this total [%(default)s]")
    a = ap.parse_args()

    if a.mode == "start":
        m = mode_start(a.workspace, a.pipeline)
        note = f", {m['dirty_files']} uncommitted file(s)" if m["dirty_files"] else ""
        print(f"[provenance] {m['describe']} on {m['branch']}{note}")
        return 0

    out, failed = mode_finish(a.workspace, a.pipeline, a.config, a.status,
                              a.cap_kb, a.snakemake_log)
    print(f"[provenance] bundle: {os.path.abspath(out)}")
    if failed:
        print(f"[provenance] {len(failed)} log(s) carry an error signature; "
              f"manifest.txt names them")
    return 0


if __name__ == "__main__":
    sys.exit(main())
