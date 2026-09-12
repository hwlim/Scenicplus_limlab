# =============================================================================
# SCENIC+ workflow: orchestration only.
#
# This file resolves the config, includes the rule files, and declares the
# default target. Rules live in rules/*.smk; helpers live in rules/common.smk.
# If this file grows a `rule`, it is in the wrong place.
#
# STATUS: increment I4 of SnakemakePlan.md. All 20 steps are rules. The bash
# driver (scripts/scenicplus_run_pipeline.sh) remains a working entry point
# until I8, and both read the SAME config/config.yaml, so one workspace can be
# driven by either -- which is how they were compared.
#
# Run it through scripts/scenicplus.run.sh rather than calling snakemake
# directly: the environment preflight has to happen before a DAG is built, not
# inside one. SnakemakePlan.md, Decision 5.
# =============================================================================
import os

from snakemake.utils import min_version, validate

min_version("8.0")


# --- Config ------------------------------------------------------------------
# Relative to the WORKING directory, which is the analysis directory, not the
# pipeline directory. `scenicplus_init.sh` puts it there.
configfile: "config/config.yaml"

# The schema lives with the pipeline, so it is resolved against this file's
# directory rather than the workspace's.
validate(config, os.path.join(workflow.basedir, "schemas", "config.schema.yaml"))

# The step scripts take `--config` and read the file themselves. Absolute,
# because a cluster job's working directory is not guaranteed to be this one.
CONFIG_FILE = os.path.abspath("config/config.yaml")

# `Pipeline: "ScenicPlus"` is also a schema constraint, so reaching here means
# it matched. The bash driver checks the same line for the same reason: another
# pipeline's config must not be runnable here by accident.

include: "rules/common.smk"

# Included first, because the prefix is built from the config: pipefail, plus
# the hash seed and BLAS thread pinning that make a step's output reproducible
# across hosts. shell_prefix() in common.smk carries the measurements behind
# each part.
shell.prefix(shell_prefix())

include: "rules/prepare.smk"
include: "rules/genome.smk"
include: "rules/grn.smk"
include: "rules/report.smk"

SPECIES_INFO = species_info(config)


# --- Targets -----------------------------------------------------------------
# Decision 4 of SnakemakePlan.md, LANDED AT I6: `rule all` targets report.html,
# so an ordinary run is not finished until the run is READABLE. That makes
# looking at the output non-optional by construction rather than by discipline
# -- the failure it prevents is a green run whose outputs nobody opened, which
# this pipeline has already produced once (the 301-megapixel RSS figure was
# "successful" for a day).
#
# R19's and R20's outputs are still named even though R21 depends on both. That
# is not redundancy: `rules.R21_report.output` alone would let a future edit to
# R21's inputs silently narrow what a default run builds, and the analysis stage
# is the product -- the report is how it is read. Naming both means the target
# list states the intent rather than inheriting it.
#
# The assembly record is named for a different reason: nothing consumes it, so
# without it here the genome checks would be skipped whenever their two files
# happened to be current.
#
# An empty target list is a silent no-op, the failure mode this workflow exists
# to remove, so onstart still says so if it ever becomes one.
TARGETS = (rules.R21_report.output
           + rules.R19_postprocess_tsv.output
           + rules.R20_visualize.output
           + [stage_path("qc", "assembly.json")])


rule all:
    input:
        TARGETS


# --- Run banner, provenance, and outcome --------------------------------------
# I7: the handlers below bracket the run with a provenance bundle.
#
# `--mode start` is NOT bookkeeping that could be folded into the finish pass.
# It captures the pipeline's COMMIT while the run is starting, because the
# sibling repo's `provenance.sh` calls `git rev-parse HEAD` from its finish
# handler and therefore records the commit someone switched TO mid-run -- one
# that produced none of the outputs. An artifact whose whole job is to say what
# ran, stating something false, confidently. Captured early, reported verbatim.
#
# Failures inside a handler abort the run, and provenance must never be the
# reason a finished run is reported as failed -- so both calls are wrapped and
# a broken bundle degrades to a printed warning.
def _provenance(mode, **kw):
    import subprocess
    cmd = ["python", os.path.join(workflow.basedir, "scripts",
                                  "scenicplus_provenance.py"),
           "--mode", mode, "--workspace", ".",
           "--pipeline", workflow.basedir]
    for k, v in kw.items():
        if v:
            cmd += [f"--{k.replace('_', '-')}", str(v)]
    try:
        subprocess.run(cmd, check=True)
    except Exception as e:
        print(f"[scenicplus] provenance {mode} failed ({e}); the run itself is "
              f"unaffected and the logs are still in {stage_dir('logs')}/")


onstart:
    species = config["input"]["species"]
    print(f"[scenicplus] species={species} -> {SPECIES_INFO['assembly']}, "
          f"pycistarget name {SPECIES_INFO['pycistarget']}")
    if not TARGETS:
        print("[scenicplus] NOTHING IS TARGETED, which means a green run here "
              "would prove nothing.")
    else:
        print("[scenicplus] all 20 steps + report.html (increment I6).")
    _provenance("start")


onsuccess:
    # Name the report, not the directory. The whole point of I6 is that the run
    # is not finished until someone can read it, and a path they have to
    # assemble themselves is one they do not open.
    _provenance("finish", config=CONFIG_FILE, status="success",
                snakemake_log=log)
    print(f"[scenicplus] done. Open {os.path.abspath('report.html')}")
    print(f"[scenicplus] logs in {stage_dir('logs')}/")


onerror:
    print(f"[scenicplus] FAILED. The failing rule's own output is in "
          f"{stage_dir('logs')}/<rule>.log,")
    print("[scenicplus] which is more specific than the snakemake log above it.")
    # Said out loud because the report is now the default target, and its
    # absence after a failure looks like a second problem rather than the
    # expected consequence of the first. Snakemake will not build a target whose
    # inputs failed; the logs are what a partial run leaves behind (and, at I7,
    # the provenance bundle).
    print("[scenicplus] NO report.html: snakemake does not build a target "
          "whose inputs failed. The logs above are the record of this run.")
    # THE BUNDLE IS THE POINT ON THIS PATH. A successful run is already
    # described by report.html; this is the run that has nothing else, so the
    # bundle carries the scoped logs, names the ones with an error signature,
    # and records the commit captured at onstart.
    _provenance("finish", config=CONFIG_FILE, status="error",
                snakemake_log=log)
