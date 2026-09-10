# =============================================================================
# SCENIC+ workflow: orchestration only.
#
# This file resolves the config, includes the rule files, and declares the
# default target. Rules live in rules/*.smk; helpers live in rules/common.smk.
# If this file grows a `rule`, it is in the wrong place.
#
# STATUS: increment I0 of SnakemakePlan.md. The skeleton and the config
# contract exist; no step rules yet. The bash driver
# (scripts/scenicplus_run_pipeline.sh) remains the working entry point until
# I8, and both read the SAME config/config.yaml, so a workspace can be driven
# by either.
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

# `set -o pipefail` and nothing else. Rule bodies end `2>&1 | tee {log}`, and
# WITHOUT pipefail the exit status of that pipeline is tee's, which is 0 -- so a
# failing script reports success and the run continues on a missing output. That
# is 22 rules' worth of silent failure in the sibling pipeline, found only by
# deliberately breaking one.
#
# Not `-e`: rule bodies chain with `&&` and would change meaning. Not `-u`: site
# `module` functions dereference unset variables.
shell.prefix("set -o pipefail; ")

# `Pipeline: "ScenicPlus"` is also a schema constraint, so reaching here means
# it matched. The bash driver checks the same line for the same reason: another
# pipeline's config must not be runnable here by accident.

include: "rules/common.smk"
include: "rules/prepare.smk"
include: "rules/genome.smk"

SPECIES_INFO = species_info(config)


# --- Targets -----------------------------------------------------------------
# Decision 4 of SnakemakePlan.md: `rule all` will target report.html once I6
# lands, so that an ordinary run is not finished until the run is readable.
# Until then it targets the furthest stage that exists, which advances one
# increment at a time.
#
# I1: the region sets, where R01-R05 end and I3's GRN rules will pick up.
# I2: the checked genome pair, when the config supplies one. Conditional because
# nothing consumes it until I3 -- rules/genome.smk explains why that is a
# warning now and becomes a hard requirement then.
#
# An empty target list is a silent no-op -- the failure mode this workflow exists
# to remove -- so onstart still says so if it ever becomes one.
TARGETS = [stage_path("cistopic", "region_sets")]
if GENOME_SUPPLIED:
    TARGETS += GENOME_FILES


rule all:
    input:
        TARGETS


# --- Run banner and outcome --------------------------------------------------
# I7 replaces the two handlers below with a provenance bundle: config.used,
# git description, the logs, lsf_jobs.tsv and assembly.json. Until then they
# say where to look, which is the part people actually need when a run stops.
onstart:
    species = config["input"]["species"]
    print(f"[scenicplus] species={species} -> {SPECIES_INFO['assembly']}, "
          f"pycistarget name {SPECIES_INFO['pycistarget']}")
    if not TARGETS:
        print("[scenicplus] NOTHING IS TARGETED, which means a green run here "
              "would prove nothing.")
    else:
        print(f"[scenicplus] steps 1-5 only (increment I1). "
              f"Steps 6-20 still belong to scripts/scenicplus_run_pipeline.sh.")


onsuccess:
    print(f"[scenicplus] done. Logs in {stage_dir('logs')}/")


onerror:
    print(f"[scenicplus] FAILED. The failing rule's own output is in "
          f"{stage_dir('logs')}/<rule>.log,")
    print("[scenicplus] which is more specific than the snakemake log above it.")
