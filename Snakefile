# =============================================================================
# SCENIC+ workflow: orchestration only.
#
# This file resolves the config, includes the rule files, and declares the
# default target. Rules live in rules/*.smk; helpers live in rules/common.smk.
# If this file grows a `rule`, it is in the wrong place.
#
# STATUS: increment I3 of SnakemakePlan.md. Steps 1-18 are rules; 19-20 are
# not. The bash driver (scripts/scenicplus_run_pipeline.sh) remains the working
# entry point until I8, and both read the SAME config/config.yaml, so one
# workspace can be driven by either.
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

SPECIES_INFO = species_info(config)


# --- Targets -----------------------------------------------------------------
# Decision 4 of SnakemakePlan.md: `rule all` will target report.html once I6
# lands, so that an ordinary run is not finished until the run is readable.
# Until then it targets the furthest stage that exists, which advances one
# increment at a time.
#
# I3 targets the GRN result, plus the assembly record -- which nothing consumes,
# so without naming it here the genome checks would be skipped whenever their
# two files happened to be current.
#
# An empty target list is a silent no-op, the failure mode this workflow exists
# to remove, so onstart still says so if it ever becomes one.
TARGETS = [stage_path("grn", "scplusmdata.h5mu"), stage_path("qc", "assembly.json")]


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
        print("[scenicplus] steps 1-18 (increment I3). Steps 19-20 still belong "
              "to scripts/scenicplus_run_pipeline.sh.")


onsuccess:
    print(f"[scenicplus] done. Logs in {stage_dir('logs')}/")


onerror:
    print(f"[scenicplus] FAILED. The failing rule's own output is in "
          f"{stage_dir('logs')}/<rule>.log,")
    print("[scenicplus] which is more specific than the snakemake log above it.")
