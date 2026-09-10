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

# `Pipeline: "ScenicPlus"` is also a schema constraint, so reaching here means
# it matched. The bash driver checks the same line for the same reason: another
# pipeline's config must not be runnable here by accident.

include: "rules/common.smk"

SPECIES_INFO = species_info(config)


# --- Targets -----------------------------------------------------------------
# Decision 4 of SnakemakePlan.md: `rule all` will target report.html once I6
# lands, so that an ordinary run is not finished until the run is readable.
# Until then it targets the analysis stage, and until I1 lands there is nothing
# to target at all.
#
# An empty target list is a silent no-op, which is the failure mode this whole
# workflow exists to remove, so onstart says it out loud instead.
TARGETS = []


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
        print("[scenicplus] NOTHING IS TARGETED. This is the I0 skeleton: the "
              "config contract and the\n"
              "[scenicplus] layout exist, no step rules do. Use "
              "scripts/scenicplus_run_pipeline.sh for a\n"
              "[scenicplus] real run until I1 lands.")


onsuccess:
    print(f"[scenicplus] done. Logs in {stage_dir('logs')}/")


onerror:
    print(f"[scenicplus] FAILED. The failing rule's own output is in "
          f"{stage_dir('logs')}/<rule>.log,")
    print("[scenicplus] which is more specific than the snakemake log above it.")
