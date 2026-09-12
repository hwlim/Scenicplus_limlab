#!/bin/bash
# -----------------------------------------------------------------------------
# Keep the ORCHESTRATOR alive as its own LSF job.
#
#   cp $SCENICPLUS_PATH/scripts/scenicplus.bsub.sh .
#   $EDITOR scenicplus.bsub.sh        # the EDIT ME block below
#   bsub < scenicplus.bsub.sh
#
# WHY THIS FILE EXISTS. `scenicplus.run.sh --lsf` does NOT submit itself. It
# runs snakemake in the shell you typed it in, and snakemake submits one bsub
# per rule from there. That process has to stay alive for the WHOLE run,
# because it is the only thing that polls LSF, notices a rule finishing and
# submits the next one. Lose the shell -- ssh drops, the login node is
# rebooted, you close the laptop, someone kills long-running processes -- and
# the jobs already queued keep going while nothing schedules what comes after.
# The run does not fail; it stops, halfway, with no error anywhere.
#
# So on a cluster, submit THIS, which is a small job whose only work is to run
# the orchestrator. `tmux`/`screen`/`nohup` on a login node solve the same
# problem and are fine for a short run; a batch job is what survives the node
# going away.
#
# Read `bsub < file`, not `bsub file`: the script is fed on STDIN so LSF reads
# the #BSUB lines. `bsub scenicplus.bsub.sh` would treat it as a command and
# ignore every directive below, giving the queue's defaults instead.
#
# STATUS: this pattern is IN USE on the CCHMC cluster and is how the operator
# has been running the pipeline. The file is that working script generalised --
# the site values below are the real ones, kept as the example rather than
# replaced with placeholders, because CCHMC is the only site this has run at.
#
# THE DRIVER JOB IS NOT A COMPUTE JOB. Two cores and 16 GB is generous for a
# process whose work is polling; the real memory and cores are reserved by the
# rules it submits, per `RULE_TIERS` in rules/common.smk. What it does need is
# WALLTIME: `-W` here must cover the whole pipeline INCLUDING queue waits, not
# the longest single rule. A full run of all 21 rules measured about two hours
# of compute (RUNBOOK section 5); five hours leaves room for the queue. If this
# job hits its runlimit the same silent halt happens, so err high -- an idle
# orchestrator costs almost nothing.
# -----------------------------------------------------------------------------
#BSUB -n 2
#BSUB -W 5:00
#BSUB -M 16000
#BSUB -R "span[hosts=1]"
#BSUB -J scenicplus_snake
#BSUB -oo scenicplus_snake.%J.out
#BSUB -eo scenicplus_snake.%J.err
#
# `-oo`/`-eo` OVERWRITE. The plain `-o`/`-e` APPEND, so a second submission
# would leave one file holding two runs' output, oldest first -- the same trap
# the rule logs have (see logs/lsf/ and the report's Logs section).

set -euo pipefail

# --- EDIT ME -----------------------------------------------------------------
# These four are site- and user-specific. `install_cchmc.sh` prints the first
# three at the end of a successful install; copy them from there rather than
# guessing. The values below are the CCHMC ones, which is the only site this
# pipeline has run at.
CONDA_MODULE="anaconda3"
CONDA_PREFIX_DIR="/data/limlab/Inhee/conda-envs/scenicplus"
PIPELINE_DIR="$HOME/Scenicplus_limlab"
ANALYSIS_DIR="$PWD"           # where config/config.yaml lives

# What to ask for. `--lsf` alone already means 20 concurrent cluster jobs
# (`jobs: 20` in profiles/lsf/config.yaml); the runner's own `-j N` overrides
# it. Anything after `--` goes straight to snakemake.
RUN_ARGS=(--lsf)
# Examples, one at a time:
#   RUN_ARGS=(--lsf -j 20)                        # cap concurrent cluster jobs
#   RUN_ARGS=(--lsf -f 20)                        # redraw figures + report
#   RUN_ARGS=(--lsf -- --forcerun R21_report)     # rebuild just the report
#   RUN_ARGS=(--lsf -- --forceall)                # everything, from scratch
# --- END EDIT ME -------------------------------------------------------------

module load "$CONDA_MODULE"
eval "$(conda shell.bash hook)"
conda activate "$CONDA_PREFIX_DIR"

export SCENICPLUS_PATH="$PIPELINE_DIR"
export PATH="$SCENICPLUS_PATH/scripts:$CONDA_PREFIX_DIR/bin:$PATH"
export SCRNA_CONDA_ENV="$CONDA_PREFIX_DIR"      # scrna.tool.sh-style prefix

# NOT optional hygiene, and install_cchmc.sh says the same thing at more
# length: ~/.local/lib/pythonX.Y/site-packages comes BEFORE the env's on
# sys.path, so a copy of scenicplus / pycisTopic / pycistarget left there by an
# earlier `pip install --user` silently wins, at run time, in every step, with
# no message.
export PYTHONNOUSERSITE=1

# The env's own libstdc++, ahead of the system one. scenicplus_check.sh
# explains which import fails without it and why the failure names the wrong
# thing.
export LD_LIBRARY_PATH="$CONDA_PREFIX_DIR/lib:${LD_LIBRARY_PATH-}"

# A batch job starts wherever it was submitted from, but only if that path
# exists on the execution host. Be explicit -- the runner reads
# ./config/config.yaml and writes every stage directory relative to $PWD, so a
# wrong directory here starts a SECOND workspace rather than failing.
cd "$ANALYSIS_DIR"

echo "[scenicplus.bsub] host      : $(hostname)"
echo "[scenicplus.bsub] job       : ${LSB_JOBID:-<not under LSF>}"
echo "[scenicplus.bsub] workspace : $PWD"
echo "[scenicplus.bsub] pipeline  : $SCENICPLUS_PATH"
echo "[scenicplus.bsub] args      : ${RUN_ARGS[*]}"

# One orchestrator per workspace. Snakemake takes a lock on .snakemake/, so a
# second submission against the same directory refuses rather than racing --
# but it refuses AFTER queueing, so check `bjobs -J scenicplus_snake` first.
exec "$SCENICPLUS_PATH/scripts/scenicplus.run.sh" "${RUN_ARGS[@]}"
