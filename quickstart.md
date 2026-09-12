# Quickstart: SCENIC+ end to end

For a lab member running this pipeline, either on the CCHMC HPC or on their own
machine. Follow it in order. [RUNBOOK.md](RUNBOOK.md) carries the reasoning
behind every step, and each entry below names the section to open when something
does not go as described.

Written in plain ASCII on purpose. This file gets copied into terminals and
editors, and typographic dashes and section signs do not always survive that.

---

## What you need first

One Seurat object holding both modalities:

- RNA counts in the `RNA` assay
- ATAC peak counts in the `peaks` assay
- a categorical cell-type column

The pipeline will not start without the last of those. Fragment paths inside the
object may be dangling and that is fine, because no step reads them, so copying
the .rds alone is enough.

Note the name of the reduction you want every figure drawn on, spelled exactly
as `Reductions(obj)` prints it, for example `wnn.umap`. Step E needs it.

RUNBOOK section 0 covers building an object that qualifies, including what to do
when yours has no cell-type labels yet.

## Which path you are on

| | CCHMC HPC | your own machine |
|---|---|---|
| A. environment | install once, on a compute node | install once |
| B. reference data | already on `/data/limlab`, nothing to download | download 45.7 GB yourself |
| F. running | `scenicplus.run.sh --lsf`, through a small `run.sh` | `scenicplus.run.sh -j <cores>` |
| what to expect | the tested path | steps 1 to 3 are comfortable; steps 9 onward need the 45.7 GB of databases and far more memory than a laptop has |

**This has run end to end at CCHMC and nowhere else.** The environment recipe
also builds on a WSL2 workstation, but building is not running. The launchers,
the queue and module names, the memory unit LSF reads, the login-versus-compute
glibc split that decides which wheels pip picks, and the assumption that
`pip.conf` may carry `user = true` are all shaped by that site. On your own
machine, expect to adjust before anything runs. RUNBOOK section 1 says what to
look at first, and section 6 lists what is already known to differ.

---

## A. Install the environment, once, about 20 minutes

RUNBOOK section 1.

    git clone -b main https://github.com/hwlim/Scenicplus_limlab.git
    cd Scenicplus_limlab
    ./install_cchmc.sh /path/to/scenicplus_env
    /path/to/scenicplus_env/bin/scenicplus --help      # must print usage

If the lab already keeps a built environment on shared storage, use that rather
than building a second copy. Ask first.

Install on the node class you will run jobs on, because pip picks wheels for the
glibc of the host it runs on, and login and compute nodes can differ. If that
last file is missing, pip installed outside the prefix, which is almost always
the user site: `SCP_FROM=2 ./install_cchmc.sh <prefix>` repairs it in about two
minutes without redoing the conda layer.

## B. Reference data, once per assembly

RUNBOOK section 2, and section 2b for mouse.

Nothing before step 9 needs the cisTarget databases, so on your own machine you
can start the download and continue with the steps below meanwhile.

**On the CCHMC HPC these already exist.** As of 2026-09-09:

    /data/limlab/Resource/Scenicplus_db/cistarget/
    /data/limlab/Resource/Scenicplus_db/genome_files/

The genome files there were generated with
`$SCENICPLUS_PATH/scripts/scenicplus_make_genome_files.R`, which needs an R
environment carrying EnsDb and BSgenome. That is the `scRNA_LimLab_Snake`
environment, not the SCENIC+ one. To build a pair for another assembly:

    Rscript $SCENICPLUS_PATH/scripts/scenicplus_make_genome_files.R \
        --species mouse --out-dir /path/to/genome     # or --species human

Supplying these by hand is the normal path, not a workaround: step 7 cannot
produce chromosome sizes for anyone, because it queries an NCBI database that
has been retired. Read the `chr1 = ...` line the script prints. For mouse,
195,471,971 is mm10 and 195,154,279 is GRCm39, so that single number tells you
which assembly you are about to analyse.

## C. Set the environment variables

RUNBOOK section 3.

    export SCENICPLUS_PATH=/path/to/Scenicplus_limlab
    export PATH=$SCENICPLUS_PATH/scripts:/path/to/scenicplus_env/bin:$PATH
    export PYTHONNOUSERSITE=1
    export LD_LIBRARY_PATH=/path/to/scenicplus_env/lib:$LD_LIBRARY_PATH

`PYTHONNOUSERSITE` is not optional on a shared cluster. Your `~/.local`
packages otherwise precede the environment on the import path and can shadow it
silently.

## D. Create the analysis directory

RUNBOOK section 3.

    mkdir -p /path/to/analysis && cd /path/to/analysis
    scenicplus_init.sh                       # writes config/config.yaml

Everything from here runs from inside that directory. You edit only its
`config/config.yaml`; the pipeline itself stays where `SCENICPLUS_PATH` points.

## E. Fill in `config/config.yaml`

RUNBOOK section 3, and section 2b for mouse.

| key | value |
|---|---|
| `input.seurat_rds` | absolute path to the object |
| `input.celltype_column` | the categorical column to group by |
| `input.reduction` | the reduction every figure is drawn on, e.g. `wnn.umap` |
| `input.species` | `hsapiens` or `mmusculus` |
| `input.genome_annotation`, `input.chromsizes` | the two files from B |
| `input.ctx_db`, `input.dem_db`, `input.motif_annotations` | from B |
| `resources.n_cpu` | cores the job will have |

Leaving `input.reduction` empty is a choice rather than a default: a layout gets
computed instead, with no batch correction, which for a multi-sample integrated
object is not the picture you clustered on. A name the object does not have
fails at step 1, in the first minute, with the available names listed.

Species and assembly must agree with the genome files from B. They are not
checked against each other.

Everything else can stay at its shipped value for a first run. If runtime
matters, read RUNBOOK section 5b before starting rather than after: it lists the
four parameters that dominate it, with a cheaper value for each and what that
value costs in resolution rather than in seconds. The largest savings are in how
many topic models step 4 fits and how wide the region-to-gene threshold grid is.
That section is also where the settings for a deliberately cheap plumbing test
are spelled out, for when you want to prove the pipeline runs rather than to
believe its numbers.

## F. Plan, then run

RUNBOOK section 4.

    scenicplus.run.sh -n                       # prints the plan, runs nothing

**On your own machine:**

    scenicplus.run.sh -j 8                     # runs here, on 8 cores

**On the CCHMC HPC**, `--lsf` submits each rule as its own job, sized from a
measurement rather than from the heaviest step:

    scenicplus.run.sh --lsf -j 20              # at most 20 cluster jobs at once

Keep the invocation in a small `run.sh` in the analysis directory rather than in
shell history. That puts the settings a run used beside that run's outputs, and
makes a re-run one command. A skeleton to adapt:

```bash
#!/usr/bin/env bash
set -euo pipefail
module purge
module load anaconda3 R/4.4.0-R0            # the preflight runs Rscript here
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /path/to/scenicplus_env

export SCENICPLUS_PATH=/path/to/Scenicplus_limlab
export PATH=$SCENICPLUS_PATH/scripts:$PATH
export PYTHONNOUSERSITE=1

scenicplus.run.sh --lsf -j 20 "$@"          # -n, -f 7, -- --forceall, ...
```

Then `./run.sh` submits, and `./run.sh -f 7` forwards the flag.

**You no longer set the reproducibility variables yourself.** The workflow pins
`PYTHONHASHSEED` and the four `*_NUM_THREADS` for every rule, deriving the
thread count from `resources.n_cpu`, so the two cannot disagree. That matters:
without them two runs of the same data on the same cluster do not agree, and
both effects are measured. Python randomises string hashing per process and
SCENIC+ builds lists from sets of names, so a permuted order permutes ties and a
top-N cut keeps different rows; and the maths libraries take their thread count
from the host's cores, with reduction order following it. A 48-core and a
64-core node gave correlations differing in the last bit.

**You no longer set memory or cores either.** Each rule carries its own
reservation, measured from a real run: 4 GB for the report, 96 GB for topic
modelling, 256 GB for cistarget. RUNBOOK section 5b has the table.

**Submit from the node class the environment was built for.** The preflight runs
on the host you submit from, before anything is queued, so it checks the shell
you are typing in rather than the one the jobs will get. Two consequences, both
of which have cost time here: R must be loaded to submit and not only to run,
and a login node whose glibc or modules differ from the compute nodes fails the
check with nothing queued. Start an interactive session on the right class.

### The bash driver, if you need it

`scenicplus_run_pipeline.sh` still works and prints a notice saying it is
obsolete. It is kept as a fallback and for reproducing an older run. It submits
all twenty steps as ONE job sized for the heaviest, so the plotting step holds
cistarget's cores and memory for hours, and it produces no report and no
provenance bundle. Converting a habit:

| bash driver | workflow |
|---|---|
| `scenicplus_run_pipeline.sh --dry-run` | `scenicplus.run.sh -n` |
| `--from 7` | `-f 7` |
| `--force` | `-- --forceall` |
| `scenicplus_run_workstation.sh` | `scenicplus.run.sh -j <cores>` |
| `scenicplus_run_lsf.sh` | `scenicplus.run.sh --lsf -j 20` |

## G. Watch, then read

RUNBOOK sections 4 and 5.

    bjobs                                      # HPC only
    tail -f logs/R04_topic_modeling.log        # one log per rule
    open report.html                           # the point of the run

**Start with `report.html`.** `rule all` targets it, so a run is not finished
until it is readable. It carries the run's provenance, the genome check, the
eRegulon tables, every figure, this run's LSF accounting and its logs.

Two things in it are worth knowing before you read them. The Genome section
confirms the assembly from chromosome 1's measured length, and says so plainly
when `input.assembly` is unset and therefore nothing checked it. And the
eRegulon-activity t-SNE is the ONE figure on its own layout: every other is
drawn on the reduction you named in E, so its coordinates are not comparable
with them.

`provenance/<timestamp>_<status>/` is written at the end of every run, success
or failure, and records the commit captured when the run STARTED rather than
when it ended. On a failure it names the logs carrying an error signature.

A rule is skipped when its outputs exist, its inputs and params have not
changed, and nothing upstream is newer. So a re-run continues rather than
starting over, and changing one config value re-runs the rules that read it.
Ask for work with `-f N` or `-- --forceall` rather than by deleting outputs, and
**never delete `.snakemake/`** — four of the five rerun triggers compare against
what is recorded there, and without it a config edit silently does nothing.

## When something fails

1. **`provenance/<timestamp>_error/manifest.txt`** names the logs carrying an
   error signature, so it says where to look rather than handing you a
   directory. It is written on failure as well as success.
2. `logs/<rule>.log` holds that rule's own output. On LSF, `logs/lsf/<rule>.*.out`
   holds the same text plus the accounting, and it survives the re-run that
   fixes the problem, because `bsub -o` APPENDS while the rule log is rewritten.
   Read the LAST block.
3. **There is NO report.html after a failure**, and that is expected rather than
   a second problem: snakemake does not build a target whose inputs failed. The
   logs and the bundle are what a partial run leaves.
4. RUNBOOK section 6 lists the failures already seen and what each one means.
5. RUNBOOK section 7 records what has actually been run, so you can tell whether
   you are on a tested path or the first person to try something.
