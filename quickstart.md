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
| F. running | LSF, through a small `run.sh` | `scenicplus_run_workstation.sh` |
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

    scenicplus_run_pipeline.sh --dry-run       # prints the plan, runs nothing

**On your own machine:**

    scenicplus_run_workstation.sh              # runs here, in this shell

**On the CCHMC HPC**, every run to date has gone through
`scenicplus_run_lsf.sh`, called from a small `run.sh` kept in the analysis
directory that sets up the environment and passes the LSF settings. Keeping it
there rather than in shell history puts the settings a run used beside that
run's outputs, and makes a re-run one command. A skeleton to adapt:

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

# Reproducibility. Without these, two runs of the same data on the same
# cluster do not agree. Both effects are measured, not theoretical:
#   PYTHONHASHSEED  python randomises string hashing per process, and SCENIC+
#                   builds lists from sets of names, so the order changes every
#                   invocation. Pandas sorts are stable, so a permuted order
#                   permutes ties and a top-N cut then keeps different rows.
#   *_NUM_THREADS   the maths libraries take their thread count from the host's
#                   core count, and reduction order follows thread count. A
#                   48-core and a 64-core node gave correlations differing in
#                   the last bit.
# Keep the thread number equal to resources.n_cpu, which is what the Snakemake
# workflow pins to; the two drivers are only comparable if they agree.
export PYTHONHASHSEED=0
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMEXPR_NUM_THREADS=16

export LSF_QUEUE=normal LSF_PROJECT=scenicplus
export LSF_CORES=16 LSF_MEM_MB=128000 LSF_WALLTIME=72:00
scenicplus_run_lsf.sh "$@"                  # --from 7, --only 20, --force, ...
```

Then `./run.sh` submits, and `./run.sh --from 7` forwards the flag.

**Submit from the node class the environment was built for.** The preflight runs
on the host you submit from, before anything is queued, so it checks the shell
you are typing in rather than the one the job will get. Two consequences, both
of which have cost time here: R must be loaded to submit and not only to run,
and a login node whose glibc or modules differ from the compute nodes fails the
check with nothing queued. Start an interactive session on the right class and
run the runner there.

**That same shell becomes the job's environment.** `scenicplus_run_lsf.sh`
submits `/bin/bash -c` and loads no modules inside the job, while LSF carries
the submission environment across, so whatever `run.sh` sets is what all twenty
steps run under. `LD_LIBRARY_PATH` is the exception, since the driver prepends
the environment's `lib` itself.

The memory figure above is an example, not a measured requirement. RUNBOOK
section 5b explains what drives it: step 4 dominates, and its peak is the sum
over the topic models fitted at once.

## G. Watch, then read

RUNBOOK sections 4 and 5.

    bjobs                                      # HPC only
    tail -f logs/01_seurat_to_anndata.log      # one log per step
    ls results/plots results/tables

Every UMAP figure names in its title the layout it was drawn on. If that is not
the reduction you chose, re-read E.

A step is skipped when its output exists, the config keys it reads have not
changed, and nothing upstream is newer. So a re-run continues rather than
starting over, and changing one parameter re-runs only the steps that read it.
The driver decides all of that from files under `results/`, so ask for work with
`--from N` or `--force` rather than by deleting outputs.

## When something fails

1. `logs/NN_<step>.log` holds that step's own output, and the driver names the
   step it stopped on.
2. RUNBOOK section 6 lists the failures already seen and what each one means.
3. RUNBOOK section 7 records what has actually been run, so you can tell whether
   you are on a tested path or the first person to try something.
