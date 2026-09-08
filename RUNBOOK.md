# RUNBOOK — running this pipeline end to end

Written 2026-09-07 from an actual local run of steps 1–3 on a PBMC multiome
fixture. Every number below is measured, not estimated; anything unverified says
so.

---

## 0. The input object

A Seurat object with **RNA counts in `RNA`**, **ATAC counts in `peaks`**, and a
**categorical cell-type column**. The pipeline will not start without the last
one (`input.celltype_column`).

Prepared test object:

    /home/ihlee/limlab_code/scenicplus-pbmc400/pbmc400_annotated.rds    115 MB

    1149 cells · RNA 12210 genes · peaks 61490 · cell_type: 25 levels
    largest: Classical monocytes 243, Naive CD4 T 190, Naive CD8 T 143

It is `pbmc-run-400/5.integrated/seurat.rds` (10x PBMC multiome fixture, hg38)
with SingleR/Monaco labels joined on as `cell_type`. The source object was NOT
modified — it belongs to a snakemake workspace whose timestamp contract would
break if it were.

**Fragment paths are dangling and that is fine.** The object carries 4 Fragment
objects pointing at `/home/ihlee/limlab_code/pbmc-fixture-400/...`, which will
not exist elsewhere. **No step reads them** — verified by grep across every step
script; `scenicplus_03_create_cistopic.py` says so in its own header, and
step 01 reads only the `counts` layers. Copy the .rds alone.

**Making your own input object.** Any Seurat object works if it has those three
things. If yours lacks cell types, either use a clustering column (SCENIC+ only
needs a categorical grouping) or annotate first — e.g. `singler_annotate.R` from
`scRNA_LimLab_Snake` with a celldex reference, then join `labels_per_cell.tsv`
onto a COPY of the object.

---

## 1. Environment

    git clone -b test https://github.com/hwlim/Scenicplus_limlab.git
    cd Scenicplus_limlab
    ./install_local.sh /path/to/scenicplus_env

Uses `environment.local.yml`, not `environment.yml` — the latter does not
produce a working environment (its own header says it was never run-tested).
See that file's header for the four deviations and why each is load-bearing.

Takes ~20 min; phase 2 pulls ~200 packages from PyPI. **If the cluster gates
egress or needs an index mirror, that is where it stops** — the conda layer will
already be complete.

**Install on the node class you will run on.** Login and compute nodes can run
different OS images and therefore different glibc, and pip picks wheels for the
glibc of the host it runs on. On a host below glibc 2.28, three pins
(`pybigtools==0.1.2`, `pysam==0.22.0`, `diptest==0.11.0`) have no usable wheel
and fall back to sdists — `pybigtools` is Rust and fails with *"Cargo, the Rust
package manager, is not installed or is not on PATH"*. The installer now checks
this up front and stops with the reason instead of failing inside a compiler
twenty minutes in. It also prints the glibc it sees in its `### preflight`
block, which is the fastest way to tell the two node classes apart.

Verify:

    /path/to/scenicplus_env/bin/scenicplus --help
    # usage: scenicplus [-h] {init_snakemake,prepare_data,grn_inference} ...

**If that file does not exist**, pip installed scenicplus outside the prefix —
on a shared cluster that means the user site (`~/.local/lib/python3.11/
site-packages`, with the console script in `~/.local/bin`), because
`~/.config/pip/pip.conf` carries `user = true` or `PIP_USER` is set. Repair
without redoing the conda layer:

    SCP_FROM=2 ./install_local.sh /path/to/scenicplus_env      # ~2 min

The user site precedes the env's `site-packages` on `sys.path`, so a stale copy
there also *shadows* the env at run time, silently. Export
`PYTHONNOUSERSITE=1` alongside `PATH` whenever you run the pipeline — the
installer prints it in its "To use" block for that reason.

`LD_LIBRARY_PATH` is the same shape of problem one layer down. `pandas`,
`pyranges`, `anndata` and `mudata` reach libstdc++ through manylinux **wheels**
whose `DT_RUNPATH` does not find the env's copy, so they load the **system**
one — and a soname loaded once is never searched for again, so conda's ICU
(which needs `CXXABI_1.3.15`) then fails against it on any node whose `/lib64`
predates GCC 13. Prepending the env's `lib` makes the first load resolve inside
the env. `scenicplus_run_pipeline.sh` sets this itself; export it too when
running steps by hand.

---

## 2. cisTarget databases — 45.7 GB, needed from step 9 on

Steps 1–8 do NOT need these. Download while the early steps run.

    B=https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based
    curl -O $B/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather   # 32.8 GB -> ctx_db
    curl -O $B/hg38_screen_v10_clust.regions_vs_motifs.scores.feather     # 12.9 GB -> dem_db
    curl -O https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl   # 94 MB

Note the path is `.../hg38/screen/...`, NOT `screening_regions` — the latter
404s. Mouse lives under `mus_musculus/mm10/`.

There is no smaller option: the gene-based databases are 0.3 GB but score gene
promoters, and SCENIC+ scores ATAC **regions**.

---

## 3. Configure

    export SCENICPLUS_PATH=/path/to/Scenicplus_limlab
    export PATH=$SCENICPLUS_PATH/scripts:/path/to/scenicplus_env/bin:$PATH
    export PYTHONNOUSERSITE=1          # ~/.local wins over the env otherwise
    export LD_LIBRARY_PATH=/path/to/scenicplus_env/lib:$LD_LIBRARY_PATH

    mkdir -p /path/to/analysis && cd /path/to/analysis
    scenicplus_init.sh                 # writes config/config.yaml

Edit `config/config.yaml`:

| key | value |
|---|---|
| `input.seurat_rds` | absolute path to your .rds |
| `input.celltype_column` | `cell_type` |
| `input.celltype_scope` | `[]` for all cells, or a subset of cell types |
| `input.species` | `hsapiens` |
| `input.ctx_db` / `input.dem_db` | the two feather files |
| `input.motif_annotations` | the .tbl |
| `resources.n_cpu` | cores for the job |

`Pipeline: "ScenicPlus"` must stay — the driver aborts without it, so another
pipeline's config cannot be fed in by accident.

**`input.assembly` is dead config** — no script reads it (grep across `scripts/`
returns nothing). Setting `mm10` there does NOT switch species; use
`input.species`.

---

## 4. Run

    scenicplus_run_pipeline.sh --dry-run     # plan only, runs nothing
    scenicplus_run_pipeline.sh               # everything stale
    scenicplus_run_pipeline.sh --only 4      # one step
    scenicplus_run_pipeline.sh --from 9      # force step 9 onward
    scenicplus_run_pipeline.sh --force       # everything

    scenicplus_run_workstation.sh            # thin wrapper, runs here
    scenicplus_run_lsf.sh                    # bsubs the driver as one job

A step is skipped iff its sentinel exists, its `.cfgsha` matches the hash of the
config keys **that step** reads, and no upstream sentinel is newer. Once any
step runs, every later step is force-run in the same invocation.

Change one GRN parameter and only that stage and its downstream re-run —
verified: editing `grn.tf_to_gene_importance_method` skipped 1–11 and re-ran
12–20.

---

## 5. The 20 steps

| # | step | sentinel (under `results/`) | notes |
|---|---|---|---|
| 1 | seurat_to_anndata | `interim/seurat_export/summary.txt` | R. **23 s, 1.2 GB** |
| 2 | build_anndata | `interim/rna.h5ad` | **15 s, 0.6 GB** |
| 3 | create_cistopic | `interim/cistopic_obj.pkl` | **5 s, 1.2 GB** |
| 4 | topic_modeling | `interim/cistopic_obj_with_topics.pkl` | **LDA, uses Ray.** Heaviest early step |
| 5 | region_sets | `interim/region_sets/.done` | DARs + topic regions |
| 6 | prepare_gex_acc | `scplus_out/ACC_GEX.h5mu` | first flattened GRN stage |
| 7 | genome_annot | `scplus_out/genome_annotation.tsv` | **needs network** (biomart) |
| 8 | search_space | `scplus_out/search_space.tsv` | |
| 9 | cistarget | `scplus_out/ctx_results.hdf5` | **needs ctx_db** |
| 10 | dem | `scplus_out/dem_results.hdf5` | **needs dem_db** |
| 11 | prepare_menr | `scplus_out/cistromes_direct.h5ad` | |
| 12 | tf_to_gene | `scplus_out/tf_to_gene_adj.tsv` | |
| 13 | region_to_gene | `scplus_out/region_to_gene_adj.tsv` | |
| 14 | egrn_direct | `scplus_out/eRegulons_direct.tsv` | |
| 15 | egrn_extended | `scplus_out/eRegulons_extended.tsv` | |
| 16 | aucell_direct | `scplus_out/AUCell_direct.h5mu` | |
| 17 | aucell_extended | `scplus_out/AUCell_extended.h5mu` | |
| 18 | scplus_mudata | `scplus_out/scplusmdata.h5mu` | the GRN result |
| 19 | postprocess_tsv | `tsv/eRegulons_combined.tsv` | |
| 20 | visualize | `plots/01_umap_celltype.pdf` | |

Steps 6–18 replace what used to be one opaque inner `snakemake` call. Each is
now a direct `scenicplus` CLI invocation with its own sentinel and config hash.
All 11 subcommands and every flag were checked against the installed CLI: **0
unknown flags**.

Logs: `logs/NN_<step>.log`, one per step.

---

## 6. Known issues

**Ray fails inside a restricted sandbox.** Step 4 uses `ray.init()` via
pycisTopic (`lda_models.py:154`, unconditional — `n_cpu=1` does not avoid it).
Under Claude Code's sandbox the raylet dies constructing a **Unix domain
socket** for the plasma store, and the node times out. `/dev/shm` was 3.9 GB and
free, so it is not a shared-memory limit. Expected to be fine on HPC; run
outside the sandbox locally.

**The 7 GB WSL box cannot run steps 9+.** The cisTarget databases are 45.7 GB
and SCENIC+ is not sized for that machine. Steps 1–3 run there comfortably.

**`resources.n_cpu` is read from config only** — no env override. Step 4 does
NOT hash it, so changing it will not re-run step 4 (correct: thread count should
not change results).

**An interrupted install leaves a broken package cache.** Re-running then fails
with, for every affected package,

    CondaVerificationError: The package for r-base located at
    <pkgs>/r-base-4.5.3-h502d0c9_3 appears to be corrupted. The path
    'share/man/man1/Rscript.1' specified in the package manifest cannot be found.

Nothing is really corrupted -- the directories were partially EXTRACTED when the
job died, and conda verifies them against the package manifest. The downloaded
archives are usually fine. Fix:

    SCP_CLEAN_PKGS=1 ./install_local.sh /path/to/env

which drops the unpacked directories, keeps the archives (so nothing is
re-downloaded), and resumes into the existing env folder.

**Old conda may not solve the env.** `install_local.sh` prefers micromamba, then
mamba, then conda. A 2020-vintage conda on the classic solver may be very slow
or fail on the R + Seurat + Signac layer.

---

## 7. What has actually been run, and where

| | status |
|---|---|
| steps 1–3 | **run on real data** (PBMC 400 fixture, local) |
| step 4 | **blocked locally by the sandbox** (Ray), not attempted elsewhere |
| steps 5–20 | **never executed** |
| flattened stages 6–18 | flags validated against the Snakefile AND the installed CLI; **not run** |
| driver sentinel/cascade logic | validated by simulation: all-skip baseline, cascade from exactly the edited step, and a negative control (an unhashed key changes nothing) |

Treat the first cluster run of steps 5–20 as their real test. The history of the
sister pipeline is that every component which actually ran turned up a defect
the static checks missed — this one already has: step 3 logged
`n_regions = <cell count>` until it was run and the number looked wrong.
