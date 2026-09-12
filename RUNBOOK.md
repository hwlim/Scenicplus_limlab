# RUNBOOK — running this pipeline end to end

Written 2026-09-07 from an actual local run of steps 1–3 on a PBMC multiome
fixture, and extended since through complete cluster runs: human/hg38 twice and
mouse/mm10 once. Every number below is measured, not estimated; anything
unverified says so.

**Where this has run: CCHMC, and nowhere else.** Their LSF cluster, module
stack, queue names and filesystem. The environment *recipe* also builds on a
WSL2 workstation, but building is not running. The launchers and their queue and
module names, the memory unit LSF reads, the login-versus-compute glibc split
that decides which wheels pip picks, and the assumption that `pip.conf` may
carry `user = true` are all shaped by that site. **Somewhere else, expect to
adjust before anything runs, and treat the first run as a port rather than an
install.** Section 1 says what to look at first, and section 6 lists what is
already known to differ.

---

## Quickstart

Moved to [quickstart.md](quickstart.md), which is the path to follow: install
through submission, for the CCHMC HPC or a personal machine. It points back here
for the reasoning behind each step, and this file stays the reference rather than
the tutorial.

---

## 0. The input object

A Seurat object with **RNA counts in `RNA`**, **ATAC counts in `peaks`**, and a
**categorical cell-type column**. The pipeline will not start without the last
one (`input.celltype_column`).

The object this runbook's timings and checks come from — a 10x PBMC multiome
fixture (hg38), `5.integrated/seurat.rds` from a `scRNA_LimLab_Snake` run with
SingleR/Monaco labels joined on as `cell_type`:

    115 MB · 1149 cells · RNA 12210 genes · peaks 61490 · cell_type 25 levels
    largest: Classical monocytes 243, Naive CD4 T 190, Naive CD8 T 143

The labels were joined onto a COPY; the source object was not modified, because
it belongs to a snakemake workspace whose timestamp contract would break if it
were.

**Fragment paths are dangling and that is fine.** The object carries Fragment
objects pointing at paths recorded when it was built, which will not exist on
another machine. **No step reads them** — verified by grep across every step
script; `scenicplus_03_create_cistopic.py` says so in its own header, and
step 01 reads only the `counts` layers. Copy the .rds alone.

**Making your own input object.** Any Seurat object works if it has those three
things. If yours lacks cell types, either use a clustering column (SCENIC+ only
needs a categorical grouping) or annotate first — e.g. `singler_annotate.R` from
`scRNA_LimLab_Snake` with a celldex reference, then join `labels_per_cell.tsv`
onto a COPY of the object.

---

## 1. Environment

    git clone -b main https://github.com/hwlim/Scenicplus_limlab.git
    cd Scenicplus_limlab
    ./install_cchmc.sh /path/to/scenicplus_env

Uses `environment.cchmc.yml`, not `environment.yml` — the latter does not
produce a working environment (its own header says it was never run-tested).
See that file's header for the four deviations and why each is load-bearing.

> **Both files are named `cchmc` because that is the only site they have been
> run at.** The env recipe itself solves and builds on two machines (CCHMC
> compute nodes via conda, a WSL2 workstation via micromamba), but **the
> pipeline has run end to end on CCHMC only**. The LSF launchers, the queue and
> module names, the login-vs-compute glibc split that decides which wheels pip
> picks, and the `user = true` pip.conf assumption are all CCHMC-shaped. A first
> run at another site is a port, not an install.

Takes ~20 min; phase 2 pulls ~200 packages from PyPI. **If the cluster gates
egress or needs an index mirror, that is where it stops** — the conda layer will
already be complete.

Phase 3 adds one pip package, `snakemake-executor-plugin-cluster-generic`, which
the Snakemake workflow needs to submit jobs and the bash driver does not use at
all. `SCP_NO_CLUSTER=1` skips it; `SCP_FROM=3` adds it to an environment built
before phase 3 existed. Its version is pinned deliberately: scenicplus pins
snakemake and its plugin interfaces with `==`, so a newer plugin can only be
installed by breaking the scenicplus install. `install_cchmc.sh` and
`profiles/lsf/config.yaml` both carry the detail.

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

    SCP_FROM=2 ./install_cchmc.sh /path/to/scenicplus_env      # ~2 min

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

## 2b. Mouse (mm10) — the assembly will not take care of itself

Everything above assumes human/hg38. Mouse needs four things changed, and one of
them is a trap that does **not** announce itself.

**The trap.** Step 7 asks Ensembl what assembly it is working with and gets
whatever Ensembl serves *today*, which for mouse is **GRCm39**. Against mm10 /
GRCm38 data that is the wrong assembly — and unlike the missing chromsizes, it
does not stop anything. Chromosome names still convert, the search space still
builds, and the peak–gene links come out quietly wrong. From a 2026-05 kidney
run (`Development.md:35`):

    Download gene annotation INFO   Using genome: GRCm39
    Could not find Id on ...esearch.fcgi?db=genome&term=GRCm39

Both halves of that line are problems. The second is the dead endpoint (§6); the
first is the one that would have survived to the results.

**So build the genome files yourself.** `scripts/scenicplus_make_genome_files.R`
derives both from a pinned EnsDb, so the assembly is the one you chose rather
than the one Ensembl currently ships:

    # in the scRNA_LimLab_Snake env — it needs EnsDb + BSgenome, which
    # environment.cchmc.yml deliberately does not carry
    Rscript $SCENICPLUS_PATH/scripts/scenicplus_make_genome_files.R \
        --species mouse --out-dir /path/to/genome

It prints `chr1 = 195,471,971 bp` for mm10. That single number is the check:
GRCm39's chr1 is 195,154,279. If you ever wonder which assembly a run used, look
there.

Deriving it from `EnsDb.Mmusculus.v79` — the same annotation
`scRNA_LimLab_Snake` used to build the peaks — means the peaks and the
annotation share an assembly *by construction* rather than by two services
happening to agree.

**The four config changes:**

| key | mouse value |
|---|---|
| `input.species` | `"mmusculus"` |
| `input.genome_annotation` / `input.chromsizes` | the two files just built |
| `input.ctx_db` / `input.dem_db` | `.../mus_musculus/mm10/...` |
| `input.motif_annotations` | `motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl` (mgi, not hgnc) |

`input.assembly` is **dead config** — no script reads it; setting `mm10` there
changes nothing.

**One key, two vocabularies.** `input.species` feeds step 7 (which wants the
Ensembl short form `mmusculus`) *and* steps 9–10 (which want
`mus_musculus`). They disagree, and it works only because `load_motif_annotations`
consults the species **solely when no annotation file is given**
(`pycistarget/utils.py:98`). Since `input.motif_annotations` is always set here,
the steps 9–10 spelling is inert. Leave that path unset and it fails.

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
| `input.reduction` | the Seurat reduction every figure is drawn on, e.g. `wnn.umap`. Empty means one is computed, unintegrated |
| `input.species` | `hsapiens` |
| `input.ctx_db` / `input.dem_db` | the two feather files |
| `input.motif_annotations` | the .tbl |
| `resources.n_cpu` | cores for the job |

`Pipeline: "ScenicPlus"` must stay — the driver aborts without it, so another
pipeline's config cannot be fed in by accident.

**`input.assembly` is dead config** — no script reads it (grep across `scripts/`
returns nothing). Setting `mm10` there does NOT switch species; use
`input.species`.

**Set `input.reduction` on an integrated object.** Steps 01, 02 and 20 read it,
and it names the reduction that becomes `X_umap`: the layout the cell-type
figure and every per-eRegulon figure are drawn on. Give the name exactly as
`Reductions(obj)` prints it. A name the object does not have fails at step 01,
in the first minute, with the available names listed.

Leaving it empty is a choice, not a default. Nothing then supplies a layout, so
step 02 computes a UMAP from the RNA matrix with no batch correction and step 20
falls back to a UMAP of eRegulon activity. Both are defensible for one sample;
for a multi-sample integrated object neither is the picture you looked at in
Seurat, and the figure does not say which one it is. Every figure now carries
the layout's name in its title, so this can be read off the output rather than
the log.

The prefix rule this replaced matched none of Seurat's actual reduction names —
`wnn.umap`, `rna.umap`, `umap.harmony` all fell through — so every run so far
plotted a recomputed unintegrated layout. Grep any past run for the line that
records it:

    grep -h 'X_umap\|No UMAP' logs/02_build_anndata.log

**This affects figures only.** No part of the GRN inference reads a reduction:
with `scenicplus.is_multiome: true` cells are paired by barcode rather than
through metacells, topics come from the fragment counts, and DARs and RSS come
from `celltype_column`. A finished run therefore does not need recomputing — see
§5b for redrawing its figures alone.

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
| 7 | genome_annot | `scplus_out/genome_annotation.tsv` | biomart; **chromsizes half is broken upstream** — see below |
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

## 5b. Making a fresh run faster

**Measure first.** Every step writes `logs/NN_*.log` and the driver prints a
banner per step, so the real per-step wall clock is already on disk:

    for f in logs/*.log; do
      printf "%-28s %s -> %s\n" "$(basename "$f")" \
        "$(head -1 "$f" | grep -oE '[0-9]{2}:[0-9]{2}:[0-9]{2}' | head -1)" \
        "$(tail -5 "$f" | grep -oE '[0-9]{2}:[0-9]{2}:[0-9]{2}' | tail -1)"
    done

Tune what that says is slow. The knobs below are ordered by expected leverage,
but ONLY steps 1-3 have measured times in the table above; the rest is reasoning
about what each parameter multiplies, not a benchmark.

| knob | default | cheaper | what it costs you |
|---|---|---|---|
| `cistopic.n_topics` | `[2,5,10,20,30,40,50]` | `[10,20,30,40]` | step 4 runs one LDA per value, all at once (one CPU each). Wall clock is the SLOWEST model, so dropping the largest counts helps most; peak memory is the SUM, so it helps there too. Dropping the small ones saves memory, not time. |
| `cistopic.n_iter` | 150 | 100 | Gibbs iterations, linear in step 4's time. Below ~100 the topic model may not have converged — fine for a plumbing test, not for results. |
| `scenicplus.search_space_*` | `1000 150000` | `1000 75000` | halves the region-gene pairs steps 13-15 grind through. On the PBMC fixture 150 kb gave 185k links with a median distance of 53 kb, so 75 kb keeps most of the mass. Genuinely changes the biology: distal enhancers beyond the cut are invisible. |
| `grn.gsea_n_perm` | 1000 | 250 | linear in steps 14-15. Coarsens the eRegulon p-values; use for a smoke run only. |
| `grn.quantile_thresholds_region_to_gene` + `top_n_regionTogenes_per_gene` | 3 + 3 values | 1 + 1 | each combination is a separate region set to score downstream. Going from 3x3 to 1x1 cuts that work ~9x. |
| `resources.n_cpu` | 16 | — | helps steps 9-13, which parallelise. See the step-4 memory note: there it is a MULTIPLIER on peak memory, not just on speed. |

What will NOT get faster by tuning: **steps 9 and 10** read the 32.8 GB and
12.9 GB cisTarget feathers, and that I/O dominates them. They are also
sentinel-cached, so the cost is paid once per workspace, not per re-run.

**Redrawing the figures of a finished run costs one step.** Setting
`input.reduction` changes the `.cfgsha` of steps 1, 2 and 20, and once step 1
runs the cascade forces all twenty — which is the right default, but pointless
here, since no stage between them reads a reduction. Step 20 reads the layout
straight out of step 01's `seurat_export/embedding_<name>.tsv`, which a finished
workspace already has for every reduction the object carried:

    # set input.reduction first, then
    scenicplus_run_pipeline.sh --only 20

`--only` bypasses the cascade, so this touches nothing but `results/plots/`.
Check the new titles: each figure names the layout it was drawn on.

**For a plumbing test rather than a result**, the combination that changes the
least science per second saved is `n_topics: [10,20,30]` + `gsea_n_perm: 250` +
a single quantile/top-n value. Put those in a separate analysis directory --
changing them in place re-runs from whichever step reads them onward, and the
`.cfgsha` cascade means step 4's list re-runs everything after it.

**Step 7 cannot produce chromsizes, for anyone.** It derives them from an NCBI
E-utilities lookup against `db=genome`, a database NCBI has retired: it answers
HTTP 200 with `<Count>0</Count>` for every term, including plain
`Homo sapiens`. Verified from two unrelated networks — this is not a firewall.
SCENIC+ catches the miss, logs *"Chromosome sizes will not be returned"*, and
**exits 0**, so the driver marks step 7 done.

The same branch also performs the **UCSC chromosome-name conversion**, so the
annotation it does write is Ensembl-style (`1`, `2`, `X`, `MT`) and will overlap
nothing against `chr`-prefixed ATAC regions. One dead endpoint, two failures,
two stages apart.

Set `input.genome_annotation` and `input.chromsizes` (§3) and step 7 makes no
network call at all. Build them once:

    curl -O https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
    awk 'BEGIN{OFS="\t"; print "Chromosome","Start","End"} {print $1,0,$2}' \
        hg38.chrom.sizes > chromsizes.tsv
    # and convert the annotation biomart returned:
    awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print;next}
      {if($1=="MT")$1="chrM"; else if($1!~/^chr/)$1="chr"$1; print}' \
      genome_annotation.tsv > genome_annotation.ucsc.tsv

Step 7 refuses if the two disagree in naming; step 8 reports all three sources
if they still do.

**Barcodes must match between steps 3 and 6.** `create_cistopic_object` defaults
to `tag_cells=True`, which appends `___<project>` to every ATAC cell name, while
the RNA AnnData keeps the plain barcode — so step 6 intersects two disjoint sets
and dies with *"No cells found which are present in both assays"*, even though
both sides came from the same Seurat object. Step 3 now passes
`tag_cells=False`. A workspace whose `cistopic_obj.pkl` predates that fix has
tagged names already; either re-run `--from 3` (redoes topic modeling) or set

    scenicplus:
      bc_transform_func: 'lambda x: x + "___scenicplus_run"'

which re-runs step 6 onward only. The two must change together — the lambda is
applied to the RNA barcodes to map them onto the ATAC ones, so it is wrong
against an untagged cisTopic object. Step 6 prints both name shapes and the
exact lambda when the intersection is empty.

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

    SCP_CLEAN_PKGS=1 ./install_cchmc.sh /path/to/env

which drops the unpacked directories, keeps the archives (so nothing is
re-downloaded), and resumes into the existing env folder.

**Old conda may not solve the env.** `install_cchmc.sh` prefers micromamba, then
mamba, then conda. A 2020-vintage conda on the classic solver may be very slow
or fail on the R + Seurat + Signac layer.

---

## 7. What has actually been run, and where

Updated 2026-09-09.

| | status |
|---|---|
| **steps 1–20, CCHMC cluster, human/hg38** | **run to completion, twice** |
| — first pass | steps 1–19 with a pre-fix step 3 (tagged barcodes) plus a `bc_transform_func` workaround; step 20 after the plotnine fix |
| — second pass | **from scratch on the current code**, stopping once at step 7 to build the genome files, which were supplied through `input.genome_annotation` / `input.chromsizes` |
| steps 1–3, local WSL workstation | run on the PBMC-400 fixture |
| step 4, local | **blocked by the sandbox** (Ray's plasma socket); never attempted locally since |
| **mouse / mm10** | **run to completion 2026-09-09**, reported by the operator; the artifacts have not been examined here. It is what turned up the `input.reduction` defect below |
| **Snakemake workflow, steps 1-18** | **validated on the cluster 2026-09-11**: run against the bash driver on the same data, with reproducibility pinning. Every text output matches by md5; every figure PNG matches by md5 and was read side by side; residual numeric differences under 1e-9. Compare PNGs, never PDFs — matplotlib stamps `/CreationDate`, so identical figures hash differently |
| **Per-rule LSF resources (I5)** | **built, NOT yet cluster-validated.** The 2026-09-11 rerun finished cleanly but ran on PRE-I5 code -- its epilogues report `Total Requested Memory: 128000.00 MB` for R19 and R20, which is the old uniform tier; under I5 they ask for 16000. It did supply R19/R20's real usage, now in `tests/measured_resources.tsv`, because what a rule USES is independent of what it reserved. The tiers themselves still want a run on this code |
| **`report.html` (I6)** | **built 2026-09-11; rendered once for real, one defect found and fixed.** `rule all` targets it, so a run is not finished until it is readable. The real render flagged the report's OWN log as empty -- a false alarm on every run, since `tee` creates it before this script prints; the row is now annotated and excluded. Otherwise gated only against a synthetic workspace (`tests/report_render.sh`, deliberately PARTIAL, asserting what the page says is MISSING) |
| **Provenance bundle (I7)** | **built 2026-09-11, local gates only.** `provenance/<ts>_<status>/` per run, from onsuccess AND onerror: manifest, config.used.yaml, this run's logs, lsf_jobs.tsv, assembly.json, snakemake.log. The commit is captured at ONSTART, so a checkout mid-run cannot make it name a commit that produced nothing. No bundle from a real GRN run yet |
| **Retiring the bash driver (I8)** | **HELD 2026-09-12, deliberately.** Nothing built since I3 has completed a cluster run on the current code, so deleting the driver would remove the fallback before the replacement is proven. One end-to-end run settles I5's tiers, R21's tier, I6's report, I7's bundle and the three new figures at once. Before it: confirm a node has the 256000 MB R09 reserves, or the job pends rather than fails |
| reproducibility | `PYTHONHASHSEED` and the BLAS thread variables must be pinned, or two runs of the same data DISAGREE -- by eRegulon membership, not just by bytes. The workflow pins them for every rule; a `run.sh` driving the bash pipeline must export them itself (Quickstart G) |
| `input.reduction` | **confirmed on the cluster 2026-09-09**: a named reduction produces the figure it names. The `--only 20` redraw of a finished run is the exercised path |
| submission | every cluster run has used `scenicplus_run_lsf.sh`, driven by a per-workspace `run.sh` ([quickstart.md](quickstart.md) section F). `scenicplus_run_lsf_cchmc.sh` has never been the path in use |
| driver sentinel/cascade logic | validated by simulation, then in practice — a `grn.*` edit re-ran 12–20 and skipped 1–11 |

**Outputs that have been looked at,** as opposed to merely produced:
`search_space.tsv` (185,070 links, 56,032 of 61,490 peaks, median TSS distance
53 kb, zero non-standard contigs) and `03_rss_per_celltype` (recovers SPIB and
BCL11A in naive B, LEF1 in naive CD4 T, KLF4 in classical monocytes, CEBPA and
MAFB in intermediate monocytes, TBX21 in effector/MAIT — consistent with the
FigR result on the same data). **Not** looked at: the eRegulon tables, and the
topic count LDA settled on.

**`scenicplus_make_genome_files.R`** has been run for both species but never
consumed by a pipeline run. Its human output was cross-checked against the
chromsizes that produced the working PBMC run: all 25 shared chromosomes have
identical lengths.

### The prediction in the old version of this section held

It said to treat the first cluster run as the real test, because every component
that actually ran turned up a defect the static checks missed. It did — **eight
failures that stopped a run**, and two more found only by looking at an output:

| # | stopped at | the two things that had to agree |
|---|---|---|
| 1 | install "succeeded" | user site vs env `site-packages` |
| 2 | pip picked an sdist | login-node glibc vs compute-node glibc |
| 3 | OOM at 128 GB | `rusage[mem]` reserved vs `-M` enforced |
| 4 | step 6 import | a wheel's libstdc++ vs conda's |
| 5 | step 6 cell overlap | cisTopic's tagged barcodes vs the RNA AnnData's |
| 6 | step 7 exited 0 | what it wrote vs what step 8 needed |
| 7 | step 8 `KeyError` | Ensembl vs UCSC chromosome names |
| 8 | step 20 `savefig` | matplotlib's API vs plotnine's |
| 9 | *nothing* | a 301-megapixel figure, found by opening it |
| 10 | *nothing* | Seurat's reduction NAMES vs the `umap` prefix step 2 tested for |

Every one named neither side. All are now fixed at the source with a check that
says which disagreed.

**So the same warning, pointed at what is still untested:** mouse is the next
real test, and the parameters have never been examined at all — the topic-count
sweep, the DAR thresholds and the search-space width are as shipped, chosen by
nobody. A green run says the plumbing holds, not that the numbers mean anything.
