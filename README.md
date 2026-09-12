# SCENIC+ pipeline (Seurat -> TF-gene network)

A 20-step pipeline that runs SCENIC+ end-to-end on a single Seurat `.rds`
multiome object (RNA + peaks/ATAC) and produces:

- `5.analysis/tsv/` — comprehensive TSVs of TF-region-gene links, AUC, RSS, TF summary
- `5.analysis/plots/` — interpretable PDF + PNG plots (one per file)
- `4.grn/scplusmdata.h5mu` — the canonical SCENIC+ MuData
- `report.html` — the run, readable: provenance, genome check, tables, every
  figure, this run's LSF accounting and logs
- `provenance/<timestamp>_<status>/` — a bundle written on success AND failure

**A Snakemake workflow runs it, and `scenicplus.run.sh` is the entry point.**
Each of the 21 rules is scheduled separately with its own measured LSF
reservation, so the plotting step no longer holds the motif-enrichment step's
sixteen cores and 256 GB. `rule all` targets `report.html`, so a run is not
finished until it is readable.

SCENIC+'s own GRN inference (steps 06-18) was flattened out of its opaque inner
snakemake into discrete `scenicplus` CLI stages, one rule each, dispatched by
`scenicplus_06_grn_stage.py`.

`scenicplus_run_pipeline.sh`, the original bash driver, still works and prints
an obsolescence notice. It is kept as a fallback and for reproducing older
runs. Do not start there: it submits all twenty steps as ONE job sized for the
heaviest, cannot run independent stages at once, hashes config rather than
code, and produces neither a report nor a provenance bundle.

## Layout

The pipeline is meant to live in **one central location** (referenced via the
`SCENICPLUS_PATH` environment variable) and is invoked from per-analysis
directories. Each analysis directory contains its own `config/config.yaml`
and outputs; the pipeline code itself is shared and never edited per-run.

All executables share a `scenicplus_` prefix to avoid collisions on `PATH`.

```
$SCENICPLUS_PATH/                          # central pipeline (set by the user)
├── Snakefile                              # orchestration; includes rules/
├── rules/                                 # common.smk + one file per stage group
├── schemas/config.schema.yaml             # rejects an unknown or mistyped key
├── profiles/lsf/                          # one bsub per rule
├── config/config.yaml                     # template — copied into each analysis dir
├── tests/                                 # local gates; none needs a cluster
└── scripts/                               # all executables ($PATH)
    ├── scenicplus.run.sh                  # THE ENTRY POINT
    ├── scenicplus_init.sh                 # bootstrap a new analysis directory
    ├── scenicplus_check.sh                # preflight: required Python + R packages
    ├── scenicplus_genome_prepare.py       # R07: validates the genome pair
    ├── scenicplus_09_report.py            # R21: builds report.html
    ├── scenicplus_provenance.py           # onstart/onsuccess/onerror bundles
    ├── scenicplus_helper.py               # YAML-slice + sha256 helper (driver-era)
    ├── scenicplus_run_pipeline.sh         # OBSOLETE 20-step driver, kept
    ├── scenicplus_run_workstation.sh      # OBSOLETE thin wrapper
    ├── scenicplus_run_lsf.sh              # OBSOLETE bsub wrapper
    ├── scenicplus_01_seurat_to_anndata.R  # step scripts (called by the rules
    ├── scenicplus_02_build_anndata.py     #  — not invoked by users)
    ├── scenicplus_03_create_cistopic.py
    ├── scenicplus_04_topic_modeling.py
    ├── scenicplus_05_region_sets.py
    ├── scenicplus_06_grn_stage.py         # dispatches GRN stages 06-18
    ├── scenicplus_07_postprocess_tsv.py   # step 19
    └── scenicplus_08_visualize.py         # step 20

<analysis-dir>/                            # one of these per dataset
├── config/config.yaml                     # the only file the user edits
├── 0.input/  1.export/  2.anndata/  3.cistopic/  4.grn/  5.analysis/
├── QC/                                    # model selection, per-stage counts
├── logs/                                  # one per rule, plus logs/lsf/
├── provenance/<timestamp>_<status>/       # one bundle per run
└── report.html                            # what `rule all` targets
```

The numbered directories are pipeline STAGES, in order. The mapping from the
old driver's `results/` and `interim/` layout is in RUNBOOK section 5.

## Setup

1. Place this repository at the location you want to be the central
   `SCENICPLUS_PATH`, then export it and add the launcher scripts to `PATH`.
   Add this to your shell rc (`~/.bashrc`, `~/.zshrc`):

   ```bash
   export SCENICPLUS_PATH=/abs/path/to/ScenicPlus
   export PATH="$SCENICPLUS_PATH/scripts:$PATH"
   ```

2. Create the conda environment (see section #4 for macOS):

   ```bash
   ./install_cchmc.sh /path/to/scenicplus_env
   ```

   **Use `install_cchmc.sh`, not `conda env create -f environment.yml`.**
   `environment.yml` is the original recipe and its own header says the
   pipeline was never run-tested against it; it does not produce a working
   environment. `install_cchmc.sh` builds from `environment.cchmc.yml`, whose
   header lists the four deviations and why each is load-bearing. Takes about
   20 minutes. RUNBOOK section 1 covers the rest, including what to do when the
   cluster gates outbound PyPI.

   Both `cchmc` files are named for the only site they have been run at. A
   first run elsewhere is a port, not an install.

   The environment provides Python 3.11.8, R with Seurat / Matrix / optparse,
   and `scenicplus` (which pulls in `pycisTopic`, `pycistarget`, scanpy,
   anndata, mudata, pyranges, pybedtools and the rest).

   You can also use a pre-existing env — just point the launcher at it via
   `SCENICPLUS_ENV=<env_name>` and it will `conda activate` for you.

   Verify the env at any time:

   ```bash
   scenicplus_check.sh                              # check current shell
   SCENICPLUS_ENV=scenicplus_limlab scenicplus_check.sh    # check a specific env
   ```

   `scenicplus.run.sh` runs this preflight before handing off to snakemake --
   deliberately in the runner rather than as a rule, so a bad environment is
   found once, up front, instead of once per rule after twenty queue waits.
   Set `SCENICPLUS_SKIP_CHECK=1` to skip it (useful when the compute node has a
   different env from the submitting host).

3. Download external resources (one-time):
   - **cisTarget databases (rankings + scores)** — **required**. These large
     feather files are not bundled in any pip package; you must download and
     point `input.ctx_db` / `input.dem_db` at them.
     <https://resources.aertslab.org/cistarget/databases/>
   - **Motif-to-TF annotations** — auto-downloaded by pycistarget at runtime
     for `homo_sapiens`, `mus_musculus`, and `drosophila_melanogaster` if
     `input.motif_annotations` is left blank. Pre-fetching is still
     recommended (compute nodes often lack outbound HTTPS, and required for
     any other species).
     <https://resources.aertslab.org/cistarget/motif2tf/>

4. For mac os, the currrent environment.yml may not work straight due to limited
   availability of the precompiled packages. For example, scenicplus installation
   via pip would gives error regarding pybedtools and setuptools. In such case,
   step-by-step installation as shown below is advised.

   ```bash
   conda env create -f "$SCENICPLUS_PATH/environment_macos.yml"
   conda activate scenicplus_limlab
   pip install --no-build-isolation pybedtools==0.9.1
   pip install "scenicplus @ git+https://github.com/aertslab/scenicplus.git"
   ```

## Running an analysis

```bash
# 1. Create / enter the analysis directory.
mkdir -p ~/projects/cerebellum_scenicplus
cd       ~/projects/cerebellum_scenicplus

# 2. Copy the config template into this directory (creates ./config/config.yaml).
scenicplus_init.sh
#   scenicplus_init.sh -f       # to overwrite an existing config.yaml
#   scenicplus_init.sh /path    # to initialize a different directory

# 3. Edit the config: paths to RDS, ctx/dem DBs, motif annotations, species.
$EDITOR config/config.yaml

# 4. Look at the plan before spending anything.
scenicplus.run.sh -n

# 5a. Workstation
SCENICPLUS_ENV=scenicplus_limlab scenicplus.run.sh -j 8

# 5b. LSF — one job per rule, at most 20 in flight
SCENICPLUS_ENV=scenicplus_limlab scenicplus.run.sh --lsf -j 20
```

`scenicplus.run.sh` runs from the **current working directory** (your analysis
dir) and reads `./config/config.yaml`. It runs the preflight first, then hands
off to snakemake; anything after `--` goes straight to snakemake.

A run ends with `report.html`. Open that first.

### Re-running

Snakemake decides what is stale, from the file times, the code and the config
values each rule declares. Editing a config key re-runs the rules that read it
and everything downstream of those, and nothing else — `cistopic.dar_log2fc_thr`
re-runs step 05 onward while step 04's LDA output stays put.

```bash
scenicplus.run.sh -n              # show what would run
scenicplus.run.sh -f 4            # force rule R04_* and its DEPENDENTS
scenicplus.run.sh -- --forceall   # force everything
```

**`-f N` names a rule, not a range.** It becomes `--forcerun R<NN>_*`, so it
re-runs that rule and what depends on it. The DAG forks, so `-f 9` leaves R10
and R13 alone and `-f 14` leaves R15 and R17 alone. The obsolete driver's
`--from N` did mean "every step numbered N or higher"; the two are not
equivalent. RUNBOOK section 4 has the table.

**Never delete `.snakemake/`.** It holds the per-output record that four of the
five rerun triggers compare against, and without it a config change silently
does nothing. Recover with `-- --forceall`, never with `--touch`.

The workflow validates `config.yaml` against `schemas/config.schema.yaml` at
DAG-build time, before any job is submitted: a mistyped key is named and
refused rather than silently falling back to a default. It also requires
`Pipeline: "ScenicPlus"`, so an unrelated pipeline's config cannot be fed in by
accident.

## Inputs (set in `config/config.yaml`)

1. `input.seurat_rds` — Seurat `.rds` with:
   - `RNA` assay containing **normalized** expression in `data` and raw counts in `counts`,
   - `peaks` assay containing peak counts (regions named `chr-start-end` or `chr:start-end`),
   - cell-type column (`input.celltype_column`) in `@meta.data`.
   - (optional) UMAP/PCA reductions — every one is exported;
     `input.reduction` decides which the figures use.
2. `input.reduction` (optional, but set it on an integrated object) — the
   Seurat reduction every figure is drawn on, named exactly as
   `Reductions(obj)` prints it (`wnn.umap`, `umap.harmony`, ...). Empty means
   a layout is computed instead, with no batch correction, which for a
   multi-sample integrated object is not the one you clustered on. A name the
   object does not have fails at step 01, with the available names listed.
   Figures only: no GRN stage reads a reduction.
3. `input.celltype_scope` (optional) — a list of cell-type values from
   `celltype_column`. When set, step 01 drops all other cells before exporting,
   so topic modeling, DARs, metacells, and RSS all operate on the focused
   universe (good for sensitivity when you care about a specific lineage).
   Leave empty / omit to use every cell.
4. `input.ctx_db`, `input.dem_db` — cisTarget feather databases.
5. `input.motif_annotations` — motif-to-TF table.

## Outputs

| Path | Description |
| ---- | ----------- |
| `5.analysis/tsv/eRegulons_direct.tsv`     | TF-region-gene with importance, rho, triplet rank (direct annotations) |
| `5.analysis/tsv/eRegulons_extended.tsv`   | same, extended (orthology / motif-similarity) annotations |
| `5.analysis/tsv/eRegulons_combined.tsv`   | union of direct + extended |
| `5.analysis/tsv/TF_summary.tsv`           | per-TF target counts and mean importance |
| `5.analysis/tsv/AUC_gene_per_cell.tsv`    | eRegulon (gene-based) AUC per cell |
| `5.analysis/tsv/AUC_region_per_cell.tsv`  | eRegulon (region-based) AUC per cell |
| `5.analysis/tsv/RSS_per_celltype.tsv`     | regulon specificity scores |
| `5.analysis/tsv/eRegulons_per_celltype.tsv` | top-N eRegulons per cell type |
| `5.analysis/plots/01_umap_celltype.{pdf,png}`         | UMAP coloured by cell type |
| `5.analysis/plots/02_umap_eRegulon_<TF>.{pdf,png}`    | UMAP per top eRegulon AUC |
| `5.analysis/plots/03_rss_per_celltype.{pdf,png}`      | RSS rank plot |
| `5.analysis/plots/04_heatmap_dotplot_direct.{pdf,png}` | gene/region AUC heatmap-dotplot |
| `5.analysis/plots/05_heatmap_dotplot_extended.{pdf,png}` | same, extended cistromes |
| `5.analysis/plots/06_TF_target_count.{pdf,png}`       | top TFs by # target genes |
| `5.analysis/plots/07_TF_importance_distribution.{pdf,png}` | TF→gene importance density |
| `5.analysis/plots/08_eGRN_network_top<N>.{pdf,png}`    | TF-target network |

## Pipeline stages

| # | Step              | Inputs                              | Outputs |
| - | ----------------- | ----------------------------------- | ------- |
| 1 | `seurat_to_anndata` | `input.seurat_rds`                 | mtx/tsv artifacts |
| 2 | `build_anndata`     | mtx artifacts                      | `2.anndata/rna.h5ad` |
| 3 | `create_cistopic`   | ATAC mtx                           | `3.cistopic/cistopic_obj.pkl` |
| 4 | `topic_modeling`    | cistopic_obj.pkl                   | with-topics pkl |
| 5 | `region_sets`       | with-topics pkl                    | `3.cistopic/region_sets/` |
| 6 | `prepare_gex_acc`   | adata + cistopic                   | `4.grn/ACC_GEX.h5mu` |
| 7 | `genome_annot`      | biomart                            | `4.grn/genome_annotation.tsv`, `chromsizes.tsv` |
| 8 | `search_space`      | ACC_GEX + genome_annot             | `4.grn/search_space.tsv` |
| 9 | `cistarget`         | region_sets + ctx_db               | `4.grn/ctx_results.hdf5` |
| 10 | `dem`              | region_sets + dem_db               | `4.grn/dem_results.hdf5` |
| 11 | `prepare_menr`     | cistarget + dem + ACC_GEX          | `4.grn/tf_names.txt`, `cistromes_{direct,extended}.h5ad` |
| 12 | `tf_to_gene`       | ACC_GEX + tf_names                 | `4.grn/tf_to_gene_adj.tsv` |
| 13 | `region_to_gene`   | ACC_GEX + search_space             | `4.grn/region_to_gene_adj.tsv` |
| 14 | `egrn_direct`      | adjacencies + cistromes_direct     | `4.grn/eRegulons_direct.tsv` |
| 15 | `egrn_extended`    | adjacencies + cistromes_extended   | `4.grn/eRegulons_extended.tsv` |
| 16 | `aucell_direct`    | eRegulons_direct + ACC_GEX         | `4.grn/AUCell_direct.h5mu` |
| 17 | `aucell_extended`  | eRegulons_extended + ACC_GEX       | `4.grn/AUCell_extended.h5mu` |
| 18 | `scplus_mudata`    | AUCell + eRegulons + ACC_GEX       | `4.grn/scplusmdata.h5mu` |
| 19 | `postprocess_tsv`  | scplusmdata.h5mu                   | TSVs in `5.analysis/tsv/` |
| 20 | `visualize`        | scplusmdata.h5mu                   | plots in `5.analysis/plots/` |

Each step above is one rule, named `R01_*` through `R20_*`, plus `R21_report`.
Steps 6-18 are the SCENIC+ GRN inference DAG, each a discrete `scenicplus` CLI
call (via `scenicplus_06_grn_stage.py`) that formerly ran inside one opaque
inner snakemake. Flattening it also recovered the parallelism it had: cistarget
and dem run together, as do tf_to_gene and region_to_gene, the two eGRN
branches and the two AUCell branches.

Each rule declares the config keys its script reads, so editing a GRN parameter
re-runs the rules that read it and their dependents — changing
`grn.rho_threshold` re-runs from `egrn_direct` (14), leaving the expensive
`cistarget`/`dem` motif enrichment alone.

`workflow.md` describes the same steps as the bash driver ran them; its
per-step inputs, options and API calls still hold, its orchestration no longer
does. RUNBOOK section 5 has the reservations and the measured runtimes.

## Notes

- The Seurat object is assumed to be **already preprocessed** (normalized RNA,
  cell type annotation, ideally a UMAP reduction). The pipeline reuses the
  Seurat UMAP for downstream plots when present.
- Peak counts only — fragment files are not required. Topic modeling uses the
  matrix-only path through pycisTopic.
- DARs are computed 1-vs-rest per `celltype_column`. Topic-binarized region
  sets (Otsu + top-3k) are also generated and used for motif enrichment.
