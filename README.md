# SCENIC+ pipeline (Seurat -> TF-gene network)

A 9-step pipeline that runs SCENIC+ end-to-end on a single Seurat `.rds`
multiome object (RNA + peaks/ATAC) and produces:

- `results/tables/` — comprehensive TSVs of TF-region-gene links, AUC, RSS, TF summary
- `results/plots/` — interpretable PDF + PNG plots (one per file)
- `results/scplus_out/scplusmdata.h5mu` — the canonical SCENIC+ MuData

The pipeline is driven by a plain bash master script (`scenicplus_run_pipeline.sh`) that
walks the steps sequentially, skipping any step whose output is up-to-date
relative to the input artifacts and the relevant slice of `config.yaml`. (The
SCENIC+ inner pipeline at step 07 still uses its own snakemake under the hood.)

## Layout

The pipeline is meant to live in **one central location** (referenced via the
`SCENICPLUS_PATH` environment variable) and is invoked from per-analysis
directories. Each analysis directory contains its own `config/config.yaml`
and outputs; the pipeline code itself is shared and never edited per-run.

All executables share a `scenicplus_` prefix to avoid collisions on `PATH`.

```
$SCENICPLUS_PATH/                          # central pipeline (set by the user)
├── config/config.yaml                     # template — copied into each analysis dir
└── scripts/                               # all executables ($PATH)
    ├── scenicplus_init.sh                 # bootstrap a new analysis directory
    ├── scenicplus_check.sh                # preflight: required Python + R packages
    ├── scenicplus_helper.py               # YAML-slice + sha256 helper
    ├── scenicplus_run_pipeline.sh         # the 9-step master driver
    ├── scenicplus_run_workstation.sh      # thin wrapper for workstations
    ├── scenicplus_run_lsf.sh              # bsub wrapper for LSF clusters
    ├── scenicplus_01_seurat_to_anndata.R  # step scripts (called by the
    ├── scenicplus_02_build_anndata.py     #  driver — not invoked by users)
    ├── scenicplus_03_create_cistopic.py
    ├── scenicplus_04_topic_modeling.py
    ├── scenicplus_05_region_sets.py
    ├── scenicplus_06_init_inner.py
    ├── scenicplus_07_postprocess_tsv.py
    └── scenicplus_08_visualize.py

<analysis-dir>/                            # one of these per dataset
├── config/config.yaml                     # the only file the user edits
├── results/                               # outputs (created by pipeline)
└── logs/                                  # per-step logs
```

## Setup

1. Place this repository at the location you want to be the central
   `SCENICPLUS_PATH`, then export it and add the launcher scripts to `PATH`.
   Add this to your shell rc (`~/.bashrc`, `~/.zshrc`):

   ```bash
   export SCENICPLUS_PATH=/abs/path/to/ScenicPlus
   export PATH="$SCENICPLUS_PATH/scripts:$PATH"
   ```

2. Create the conda environment (See section #4 for macos). `environment.yml` at the repo root pins
   Python 3.11.8, R + Seurat / Matrix / optparse, and a single pip line for
   `scenicplus` (which pulls `pycisTopic` and `pycistarget` in transitively
   along with the rest of the Python stack — scanpy / anndata / mudata /
   pyranges / pybedtools / …):

   ```bash
   conda env create -f "$SCENICPLUS_PATH/environment.yml"
   conda activate scenicplus_limlab
   ```

   If your cluster has no outbound internet, clone the SCENIC+ repo and
   `pip install -e /path/to/scenicplus` from a local checkout instead.

   You can also use a pre-existing env — just point the launcher at it via
   `SCENICPLUS_ENV=<env_name>` and it will `conda activate` for you.

   Verify the env at any time:

   ```bash
   scenicplus_check.sh                              # check current shell
   SCENICPLUS_ENV=scenicplus scenicplus_check.sh    # check a specific env
   ```

   The launcher scripts run this preflight automatically before kicking off
   Snakemake. Set `SCENICPLUS_SKIP_CHECK=1` to skip it (useful when the
   cluster compute node has a different env from the submitting host).

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

# 4a. Workstation
SCENICPLUS_ENV=scenicplus scenicplus_run_workstation.sh

# 4b. LSF (single bsub'd driver job)
SCENICPLUS_ENV=scenicplus LSF_QUEUE=long LSF_CORES=16 \
    scenicplus_run_lsf.sh
```

Both launchers ultimately call `scenicplus_run_pipeline.sh`, which always runs from the
**current working directory** (your analysis dir) and reads `./config/config.yaml`.

### Re-running and staleness

Each step writes a sentinel output plus a `<sentinel>.cfgsha` sidecar
containing sha256 of the config keys that step actually reads. A step is
**skipped** iff:

1. Its sentinel exists, **and**
2. Its `.cfgsha` matches the current config slice, **and**
3. No upstream sentinel is newer than its own.

Once a step has to run, every later step is force-run in the same invocation.
Outputs are written to `<sentinel>.partial` first and atomically renamed, so
crashes never leave a half-written file that looks fresh.

CLI flags (forwarded through `scenicplus_run_workstation.sh` / `scenicplus_run_lsf.sh`):

```bash
scenicplus_run_pipeline.sh             # run everything that's stale
scenicplus_run_pipeline.sh --dry-run   # show which steps would run
scenicplus_run_pipeline.sh --from 4    # force re-run from step 4 onward
scenicplus_run_pipeline.sh --only 4    # run only step 4
scenicplus_run_pipeline.sh --force     # force re-run everything
```

Because the staleness check is per-step on a config *slice*, editing
`cistopic.dar_log2fc_thr` only invalidates step 05 onward; LDA outputs from
step 04 stay cached.

`scenicplus_run_pipeline.sh` validates that `config.yaml` has `Pipeline: "ScenicPlus"`;
if that line is missing it aborts immediately, so an unrelated pipeline's
config cannot be fed in by accident.

## Inputs (set in `config/config.yaml`)

1. `input.seurat_rds` — Seurat `.rds` with:
   - `RNA` assay containing **normalized** expression in `data` and raw counts in `counts`,
   - `peaks` assay containing peak counts (regions named `chr-start-end` or `chr:start-end`),
   - cell-type column (`input.celltype_column`) in `@meta.data`.
   - (optional) UMAP/PCA reductions — these are reused if present.
2. `input.celltype_scope` (optional) — a list of cell-type values from
   `celltype_column`. When set, step 01 drops all other cells before exporting,
   so topic modeling, DARs, metacells, and RSS all operate on the focused
   universe (good for sensitivity when you care about a specific lineage).
   Leave empty / omit to use every cell.
3. `input.ctx_db`, `input.dem_db` — cisTarget feather databases.
4. `input.motif_annotations` — motif-to-TF table.

## Outputs

| Path | Description |
| ---- | ----------- |
| `results/tables/eRegulons_direct.tsv`     | TF-region-gene with importance, rho, triplet rank (direct annotations) |
| `results/tables/eRegulons_extended.tsv`   | same, extended (orthology / motif-similarity) annotations |
| `results/tables/eRegulons_combined.tsv`   | union of direct + extended |
| `results/tables/TF_summary.tsv`           | per-TF target counts and mean importance |
| `results/tables/AUC_gene_per_cell.tsv`    | eRegulon (gene-based) AUC per cell |
| `results/tables/AUC_region_per_cell.tsv`  | eRegulon (region-based) AUC per cell |
| `results/tables/RSS_per_celltype.tsv`     | regulon specificity scores |
| `results/tables/eRegulons_per_celltype.tsv` | top-N eRegulons per cell type |
| `results/plots/01_umap_celltype.{pdf,png}`         | UMAP coloured by cell type |
| `results/plots/02_umap_eRegulon_<TF>.{pdf,png}`    | UMAP per top eRegulon AUC |
| `results/plots/03_rss_per_celltype.{pdf,png}`      | RSS rank plot |
| `results/plots/04_heatmap_dotplot_direct.{pdf,png}` | gene/region AUC heatmap-dotplot |
| `results/plots/05_heatmap_dotplot_extended.{pdf,png}` | same, extended cistromes |
| `results/plots/06_TF_target_count.{pdf,png}`       | top TFs by # target genes |
| `results/plots/07_TF_importance_distribution.{pdf,png}` | TF→gene importance density |
| `results/plots/08_eGRN_network_top<N>.{pdf,png}`    | TF-target network |

## Pipeline stages

| # | Step              | Inputs                              | Outputs |
| - | ----------------- | ----------------------------------- | ------- |
| 1 | `seurat_to_anndata` | `input.seurat_rds`                 | mtx/tsv artifacts |
| 2 | `build_anndata`     | mtx artifacts                      | `interim/rna.h5ad` |
| 3 | `create_cistopic`   | ATAC mtx                           | `interim/cistopic_obj.pkl` |
| 4 | `topic_modeling`    | cistopic_obj.pkl                   | with-topics pkl |
| 5 | `region_sets`       | with-topics pkl                    | `interim/region_sets/` |
| 6 | `init_scenicplus`   | adata + cistopic + region_sets     | `scplus_pipeline/Snakemake/config/config.yaml` |
| 7 | `run_scenicplus`    | scenicplus snakemake (inner)       | `scplus_out/scplusmdata.h5mu` |
| 8 | `postprocess_tsv`   | scplusmdata.h5mu                   | TSVs in `results/tables/` |
| 9 | `visualize`         | scplusmdata.h5mu                   | plots in `results/plots/` |

See `workflow.md` for per-step inputs/outputs, key options, and the core
SCENIC+ / pycisTopic API calls each step makes.

## Notes

- The Seurat object is assumed to be **already preprocessed** (normalized RNA,
  cell type annotation, ideally a UMAP reduction). The pipeline reuses the
  Seurat UMAP for downstream plots when present.
- Peak counts only — fragment files are not required. Topic modeling uses the
  matrix-only path through pycisTopic.
- DARs are computed 1-vs-rest per `celltype_column`. Topic-binarized region
  sets (Otsu + top-3k) are also generated and used for motif enrichment.
