# R Command-Line Script

## Goal
Build a pipeline to run a SCENIC+ to infer TF-gene networks for a given single cell multiome data.

## Note
Markdown folder contains following documents
- Markdown/API.txt: API reference
- Markdown/human_cerebellum_scRNA_pp.md: Preprocessing the scRNA-seq data
- Markdown/human_cerebellum.md: Running SCENIC+
- Markdown/human_cerebellum_ctx_db.md: Creating custom cistarget database
- Markdown/Perturbation_simulation.md: Perturbation simulation
Refer to the first three files only for technical detail and scenic+ API usages; No need to look at the custom database and perturbation simulation at the moment.

## Input
A seurat object file (rds) containing cell type annotation and normalized gene expression where
RNA assay is stored in RNA assay
ATAC assay is stored in peaks assay

Optionally, the user can restrict the analysis to a focused subset of cell
types by setting `input.celltype_scope` in `config.yaml` to a list of values
from `input.celltype_column`. Step 01 drops every other cell before exporting,
so topic modeling, DARs, metacells, and RSS all operate on the focused
universe (better sensitivity for lineage-specific analyses). Empty / omitted
list = use every cell.

## Output
- Comprehensive list of predicted TF-target gene including their relevant quantitative measures and statistics
- Comprehensive and interpretable visualizations of the results including intermediate results if helpful
- Plots both in pdf and png format. One plot per file.
- Comprehensive data sheet in tsv format matching the visualziation results

## Implementation note
- Sequential bash-driven implementation starting from a seurat object. The
  outer orchestrator is `scripts/scenicplus_run_pipeline.sh` (the master
  driver) which walks 9 step scripts (`scenicplus_01_*.R`,
  `scenicplus_02_*.py`, …, `scenicplus_08_*.py`) in the same `scripts/`
  folder. All shipped executables share a `scenicplus_` prefix to avoid
  PATH collisions with other pipelines. The SCENIC+ inner pipeline at step 07
  still uses its own snakemake under the hood — that's unavoidable and is
  what `snakemake` in `environment.yml` is for.
- Relevant settings parameters defined in a yaml file (`config/config.yaml`).
- The master driver checks intermediate results before running each step. A
  step is skipped iff (a) its sentinel output exists, (b) the sha256 of the
  config keys it actually reads matches the `<sentinel>.cfgsha` sidecar
  written on the previous run, and (c) no upstream sentinel is newer than its
  own. Once a step has to run, every later step is force-run in the same
  invocation. Per-step outputs are written to `<sentinel>.partial` and
  atomically renamed so a crash never leaves a half-written file that looks
  fresh.
- CLI: `scenicplus_run_pipeline.sh [--dry-run|--from N|--only N|--force]`.
- Bash launchers are thin wrappers around `scenicplus_run_pipeline.sh`:
  - `scenicplus_run_workstation.sh` — `exec`s the driver in the current shell.
  - `scenicplus_run_lsf.sh` — bsubs the driver as a single LSF job
    (cores/mem/walltime via env vars). The bsub'd job is sized for step 07,
    which runs its own internal snakemake with that many cores.
- Assume that the seurat object already contains cluster, cell type
  annotation, and UMAP projection. Incorporate those as much as possible.
- Separate implementation and per-execution configuration:
  - The pipeline lives in a central location pointed to by the env var
    `SCENICPLUS_PATH`. `config/` (template) and `scripts/` (all executables —
    launchers + step scripts) are in that central location.
  - `config/` in the central location holds a pre-filled template.
  - `scenicplus_init.sh` initiates an analysis directory by copying
    `config.yaml` from the central location into the analysis folder.
  - Users only access and edit `config.yaml` under the analysis directory.
  - `config.yaml` must contain `Pipeline: "ScenicPlus"`.
    `scenicplus_run_pipeline.sh` checks this and aborts if it's missing or
    wrong, so an unrelated pipeline's config cannot be fed in by accident.
  - Assume `SCENICPLUS_PATH` is already defined and `$SCENICPLUS_PATH/scripts`
    is on `PATH`.
  - Keep `README.md` in sync with how to set up and run the pipeline.
  - Keep the layout flat — no empty folders.