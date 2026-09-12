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
- **A Snakemake workflow is the pipeline.** `Snakefile` + `rules/` define 21
  rules (`R01_*`…`R20_*` plus `R21_report`), each with a per-rule LSF
  reservation measured from a real run, and `rule all` targets `report.html`.
  The entry point is `scripts/scenicplus.run.sh`. `SnakemakePlan.md` carries
  every increment and the defects each one turned up; resume from there.
- SCENIC+'s former inner snakemake (old step 07) is flattened into native
  stages 06-18, one `scenicplus` CLI call each, dispatched by
  `scenicplus_06_grn_stage.py`. Steps 01-05 preprocess (R/pycisTopic), 06-18
  are the GRN inference DAG, and 19-20 postprocess/visualize
  (`scenicplus_07_*.py`, `scenicplus_08_*.py` — the file-name prefixes are
  historical and no longer equal the step number). Declaring each stage's real
  inputs recovered the intra-DAG parallelism the flattening had cost:
  cistarget ∥ dem, tf_to_gene ∥ region_to_gene, and both branch pairs after.
- All shipped executables share a `scenicplus_` prefix to avoid PATH collisions
  with other pipelines. `scenicplus.run.sh` is the one exception, named to read
  as the verb it is.
- **`scripts/scenicplus_run_pipeline.sh` is OBSOLETE and deliberately KEPT** —
  as a fallback and for reproducing older runs. It walks the 20 steps
  sequentially as ONE job sized for the heaviest, resumes on
  sentinel/`.cfgsha`/cascade, and produces neither a report nor a provenance
  bundle. Do not extend it, and do not point a new user at it.
- Relevant settings parameters defined in a yaml file (`config/config.yaml`).
- **Resume, in the workflow:** snakemake decides, from file times, the code, and
  the config values each rule declares via `cfg_params()` in `rules/common.smk`.
  The step scripts read the config file themselves, so without those
  declarations a config edit would change nothing — that is the whole reason
  the helper exists. **Never delete `.snakemake/`**; recover with `--forceall`,
  never `--touch`.
- **Resume, in the obsolete driver** (below, for reference only): it checks
  intermediate results before running each step. A
  step is skipped iff (a) its sentinel output exists, (b) the sha256 of the
  config keys it actually reads matches the `<sentinel>.cfgsha` sidecar
  written on the previous run, and (c) no upstream sentinel is newer than its
  own. Once a step has to run, every later step is force-run in the same
  invocation. Per-step outputs are written to `<sentinel>.partial` and
  atomically renamed so a crash never leaves a half-written file that looks
  fresh.
- CLI: `scenicplus.run.sh [-n] [-j N] [--lsf] [-p] [-f N] [-- <snakemake args>]`.
  **`-f N` selects a RULE and its dependents, not a range of steps** — the DAG
  forks, so `-f 9` leaves R10 and R13 alone. The obsolete driver's `--from N`
  did mean every step numbered N or higher; the two are not equivalent.
- The obsolete driver and its two launchers are still present and still work:
  `scenicplus_run_pipeline.sh [--dry-run|--from N|--only N|--force]`,
  `scenicplus_run_workstation.sh` (execs it here),
  `scenicplus_run_lsf.sh` (bsubs it as ONE job sized for the heaviest GRN
  stage). That single oversized job is the reason they were retired — the
  workflow right-sizes each rule instead.
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
    Both the workflow and the driver check this and abort if it's missing or
    wrong, so an unrelated pipeline's config cannot be fed in by accident. The
    workflow additionally validates the whole config against
    `schemas/config.schema.yaml` at DAG-build time, so a mistyped key is named
    and refused instead of silently falling back to a default.
  - Assume `SCENICPLUS_PATH` is already defined and `$SCENICPLUS_PATH/scripts`
    is on `PATH`.
  - Keep `README.md` in sync with how to set up and run the pipeline.
  - Keep the layout flat — no empty folders.