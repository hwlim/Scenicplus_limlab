#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Master driver for the SCENIC+ pipeline.
#
# Replaces snakemake. Walks 9 sequential steps and skips a step iff:
#   1) its sentinel output exists, AND
#   2) its sentinel.cfgsha matches sha256 of the relevant config slice, AND
#   3) the upstream sentinel is not newer than this step's sentinel.
# Once any step runs, every later step is force-run in the same invocation
# (cascade). Each step writes to <sentinel>.partial first, then renames.
#
# Usage:
#   scenicplus_run_pipeline.sh                    # run everything that's stale
#   scenicplus_run_pipeline.sh --dry-run          # show plan without executing
#   scenicplus_run_pipeline.sh --from N           # force re-run from step N onward
#   scenicplus_run_pipeline.sh --only N           # run only step N
#   scenicplus_run_pipeline.sh --force            # force re-run everything
# -----------------------------------------------------------------------------
set -euo pipefail

if [[ -z "${SCENICPLUS_PATH:-}" ]]; then
    echo "[scenicplus_run_pipeline] ERROR: SCENICPLUS_PATH is not set." >&2
    exit 1
fi

CONFIG="$PWD/config/config.yaml"
SCRIPT_DIR="$SCENICPLUS_PATH/scripts"
HELPER="$SCRIPT_DIR/scenicplus_helper.py"

if [[ ! -f "$CONFIG" ]]; then
    echo "[scenicplus_run_pipeline] ERROR: $CONFIG not found." >&2
    echo "  Run scenicplus_init.sh in this directory first." >&2
    exit 1
fi

# ---- CLI -------------------------------------------------------------------
DRY_RUN=0; FORCE_FROM=0; ONLY_STEP=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) DRY_RUN=1; shift;;
        --force)   FORCE_FROM=1; shift;;
        --from)    FORCE_FROM="$2"; shift 2;;
        --only)    ONLY_STEP="$2"; shift 2;;
        -h|--help) sed -n '2,20p' "$0"; exit 0;;
        *) echo "[scenicplus_run_pipeline] unknown arg: $1" >&2; exit 2;;
    esac
done

# ---- Optional conda activation --------------------------------------------
if [[ -n "${SCENICPLUS_ENV:-}" ]]; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$SCENICPLUS_ENV"
fi

if [[ "${SCENICPLUS_SKIP_CHECK:-0}" != "1" ]]; then
    "$SCENICPLUS_PATH/scripts/scenicplus_check.sh"
fi

mkdir -p logs results tmp

# ---- Validate Pipeline identity --------------------------------------------
PIPELINE_ID=$(python "$HELPER" get "$CONFIG" Pipeline || true)
if [[ "$PIPELINE_ID" != "ScenicPlus" ]]; then
    echo "[scenicplus_run_pipeline] ERROR: config.yaml has Pipeline=$PIPELINE_ID (expected ScenicPlus)." >&2
    echo "  Re-run scenicplus_init.sh to refresh the template." >&2
    exit 1
fi

# ---- Resolve paths from config ---------------------------------------------
ROOT="$(python "$HELPER" get "$CONFIG" output.root)"
TMP_DIR="$(python "$HELPER" get "$CONFIG" output.tmp)"
CT_COL="$(python "$HELPER" get "$CONFIG" input.celltype_column)"
N_CPU="$(python "$HELPER" get "$CONFIG" resources.n_cpu)"
[[ "$ROOT" = /* ]] || ROOT="$PWD/$ROOT"
[[ "$TMP_DIR" = /* ]] || TMP_DIR="$PWD/$TMP_DIR"

INTERIM="$ROOT/interim"
SCPLUS_PIPELINE="$ROOT/scplus_pipeline"
SCPLUS_OUT="$ROOT/scplus_out"
TSV_DIR="$ROOT/tables"
PLOT_DIR="$ROOT/plots"

# ---- Skip-rule logic -------------------------------------------------------
CASCADE=0
declare -a UPSTREAM_SENTINELS=()

step_is_fresh() {
    # $1 = sentinel path, $2 = expected hash
    local sentinel="$1" want="$2" sidecar
    [[ -e "$sentinel" ]] || return 1
    sidecar="${sentinel}.cfgsha"
    [[ -f "$sidecar" ]] || return 1
    [[ "$(cat "$sidecar")" == "$want" ]] || return 1
    # Re-run if any upstream sentinel is newer than this one.
    local up
    for up in "${UPSTREAM_SENTINELS[@]}"; do
        if [[ "$up" -nt "$sentinel" ]]; then
            return 1
        fi
    done
    return 0
}

run_step() {
    # $1 step_id  $2 step_name  $3 cfg_keys_csv  $4 sentinel  $5 cmd
    local id="$1" name="$2" keys="$3" sentinel="$4" cmd="$5"
    local want forced=0
    want="$(python "$HELPER" hash "$CONFIG" "$keys")"

    if [[ "$ONLY_STEP" -ne 0 && "$ONLY_STEP" -ne "$id" ]]; then
        echo "[skip-only]  $id $name"
        UPSTREAM_SENTINELS+=("$sentinel")
        return
    fi
    if [[ "$FORCE_FROM" -ne 0 && "$id" -ge "$FORCE_FROM" ]]; then forced=1; fi
    if [[ "$CASCADE" -eq 1 ]]; then forced=1; fi

    if [[ "$forced" -eq 0 ]] && step_is_fresh "$sentinel" "$want"; then
        echo "[skip]       $id $name  ($(basename "$sentinel"))"
        UPSTREAM_SENTINELS+=("$sentinel")
        return
    fi

    if [[ "$DRY_RUN" -eq 1 ]]; then
        echo "[would-run]  $id $name"
        UPSTREAM_SENTINELS+=("$sentinel")
        CASCADE=1
        return
    fi

    echo "[run]        $id $name"
    rm -f "$sentinel" "${sentinel}.cfgsha"
    mkdir -p "$(dirname "$sentinel")"
    # Subshell so a `cd` inside the cmd (e.g. step 07) doesn't leak into the
    # parent driver and skew $PWD for later steps.
    ( eval "$cmd" )
    if [[ ! -e "$sentinel" ]]; then
        echo "[scenicplus_run_pipeline] ERROR: step $id $name finished but $sentinel was not created." >&2
        exit 1
    fi
    echo -n "$want" > "${sentinel}.cfgsha"
    UPSTREAM_SENTINELS+=("$sentinel")
    CASCADE=1
}

# ---- Step definitions ------------------------------------------------------
S01="$INTERIM/seurat_export/summary.txt"
S02="$INTERIM/rna.h5ad"
S03="$INTERIM/cistopic_obj.pkl"
S04="$INTERIM/cistopic_obj_with_topics.pkl"
S05="$INTERIM/region_sets/.done"
S06="$SCPLUS_PIPELINE/Snakemake/config/config.yaml"
S07="$SCPLUS_OUT/scplusmdata.h5mu"
S08="$TSV_DIR/eRegulons_combined.tsv"
S09="$PLOT_DIR/01_umap_celltype.pdf"

echo "[scenicplus_run_pipeline] N_CPU=$N_CPU  DRY_RUN=$DRY_RUN  FORCE_FROM=$FORCE_FROM  ONLY=$ONLY_STEP"

CELLTYPE_SCOPE="$(python "$HELPER" getcsv "$CONFIG" input.celltype_scope)"

run_step 1 seurat_to_anndata \
    "input.seurat_rds,input.celltype_column,input.celltype_scope" \
    "$S01" "
        Rscript '$SCRIPT_DIR/scenicplus_01_seurat_to_anndata.R' \
            --rds '$(python "$HELPER" get "$CONFIG" input.seurat_rds)' \
            --celltype_col '$CT_COL' \
            --celltype_scope '$CELLTYPE_SCOPE' \
            --out_dir '$INTERIM/seurat_export' \
            > 'logs/01_seurat_to_anndata.log' 2>&1
    "

run_step 2 build_anndata \
    "input.celltype_column" \
    "$S02" "
        python '$SCRIPT_DIR/scenicplus_02_build_anndata.py' \
            --in_dir '$INTERIM/seurat_export' \
            --out_h5ad '$S02' \
            --celltype_col '$CT_COL' \
            > 'logs/02_build_anndata.log' 2>&1
    "

run_step 3 create_cistopic \
    "" \
    "$S03" "
        python '$SCRIPT_DIR/scenicplus_03_create_cistopic.py' \
            --in_dir '$INTERIM/seurat_export' \
            --out_pkl '$S03' \
            > 'logs/03_create_cistopic.log' 2>&1
    "

run_step 4 topic_modeling \
    "cistopic.n_topics,cistopic.n_iter,cistopic.alpha,cistopic.alpha_by_topic,cistopic.eta,cistopic.eta_by_topic,cistopic.random_state,resources.seed" \
    "$S04" "
        python '$SCRIPT_DIR/scenicplus_04_topic_modeling.py' \
            --in_pkl '$S03' \
            --out_pkl '$S04' \
            --config '$CONFIG' \
            --tmp_dir '$TMP_DIR/lda' \
            --n_cpu '$N_CPU' \
            > 'logs/04_topic_modeling.log' 2>&1
    "

run_step 5 region_sets \
    "input.celltype_column,cistopic.dar_adjpval_thr,cistopic.dar_log2fc_thr" \
    "$S05" "
        python '$SCRIPT_DIR/scenicplus_05_region_sets.py' \
            --in_pkl '$S04' \
            --out_dir '$INTERIM/region_sets' \
            --config '$CONFIG' \
            --n_cpu '$N_CPU' \
            > 'logs/05_region_sets.log' 2>&1
        touch '$S05'
    "

run_step 6 init_scenicplus \
    "input.species,input.assembly,input.biomart_host,input.ctx_db,input.dem_db,input.motif_annotations,scenicplus,grn,resources.n_cpu,resources.seed,output.tmp" \
    "$S06" "
        python '$SCRIPT_DIR/scenicplus_06_init_inner.py' \
            --config '$CONFIG' \
            --out_dir '$SCPLUS_PIPELINE' \
            --cistopic_obj '$S04' \
            --adata '$S02' \
            --region_sets '$INTERIM/region_sets' \
            --scplus_out '$SCPLUS_OUT' \
            > 'logs/06_init_scenicplus.log' 2>&1
    "

run_step 7 run_scenicplus \
    "resources.n_cpu" \
    "$S07" "
        cd '$SCPLUS_PIPELINE/Snakemake' && \
        snakemake --cores '$N_CPU' --rerun-incomplete \
            > 'logs/07_run_scenicplus.log' 2>&1
    "

run_step 8 postprocess_tsv \
    "input.celltype_column,visualization.top_n_eRegulons_per_celltype" \
    "$S08" "
        python '$SCRIPT_DIR/scenicplus_07_postprocess_tsv.py' \
            --scplus_mdata '$S07' \
            --config '$CONFIG' \
            --out_dir '$TSV_DIR' \
            > 'logs/08_postprocess_tsv.log' 2>&1
    "

run_step 9 visualize \
    "input.celltype_column,visualization" \
    "$S09" "
        python '$SCRIPT_DIR/scenicplus_08_visualize.py' \
            --scplus_mdata '$S07' \
            --config '$CONFIG' \
            --out_dir '$PLOT_DIR' \
            > 'logs/09_visualize.log' 2>&1
    "

echo "[scenicplus_run_pipeline] done."
