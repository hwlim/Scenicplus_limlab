#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Master driver for the SCENIC+ pipeline.
#
# Replaces snakemake (including SCENIC+'s former inner snakemake, now flattened
# into native stages 06-18). Walks 20 sequential steps and skips a step iff:
#   1) its sentinel output exists, AND
#   2) its sentinel.cfgsha matches sha256 of the relevant config slice, AND
#   3) the upstream sentinel is not newer than this step's sentinel.
# Once any step runs, every later step is force-run in the same invocation
# (cascade). Each step writes to <sentinel>.partial first, then renames.
#
# NOTE: flattening the inner snakemake into serial stages means the intra-DAG
# parallelism SCENIC+'s snakemake exploited (e.g. cistarget || dem) is now
# sequential; each stage is still multi-threaded internally via resources.n_cpu.
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

## ---- Optional conda activation --------------------------------------------
#if [[ -n "${SCENICPLUS_ENV:-}" ]]; then
#    # shellcheck disable=SC1091
#    source "$(conda info --base)/etc/profile.d/conda.sh"
#    conda activate "$SCENICPLUS_ENV"
#fi
#
# if [[ "${SCENICPLUS_SKIP_CHECK:-0}" != "1" ]]; then
#     "$SCENICPLUS_PATH/scripts/scenicplus_check.sh"
# fi

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
    # Subshell so anything a cmd does to shell state doesn't leak into the
    # parent driver and skew $PWD / env for later steps.
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
# Steps 06-18 are the flattened SCENIC+ GRN inference DAG (formerly one opaque
# inner snakemake). Each stage is a discrete `scenicplus` CLI call via
# scenicplus_06_grn_stage.py; its sentinel is that stage's primary output.
# Sentinels follow the DAG's topological order, so no stage depends on a later
# one and the driver's linear upstream/cascade logic stays correct (upstream
# tracking is intentionally conservative: e.g. dem is treated as downstream of
# genome_annot even in the unbalanced branch, which over-runs but never
# under-runs).
S06="$SCPLUS_OUT/ACC_GEX.h5mu"              # prepare_gex_acc
S07="$SCPLUS_OUT/genome_annotation.tsv"     # genome_annot (+chromsizes)
S08="$SCPLUS_OUT/search_space.tsv"          # search_space
S09="$SCPLUS_OUT/ctx_results.hdf5"          # cistarget (+html)
S10="$SCPLUS_OUT/dem_results.hdf5"          # dem (+html)
S11="$SCPLUS_OUT/cistromes_direct.h5ad"     # prepare_menr (+tf_names, cistromes_extended)
S12="$SCPLUS_OUT/tf_to_gene_adj.tsv"        # tf_to_gene
S13="$SCPLUS_OUT/region_to_gene_adj.tsv"    # region_to_gene
S14="$SCPLUS_OUT/eRegulons_direct.tsv"      # egrn_direct
S15="$SCPLUS_OUT/eRegulons_extended.tsv"    # egrn_extended
S16="$SCPLUS_OUT/AUCell_direct.h5mu"        # aucell_direct
S17="$SCPLUS_OUT/AUCell_extended.h5mu"      # aucell_extended
S18="$SCPLUS_OUT/scplusmdata.h5mu"          # scplus_mudata (final GRN output)
S19="$TSV_DIR/eRegulons_combined.tsv"       # postprocess_tsv
S20="$PLOT_DIR/01_umap_celltype.pdf"        # visualize
GRN_STAGE="$SCRIPT_DIR/scenicplus_06_grn_stage.py"

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

run_step 6 prepare_gex_acc \
    "scenicplus.is_multiome,scenicplus.bc_transform_func" \
    "$S06" "
        python '$GRN_STAGE' --stage prepare_gex_acc \
            --config '$CONFIG' --scplus_out '$SCPLUS_OUT' \
            --cistopic_obj '$S04' --adata '$S02' \
            > 'logs/06_prepare_gex_acc.log' 2>&1
    "

run_step 7 genome_annot \
    "input.species,input.biomart_host" \
    "$S07" "
        python '$GRN_STAGE' --stage genome_annot \
            --config '$CONFIG' --scplus_out '$SCPLUS_OUT' \
            > 'logs/07_genome_annot.log' 2>&1
    "

run_step 8 search_space \
    "scenicplus.search_space_upstream,scenicplus.search_space_downstream,scenicplus.search_space_extend_tss" \
    "$S08" "
        python '$GRN_STAGE' --stage search_space \
            --config '$CONFIG' --scplus_out '$SCPLUS_OUT' \
            > 'logs/08_search_space.log' 2>&1
    "

run_step 9 cistarget \
    "input.ctx_db,input.motif_annotations,input.species,scenicplus.fraction_overlap_w_ctx_database,scenicplus.ctx_auc_threshold,scenicplus.ctx_nes_threshold,scenicplus.ctx_rank_threshold,scenicplus.annotation_version,scenicplus.motif_similarity_fdr,scenicplus.orthologous_identity_threshold,scenicplus.annotations_to_use,resources.n_cpu,output.tmp" \
    "$S09" "
        python '$GRN_STAGE' --stage cistarget \
            --config '$CONFIG' --scplus_out '$SCPLUS_OUT' \
            --region_sets '$INTERIM/region_sets' \
            > 'logs/09_cistarget.log' 2>&1
    "

run_step 10 dem \
    "input.dem_db,input.motif_annotations,input.species,scenicplus.fraction_overlap_w_dem_database,scenicplus.dem_max_bg_regions,scenicplus.dem_balance_number_of_promoters,scenicplus.dem_promoter_space,scenicplus.dem_adj_pval_thr,scenicplus.dem_log2fc_thr,scenicplus.dem_mean_fg_thr,scenicplus.dem_motif_hit_thr,scenicplus.annotation_version,scenicplus.motif_similarity_fdr,scenicplus.orthologous_identity_threshold,scenicplus.annotations_to_use,resources.n_cpu,resources.seed,output.tmp" \
    "$S10" "
        python '$GRN_STAGE' --stage dem \
            --config '$CONFIG' --scplus_out '$SCPLUS_OUT' \
            --region_sets '$INTERIM/region_sets' \
            > 'logs/10_dem.log' 2>&1
    "

run_step 11 prepare_menr \
    "scenicplus.direct_annotation,scenicplus.extended_annotation" \
    "$S11" "
        python '$GRN_STAGE' --stage prepare_menr \
            --config '$CONFIG' --scplus_out '$SCPLUS_OUT' \
            > 'logs/11_prepare_menr.log' 2>&1
    "

run_step 12 tf_to_gene \
    "grn.tf_to_gene_importance_method,resources.n_cpu,resources.seed,output.tmp" \
    "$S12" "
        python '$GRN_STAGE' --stage tf_to_gene \
            --config '$CONFIG' --scplus_out '$SCPLUS_OUT' \
            > 'logs/12_tf_to_gene.log' 2>&1
    "

run_step 13 region_to_gene \
    "grn.region_to_gene_importance_method,grn.region_to_gene_correlation_method,resources.n_cpu,output.tmp" \
    "$S13" "
        python '$GRN_STAGE' --stage region_to_gene \
            --config '$CONFIG' --scplus_out '$SCPLUS_OUT' \
            > 'logs/13_region_to_gene.log' 2>&1
    "

run_step 14 egrn_direct \
    "grn.order_regions_to_genes_by,grn.order_TFs_to_genes_by,grn.gsea_n_perm,grn.quantile_thresholds_region_to_gene,grn.top_n_regionTogenes_per_gene,grn.top_n_regionTogenes_per_region,grn.min_regions_per_gene,grn.rho_threshold,grn.min_target_genes,input.ctx_db,resources.n_cpu,output.tmp" \
    "$S14" "
        python '$GRN_STAGE' --stage egrn_direct \
            --config '$CONFIG' --scplus_out '$SCPLUS_OUT' \
            > 'logs/14_egrn_direct.log' 2>&1
    "

run_step 15 egrn_extended \
    "grn.order_regions_to_genes_by,grn.order_TFs_to_genes_by,grn.gsea_n_perm,grn.quantile_thresholds_region_to_gene,grn.top_n_regionTogenes_per_gene,grn.top_n_regionTogenes_per_region,grn.min_regions_per_gene,grn.rho_threshold,grn.min_target_genes,input.ctx_db,resources.n_cpu,output.tmp" \
    "$S15" "
        python '$GRN_STAGE' --stage egrn_extended \
            --config '$CONFIG' --scplus_out '$SCPLUS_OUT' \
            > 'logs/15_egrn_extended.log' 2>&1
    "

run_step 16 aucell_direct \
    "resources.n_cpu" \
    "$S16" "
        python '$GRN_STAGE' --stage aucell_direct \
            --config '$CONFIG' --scplus_out '$SCPLUS_OUT' \
            > 'logs/16_aucell_direct.log' 2>&1
    "

run_step 17 aucell_extended \
    "resources.n_cpu" \
    "$S17" "
        python '$GRN_STAGE' --stage aucell_extended \
            --config '$CONFIG' --scplus_out '$SCPLUS_OUT' \
            > 'logs/17_aucell_extended.log' 2>&1
    "

run_step 18 scplus_mudata \
    "" \
    "$S18" "
        python '$GRN_STAGE' --stage scplus_mudata \
            --config '$CONFIG' --scplus_out '$SCPLUS_OUT' \
            > 'logs/18_scplus_mudata.log' 2>&1
    "

run_step 19 postprocess_tsv \
    "input.celltype_column,visualization.top_n_eRegulons_per_celltype" \
    "$S19" "
        python '$SCRIPT_DIR/scenicplus_07_postprocess_tsv.py' \
            --scplus_mdata '$S18' \
            --config '$CONFIG' \
            --out_dir '$TSV_DIR' \
            > 'logs/19_postprocess_tsv.log' 2>&1
    "

run_step 20 visualize \
    "input.celltype_column,visualization" \
    "$S20" "
        python '$SCRIPT_DIR/scenicplus_08_visualize.py' \
            --scplus_mdata '$S18' \
            --config '$CONFIG' \
            --out_dir '$PLOT_DIR' \
            > 'logs/20_visualize.log' 2>&1
    "

echo "[scenicplus_run_pipeline] done."
