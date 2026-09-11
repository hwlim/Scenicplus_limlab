# =============================================================================
# R19-R20: the GRN result, turned into tables and figures.
#
# Both read only `scplusmdata.h5mu`, so they are independent of each other and
# run together. R20 additionally reads step 01's embedding, which is why the
# layout a figure is drawn on does not have to survive nine intermediate files
# -- see scenicplus_08_visualize.py.
#
# WHICH OUTPUTS ARE DECLARED, and why it is not all of them. The plan's rule is
# to declare every output including side effects, because a stage that exits 0
# without its outputs is the failure Snakemake catches for free. Two kinds resist
# that:
#
#   * DATA-DEPENDENT NAMES. `02_umap_eRegulon_<name>` is one figure per top
#     eRegulon, and `RSS_per_celltype.tsv` depends on the cell-type column being
#     present in the MuData. Names that are not knowable at DAG-build time
#     cannot be declared.
#   * UPSTREAM-OPTIONAL. `03_rss_per_celltype` is computed inside a try/except
#     that prints and continues, so a sparse cell type can legitimately leave it
#     out. Declaring it would convert a warning into a failed run.
#
# Everything else IS declared, which is what makes the plotnine regression
# (04/05 are the heatmap-dotplots that crashed on `.savefig`) a MissingOutput
# failure rather than a run that finishes without them.
# =============================================================================

TSV = stage_path("analysis", "tsv")
PLOTS = stage_path("analysis", "plots")

# `08_eGRN_network_top<N>` carries the config value in its filename, so the
# declared name follows the config rather than being hardcoded.
_NET_TOP = config["visualization"]["network_top_n_TFs"]


def _fig(stem):
    """Both formats. One plot per file, both PDF and PNG, is the project's
    stated output contract -- see CLAUDE.md."""
    return [os.path.join(PLOTS, f"{stem}.{ext}") for ext in ("pdf", "png")]


rule R19_postprocess_tsv:
    """The eRegulon tables a collaborator can read without python.

    Four outputs are unconditional. The AUC matrices and the RSS tables are
    gated on the AUC modalities being non-empty and the cell-type column being
    found, so they are written but not declared.
    """
    input:
        mdata=stage_path("grn", "scplusmdata.h5mu"),
        script=script_path("scenicplus_07_postprocess_tsv.py"),
    output:
        direct=os.path.join(TSV, "eRegulons_direct.tsv"),
        extended=os.path.join(TSV, "eRegulons_extended.tsv"),
        combined=os.path.join(TSV, "eRegulons_combined.tsv"),
        tf_summary=os.path.join(TSV, "TF_summary.tsv"),
    log:
        log_path("R19_postprocess_tsv")
    threads: 1
    resources:
        mem_mb=mem("R19_postprocess_tsv"),
        runtime=rt("R19_postprocess_tsv"),
    params:
        cfg=cfg_params("input.celltype_column",
                       "visualization.top_n_eRegulons_per_celltype"),
        config_file=CONFIG_FILE,
        out_dir=TSV,
    shell:
        "python {input.script}"
        " --scplus_mdata {input.mdata:q}"
        " --config {params.config_file:q}"
        " --out_dir {params.out_dir:q}"
        " 2>&1 | tee {log}"


rule R20_visualize:
    """The figures.

    Reads the embedding from step 01 rather than the MuData, because the
    eRegulon object is concatenated from the AUC modalities alone and inherits
    no layout. That also means this rule can be re-run on its own to redraw a
    finished run: `--force` on one of its outputs touches nothing upstream.
    """
    input:
        mdata=stage_path("grn", "scplusmdata.h5mu"),
        # Only when the config names a reduction; otherwise the script falls
        # back to a UMAP of eRegulon activity and needs no export.
        embedding=([os.path.join(EXPORT_DIR, f"embedding_{_REDUCTION}.tsv")]
                   if _REDUCTION else []),
        script=script_path("scenicplus_08_visualize.py"),
    output:
        *_fig("01_umap_celltype"),
        *_fig("04_heatmap_dotplot_direct"),
        *_fig("05_heatmap_dotplot_extended"),
        *_fig("06_TF_target_count"),
        *_fig("07_TF_importance_distribution"),
        *_fig(f"08_eGRN_network_top{_NET_TOP}"),
    log:
        log_path("R20_visualize")
    threads: 1
    resources:
        mem_mb=mem("R20_visualize"),
        runtime=rt("R20_visualize"),
    params:
        cfg=cfg_params("input.celltype_column", "input.reduction",
                       "visualization"),
        config_file=CONFIG_FILE,
        out_dir=PLOTS,
        embedding_dir=EXPORT_DIR,
        reduction=opt_arg("--reduction", _REDUCTION),
    shell:
        "python {input.script}"
        " --scplus_mdata {input.mdata:q}"
        " --config {params.config_file:q}"
        " --out_dir {params.out_dir:q}"
        " --embedding_dir {params.embedding_dir:q}"
        " {params.reduction}"
        " 2>&1 | tee {log}"
