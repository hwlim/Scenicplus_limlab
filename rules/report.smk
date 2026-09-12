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
        # resources.seed is tracked because the eRegulon t-SNE is STOCHASTIC.
        # Without it, changing the seed would leave the figure untouched and
        # snakemake would report nothing to do -- a silent no-op of exactly the
        # kind the sibling repo just spent an increment removing.
        cfg=cfg_params("input.celltype_column", "input.reduction",
                       "visualization", "resources.seed"),
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


# =============================================================================
# R21: report.html -- the run, made readable.
#
# NOT one of the twenty steps. The bash driver has no equivalent, which is the
# point: a finished run used to leave ~20 loose PNGs and 4 TSVs with nothing
# saying which mattered or whether the run was healthy. `rule all` targets this
# file (SnakemakePlan.md Decision 4), so an ordinary run is not done until the
# run is readable -- non-optional by construction rather than by discipline.
#
# WHAT IT DEPENDS ON, AND WHAT IT DELIBERATELY DOES NOT.
#
# Declared inputs are R19's four tables and R20's figures, so the report cannot
# be built from a stale analysis stage and re-runs when either does.
#
# NOT declared: `logs/`, `logs/lsf/` and `QC/assembly.json`. Two different
# reasons, and neither is an oversight.
#
#   * The LOGS ARE WRITTEN BY THE JOBS THIS RULE WAITS FOR, and a rule's own log
#     does not exist while it runs. Declaring them would be a dependency on
#     files whose final content postdates the dependency -- the same shape as
#     the sibling repo's "an artifact written DURING a run cannot describe that
#     run completely". The report reads whatever is on disk when it runs and
#     says what is missing.
#   * `assembly.json` IS a real product of R07, but R07 is upstream of R19/R20
#     already, so depending on it adds an edge that changes nothing. It is read
#     opportunistically and reported as absent when it is not there -- which is
#     the honest outcome for a workspace whose genome pair was supplied without
#     the checks having run.
#
# THE REPORT IS THE DEFAULT TARGET, NOT THE ONLY ARTIFACT. Snakemake will not
# build it if an upstream rule failed, so a partial run leaves no report -- an
# accepted consequence recorded in the plan. The per-rule logs and (at I7) the
# provenance bundle are what a failed run leaves behind.
# =============================================================================
rule R21_report:
    """One self-contained HTML page: run, genome, tables, figures, compute."""
    input:
        tables=rules.R19_postprocess_tsv.output,
        figures=rules.R20_visualize.output,
        script=script_path("scenicplus_09_report.py"),
    output:
        html="report.html",
    log:
        log_path("R21_report")
    threads: 1
    resources:
        mem_mb=mem("R21_report"),
        runtime=rt("R21_report"),
    params:
        # Every value the page prints comes from files it reads at run time, so
        # there is no config key here whose change should rebuild it. The two
        # that shape the PAGE rather than its content are params on purpose:
        # editing either re-renders, which is cheap and correct.
        config_file=CONFIG_FILE,
        pipeline=workflow.basedir,
        head_rows=config.get("report", {}).get("head_rows", 15),
        max_embed_mb=config.get("report", {}).get("max_embed_mb", 4.0),
    shell:
        "python {input.script}"
        " --workspace ."
        " --config {params.config_file:q}"
        " --pipeline {params.pipeline:q}"
        " --out {output.html:q}"
        " --head-rows {params.head_rows}"
        " --max-embed-mb {params.max_embed_mb}"
        # {log} is this rule's own log. It exists from the moment the job
        # starts (tee creates it) and is still EMPTY when the script reads
        # logs/, because this script's output is printed afterwards.
        # Naming it lets the page annotate that row instead of raising a
        # false empty-log alarm on every single run.
        " --self-log {log}"
        " 2>&1 | tee {log}"
