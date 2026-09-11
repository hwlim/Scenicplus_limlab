# =============================================================================
# R06, R08-R18: the GRN inference DAG.
#
# Every rule is one `scenicplus_06_grn_stage.py --stage X` call, exactly as the
# bash driver makes it. That script stays the single place that knows each
# stage's flags; it took a run through all twenty steps to get right and nothing
# here second-guesses it.
#
# WHAT CHANGES IS THE SHAPE. The driver walks these thirteen stages in a line,
# because flattening the inner snakemake traded away its intra-DAG parallelism
# (`CLAUDE.md` says so). Declaring each stage's real inputs gets it back:
#
#     R09 cistarget   || R10 dem            both read only the region sets
#     R12 tf_to_gene  || R13 region_to_gene
#     R14 egrn_direct || R15 egrn_extended
#     R16 aucell_direct || R17 aucell_extended
#
# The dependencies are taken from each stage's ARGUMENTS, not from the driver's
# ordering, and one of them is not obvious: R12 needs `tf_names.txt`, which
# R11 produces. The driver satisfied that by accident of being sequential.
#
# `--scplus_out` is a DIRECTORY and every stage derives its own paths inside it
# (out_paths() in the stage script). So each rule names the file(s) that stage
# writes, and `4.grn/` accumulates.
# =============================================================================

GRN = stage_path("grn")          # what --scplus_out points at
GRN_STAGE = script_path("scenicplus_06_grn_stage.py")


def grn_out(name):
    return os.path.join(GRN, name)


# The genome pair stops being optional here: R08 cannot build a search space
# without it, and the download that would otherwise fill the gap cannot produce
# chromsizes for anyone. A parse-time refusal beats a MissingInputException
# naming an empty string.
if not GENOME_SUPPLIED:
    raise ValueError(
        "input.genome_annotation and input.chromsizes must be set: the GRN "
        "stages need them.\n"
        "  Step 7's download CANNOT produce chromsizes -- NCBI's db=genome is "
        "retired, it exits 0\n"
        "  having written only the annotation, and that annotation is "
        "Ensembl-style besides.\n"
        "  Build both once per assembly with "
        "scripts/scenicplus_make_genome_files.R, in an environment\n"
        "  that has EnsDb and BSgenome (the scRNA_LimLab_Snake one). RUNBOOK "
        "section 2b has the detail.")


rule R06_prepare_gex_acc:
    """Pair the expression and accessibility sides into one MuData.

    With `scenicplus.is_multiome: true` the pairing is by BARCODE, which is why
    step 3's cell naming has to match the AnnData's -- the failure that read as
    "no cells in both assays".
    """
    input:
        cistopic=stage_path("cistopic", "cistopic_obj_with_topics.pkl"),
        adata=stage_path("anndata", "rna.h5ad"),
        script=GRN_STAGE,
    output:
        mudata=grn_out("ACC_GEX.h5mu"),
    log:
        log_path("R06_prepare_gex_acc")
    threads: 1
    resources:
        mem_mb=mem("R06_prepare_gex_acc"),
        runtime=rt("R06_prepare_gex_acc"),
    params:
        cfg=cfg_params("scenicplus.is_multiome", "scenicplus.bc_transform_func"),
        config_file=CONFIG_FILE,
        scplus_out=GRN,
    shell:
        "python {input.script} --stage prepare_gex_acc"
        " --config {params.config_file:q} --scplus_out {params.scplus_out:q}"
        " --cistopic_obj {input.cistopic:q} --adata {input.adata:q}"
        " 2>&1 | tee {log}"


rule R08_search_space:
    """Which regions may be linked to which genes.

    Joins three files by chromosome -- the MuData's regions, the annotation and
    the chromsizes -- which is why R07 checks all three agree before this runs
    rather than leaving a pandas KeyError to explain it.
    """
    input:
        mudata=grn_out("ACC_GEX.h5mu"),
        annotation=grn_out("genome_annotation.tsv"),
        chromsizes=grn_out("chromsizes.tsv"),
        script=GRN_STAGE,
    output:
        search_space=grn_out("search_space.tsv"),
    log:
        log_path("R08_search_space")
    threads: 1
    resources:
        mem_mb=mem("R08_search_space"),
        runtime=rt("R08_search_space"),
    params:
        cfg=cfg_params("scenicplus.search_space_upstream",
                       "scenicplus.search_space_downstream",
                       "scenicplus.search_space_extend_tss"),
        config_file=CONFIG_FILE,
        scplus_out=GRN,
    shell:
        "python {input.script} --stage search_space"
        " --config {params.config_file:q} --scplus_out {params.scplus_out:q}"
        " 2>&1 | tee {log}"


rule R09_cistarget:
    """Motif enrichment by ranking. Reads the 32.8 GB ranking database, which
    is what dominates it -- tuning the parameters will not help."""
    input:
        region_sets=stage_path("cistopic", "region_sets"),
        db=config["input"]["ctx_db"],
        motifs=config["input"]["motif_annotations"],
        script=GRN_STAGE,
    output:
        result=grn_out("ctx_results.hdf5"),
        # Declared because the stage writes it, per the plan's rule about
        # outputs written as side effects. If a configuration turns out not to
        # produce it, the MissingOutputException says so and it moves out.
        html=grn_out("ctx_results.html"),
    log:
        log_path("R09_cistarget")
    threads: n_cpu()
    resources:
        mem_mb=mem("R09_cistarget"),
        runtime=rt("R09_cistarget"),
    params:
        cfg=cfg_params("input.ctx_db", "input.motif_annotations", "input.species",
                       "scenicplus.fraction_overlap_w_ctx_database",
                       "scenicplus.ctx_auc_threshold", "scenicplus.ctx_nes_threshold",
                       "scenicplus.ctx_rank_threshold", "scenicplus.annotation_version",
                       "scenicplus.motif_similarity_fdr",
                       "scenicplus.orthologous_identity_threshold",
                       "scenicplus.annotations_to_use", "resources.n_cpu",
                       "output.tmp"),
        config_file=CONFIG_FILE,
        scplus_out=GRN,
    shell:
        "python {input.script} --stage cistarget"
        " --config {params.config_file:q} --scplus_out {params.scplus_out:q}"
        " --region_sets {input.region_sets:q}"
        " 2>&1 | tee {log}"


rule R10_dem:
    """Motif enrichment by differential scoring. Reads the 12.9 GB score
    database. Independent of R09, so the two run together."""
    input:
        region_sets=stage_path("cistopic", "region_sets"),
        db=config["input"]["dem_db"],
        motifs=config["input"]["motif_annotations"],
        script=GRN_STAGE,
    output:
        result=grn_out("dem_results.hdf5"),
        html=grn_out("dem_results.html"),
    log:
        log_path("R10_dem")
    threads: n_cpu()
    resources:
        mem_mb=mem("R10_dem"),
        runtime=rt("R10_dem"),
    params:
        cfg=cfg_params("input.dem_db", "input.motif_annotations", "input.species",
                       "scenicplus.fraction_overlap_w_dem_database",
                       "scenicplus.dem_max_bg_regions",
                       "scenicplus.dem_balance_number_of_promoters",
                       "scenicplus.dem_promoter_space", "scenicplus.dem_adj_pval_thr",
                       "scenicplus.dem_log2fc_thr", "scenicplus.dem_mean_fg_thr",
                       "scenicplus.dem_motif_hit_thr", "scenicplus.annotation_version",
                       "scenicplus.motif_similarity_fdr",
                       "scenicplus.orthologous_identity_threshold",
                       "scenicplus.annotations_to_use", "resources.n_cpu",
                       "resources.seed", "output.tmp"),
        config_file=CONFIG_FILE,
        scplus_out=GRN,
    shell:
        "python {input.script} --stage dem"
        " --config {params.config_file:q} --scplus_out {params.scplus_out:q}"
        " --region_sets {input.region_sets:q}"
        " 2>&1 | tee {log}"


rule R11_prepare_menr:
    """Both enrichment results into cistromes, and the TF list everything
    downstream is keyed on."""
    input:
        ctx=grn_out("ctx_results.hdf5"),
        dem=grn_out("dem_results.hdf5"),
        mudata=grn_out("ACC_GEX.h5mu"),
        script=GRN_STAGE,
    output:
        direct=grn_out("cistromes_direct.h5ad"),
        extended=grn_out("cistromes_extended.h5ad"),
        tf_names=grn_out("tf_names.txt"),
    log:
        log_path("R11_prepare_menr")
    threads: 1
    resources:
        mem_mb=mem("R11_prepare_menr"),
        runtime=rt("R11_prepare_menr"),
    params:
        cfg=cfg_params("scenicplus.direct_annotation",
                       "scenicplus.extended_annotation"),
        config_file=CONFIG_FILE,
        scplus_out=GRN,
    shell:
        "python {input.script} --stage prepare_menr"
        " --config {params.config_file:q} --scplus_out {params.scplus_out:q}"
        " 2>&1 | tee {log}"


rule R12_tf_to_gene:
    """TF-to-gene importances by gradient boosting. Genuinely CPU-parallel.

    `tf_names.txt` comes from R11, which the driver's linear order satisfied by
    accident rather than by declaration.
    """
    input:
        mudata=grn_out("ACC_GEX.h5mu"),
        tf_names=grn_out("tf_names.txt"),
        script=GRN_STAGE,
    output:
        adj=grn_out("tf_to_gene_adj.tsv"),
    log:
        log_path("R12_tf_to_gene")
    threads: n_cpu()
    resources:
        mem_mb=mem("R12_tf_to_gene"),
        runtime=rt("R12_tf_to_gene"),
    params:
        cfg=cfg_params("grn.tf_to_gene_importance_method", "resources.n_cpu",
                       "resources.seed", "output.tmp"),
        config_file=CONFIG_FILE,
        scplus_out=GRN,
    shell:
        "python {input.script} --stage tf_to_gene"
        " --config {params.config_file:q} --scplus_out {params.scplus_out:q}"
        " 2>&1 | tee {log}"


rule R13_region_to_gene:
    """Region-to-gene importances and correlations, within the search space."""
    input:
        mudata=grn_out("ACC_GEX.h5mu"),
        search_space=grn_out("search_space.tsv"),
        script=GRN_STAGE,
    output:
        adj=grn_out("region_to_gene_adj.tsv"),
    log:
        log_path("R13_region_to_gene")
    threads: n_cpu()
    resources:
        mem_mb=mem("R13_region_to_gene"),
        runtime=rt("R13_region_to_gene"),
    params:
        cfg=cfg_params("grn.region_to_gene_importance_method",
                       "grn.region_to_gene_correlation_method",
                       "resources.n_cpu", "output.tmp"),
        config_file=CONFIG_FILE,
        scplus_out=GRN,
    shell:
        "python {input.script} --stage region_to_gene"
        " --config {params.config_file:q} --scplus_out {params.scplus_out:q}"
        " 2>&1 | tee {log}"


# R14/R15 and R16/R17 differ only in direct vs extended annotations, but they
# stay separate rules rather than becoming one wildcard rule: each takes a
# different cistrome file and writes a different output, the stage script
# already dispatches on the name, and `--forcerun R15_egrn_extended` naming one
# of them is worth more than the lines saved.
_EGRN_PARAMS = ("grn.order_regions_to_genes_by", "grn.order_TFs_to_genes_by",
                "grn.gsea_n_perm", "grn.quantile_thresholds_region_to_gene",
                "grn.top_n_regionTogenes_per_gene",
                "grn.top_n_regionTogenes_per_region", "grn.min_regions_per_gene",
                "grn.rho_threshold", "grn.min_target_genes", "input.ctx_db",
                "resources.n_cpu", "output.tmp")


rule R14_egrn_direct:
    """eRegulons from directly annotated motifs."""
    input:
        tf2g=grn_out("tf_to_gene_adj.tsv"),
        r2g=grn_out("region_to_gene_adj.tsv"),
        cistromes=grn_out("cistromes_direct.h5ad"),
        db=config["input"]["ctx_db"],
        script=GRN_STAGE,
    output:
        eregulons=grn_out("eRegulons_direct.tsv"),
    log:
        log_path("R14_egrn_direct")
    threads: n_cpu()
    resources:
        mem_mb=mem("R14_egrn_direct"),
        runtime=rt("R14_egrn_direct"),
    params:
        cfg=cfg_params(*_EGRN_PARAMS),
        config_file=CONFIG_FILE,
        scplus_out=GRN,
    shell:
        "python {input.script} --stage egrn_direct"
        " --config {params.config_file:q} --scplus_out {params.scplus_out:q}"
        " 2>&1 | tee {log}"


rule R15_egrn_extended:
    """eRegulons from orthology- and similarity-extended motif annotations."""
    input:
        tf2g=grn_out("tf_to_gene_adj.tsv"),
        r2g=grn_out("region_to_gene_adj.tsv"),
        cistromes=grn_out("cistromes_extended.h5ad"),
        db=config["input"]["ctx_db"],
        script=GRN_STAGE,
    output:
        eregulons=grn_out("eRegulons_extended.tsv"),
    log:
        log_path("R15_egrn_extended")
    threads: n_cpu()
    resources:
        mem_mb=mem("R15_egrn_extended"),
        runtime=rt("R15_egrn_extended"),
    params:
        cfg=cfg_params(*_EGRN_PARAMS),
        config_file=CONFIG_FILE,
        scplus_out=GRN,
    shell:
        "python {input.script} --stage egrn_extended"
        " --config {params.config_file:q} --scplus_out {params.scplus_out:q}"
        " 2>&1 | tee {log}"


rule R16_aucell_direct:
    """Per-cell activity of each direct eRegulon."""
    input:
        eregulons=grn_out("eRegulons_direct.tsv"),
        mudata=grn_out("ACC_GEX.h5mu"),
        script=GRN_STAGE,
    output:
        auc=grn_out("AUCell_direct.h5mu"),
    log:
        log_path("R16_aucell_direct")
    threads: n_cpu()
    resources:
        mem_mb=mem("R16_aucell_direct"),
        runtime=rt("R16_aucell_direct"),
    params:
        cfg=cfg_params("resources.n_cpu"),
        config_file=CONFIG_FILE,
        scplus_out=GRN,
    shell:
        "python {input.script} --stage aucell_direct"
        " --config {params.config_file:q} --scplus_out {params.scplus_out:q}"
        " 2>&1 | tee {log}"


rule R17_aucell_extended:
    """Per-cell activity of each extended eRegulon."""
    input:
        eregulons=grn_out("eRegulons_extended.tsv"),
        mudata=grn_out("ACC_GEX.h5mu"),
        script=GRN_STAGE,
    output:
        auc=grn_out("AUCell_extended.h5mu"),
    log:
        log_path("R17_aucell_extended")
    threads: n_cpu()
    resources:
        mem_mb=mem("R17_aucell_extended"),
        runtime=rt("R17_aucell_extended"),
    params:
        cfg=cfg_params("resources.n_cpu"),
        config_file=CONFIG_FILE,
        scplus_out=GRN,
    shell:
        "python {input.script} --stage aucell_extended"
        " --config {params.config_file:q} --scplus_out {params.scplus_out:q}"
        " 2>&1 | tee {log}"


rule R18_scplus_mudata:
    """Everything into one MuData: the GRN result the analysis stage reads."""
    input:
        mudata=grn_out("ACC_GEX.h5mu"),
        auc_direct=grn_out("AUCell_direct.h5mu"),
        auc_extended=grn_out("AUCell_extended.h5mu"),
        ereg_direct=grn_out("eRegulons_direct.tsv"),
        ereg_extended=grn_out("eRegulons_extended.tsv"),
        script=GRN_STAGE,
    output:
        mdata=grn_out("scplusmdata.h5mu"),
    log:
        log_path("R18_scplus_mudata")
    threads: 1
    resources:
        mem_mb=mem("R18_scplus_mudata"),
        runtime=rt("R18_scplus_mudata"),
    params:
        config_file=CONFIG_FILE,
        scplus_out=GRN,
    shell:
        "python {input.script} --stage scplus_mudata"
        " --config {params.config_file:q} --scplus_out {params.scplus_out:q}"
        " 2>&1 | tee {log}"
