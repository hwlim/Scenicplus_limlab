# =============================================================================
# R01-R05: from a Seurat object to the region sets the GRN inference consumes.
#
# The scripts are the bash driver's, unchanged, called with the same flags. What
# changes is the scheduling: sentinels and `.cfgsha` sidecars become Snakemake's
# input/params/code triggers.
#
# Each rule's `params` carries the config keys ITS script reads, taken from the
# driver's own slice lists. That is what makes a config edit re-run the right
# step and only that step -- see cfg_params() in common.smk for why declaring
# them is not optional.
#
# Resources are deliberately uniform. The threads are real, because which steps
# take --n_cpu is known; the memory is the single tier from common.smk, which is
# the status quo rather than a measurement. I5 replaces it with numbers from a
# run that records peak RSS per step. Guessing them here would look like
# evidence.
# =============================================================================

EXPORT_DIR = stage_path("export", "seurat_export")

# The fixed part of step 01's output. The embeddings are variable in number, one
# per reduction in the object, so they cannot all be named -- but the ONE that
# matters can be, when the config names it. Declaring it turns "the layout you
# asked for was never written" into a MissingOutputException at step 1 rather
# than a puzzle at step 20.
EXPORT_FILES = [
    os.path.join(EXPORT_DIR, f)
    for f in (
        "rna_raw_counts.mtx", "rna_norm.mtx", "rna_features.tsv",
        "rna_barcodes.tsv", "atac_counts.mtx", "atac_barcodes.tsv",
        "atac_regions.tsv", "cell_metadata.tsv", "reductions.txt",
        "summary.txt",
    )
]

_REDUCTION = config["input"].get("reduction", "").strip()
if _REDUCTION:
    EXPORT_FILES.append(os.path.join(EXPORT_DIR, f"embedding_{_REDUCTION}.tsv"))


rule R01_seurat_export:
    """Seurat .rds -> matrices, metadata and embeddings on disk."""
    input:
        rds=config["input"]["seurat_rds"],
        script=script_path("scenicplus_01_seurat_to_anndata.R"),
    output:
        EXPORT_FILES,
    log:
        log_path("R01_seurat_export"),
    threads: 1
    resources:
        mem_mb=mem("R01_seurat_export"),
        runtime=rt("R01_seurat_export"),
    params:
        cfg=cfg_params("input.seurat_rds", "input.celltype_column",
                       "input.celltype_scope", "input.reduction"),
        out_dir=EXPORT_DIR,
        celltype_col=config["input"]["celltype_column"],
        # Omitted rather than passed empty: see opt_arg() in common.smk for the
        # measurement behind that.
        scope=opt_arg("--celltype_scope",
                      ",".join(str(x) for x in
                               config["input"].get("celltype_scope", []) or [])),
        reduction=opt_arg("--reduction", _REDUCTION),
    shell:
        "Rscript {input.script}"
        " --rds {input.rds:q}"
        " --celltype_col {params.celltype_col:q}"
        " {params.scope}"
        " {params.reduction}"
        " --out_dir {params.out_dir:q}"
        " 2>&1 | tee {log}"


rule R02_build_anndata:
    """The matrices, packed into the AnnData SCENIC+ expects."""
    input:
        export=EXPORT_FILES,
        script=script_path("scenicplus_02_build_anndata.py"),
    output:
        h5ad=stage_path("anndata", "rna.h5ad"),
    log:
        log_path("R02_build_anndata"),
    threads: 1
    resources:
        mem_mb=mem("R02_build_anndata"),
        runtime=rt("R02_build_anndata"),
    params:
        cfg=cfg_params("input.celltype_column", "input.reduction"),
        in_dir=EXPORT_DIR,
        celltype_col=config["input"]["celltype_column"],
        reduction=opt_arg("--reduction", _REDUCTION),
    shell:
        "python {input.script}"
        " --in_dir {params.in_dir:q}"
        " --out_h5ad {output.h5ad:q}"
        " --celltype_col {params.celltype_col:q}"
        " {params.reduction}"
        " 2>&1 | tee {log}"


rule R03_create_cistopic:
    """The ATAC side, as a cisTopic object.

    No config params: this script reads no config at all, which is why the
    driver's slice for step 3 is empty. That emptiness is also the driver's
    blind spot -- it hashes CONFIG, not CODE, so when this script changed to
    stop tagging cell names, no workspace re-ran it and the old tagged names
    survived.

    Snakemake does NOT fix that for free, which is worth stating because the
    plan assumed it did. Its `code` trigger covers a rule's own text, not a
    script the rule shells out to. Measured: edit an external script and
    snakemake reports "Nothing to be done". Declaring the script as an INPUT is
    what closes it, which is why every rule here does.
    """
    input:
        export=EXPORT_FILES,
        script=script_path("scenicplus_03_create_cistopic.py"),
    output:
        pkl=stage_path("cistopic", "cistopic_obj.pkl"),
    log:
        log_path("R03_create_cistopic"),
    threads: 1
    resources:
        mem_mb=mem("R03_create_cistopic"),
        runtime=rt("R03_create_cistopic"),
    params:
        in_dir=EXPORT_DIR,
    shell:
        "python {input.script}"
        " --in_dir {params.in_dir:q}"
        " --out_pkl {output.pkl:q}"
        " 2>&1 | tee {log}"


rule R04_topic_modeling:
    """LDA over the cell-by-region matrix: the heaviest early step.

    Memory is the SUM over the models fitted at once, not the largest of them,
    because `cistopic.n_topics` fits one model per value concurrently. That is
    what makes this the step whose reservation matters most, and why it is the
    first place I5 should point a measurement.
    """
    input:
        pkl=stage_path("cistopic", "cistopic_obj.pkl"),
        script=script_path("scenicplus_04_topic_modeling.py"),
    output:
        pkl=stage_path("cistopic", "cistopic_obj_with_topics.pkl"),
    log:
        log_path("R04_topic_modeling"),
    threads: n_cpu()
    resources:
        mem_mb=mem("R04_topic_modeling"),
        runtime=rt("R04_topic_modeling"),
    params:
        cfg=cfg_params("cistopic.n_topics", "cistopic.n_iter", "cistopic.alpha",
                       "cistopic.alpha_by_topic", "cistopic.eta",
                       "cistopic.eta_by_topic", "cistopic.random_state",
                       "resources.seed"),
        config_file=CONFIG_FILE,
        tmp_dir=os.path.join(config.get("output", {}).get("tmp", "tmp"), "lda"),
    shell:
        "python {input.script}"
        " --in_pkl {input.pkl:q}"
        " --out_pkl {output.pkl:q}"
        " --config {params.config_file:q}"
        " --tmp_dir {params.tmp_dir:q}"
        " --n_cpu {threads}"
        " 2>&1 | tee {log}"


rule R05_region_sets:
    """Binarized topics and differentially accessible regions, as BED sets.

    A directory output, because the number of sets follows the topic count and
    the cell types present. The driver used a `.done` sentinel for the same
    reason; Snakemake can express the directory itself, so the sentinel goes.
    """
    input:
        pkl=stage_path("cistopic", "cistopic_obj_with_topics.pkl"),
        script=script_path("scenicplus_05_region_sets.py"),
    output:
        directory(stage_path("cistopic", "region_sets")),
    log:
        log_path("R05_region_sets"),
    threads: n_cpu()
    resources:
        mem_mb=mem("R05_region_sets"),
        runtime=rt("R05_region_sets"),
    params:
        cfg=cfg_params("input.celltype_column", "cistopic.dar_adjpval_thr",
                       "cistopic.dar_log2fc_thr"),
        config_file=CONFIG_FILE,
    shell:
        "python {input.script}"
        " --in_pkl {input.pkl:q}"
        " --out_dir {output:q}"
        " --config {params.config_file:q}"
        " --n_cpu {threads}"
        " 2>&1 | tee {log}"
