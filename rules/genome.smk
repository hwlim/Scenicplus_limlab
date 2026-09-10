# =============================================================================
# R07: the genome annotation and chromosome sizes.
#
# THIS RULE DOES NOT GENERATE THEM, and that is a decision rather than an
# omission. Two reasons, the first decisive:
#
#   1. It cannot. Building them needs EnsDb + BSgenome, and the SCENIC+
#      environment has NONE of that stack -- measured: EnsDb.Hsapiens.v86,
#      EnsDb.Mmusculus.v79, both BSgenome packages, ensembldb and
#      AnnotationFilter are all absent. They live in the scRNA_LimLab_Snake
#      environment, which is where scenicplus_make_genome_files.R says to run.
#   2. They are reference data, not run outputs. They depend only on species and
#      assembly, so they belong beside the 45.7 GB of cisTarget databases and are
#      built once per assembly, by hand, like those. SnakemakePlan.md's
#      Decision 1 asks for a shared reference directory; pointing
#      input.genome_annotation / input.chromsizes into one satisfies it without
#      a new config key.
#
# What the rule does instead is the part that was never done anywhere: CHECK the
# pair, against each other AND against the peaks this run actually has, before
# the expensive stages consume them. See scenicplus_genome_prepare.py for the
# four checks and why each one exists.
#
# The consequence worth knowing: this depends on R01, because the peak
# chromosomes come from the export. The driver's step 7 had no such dependency
# and no such check.
# =============================================================================

_ANNOTATION = config["input"].get("genome_annotation", "").strip()
_CHROMSIZES = config["input"].get("chromsizes", "").strip()
GENOME_SUPPLIED = bool(_ANNOTATION and _CHROMSIZES)

if bool(_ANNOTATION) != bool(_CHROMSIZES):
    raise ValueError(
        "set input.genome_annotation and input.chromsizes BOTH or NEITHER.\n"
        f"  genome_annotation: {_ANNOTATION or '(empty)'}\n"
        f"  chromsizes:        {_CHROMSIZES or '(empty)'}\n"
        "  They come out of one call, so supplying one still leaves the other to\n"
        "  the download -- and the download cannot produce chromsizes for anyone\n"
        "  (NCBI's db=genome is retired; aertslab/scenicplus#640).")

GENOME_FILES = [
    stage_path("grn", "genome_annotation.tsv"),
    stage_path("grn", "chromsizes.tsv"),
    stage_path("qc", "assembly.json"),
]

if not GENOME_SUPPLIED:
    # Not a parse error, because nothing consumes these until I3. It becomes one
    # then: the search space cannot be built without them, and the download that
    # would otherwise fill the gap exits 0 having written half of what it should.
    print("[scenicplus] input.genome_annotation / input.chromsizes are unset, so "
          "R07 is not available.\n"
          "[scenicplus]   Build them once per assembly with "
          "scripts/scenicplus_make_genome_files.R, in an\n"
          "[scenicplus]   environment that has EnsDb and BSgenome, then point the "
          "config at them.\n"
          "[scenicplus]   RUNBOOK section 2b has the detail, including why this is "
          "the normal path.")


if GENOME_SUPPLIED:

    rule R07_genome_annot:
        """Check a supplied annotation + chromsizes pair, then place it.

        Chromosome 1's length is checked against the assembly, which is the one
        number that distinguishes mm10 from GRCm39 -- the failure that reached
        a real kidney run and would not have announced itself.
        """
        input:
            annotation=_ANNOTATION,
            chromsizes=_CHROMSIZES,
            # The cross-check needs the peaks, so this waits for the export.
            regions=os.path.join(EXPORT_DIR, "atac_regions.tsv"),
            script=script_path("scenicplus_genome_prepare.py"),
        output:
            annotation=stage_path("grn", "genome_annotation.tsv"),
            chromsizes=stage_path("grn", "chromsizes.tsv"),
            record=stage_path("qc", "assembly.json"),
        log:
            log_path("R07_genome_annot"),
        threads: 1
        resources:
            mem_mb=mem("default", config),
            runtime=rt("default", config),
        params:
            cfg=cfg_params("input.genome_annotation", "input.chromsizes",
                           "input.species", "input.assembly"),
            species=config["input"]["species"],
            assembly=opt_arg("--assembly", config["input"].get("assembly", "")),
        shell:
            "python {input.script}"
            " --annotation {input.annotation:q}"
            " --chromsizes {input.chromsizes:q}"
            " --atac_regions {input.regions:q}"
            " --out_annotation {output.annotation:q}"
            " --out_chromsizes {output.chromsizes:q}"
            " --record {output.record:q}"
            " --species {params.species:q}"
            " {params.assembly}"
            " 2>&1 | tee {log}"
