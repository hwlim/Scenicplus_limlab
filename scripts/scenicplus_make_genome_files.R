#!/usr/bin/env Rscript
# Build genome_annotation.tsv and chromsizes.tsv for SCENIC+, on the assembly
# YOUR DATA IS ON.
#
# WHY THIS EXISTS. `scenicplus prepare_data download_genome_annotations` cannot
# produce chromsizes at all — it queries NCBI's retired `db=genome`, which
# answers HTTP 200 with <Count>0</Count> for every term (aertslab/scenicplus
# #640). For human that is merely an obstacle. FOR MOUSE IT HIDES A WORSE ONE:
# the download takes its assembly from whatever Ensembl currently serves, which
# for mouse is GRCm39. A run against mm10/GRCm38 data then gets an annotation on
# the WRONG ASSEMBLY, and that failure does not announce itself — chromosome
# names still convert, the search space still builds, and the peak-gene links
# are quietly wrong. Recorded in Development.md:35 from a 2026-05 kidney run:
#
#     Download gene annotation INFO   Using genome: GRCm39
#
# So for mouse, supplying these files is not a workaround for a broken
# endpoint. It is the only way to be sure which assembly you analysed.
#
#   scenicplus_make_genome_files.R --species mouse --out-dir .
#
# WHICH ENVIRONMENT. This needs EnsDb + BSgenome, which are in the
# scRNA_LimLab_Snake pipeline env, NOT in environment.cchmc.yml — the SCENIC+
# env deliberately carries only Seurat/Signac on the R side. Run it there, then
# point config.yaml at the two files it writes. Deriving the annotation from the
# same EnsDb the upstream pipeline used is the point: the peaks and the
# annotation then come from one assembly by construction, rather than by
# agreement between two services.
#
#   EnsDb.Mmusculus.v79  = GRCm38 = mm10
#   EnsDb.Hsapiens.v86   = GRCh38 = hg38

suppressPackageStartupMessages({
  library(optparse)
})

opt <- parse_args(OptionParser(
  option_list = list(
    make_option("--species", type = "character", default = NULL, metavar = "S",
                help = "REQUIRED. mouse | human."),
    make_option("--out-dir", type = "character", default = ".", metavar = "DIR",
                help = "Where to write the two files [default: %default]"),
    make_option("--prefix", type = "character", default = "", metavar = "P",
                help = paste("Optional filename prefix, e.g. `mm10.` -- useful",
                             "when several assemblies live side by side.")),
    make_option("--keep-scaffolds", action = "store_true", default = FALSE,
                help = paste("Keep non-standard contigs. Off by default: the",
                             "download path subsets to assembled molecules and",
                             "a search space over chrUn_* helps nobody."))),
  usage = "%prog --species mouse --out-dir .",
  description = paste("\nWrite genome_annotation.tsv and chromsizes.tsv for",
                      "SCENIC+ from a pinned EnsDb,\nso the assembly is the one",
                      "you chose rather than the one Ensembl serves today.\n")),
  convert_hyphens_to_underscores = TRUE)

# optparse omits NULL-default options, and `$` partial-matches on lists, so a
# name absent from the result can silently resolve to a longer one. Give every
# declared name an exact match to find. (Same fix as figr_grn.R.)
options(warnPartialMatchDollar = TRUE)

if (is.null(opt$species)) stop("--species is required (mouse | human)")
species <- tolower(opt$species)
if (!species %in% c("mouse", "human"))
  stop("--species must be mouse or human (got: ", opt$species, ")")

refs <- switch(species,
  mouse = list(ensdb = "EnsDb.Mmusculus.v79",  bsg = "BSgenome.Mmusculus.UCSC.mm10",
               build = "mm10", ensembl = "GRCm38"),
  human = list(ensdb = "EnsDb.Hsapiens.v86",   bsg = "BSgenome.Hsapiens.UCSC.hg38",
               build = "hg38", ensembl = "GRCh38"))

for (p in c(refs$ensdb, refs$bsg, "GenomicRanges", "ensembldb", "AnnotationFilter"))
  if (!requireNamespace(p, quietly = TRUE))
    stop(p, " is not installed. This script needs the EnsDb/BSgenome stack, ",
         "which lives in the scRNA_LimLab_Snake pipeline env rather than the ",
         "SCENIC+ env -- run it there.")

suppressPackageStartupMessages({
  library(refs$ensdb, character.only = TRUE)
  library(refs$bsg,   character.only = TRUE)
  library(GenomicRanges)
})
ensdb <- get(refs$ensdb)
bsg   <- get(refs$bsg)

message(sprintf("%s: %s (%s = %s)", species, refs$ensdb, refs$ensembl, refs$build))

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(f) file.path(opt$out_dir, paste0(opt$prefix, f))

# ---------------------------------------------------------------------------
# chromsizes.tsv -- Chromosome / Start / End, header row, tab separated.
#
# search_space reads it with pandas.read_table, so the header and the Start
# column are both load-bearing: a raw UCSC .chrom.sizes has neither and pandas
# silently takes the first chromosome as the header.
# ---------------------------------------------------------------------------
sl <- GenomeInfoDb::seqlengths(bsg)
if (!opt$keep_scaffolds)
  sl <- sl[names(sl) %in% GenomeInfoDb::standardChromosomes(bsg)]
chromsizes <- data.frame(Chromosome = names(sl), Start = 0L,
                         End = as.integer(sl), row.names = NULL)
chromsizes <- chromsizes[order(chromsizes$Chromosome), ]
write.table(chromsizes, out("chromsizes.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
message(sprintf("  chromsizes.tsv: %d chromosome(s)", nrow(chromsizes)))

# ---------------------------------------------------------------------------
# genome_annotation.tsv -- the columns get_search_space requires:
#   Chromosome, Start, End, Strand, Gene, Transcription_Start_Site,
#   Transcript_type
# One row per TRANSCRIPT (not per gene): the TSS is what the search space is
# built around, and a gene with several TSSs contributes each of them, which is
# what the download path produces too.
# ---------------------------------------------------------------------------
# TxBiotypeFilter lives in AnnotationFilter, not ensembldb -- ensembldb 2.34
# does not export it, and there is no formula interface on transcripts() here
# either. Checked against the installed version rather than assumed.
tx <- ensembldb::transcripts(
  ensdb, columns = c("gene_name", "tx_biotype", "seq_name"),
  filter = AnnotationFilter::TxBiotypeFilter("protein_coding"))

ann <- data.frame(
  Chromosome = paste0("chr", as.character(GenomicRanges::seqnames(tx))),
  Start      = GenomicRanges::start(tx),
  End        = GenomicRanges::end(tx),
  Strand     = ifelse(as.character(GenomicRanges::strand(tx)) == "+", "+", "-"),
  Gene       = as.character(GenomicRanges::mcols(tx)$gene_name),
  Transcription_Start_Site = ifelse(
    as.character(GenomicRanges::strand(tx)) == "+",
    GenomicRanges::start(tx), GenomicRanges::end(tx)),
  Transcript_type = as.character(GenomicRanges::mcols(tx)$tx_biotype),
  stringsAsFactors = FALSE, row.names = NULL)

# Ensembl spells the mitochondrion MT; UCSC spells it chrM. Everything else is
# just the chr prefix added above.
ann$Chromosome[ann$Chromosome == "chrMT"] <- "chrM"
if (!opt$keep_scaffolds)
  ann <- ann[ann$Chromosome %in% GenomeInfoDb::standardChromosomes(bsg), ]
ann <- ann[nzchar(ann$Gene) & !is.na(ann$Gene), ]
ann <- ann[order(ann$Chromosome, ann$Start), ]

write.table(ann, out("genome_annotation.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
message(sprintf("  genome_annotation.tsv: %d transcript(s), %d gene(s), %d chromosome(s)",
                nrow(ann), length(unique(ann$Gene)), length(unique(ann$Chromosome))))

# A sanity line worth printing every time: chromosome 1's length IS the
# assembly. mm10 chr1 is 195,471,971 and GRCm39 chr1 is 195,154,279, so this
# one number distinguishes the assembly this whole script exists to pin down.
c1 <- chromsizes$End[chromsizes$Chromosome == "chr1"]
if (length(c1) == 1)
  message(sprintf("  chr1 = %s bp  <- this is what identifies the assembly as %s",
                  format(c1, big.mark = ","), refs$build))

message("\nPoint config.yaml at these:")
message("  input:")
message("    genome_annotation: ", normalizePath(out("genome_annotation.tsv")))
message("    chromsizes: ",        normalizePath(out("chromsizes.tsv")))
message("    species: ", if (species == "mouse") "\"mmusculus\"" else "\"hsapiens\"")
