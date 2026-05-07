#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# Split a Seurat .rds object into the inputs SCENIC+ needs:
#   - rna.h5ad           : AnnData with normalized log-counts and raw counts in .raw
#   - atac_counts.mtx    : sparse peak x cell ATAC count matrix
#   - atac_barcodes.tsv  : cell barcodes (matching RNA, in same order)
#   - atac_regions.tsv   : peak coordinates as "chr:start-end"
#   - cell_metadata.tsv  : cell metadata (cell-type column required)
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Matrix)
})

option_list <- list(
  make_option(c("--rds"),            type = "character"),
  make_option(c("--celltype_col"),   type = "character"),
  make_option(c("--celltype_scope"), type = "character", default = "",
              help = "Comma-separated cell-type values to keep. Empty = use all cells."),
  make_option(c("--out_dir"),        type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(file.exists(opt$rds))
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

message("[seurat_to_anndata] Reading: ", opt$rds)
obj <- readRDS(opt$rds)

if (!"RNA" %in% Assays(obj))    stop("Seurat object missing 'RNA' assay")
if (!"peaks" %in% Assays(obj))  stop("Seurat object missing 'peaks' assay")
if (!opt$celltype_col %in% colnames(obj@meta.data)) {
  stop("celltype_column '", opt$celltype_col, "' not found in metadata")
}

# --- Optional cell-type subsetting (drives sensitivity for focused analyses)
scope <- trimws(strsplit(opt$celltype_scope, ",", fixed = TRUE)[[1]])
scope <- scope[nchar(scope) > 0]
if (length(scope) > 0) {
  available <- unique(as.character(obj@meta.data[[opt$celltype_col]]))
  missing <- setdiff(scope, available)
  if (length(missing) > 0) {
    stop("celltype_scope contains values not present in '", opt$celltype_col,
         "': ", paste(missing, collapse = ", "))
  }
  keep <- which(as.character(obj@meta.data[[opt$celltype_col]]) %in% scope)
  if (length(keep) == 0) stop("celltype_scope filter removed every cell.")
  message("[seurat_to_anndata] celltype_scope: keeping ", length(keep),
          "/", ncol(obj), " cells across ", length(scope), " cell types")
  obj <- obj[, keep]
}

# --- RNA: write h5ad-friendly artifacts. We avoid the SeuratDisk dependency by
#     writing a 10X-style mtx that the python step ingests via scanpy. The raw
#     (un-normalized) counts go to rna_raw_counts.mtx; the normalized log-counts
#     go to rna_norm.mtx so both can be packed into AnnData (.X and .raw).
DefaultAssay(obj) <- "RNA"
rna_raw  <- GetAssayData(obj, assay = "RNA", slot = "counts")
rna_norm <- GetAssayData(obj, assay = "RNA", slot = "data")
if (length(rna_norm@x) == 0) {
  stop("RNA 'data' slot is empty — normalize the RNA assay before running.")
}

writeMM(rna_raw,  file.path(opt$out_dir, "rna_raw_counts.mtx"))
writeMM(rna_norm, file.path(opt$out_dir, "rna_norm.mtx"))
write.table(rownames(rna_raw),
            file.path(opt$out_dir, "rna_features.tsv"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(colnames(rna_raw),
            file.path(opt$out_dir, "rna_barcodes.tsv"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# --- ATAC counts (peaks x cells, raw)
atac_raw <- GetAssayData(obj, assay = "peaks", slot = "counts")
writeMM(atac_raw, file.path(opt$out_dir, "atac_counts.mtx"))
write.table(colnames(atac_raw),
            file.path(opt$out_dir, "atac_barcodes.tsv"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Normalize peak names to "chr:start-end" used by pycisTopic/SCENIC+.
peak_names <- rownames(atac_raw)
peak_names <- gsub("[-_]", ":", peak_names, fixed = FALSE)
peak_names <- sub("^([^:]+):([^:]+):([^:]+)$", "\\1:\\2-\\3", peak_names)
write.table(peak_names,
            file.path(opt$out_dir, "atac_regions.tsv"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# --- Cell metadata
meta <- obj@meta.data
meta$barcode <- rownames(meta)
write.table(meta,
            file.path(opt$out_dir, "cell_metadata.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- Pre-computed dimensional reductions (UMAP / PCA) — keep them so we can
#     reuse the curated layout instead of recomputing one downstream.
reductions <- Reductions(obj)
write_embedding <- function(name) {
  emb <- Embeddings(obj, reduction = name)
  emb_df <- as.data.frame(emb)
  emb_df$barcode <- rownames(emb_df)
  emb_df <- emb_df[, c("barcode", setdiff(colnames(emb_df), "barcode"))]
  write.table(emb_df,
              file.path(opt$out_dir, paste0("embedding_", name, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
for (rname in reductions) {
  tryCatch(write_embedding(rname),
           error = function(e) message("Skipping reduction '", rname, "': ", e$message))
}
writeLines(reductions, file.path(opt$out_dir, "reductions.txt"))

# --- Provenance
writeLines(
  c(
    paste0("celltype_column: ", opt$celltype_col),
    paste0("celltype_scope: ",
           if (length(scope) > 0) paste(scope, collapse = ",") else "(all)"),
    paste0("n_cells: ", ncol(rna_raw)),
    paste0("n_genes: ", nrow(rna_raw)),
    paste0("n_peaks: ", nrow(atac_raw)),
    paste0("celltypes: ",
           paste(sort(unique(as.character(meta[[opt$celltype_col]]))), collapse = ","))
  ),
  file.path(opt$out_dir, "summary.txt")
)

message("[seurat_to_anndata] Done -> ", opt$out_dir)
