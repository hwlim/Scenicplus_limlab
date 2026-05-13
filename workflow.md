# SCENIC+ pipeline — step-by-step

End-to-end flow from a Seurat `.rds` (RNA + peaks) to a TF-gene regulatory
network. Each step is a stage in `scripts/scenicplus_run_pipeline.sh`; the
per-step script it invokes lives in the same `scripts/` folder
(`scenicplus_01_*.R`, `scenicplus_02_*.py`, …). Configurable knobs come from
`config/config.yaml`. Steps skip themselves when their sentinel output exists
and its `.cfgsha` sidecar still matches the relevant config slice.

```
Seurat .rds
   │  (01) scenicplus_01_seurat_to_anndata.R   — Seurat → mtx artifacts
   ▼
interim/seurat_export/
   │  (02) scenicplus_02_build_anndata.py      — mtx → AnnData
   │  (03) scenicplus_03_create_cistopic.py    — peak mtx → CistopicObject
   ▼                                ▼
interim/rna.h5ad        interim/cistopic_obj.pkl
                                   │  (04) scenicplus_04_topic_modeling.py — LDA + best-model
                                   ▼
                       interim/cistopic_obj_with_topics.pkl
                                   │  (05) scenicplus_05_region_sets.py    — DARs + topic beds
                                   ▼
                       interim/region_sets/{DARs_*,Topics_otsu,Topics_top_3k}/*.bed
   │                               │
   └──────────────┬────────────────┘
                  │  (06) scenicplus_06_init_inner.sh        — scaffold + patch SCENIC+ config
                  ▼
       results/scplus_pipeline/Snakemake/config/config.yaml
                  │  (07) snakemake (SCENIC+ inner pipeline, invoked by driver)
                  ▼
       results/scplus_out/scplusmdata.h5mu
                  │  (08) scenicplus_07_postprocess_tsv.py   — TSV deliverables
                  │  (09) scenicplus_08_visualize.py         — PDF + PNG plots
                  ▼
       results/tables/*.tsv,  results/plots/*.{pdf,png}
```

---

## 01 · `seurat_to_anndata` (R)

Split a Seurat object into Matrix-Market artifacts that the Python steps can
ingest without `SeuratDisk`. Also exports any pre-computed reductions (UMAP /
PCA) so the curated layout flows through.

- **In:** `input.seurat_rds` (must have `RNA` + `peaks` assays and
  `input.celltype_column` in `@meta.data`)
- **Out:** `interim/seurat_export/{rna_raw_counts.mtx, rna_norm.mtx,
  rna_features.tsv, rna_barcodes.tsv, atac_counts.mtx, atac_barcodes.tsv,
  atac_regions.tsv, cell_metadata.tsv, embedding_<reduction>.tsv,
  reductions.txt, summary.txt}`
- **Key options:** `--celltype_col`, `--celltype_scope` (comma-separated;
  if non-empty, the Seurat object is subset to those cell types before export
  — every downstream step then operates on the focused universe)

```r
obj <- readRDS(opt$rds)

# Optional cell-type subsetting (boosts sensitivity for focused analyses).
scope <- trimws(strsplit(opt$celltype_scope, ",", fixed = TRUE)[[1]])
scope <- scope[nchar(scope) > 0]
if (length(scope) > 0) {
  keep <- which(as.character(obj@meta.data[[opt$celltype_col]]) %in% scope)
  obj  <- obj[, keep]
}

rna_raw  <- GetAssayData(obj, assay = "RNA",   slot = "counts")
rna_norm <- GetAssayData(obj, assay = "RNA",   slot = "data")
atac_raw <- GetAssayData(obj, assay = "peaks", slot = "counts")

writeMM(rna_raw,  ".../rna_raw_counts.mtx")
writeMM(rna_norm, ".../rna_norm.mtx")
writeMM(atac_raw, ".../atac_counts.mtx")

# normalize peak names to "chr:start-end"
peak_names <- gsub("[-_]", ":", rownames(atac_raw))
peak_names <- sub("^([^:]+):([^:]+):([^:]+)$", "\\1:\\2-\\3", peak_names)
```

---

## 02 · `build_anndata` (Python / scanpy)

Pack the RNA mtx into an `AnnData`. Normalized log-counts go to `.X`, raw
counts to `.raw`. Seurat reductions are re-attached to `.obsm["X_umap"]` etc.;
if no UMAP exists, a quick scanpy UMAP is computed for plotting.

- **In:** `interim/seurat_export/`
- **Out:** `interim/rna.h5ad`
- **Key options:** `--celltype_col`

```python
adata     = ad.AnnData(X=norm.T.tocsr(), obs=meta, var=pd.DataFrame(index=genes))
adata_raw = ad.AnnData(X=raw.T.tocsr(),  obs=meta, var=pd.DataFrame(index=genes))
adata.var_names_make_unique();  adata_raw.var_names_make_unique()
adata.raw = adata_raw

# Re-import Seurat embeddings if present, else compute a fallback UMAP
adata.obsm["X_umap"] = emb.loc[adata.obs_names].to_numpy()
# fallback:
sc.pp.highly_variable_genes(adata_tmp, flavor="seurat", n_top_genes=2000)
sc.pp.scale(adata_tmp, max_value=10, zero_center=False)
sc.tl.pca(adata_tmp); sc.pp.neighbors(adata_tmp); sc.tl.umap(adata_tmp)

adata.write_h5ad(args.out_h5ad)
```

---

## 03 · `create_cistopic` (pycisTopic)

Build a `CistopicObject` from the peak-count matrix (matrix-only path —
fragment files are not required). Cell metadata is attached so DARs can be
called per cell type later.

- **In:** `atac_counts.mtx`, `atac_barcodes.tsv`, `atac_regions.tsv`,
  `cell_metadata.tsv`
- **Out:** `interim/cistopic_obj.pkl`
- **Key options:** `--project` (label baked into the object)

```python
from pycisTopic.cistopic_class import create_cistopic_object

counts_df = pd.DataFrame.sparse.from_spmatrix(counts, index=regions, columns=barcodes)

cto = create_cistopic_object(
    fragment_matrix = counts_df,
    cell_names      = barcodes,
    region_names    = regions,
    project         = "scenicplus_run",
)
cto.add_cell_data(meta)
pickle.dump(cto, open(args.out_pkl, "wb"))
```

---

## 04 · `topic_modeling` (pycisTopic LDA)

Fit collapsed-Gibbs LDA across a sweep of topic counts, then auto-select the
best model by pycisTopic's combined score. The selected model is attached to
the object via `add_LDA_model()`.

- **In:** `interim/cistopic_obj.pkl`
- **Out:** `interim/cistopic_obj_with_topics.pkl`
- **Key options** (from `config.cistopic`):
  `n_topics`, `n_iter`, `alpha`, `alpha_by_topic`, `eta`, `eta_by_topic`,
  `random_state`; CLI: `--n_cpu`, `--tmp_dir`

```python
from pycisTopic.lda_models import run_cgs_models, evaluate_models

models = run_cgs_models(
    cto,
    n_topics       = ct_cfg["n_topics"],          # e.g. [2,5,10,20,30,40,50]
    n_cpu          = args.n_cpu,
    n_iter         = ct_cfg["n_iter"],            # 150
    random_state   = ct_cfg["random_state"],
    alpha          = ct_cfg["alpha"],             # 50
    alpha_by_topic = ct_cfg["alpha_by_topic"],    # True
    eta            = ct_cfg["eta"],               # 0.1
    eta_by_topic   = ct_cfg["eta_by_topic"],      # False
    save_path      = args.tmp_dir,
)
best = evaluate_models(models, select_model=None, return_model=True, plot=False)
cto.add_LDA_model(best)
```

---

## 05 · `region_sets` (DARs + topic binarization)

Produce the BED folder SCENIC+ motif-enrichment consumes. Two flavours of
topic-based region sets and one DAR set per cell type.

- **In:** `interim/cistopic_obj_with_topics.pkl`
- **Out:** `interim/region_sets/{Topics_otsu,Topics_top_3k,DARs_<celltype_col>}/*.bed`
- **Key options** (from `config.cistopic`): `dar_adjpval_thr`, `dar_log2fc_thr`;
  CLI: `--n_cpu`

```python
from pycisTopic.topic_binarization import binarize_topics
from pycisTopic.diff_features    import (
    impute_accessibility, normalize_scores, find_diff_features,
)

# Topic-based regions
region_bin_otsu = binarize_topics(cto, method="otsu")
region_bin_top  = binarize_topics(cto, method="ntop", ntop=3000)

# DARs per cell type (1-vs-rest)
imputed  = impute_accessibility(cto, selected_cells=None, selected_regions=None,
                                scale_factor=10**6, chunk_size=20000)
norm_imp = normalize_scores(imputed, scale_factor=10**4)
markers  = find_diff_features(
    cto, norm_imp, variable=celltype_col,
    adjpval_thr = ct_cfg["dar_adjpval_thr"],   # 0.05
    log2fc_thr  = ct_cfg["dar_log2fc_thr"],    # 0.5
    n_cpu       = args.n_cpu,
    split_pattern="-",
)
# each set written as chrom\tstart\tend\tname BED file
```

---

## 06 · `init_scenicplus` (scaffold + patch)

Run `scenicplus init_snakemake` to scaffold the inner SCENIC+ snakemake
project, then overwrite its `config.yaml` with paths and parameters from our
top-level config so the inner pipeline points at our cistopic object, AnnData,
region sets, and cisTarget databases.

- **In:** `interim/cistopic_obj_with_topics.pkl`, `interim/rna.h5ad`,
  `interim/region_sets/`, `config/config.yaml`
- **Out:** `results/scplus_pipeline/Snakemake/config/config.yaml`
- **Key options:** all of `config.scenicplus.*`, `config.grn.*`,
  `input.{species,assembly,biomart_host,ctx_db,dem_db,motif_annotations}`

```python
subprocess.run(["scenicplus", "init_snakemake", "--out_dir", out_dir], check=True)

sc_cfg["input_data"] = {
    "cisTopic_obj_fname":         cistopic_obj,
    "GEX_anndata_fname":          adata,
    "region_set_folder":          region_sets_dir,
    "ctx_db_fname":               cfg["input"]["ctx_db"],
    "dem_db_fname":               cfg["input"]["dem_db"],
    "path_to_motif_annotations":  cfg["input"]["motif_annotations"],
}
sc_cfg["params_general"]          = {"temp_dir": ..., "n_cpu": ..., "seed": ...}
sc_cfg["params_data_preparation"] = {  # search_space_*, biomart_host, species, ...
    "is_multiome": True, "search_space_upstream": "1000 150000", ...
}
sc_cfg["params_motif_enrichment"] = {  # ctx/dem thresholds, motif-similarity FDR
    "ctx_auc_threshold": 0.005, "ctx_nes_threshold": 3.0, ...
}
sc_cfg["params_inference"]        = {  # GBM/RF, rho_threshold, min_target_genes, ...
    "tf_to_gene_importance_method": "GBM", "min_target_genes": 10, ...
}
yaml.safe_dump(sc_cfg, open(scplus_cfg_path, "w"), sort_keys=False)
```

---

## 07 · `run_scenicplus` (the SCENIC+ inner pipeline)

Hand off to the scaffolded SCENIC+ snakemake project. This is where the heavy
lifting happens: metacells → motif enrichment (cisTarget + DEM) → cistromes →
TF-to-gene + region-to-gene adjacencies → eRegulon construction → AUCell.

- **In:** `results/scplus_pipeline/Snakemake/config/config.yaml`
- **Out:** `results/scplus_out/scplusmdata.h5mu` (plus per-stage intermediates:
  `tf_to_gene_adj.tsv`, `region_to_gene_adj.tsv`, `eRegulons_direct.tsv`,
  `eRegulons_extended.tsv`, `AUCell_{direct,extended}.h5mu`,
  `cistromes_{direct,extended}.h5ad`, `dem_results.hdf5`, `ctx_results.hdf5`)
- **Key options:** `threads` = `config.resources.n_cpu`

```bash
cd results/scplus_pipeline/Snakemake
snakemake --cores ${N_CPU} --rerun-incomplete
```

---

## 08 · `postprocess_tsv`

Read the final `scplusmdata.h5mu` and dump comprehensive TSVs. `RSS` is
computed via SCENIC+'s `regulon_specificity_scores`.

- **In:** `results/scplus_out/scplusmdata.h5mu`
- **Out:** `results/tables/{eRegulons_direct, eRegulons_extended,
  eRegulons_combined, TF_summary, AUC_gene_per_cell, AUC_region_per_cell,
  RSS_per_celltype, eRegulons_per_celltype}.tsv`
- **Key options:** `config.visualization.top_n_eRegulons_per_celltype`

```python
import mudata as mu
from scenicplus.RSS import regulon_specificity_scores

md       = mu.read(args.scplus_mdata)
direct   = md.uns["direct_e_regulon_metadata"]      # TF, Region, Gene, importance, rho, triplet_rank, ...
extended = md.uns["extended_e_regulon_metadata"]
combined = pd.concat([direct.assign(annotation_source="direct"),
                      extended.assign(annotation_source="extended")])

# per-cell AUC matrices (cells × eRegulons)
for k in ["direct_gene_based_AUC", "extended_gene_based_AUC",
          "direct_region_based_AUC", "extended_region_based_AUC"]:
    a = md[k]                                       # AnnData
    pd.DataFrame(a.X.toarray(), index=a.obs_names, columns=a.var_names) \
      .to_csv(...)

# RSS per cell type
rss = regulon_specificity_scores(
    scplus_mudata = md,
    variable      = celltype_col,
    modalities    = ["direct_gene_based_AUC", "extended_gene_based_AUC"],
)
```

---

## 09 · `visualize`

Render one-plot-per-file PDFs *and* PNGs. UMAP, RSS rank plot, and SCENIC+'s
`heatmap_dotplot` come from the SCENIC+ API; bar/density/network plots are
matplotlib + networkx.

- **In:** `results/scplus_out/scplusmdata.h5mu`
- **Out:** `results/plots/*.pdf` + `*.png` (UMAPs, RSS, heatmap-dotplots,
  TF target-count, importance density, eGRN network)
- **Key options:** `config.visualization.{top_n_eRegulons_per_celltype,
  network_top_n_TFs, network_top_n_targets_per_TF}`

```python
from scenicplus.RSS              import regulon_specificity_scores, plot_rss
from scenicplus.plotting.dotplot import heatmap_dotplot

# RSS rank plot
rss = regulon_specificity_scores(scplus_mudata=md, variable=celltype_col,
                                 modalities=["direct_gene_based_AUC",
                                             "extended_gene_based_AUC"])
plot_rss(data_matrix=rss, top_n=viz["top_n_eRegulons_per_celltype"], num_columns=cols)

# eRegulon AUC heatmap-dotplot
heatmap_dotplot(
    scplus_mudata        = md,
    color_modality       = "direct_gene_based_AUC",
    size_modality        = "direct_region_based_AUC",
    group_variable       = celltype_col,
    eRegulon_metadata_key= "direct_e_regulon_metadata",
    color_feature_key    = "Gene_signature_name",
    size_feature_key     = "Region_signature_name",
    feature_name_key     = "eRegulon_name",
    sort_data_by         = "direct_gene_based_AUC",
    orientation          = "horizontal",
)

# UMAP per top eRegulon (eRegulon AUC values are the colour)
sc.pl.umap(er_adata, color=eRegulon_name, cmap="viridis", ...)
```

---

## Configuration cheat-sheet

| Section                 | Drives                                                  |
| ----------------------- | ------------------------------------------------------- |
| `input.*`               | RDS path, cell-type column, optional `celltype_scope`, species/assembly, cisTarget DBs, motif annotations |
| `output.{root,tmp}`     | output root and scratch dirs                            |
| `cistopic.*`            | LDA hyperparameters (step 04) + DAR thresholds (step 05) |
| `scenicplus.*`          | search-space, motif-enrichment, ctx/DEM thresholds (step 06/07) |
| `grn.*`                 | TF/region-to-gene importance method, `rho_threshold`, `min_target_genes` (step 07) |
| `resources.{n_cpu,seed}`| parallelism + reproducibility                           |
| `visualization.*`       | top-N controls for tables and plots (step 08/09)        |
