#!/usr/bin/env python
"""
Assemble the AnnData (.h5ad) SCENIC+ expects from the artifacts written by
the R step. Stores raw counts in .raw and normalized log counts in .X, and
copies any pre-computed Seurat reductions (UMAP / PCA / etc.) into .obsm so
the curated cell-type layout is preserved downstream.
"""
import argparse
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp


def load_mtx(path: Path) -> sp.csr_matrix:
    return scipy.io.mmread(str(path)).tocsr()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--in_dir", required=True)
    p.add_argument("--out_h5ad", required=True)
    p.add_argument("--celltype_col", required=True)
    args = p.parse_args()

    in_dir = Path(args.in_dir)

    raw = load_mtx(in_dir / "rna_raw_counts.mtx")           # genes x cells
    norm = load_mtx(in_dir / "rna_norm.mtx")                # genes x cells
    genes = pd.read_csv(in_dir / "rna_features.tsv", header=None)[0].tolist()
    cells = pd.read_csv(in_dir / "rna_barcodes.tsv", header=None)[0].tolist()
    meta = pd.read_csv(in_dir / "cell_metadata.tsv", sep="\t")
    meta = meta.set_index("barcode").loc[cells]

    # SCENIC+ expects cells x genes
    adata_raw = ad.AnnData(
        X=raw.T.tocsr(),
        obs=meta.copy(),
        var=pd.DataFrame(index=genes),
    )
    adata = ad.AnnData(
        X=norm.T.tocsr(),
        obs=meta.copy(),
        var=pd.DataFrame(index=genes),
    )
    adata.var_names_make_unique()
    adata_raw.var_names_make_unique()
    adata.raw = adata_raw

    # Re-attach Seurat reductions if exported by the R step.
    reductions_file = in_dir / "reductions.txt"
    if reductions_file.exists():
        reductions = [r for r in reductions_file.read_text().splitlines() if r]
        for rname in reductions:
            emb_path = in_dir / f"embedding_{rname}.tsv"
            if not emb_path.exists():
                continue
            emb = pd.read_csv(emb_path, sep="\t").set_index("barcode")
            emb = emb.loc[adata.obs_names]
            key = "X_umap" if rname.lower().startswith("umap") else f"X_{rname.lower()}"
            adata.obsm[key] = emb.to_numpy()
            print(f"[build_anndata] Imported reduction '{rname}' -> obsm['{key}']")

    # If no UMAP came from Seurat, compute a quick one so visualization works.
    if "X_umap" not in adata.obsm:
        print("[build_anndata] No UMAP in Seurat object — computing one for plots.")
        try:
            adata_tmp = adata.copy()
            sc.pp.highly_variable_genes(adata_tmp, flavor="seurat", n_top_genes=2000)
            sc.pp.scale(adata_tmp, max_value=10, zero_center=False)
            sc.tl.pca(adata_tmp, n_comps=min(50, adata_tmp.n_vars - 1))
            sc.pp.neighbors(adata_tmp)
            sc.tl.umap(adata_tmp)
            adata.obsm["X_umap"] = adata_tmp.obsm["X_umap"]
        except Exception as e:
            print(f"[build_anndata] UMAP fallback failed: {e}")

    adata.obs[args.celltype_col] = adata.obs[args.celltype_col].astype("category")

    out = Path(args.out_h5ad)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(out.suffix + ".partial")
    adata.write_h5ad(tmp)
    tmp.replace(out)
    print(f"[build_anndata] Wrote: {out}  (n_obs={adata.n_obs}, n_vars={adata.n_vars})")


if __name__ == "__main__":
    main()
