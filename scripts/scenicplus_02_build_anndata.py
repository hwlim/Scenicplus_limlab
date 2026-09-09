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
    p.add_argument("--reduction", default="",
                   help="Seurat reduction to expose as obsm['X_umap']. Empty "
                        "means compute one here, WITHOUT batch correction.")
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
    #
    # Every reduction is imported under obsm["X_<name>"]; exactly one of them
    # becomes X_umap, the key every plotting call actually reads. That choice
    # used to be made by name prefix -- anything starting with "umap" -- which
    # matches NONE of the names Seurat produces in practice: wnn.umap,
    # rna.umap, umap.harmony. An object carrying three curated layouts
    # therefore fell through to the fallback below and got a fresh PCA-UMAP
    # with no batch correction, which on a multi-sample integrated object is
    # the wrong layout and reports itself only as one line in this log.
    chosen = args.reduction.strip()
    imported = []
    reductions_file = in_dir / "reductions.txt"
    if reductions_file.exists():
        reductions = [r for r in reductions_file.read_text().splitlines() if r]
        for rname in reductions:
            emb_path = in_dir / f"embedding_{rname}.tsv"
            if not emb_path.exists():
                continue
            emb = pd.read_csv(emb_path, sep="\t").set_index("barcode")
            emb = emb.loc[adata.obs_names]
            key = f"X_{rname.lower()}"
            adata.obsm[key] = emb.to_numpy()
            imported.append(rname)
            print(f"[build_anndata] Imported reduction '{rname}' -> obsm['{key}']")

    if chosen:
        # Assigned AFTER the loop, so a reduction literally named "umap" (which
        # lands on X_umap by the rule above) cannot outrank the chosen one.
        key = f"X_{chosen.lower()}"
        if key not in adata.obsm:
            raise SystemExit(
                f"[build_anndata] input.reduction {chosen!r} was not exported "
                f"by step 01.\n"
                f"  exported: {', '.join(imported) if imported else '(none)'}\n"
                f"  Set input.reduction to one of those, or leave it empty to "
                f"compute a UMAP here.")
        adata.obsm["X_umap"] = adata.obsm[key]
        print(f"[build_anndata] X_umap <- '{chosen}'  (input.reduction)")

    # Last resort: nothing was chosen and nothing already occupies X_umap.
    # Say plainly what this produces, because the layout is an unintegrated
    # PCA-UMAP of the RNA matrix -- reasonable for one sample, misleading for
    # anything integrated, and indistinguishable from a curated layout once it
    # is sitting in X_umap.
    if "X_umap" not in adata.obsm:
        print("[build_anndata] WARNING: no reduction chosen and no X_umap — "
              "computing one from the RNA matrix, with NO batch correction.")
        if imported:
            print("[build_anndata]          this object HAS reductions: "
                  + ", ".join(imported))
            print("[build_anndata]          set input.reduction to the one you "
                  "want plotted.")
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
