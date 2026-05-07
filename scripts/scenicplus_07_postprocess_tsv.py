#!/usr/bin/env python
"""
Read the final scplusmdata.h5mu and write a comprehensive set of TSVs:

  eRegulons_direct.tsv          TF-region-gene with importance / rho / triplet rank
  eRegulons_extended.tsv        same, extended motif-to-TF annotations
  eRegulons_combined.tsv        union of direct + extended
  TF_summary.tsv                per-TF: n_targets, n_regions, mean importance, ...
  AUC_gene_per_cell.tsv         eRegulon target-gene AUC per cell
  AUC_region_per_cell.tsv       eRegulon target-region AUC per cell
  RSS_per_<celltype>.tsv        eRegulon specificity score per cell type
  eRegulons_per_celltype.tsv    top-N eRegulons per cell type by RSS
"""
import argparse
from pathlib import Path

import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd
import yaml


def per_tf_summary(meta: pd.DataFrame) -> pd.DataFrame:
    if meta.empty:
        return pd.DataFrame()
    g = meta.groupby("TF", observed=True)
    out = pd.DataFrame({
        "n_target_genes":   g["Gene"].nunique(),
        "n_target_regions": g["Region"].nunique(),
        "n_triplets":       g.size(),
        "mean_importance_TF2G": g["importance_TF2G"].mean(),
        "mean_importance_R2G": g["importance_R2G"].mean(),
        "mean_rho_TF2G":    g["rho_TF2G"].mean(),
        "mean_rho_R2G":     g["rho_R2G"].mean(),
        "mean_triplet_rank": g["triplet_rank"].mean(),
        "n_eRegulons":      g["eRegulon_name"].nunique(),
    }).reset_index()
    out = out.sort_values("n_target_genes", ascending=False)
    return out


def concat_auc(scplus_mdata, gene_keys, region_keys):
    def _stack(keys):
        frames = []
        for k in keys:
            if k in scplus_mdata.mod:
                a = scplus_mdata[k]
                df = pd.DataFrame(
                    a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X),
                    index=a.obs_names, columns=a.var_names,
                )
                frames.append(df)
        if not frames:
            return pd.DataFrame()
        return pd.concat(frames, axis=1)
    return _stack(gene_keys), _stack(region_keys)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--scplus_mdata", required=True)
    p.add_argument("--config", required=True)
    p.add_argument("--out_dir", required=True)
    args = p.parse_args()

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)
    celltype_col = cfg["input"]["celltype_column"]
    top_n = cfg["visualization"]["top_n_eRegulons_per_celltype"]

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    md = mu.read(args.scplus_mdata)

    # eRegulon metadata
    direct = md.uns.get("direct_e_regulon_metadata", pd.DataFrame())
    extended = md.uns.get("extended_e_regulon_metadata", pd.DataFrame())
    if isinstance(direct, dict):
        direct = pd.DataFrame(direct)
    if isinstance(extended, dict):
        extended = pd.DataFrame(extended)

    direct.to_csv(out_dir / "eRegulons_direct.tsv", sep="\t", index=False)
    extended.to_csv(out_dir / "eRegulons_extended.tsv", sep="\t", index=False)

    combined = pd.concat([direct.assign(annotation_source="direct"),
                          extended.assign(annotation_source="extended")],
                         ignore_index=True)
    combined.to_csv(out_dir / "eRegulons_combined.tsv", sep="\t", index=False)
    per_tf_summary(combined).to_csv(out_dir / "TF_summary.tsv", sep="\t", index=False)

    # AUC matrices (cells x eRegulons)
    auc_gene, auc_region = concat_auc(
        md,
        gene_keys=["direct_gene_based_AUC", "extended_gene_based_AUC"],
        region_keys=["direct_region_based_AUC", "extended_region_based_AUC"],
    )
    if not auc_gene.empty:
        auc_gene.to_csv(out_dir / "AUC_gene_per_cell.tsv", sep="\t")
    if not auc_region.empty:
        auc_region.to_csv(out_dir / "AUC_region_per_cell.tsv", sep="\t")

    # RSS
    obs = md.obs.copy()
    candidate_keys = [f"scRNA_counts:{celltype_col}", celltype_col]
    rss_var = next((k for k in candidate_keys if k in obs.columns), None)

    if rss_var is not None and not auc_gene.empty:
        from scenicplus.RSS import regulon_specificity_scores
        try:
            rss = regulon_specificity_scores(
                scplus_mudata=md,
                variable=rss_var,
                modalities=["direct_gene_based_AUC", "extended_gene_based_AUC"],
            )
            rss.to_csv(out_dir / "RSS_per_celltype.tsv", sep="\t")

            top_per_ct = []
            for ct in rss.index:
                top = rss.loc[ct].sort_values(ascending=False).head(top_n)
                for rank, (er, score) in enumerate(top.items(), start=1):
                    top_per_ct.append({
                        "celltype": ct, "rank": rank,
                        "eRegulon": er, "RSS": score,
                    })
            pd.DataFrame(top_per_ct).to_csv(
                out_dir / "eRegulons_per_celltype.tsv", sep="\t", index=False)
        except Exception as e:
            print(f"[postprocess] RSS failed: {e}")
    else:
        print(f"[postprocess] No cell-type column matched ({candidate_keys}) — skipping RSS")

    print(f"[postprocess] TSVs written to: {out_dir}")


if __name__ == "__main__":
    main()
