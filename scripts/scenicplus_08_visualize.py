#!/usr/bin/env python
"""
Generate one-plot-per-file PDFs *and* PNGs from the SCENIC+ MuData.

Plots:
  01_umap_celltype.{pdf,png}            UMAP coloured by cell type
  02_umap_eRegulon_<name>.{pdf,png}     UMAP per top eRegulon AUC
  03_rss_per_celltype.{pdf,png}         RSS rank plot
  04_heatmap_dotplot_direct.{pdf,png}   eRegulon AUC heatmap-dotplot (direct)
  05_heatmap_dotplot_extended.{pdf,png}                          (extended)
  06_TF_target_count.{pdf,png}          n_target_genes per TF (bar)
  07_TF_importance_distribution.{pdf,png}  importance density per TF (top N)
  08_eGRN_network_top<N>.{pdf,png}      TF-target network for top-N TFs
"""
import argparse
import warnings
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mudata as mu
import numpy as np
import pandas as pd
import scanpy as sc
import yaml

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0


def save(fig, out_dir: Path, name: str):
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(out_dir / f"{name}.png", bbox_inches="tight", dpi=200)
    plt.close(fig)


def make_eRegulon_adata(md):
    parts = []
    for k in ["direct_gene_based_AUC", "extended_gene_based_AUC"]:
        if k in md.mod:
            parts.append(md[k])
    if not parts:
        return None
    a = ad.concat(parts, axis=1, merge="unique")
    a.obs = md.obs.loc[a.obs_names]
    return a


def umap_celltype(adata, ct_col, out_dir):
    if "X_umap" not in adata.obsm:
        sc.pp.neighbors(adata, use_rep="X")
        sc.tl.umap(adata)
    fig, ax = plt.subplots(figsize=(7, 6))
    sc.pl.umap(adata, color=ct_col, ax=ax, show=False, frameon=False, legend_loc="on data")
    save(fig, out_dir, "01_umap_celltype")


def umap_eRegulons(adata, eRegulons, out_dir):
    for er in eRegulons:
        if er not in adata.var_names:
            continue
        fig, ax = plt.subplots(figsize=(6, 5))
        sc.pl.umap(adata, color=er, ax=ax, show=False, frameon=False, cmap="viridis")
        save(fig, out_dir, f"02_umap_eRegulon_{er.replace('/', '_').replace(' ', '_')}")


def rss_plot(rss, out_dir, top_n):
    from scenicplus.RSS import plot_rss
    n_groups = rss.shape[0]
    cols = min(5, max(1, n_groups))
    fig = plot_rss(data_matrix=rss, top_n=top_n, num_columns=cols,
                   figsize=(4 * cols, 3 * int(np.ceil(n_groups / cols))))
    if fig is None:
        fig = plt.gcf()
    save(fig, out_dir, "03_rss_per_celltype")


def heatmap_dotplot(md, ct_col, out_dir, source):
    from scenicplus.plotting.dotplot import heatmap_dotplot
    color_mod = f"{source}_gene_based_AUC"
    size_mod  = f"{source}_region_based_AUC"
    meta_key  = f"{source}_e_regulon_metadata"
    if color_mod not in md.mod or size_mod not in md.mod or meta_key not in md.uns:
        print(f"[viz] heatmap_dotplot skipped for {source}")
        return
    fig = heatmap_dotplot(
        scplus_mudata=md,
        color_modality=color_mod, size_modality=size_mod,
        group_variable=ct_col,
        eRegulon_metadata_key=meta_key,
        color_feature_key="Gene_signature_name",
        size_feature_key="Region_signature_name",
        feature_name_key="eRegulon_name",
        sort_data_by=color_mod,
        orientation="horizontal",
        figsize=(16, 6),
    )
    if fig is None:
        fig = plt.gcf()
    save(fig, out_dir, f"04_heatmap_dotplot_{source}" if source == "direct"
         else f"05_heatmap_dotplot_{source}")


def tf_target_count(meta, out_dir, top_n=40):
    if meta.empty:
        return
    counts = meta.groupby("TF", observed=True)["Gene"].nunique() \
                 .sort_values(ascending=False).head(top_n)
    fig, ax = plt.subplots(figsize=(8, max(4, 0.25 * len(counts))))
    counts[::-1].plot(kind="barh", ax=ax, color="steelblue")
    ax.set_xlabel("# target genes")
    ax.set_ylabel("")
    ax.set_title(f"Top {len(counts)} TFs by # target genes")
    save(fig, out_dir, "06_TF_target_count")


def tf_importance_density(meta, out_dir, top_n=10):
    if meta.empty:
        return
    top_tfs = meta.groupby("TF", observed=True)["Gene"].nunique() \
                  .sort_values(ascending=False).head(top_n).index
    fig, ax = plt.subplots(figsize=(8, 5))
    for tf in top_tfs:
        sub = meta.loc[meta["TF"] == tf, "importance_TF2G"].dropna()
        if sub.empty:
            continue
        sub.plot(kind="kde", ax=ax, label=tf)
    ax.set_xlabel("importance_TF2G")
    ax.set_title(f"TF→gene importance distribution (top {top_n} TFs)")
    ax.legend(fontsize=8, ncol=2)
    save(fig, out_dir, "07_TF_importance_distribution")


def network_plot(meta, out_dir, top_tfs, top_targets):
    if meta.empty:
        return
    try:
        import networkx as nx
    except ImportError:
        print("[viz] networkx missing — skipping network plot")
        return

    tf_size = meta.groupby("TF", observed=True)["Gene"].nunique() \
                  .sort_values(ascending=False).head(top_tfs).index.tolist()
    sub = meta[meta["TF"].isin(tf_size)].copy()
    sub["abs_imp"] = sub["importance_TF2G"].abs()
    sub = sub.sort_values("abs_imp", ascending=False) \
             .groupby("TF", observed=True).head(top_targets)

    g = nx.DiGraph()
    for _, r in sub.iterrows():
        g.add_edge(r["TF"], r["Gene"], weight=float(r["importance_TF2G"]))

    pos = nx.spring_layout(g, seed=0, k=0.5)
    fig, ax = plt.subplots(figsize=(12, 10))
    tf_nodes = [n for n in g.nodes if n in tf_size]
    target_nodes = [n for n in g.nodes if n not in tf_size]
    nx.draw_networkx_nodes(g, pos, nodelist=tf_nodes, node_size=700,
                           node_color="tomato", ax=ax, label="TF")
    nx.draw_networkx_nodes(g, pos, nodelist=target_nodes, node_size=80,
                           node_color="lightblue", ax=ax, label="Target")
    nx.draw_networkx_edges(g, pos, alpha=0.3, arrowsize=8, ax=ax)
    nx.draw_networkx_labels(g, pos,
                            labels={n: n for n in tf_nodes},
                            font_size=10, ax=ax)
    ax.set_axis_off()
    ax.legend(loc="upper right")
    ax.set_title(f"eGRN: top {top_tfs} TFs (top {top_targets} targets each)")
    save(fig, out_dir, f"08_eGRN_network_top{top_tfs}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--scplus_mdata", required=True)
    p.add_argument("--config", required=True)
    p.add_argument("--out_dir", required=True)
    args = p.parse_args()

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)
    celltype_col = cfg["input"]["celltype_column"]
    viz = cfg["visualization"]
    out_dir = Path(args.out_dir)

    md = mu.read(args.scplus_mdata)

    candidate_keys = [f"scRNA_counts:{celltype_col}", celltype_col]
    rss_var = next((k for k in candidate_keys if k in md.obs.columns), None)

    # Build eRegulon AnnData for UMAP/RSS
    er_adata = make_eRegulon_adata(md)
    if er_adata is not None and rss_var is not None:
        er_adata.obs[rss_var] = er_adata.obs[rss_var].astype("category")
        umap_celltype(er_adata, rss_var, out_dir)

    # UMAP per top eRegulon (chosen by RSS if possible, otherwise by AUC variance)
    top_eRegulons = []
    if er_adata is not None:
        if rss_var is not None:
            try:
                from scenicplus.RSS import regulon_specificity_scores
                rss = regulon_specificity_scores(
                    scplus_mudata=md, variable=rss_var,
                    modalities=["direct_gene_based_AUC", "extended_gene_based_AUC"],
                )
                top_eRegulons = sorted({
                    er for ct in rss.index
                    for er in rss.loc[ct].sort_values(ascending=False)
                                 .head(viz["top_n_eRegulons_per_celltype"]).index
                })
                rss_plot(rss, out_dir, top_n=viz["top_n_eRegulons_per_celltype"])
            except Exception as e:
                print(f"[viz] RSS computation failed: {e}")

        if not top_eRegulons:
            var = pd.DataFrame(er_adata.X, index=er_adata.obs_names,
                               columns=er_adata.var_names).var()
            top_eRegulons = var.sort_values(ascending=False).head(20).index.tolist()
        umap_eRegulons(er_adata, top_eRegulons, out_dir)

    # heatmap-dotplots
    if rss_var is not None:
        heatmap_dotplot(md, rss_var, out_dir, source="direct")
        heatmap_dotplot(md, rss_var, out_dir, source="extended")

    # eRegulon metadata-derived plots
    direct = md.uns.get("direct_e_regulon_metadata", pd.DataFrame())
    extended = md.uns.get("extended_e_regulon_metadata", pd.DataFrame())
    combined = pd.concat([direct, extended], ignore_index=True) \
        if (len(direct) + len(extended)) else pd.DataFrame()
    tf_target_count(combined, out_dir, top_n=40)
    tf_importance_density(combined, out_dir, top_n=viz["network_top_n_TFs"])
    network_plot(combined, out_dir,
                 top_tfs=viz["network_top_n_TFs"],
                 top_targets=viz["network_top_n_targets_per_TF"])

    print(f"[viz] Plots written to: {out_dir}")


if __name__ == "__main__":
    main()
