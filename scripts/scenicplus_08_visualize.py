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

The layout the UMAP panels are drawn on is `input.reduction`, read from step
01's seurat_export/embedding_<name>.tsv. With no reduction set they fall back
to a UMAP of eRegulon activity, which is a different picture from the Seurat
one and not a worse one -- so either way the title names which it is.
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
    """Write one figure as both .pdf and .png, whichever library made it.

    NOT every SCENIC+ plotting function returns a matplotlib Figure.
    scenicplus.plotting.dotplot.heatmap_dotplot builds a PLOTNINE ggplot and
    returns it, and plotnine's saver is .save(), not .savefig():

        AttributeError: 'ggplot' object has no attribute 'savefig'

    which is what step 20 died with, after the whole GRN had been computed.
    scenicplus.RSS.plot_rss by contrast returns None and draws through pyplot,
    so its caller's `plt.gcf()` fallback is right for that one -- the two
    conventions sit side by side in the same package.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    if hasattr(fig, "savefig"):                      # matplotlib Figure
        fig.savefig(out_dir / f"{name}.pdf", bbox_inches="tight")
        fig.savefig(out_dir / f"{name}.png", bbox_inches="tight", dpi=200)
        plt.close(fig)
    elif hasattr(fig, "save"):                       # plotnine ggplot
        # verbose=False silences "Saving 16 x 6 in image"; limitsize=False
        # because heatmap_dotplot's figsize scales with the number of
        # eRegulons and can exceed plotnine's 25-inch guard on a real run.
        fig.save(out_dir / f"{name}.pdf", verbose=False, limitsize=False)
        fig.save(out_dir / f"{name}.png", dpi=200, verbose=False, limitsize=False)
    else:
        raise TypeError(
            f"cannot save a {type(fig).__module__}.{type(fig).__name__}: it has "
            f"neither .savefig() (matplotlib) nor .save() (plotnine)")


def region_overlap(meta, out_dir, source, top_n):
    """Pairwise Jaccard overlap of eRegulons' TARGET REGION sets.

    What it answers: how much do these regulons actually differ? Two eRegulons
    sharing most of their regions are one finding reported twice -- co-binding
    TFs, or a motif matched by a family -- and the eRegulon tables alone cannot
    show that, because each row is a triplet and the redundancy only appears
    when the sets are compared.

    BUILT HERE RATHER THAN CALLED. `scenicplus.plotting.correlation_plot`
    ships `jaccard_heatmap`, but it takes the LEGACY `SCENICPLUS` class while
    this pipeline produces MuData -- which is why `heatmap_dotplot` is called
    with `scplus_mudata=`. Constructing a legacy object to reach one plotting
    function would put a second representation of the run in the codebase, to
    be kept in step forever, for a Jaccard matrix that is a groupby and a
    pairwise loop. The region sets are already in the metadata this rule reads.

    TOP-N BY REGION COUNT, and the cut is stated on the figure. All eRegulons
    can be several hundred, which is an unreadable heatmap and a slow one --
    the matrix is quadratic. Ranking by set size keeps the regulons with
    something to overlap; a random or alphabetical cut would not.
    """
    if meta is None or len(meta) == 0 or "Region" not in meta.columns:
        print(f"[viz] region_overlap skipped for {source}: no region metadata")
        return
    sets = (meta.groupby("eRegulon_name", observed=True)["Region"]
                .agg(lambda x: frozenset(x)))
    sets = sets[sets.map(len) > 0]
    if len(sets) < 2:
        print(f"[viz] region_overlap skipped for {source}: "
              f"{len(sets)} eRegulon(s), nothing to compare")
        return
    order = sets.map(len).sort_values(ascending=False)
    keep = order.head(top_n).index.tolist()
    sets = sets[keep]

    import numpy as np
    names = list(sets.index)
    n = len(names)
    m = np.eye(n)
    for i in range(n):
        a = sets.iloc[i]
        for j in range(i + 1, n):
            b = sets.iloc[j]
            inter = len(a & b)
            # Jaccard, and 0/0 cannot arise: empty sets were dropped above.
            m[i, j] = m[j, i] = inter / len(a | b)
    df = pd.DataFrame(m, index=names, columns=names)

    import seaborn as sns
    # Clustered, because the block structure IS the finding -- an alphabetical
    # order hides exactly the groups of redundant regulons this exists to show.
    # A degenerate matrix (every pair identical) makes linkage complain rather
    # than fail, so fall back to the given order instead of losing the figure.
    try:
        g = sns.clustermap(df, cmap="rocket_r", vmin=0, vmax=1,
                           figsize=(min(2 + 0.28 * n, 24), min(2 + 0.28 * n, 24)),
                                                      xticklabels=True, yticklabels=True)
    except Exception as e:
        print(f"[viz] region_overlap: clustering failed ({e}); unclustered")
        g = sns.clustermap(df, cmap="rocket_r", vmin=0, vmax=1,
                           row_cluster=False, col_cluster=False,
                           figsize=(min(2 + 0.28 * n, 24), min(2 + 0.28 * n, 24)),
                                                      xticklabels=True, yticklabels=True)
    g.ax_heatmap.tick_params(labelsize=6)
    # seaborn parks the colourbar top-LEFT by default, where a suptitle lands on
    # top of it and its rotated label runs through the row dendrogram. Both were
    # visible only by opening the emitted PNG -- the same lesson as the
    # 301-megapixel RSS figure that was "successful" for a day. Moved to the
    # bottom-left gutter, under the dendrogram, where nothing else is drawn.
    g.cax.set_position([0.02, 0.06, 0.02, 0.14])
    g.cax.tick_params(labelsize=7)
    g.cax.set_ylabel("Jaccard", fontsize=8)
    total = len(order)
    # Title on the heatmap axes, not the figure: a suptitle sits above the
    # dendrogram and drifts with figsize, which here scales with the eRegulon
    # count.
    g.ax_heatmap.set_title(
        f"{source} eRegulons: target-region overlap "
        f"({n} of {total} shown, largest region sets first)",
        fontsize=11, pad=12)
    stem = "09_region_overlap_direct" if source == "direct" \
        else "10_region_overlap_extended"
    save(g.figure, out_dir, stem)

    # The matrix itself, because a heatmap is not a number anyone can act on
    # and the project's contract pairs every figure with a data sheet.
    #
    # BESIDE THE FIGURE, not in `tsv/`. That directory is R19's output and R19
    # is a different rule; two rules writing one directory is the ownership
    # problem the sibling repo states as an invariant, and snakemake would not
    # know about it. A .tsv in plots/ reads slightly oddly and is correct.
    df.to_csv(out_dir / f"09_region_overlap_{source}.tsv", sep="\t")


def make_eRegulon_adata(md, kind="gene"):
    """Cells x eRegulons AUC, direct and extended side by side.

    `kind` selects the feature space: "gene" for target-gene enrichment,
    "region" for target-region. Defaulted to "gene" so the existing per-eRegulon
    UMAP callers are unchanged.
    """
    parts = []
    for k in [f"direct_{kind}_based_AUC", f"extended_{kind}_based_AUC"]:
        if k in md.mod:
            parts.append(md[k])
    if not parts:
        return None
    a = ad.concat(parts, axis=1, merge="unique")
    a.obs = md.obs.loc[a.obs_names]
    return a


def eregulon_tsne(md, ct_col, out_dir, kind, seed):
    """t-SNE of CELLS in eRegulon-activity space, coloured by cell type.

    The paper's Fig 2 view, and the only figure here drawn on a layout of its
    own: every other one uses the Seurat reduction that became `X_umap`. That
    is the point of it -- if cell types separate on eRegulon activity alone,
    the regulons carry the identity, which is a claim the UMAP cannot make
    because the UMAP was computed from expression.

    IT IS ALSO WHY THE CAPTION MUST SAY SO. A reader who has scrolled past six
    figures on one layout will read a seventh as the same coordinates unless
    told otherwise, and `input.reduction` has already cost this project a run
    for exactly that class of confusion.

    SEEDED FROM `resources.seed`. t-SNE is stochastic, and an unseeded one here
    would quietly undo the reproducibility the pinning work established: two
    runs of the same data would differ in a figure while every table matched.
    The rule tracks that key so changing it redraws.
    """
    a = make_eRegulon_adata(md, kind=kind)
    if a is None or a.n_obs < 10 or a.n_vars < 2:
        print(f"[viz] {kind}-based t-SNE skipped: "
              f"{'no AUC modality' if a is None else f'{a.n_obs} cells x {a.n_vars} eRegulons'}")
        return
    if ct_col not in a.obs.columns:
        print(f"[viz] {kind}-based t-SNE skipped: no '{ct_col}' in cell metadata")
        return

    import numpy as np
    import matplotlib.patheffects as pe
    import scanpy as sc
    # sklearn REQUIRES perplexity < n_samples, and raises rather than adjusting.
    # A focused celltype_scope can leave few cells, so clamp instead of failing
    # a whole rule over a plotting parameter.
    perp = max(5.0, min(30.0, (a.n_obs - 1) / 3.0))
    b = a.copy()
    sc.pp.pca(b, n_comps=min(50, b.n_vars - 1, b.n_obs - 1))
    sc.tl.tsne(b, use_rep="X_pca", perplexity=perp, random_state=seed)

    xy = b.obsm["X_tsne"]
    labels = b.obs[ct_col].astype(str)
    order = sorted(labels.unique())
    cmap = plt.get_cmap("tab20")
    fig, ax = plt.subplots(figsize=(8, 7))
    for i, lab in enumerate(order):
        m = (labels == lab).to_numpy()
        ax.scatter(xy[m, 0], xy[m, 1], s=4, linewidths=0,
                   color=cmap(i % 20), label=lab)
    # LABEL EACH CLUSTER IN PLACE, not only in the legend. tab20 repeats after
    # 20 and its adjacent hues are already hard to tell apart -- the sibling
    # repo hit exactly this at 27 labels, where the legend became "the lookup it
    # exists to save". A name at the cluster's median position is readable
    # without matching colours at all, so the legend becomes a fallback for
    # clusters too small or too overlapped to carry text.
    for i, lab in enumerate(order):
        m = (labels == lab).to_numpy()
        cx, cy = np.median(xy[m, 0]), np.median(xy[m, 1])
        ax.text(cx, cy, lab, fontsize=7, ha="center", va="center",
                color="black", zorder=5,
                path_effects=[pe.withStroke(linewidth=2.2, foreground="white")])
    ax.set_xlabel("t-SNE 1")
    ax.set_ylabel("t-SNE 2")
    ax.set_title(f"Cells in {kind}-based eRegulon activity space\n"
                 f"{a.n_obs:,} cells x {a.n_vars} eRegulons, perplexity {perp:g}, "
                 f"seed {seed}", fontsize=10)
    # Legend outside: with 25 cell types it otherwise covers the cloud it
    # describes -- measured on the PBMC fixture, which has exactly that many.
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False,
              markerscale=3, fontsize=7,
              ncol=1 if len(order) <= 22 else 2)
    stem = "11_tsne_eRegulon_gene_based" if kind == "gene" \
        else "12_tsne_eRegulon_region_based"
    save(fig, out_dir, stem)


def attach_embedding(adata, emb_dir, reduction):
    """Put a named layout on the eRegulon object and return what it is.

    The eRegulon object is concatenated from the AUC modalities alone, so it
    inherits no embedding from the Seurat object however carefully step 02
    chose one. Read the chosen layout from step 01's export instead of hoping
    it propagated through nine intermediate files -- which also means this
    step can be re-run on its own (`--only 20`) to redraw the figures of a
    finished run without recomputing any of the GRN.

    Returns a short label naming the layout, which goes on every figure.
    """
    reduction = (reduction or "").strip()
    if not reduction:
        if "X_umap" in adata.obsm:
            return "X_umap from the MuData"
        # No layout anywhere: a UMAP of regulon activity is a real answer here,
        # it is just a different one from the Seurat layout. Name it so nobody
        # reads these panels as the cell-type UMAP they know.
        sc.pp.neighbors(adata, use_rep="X")
        sc.tl.umap(adata)
        return "UMAP of eRegulon AUC (computed here)"

    path = Path(emb_dir) / f"embedding_{reduction}.tsv"
    if not path.exists():
        raise SystemExit(
            f"[viz] input.reduction {reduction!r}: {path} does not exist.\n"
            f"  Step 01 writes one embedding_<name>.tsv per reduction. Check "
            f"the name against that directory,\n  or clear input.reduction to "
            f"plot a UMAP of eRegulon activity instead.")

    emb = pd.read_csv(path, sep="\t").set_index("barcode")
    hit = adata.obs_names.intersection(emb.index)
    if len(hit) != adata.n_obs:
        # Barcodes disagreeing between the two sides is the failure this
        # pipeline has already hit once, at step 6. Say which names are on
        # each side rather than plotting an all-NaN layout.
        raise SystemExit(
            f"[viz] {len(hit)}/{adata.n_obs} cells found in {path.name}.\n"
            f"  SCENIC+ cell:  {adata.obs_names[0]!r}\n"
            f"  Seurat cell:   {emb.index[0]!r}\n"
            f"  The two are named differently, so no layout can be attached.")
    adata.obsm["X_umap"] = emb.loc[adata.obs_names].to_numpy()[:, :2]
    return reduction


def umap_celltype(adata, ct_col, out_dir, label):
    fig, ax = plt.subplots(figsize=(7, 6))
    sc.pl.umap(adata, color=ct_col, ax=ax, show=False, frameon=False, legend_loc="on data")
    ax.set_title(f"{ct_col}  [{label}]")
    save(fig, out_dir, "01_umap_celltype")


def umap_eRegulons(adata, eRegulons, out_dir, label):
    for er in eRegulons:
        if er not in adata.var_names:
            continue
        fig, ax = plt.subplots(figsize=(6, 5))
        sc.pl.umap(adata, color=er, ax=ax, show=False, frameon=False, cmap="viridis")
        ax.set_title(f"{er}  [{label}]")
        save(fig, out_dir, f"02_umap_eRegulon_{er.replace('/', '_').replace(' ', '_')}")


def rss_plot(rss, out_dir, top_n):
    from scenicplus.RSS import plot_rss
    n_groups = rss.shape[0]
    cols = min(5, max(1, n_groups))
    # figsize is PER SUBPLOT, not overall. plot_rss multiplies it itself --
    #   figsize = (figsize[0] * num_columns, figsize[1] * num_rows)   RSS.py:37
    # -- while its docstring calls the argument "the overall size of the
    # figure". Passing an already-multiplied size squares it: with 25 cell
    # types this asked for (4*5, 3*5) and got (100, 75) inches, i.e. a
    # 20000 x 15000 px PNG that PIL refuses to open without raising
    # DecompressionBombError. Measured both ways on a 25-group matrix.
    fig = plot_rss(data_matrix=rss, top_n=top_n, num_columns=cols,
                   figsize=(4, 3))
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
    p.add_argument("--embedding_dir", default="",
                   help="step 01's seurat_export/, holding embedding_<name>.tsv")
    p.add_argument("--reduction", default="",
                   help="input.reduction: which of those to plot on. Empty = "
                        "plot on a UMAP of eRegulon activity.")
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
    emb_label = None
    if er_adata is not None:
        # Before any plotting: one layout, named, used by every panel below.
        # Attaching it here rather than inside the first plot also fixes a
        # latent crash -- umap_eRegulons runs even when no cell-type column is
        # found, and used to reach sc.pl.umap with nothing in obsm.
        emb_label = attach_embedding(er_adata, args.embedding_dir, args.reduction)
        print(f"[viz] layout: {emb_label}", flush=True)
    if er_adata is not None and rss_var is not None:
        er_adata.obs[rss_var] = er_adata.obs[rss_var].astype("category")
        umap_celltype(er_adata, rss_var, out_dir, emb_label)

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
        umap_eRegulons(er_adata, top_eRegulons, out_dir, emb_label)

    # heatmap-dotplots
    if rss_var is not None:
        heatmap_dotplot(md, rss_var, out_dir, source="direct")
        heatmap_dotplot(md, rss_var, out_dir, source="extended")

    # eRegulon metadata-derived plots
    direct = md.uns.get("direct_e_regulon_metadata", pd.DataFrame())
    extended = md.uns.get("extended_e_regulon_metadata", pd.DataFrame())
    combined = pd.concat([direct, extended], ignore_index=True) \
        if (len(direct) + len(extended)) else pd.DataFrame()
    seed = int(cfg.get("resources", {}).get("seed", 555))
    eregulon_tsne(md, celltype_col, out_dir, "gene", seed)
    eregulon_tsne(md, celltype_col, out_dir, "region", seed)
    region_overlap(direct, out_dir, "direct", viz.get("overlap_top_n", 40))
    region_overlap(extended, out_dir, "extended", viz.get("overlap_top_n", 40))
    tf_target_count(combined, out_dir, top_n=40)
    tf_importance_density(combined, out_dir, top_n=viz["network_top_n_TFs"])
    network_plot(combined, out_dir,
                 top_tfs=viz["network_top_n_TFs"],
                 top_targets=viz["network_top_n_targets_per_TF"])

    print(f"[viz] Plots written to: {out_dir}")


if __name__ == "__main__":
    main()
