#!/usr/bin/env python
"""
Generate the region-set folder consumed by SCENIC+ motif enrichment:

  region_sets/
    DARs_<celltype_col>/        one bed per cell type (DARs, 1-vs-rest)
    Topics_otsu/                one bed per binarized topic (Otsu cutoff)
    Topics_top_3k/              one bed per topic, top 3k regions
"""
import argparse
import os
import pickle
from pathlib import Path

import pandas as pd
import yaml


def write_bed(regions, path: Path):
    rows = []
    for r in regions:
        try:
            chrom, rest = r.split(":")
            start, end = rest.split("-")
            rows.append((chrom, int(start), int(end), r))
        except ValueError:
            continue
    if not rows:
        return
    df = pd.DataFrame(rows, columns=["chrom", "start", "end", "name"])
    df.to_csv(path, sep="\t", header=False, index=False)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--in_pkl", required=True, help="CistopicObject with LDA model")
    p.add_argument("--out_dir", required=True)
    p.add_argument("--config", required=True)
    p.add_argument("--n_cpu", type=int, default=8)
    args = p.parse_args()

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)
    celltype_col = cfg["input"]["celltype_column"]
    ct_cfg = cfg["cistopic"]

    with open(args.in_pkl, "rb") as fh:
        cto = pickle.load(fh)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    from pycisTopic.diff_features import (
        find_diff_features,
        impute_accessibility,
        normalize_scores,
    )
    from pycisTopic.topic_binarization import binarize_topics

    # ----- Topic-based region sets -----
    region_bin_otsu = binarize_topics(cto, method="otsu")
    region_bin_top  = binarize_topics(cto, method="ntop", ntop=3000)

    otsu_dir = out_dir / "Topics_otsu"
    top_dir  = out_dir / "Topics_top_3k"
    otsu_dir.mkdir(exist_ok=True)
    top_dir.mkdir(exist_ok=True)
    for topic, df in region_bin_otsu.items():
        write_bed(df.index.tolist(), otsu_dir / f"{topic}.bed")
    for topic, df in region_bin_top.items():
        write_bed(df.index.tolist(), top_dir / f"{topic}.bed")

    # ----- DAR region sets -----
    imputed = impute_accessibility(
        cto, selected_cells=None, selected_regions=None,
        scale_factor=10**6, chunk_size=20000,
    )
    norm_imp = normalize_scores(imputed, scale_factor=10**4)

    markers_dict = find_diff_features(
        cto,
        norm_imp,
        variable=celltype_col,
        var_features=None,
        contrasts=None,
        adjpval_thr=ct_cfg["dar_adjpval_thr"],
        log2fc_thr=ct_cfg["dar_log2fc_thr"],
        n_cpu=args.n_cpu,
        split_pattern="-",
    )

    dar_dir = out_dir / f"DARs_{celltype_col}"
    dar_dir.mkdir(exist_ok=True)
    for ct_name, df in markers_dict.items():
        if df is None or df.shape[0] == 0:
            continue
        safe = "".join(c if c.isalnum() or c in "._-" else "_" for c in str(ct_name))
        write_bed(df.index.tolist(), dar_dir / f"{safe}.bed")

    print(f"[region_sets] Wrote region sets under: {out_dir}")
    for sub in sorted(os.listdir(out_dir)):
        sub_dir = out_dir / sub
        if sub_dir.is_dir():
            n = len(list(sub_dir.glob("*.bed")))
            print(f"  {sub}: {n} bed files")


if __name__ == "__main__":
    main()
