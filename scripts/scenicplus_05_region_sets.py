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


def write_bed(regions, path: Path, min_regions: int = 0) -> bool:
    """Write one region set. -> True if written, False if skipped.

    MINIMUM SIZE, and the guard lives HERE rather than at one call site.

    A degenerate region set aborts the whole motif-enrichment stage. Measured
    on a real run: Otsu binarization produced a 29-region topic out of 40, none
    of whose regions overlapped the cisTarget database, and pycistarget's
    ValueError propagated through joblib's loky backend and tore down the pool
    -- so both R09_cistarget and R10_dem died and the 86 healthy sets were
    discarded, after many hours of upstream stages. See issue #1.

    IT IS A PROXY, NOT A CURE. What actually fails is "zero DATABASE regions
    mapped", and a region count only correlates with that: on the same run a
    49-region topic mapped 34 database regions and survived, so a larger set
    that happens to miss the cCREs would still crash. This lowers the
    probability; the fix that removes the failure mode is per-set error
    handling inside SCENIC+, which is filed upstream.

    ALL THREE WRITERS get it, which is why it is not in the Otsu loop. The
    original report saw only Otsu produce a tiny set, but DARs are 1-vs-rest
    per cell type and `find_diff_features` skips only a completely empty
    result -- a rare cell type is the same crash through a different door.

    EVERY SKIP IS LOGGED, including the two that were already silent. A thin
    region set and a healthy one must not look alike in the output, and
    `write_bed` previously returned quietly both for an all-unparseable region
    list and for an empty one.
    """
    rows = []
    malformed = 0
    for r in regions:
        try:
            chrom, rest = r.split(":")
            start, end = rest.split("-")
            rows.append((chrom, int(start), int(end), r))
        except ValueError:
            malformed += 1
    if malformed:
        print(f"  [region_sets] {path.name}: {malformed} region name(s) are not "
              f"chrom:start-end and were dropped")
    if not rows:
        print(f"  [region_sets] SKIP {path.name}: no usable regions "
              f"({len(regions)} in, {malformed} malformed)")
        return False
    if len(rows) < min_regions:
        print(f"  [region_sets] SKIP {path.name}: {len(rows)} regions "
              f"(< min_regions_per_set = {min_regions})")
        return False
    df = pd.DataFrame(rows, columns=["chrom", "start", "end", "name"])
    df.to_csv(path, sep="\t", header=False, index=False)
    return True


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
    # Config-driven, because the right value depends on the cisTarget database
    # and on fraction_overlap -- a SCREEN cCRE database is a curated subset of
    # the genome rather than a tiling of it, so a small set maps far worse
    # against it than against a tiling one.
    min_regions = int(ct_cfg.get("min_regions_per_set", 500))

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
    n_skipped = 0
    for topic, df in region_bin_otsu.items():
        n_skipped += not write_bed(df.index.tolist(),
                                   otsu_dir / f"{topic}.bed", min_regions)
    for topic, df in region_bin_top.items():
        n_skipped += not write_bed(df.index.tolist(),
                                   top_dir / f"{topic}.bed", min_regions)

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
        n_skipped += not write_bed(df.index.tolist(),
                                   dar_dir / f"{safe}.bed", min_regions)

    # The count, stated once, so a run that dropped sets says so in its last
    # line rather than only in the middle of the log.
    print(f"[region_sets] Wrote region sets under: {out_dir}"
          + (f"  ({n_skipped} set(s) skipped, min_regions_per_set={min_regions})"
             if n_skipped else ""))
    for sub in sorted(os.listdir(out_dir)):
        sub_dir = out_dir / sub
        if sub_dir.is_dir():
            n = len(list(sub_dir.glob("*.bed")))
            print(f"  {sub}: {n} bed files")


if __name__ == "__main__":
    main()
