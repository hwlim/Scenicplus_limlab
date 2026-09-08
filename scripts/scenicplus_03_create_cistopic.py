#!/usr/bin/env python
"""
Build a CistopicObject from the peak-count matrix exported by the R step.
The Seurat object usually does not carry fragment files, so we go through
the matrix-only path. Cell metadata (incl. cell type) is attached so DARs
can be computed per cell type later on.
"""
import argparse
import pickle
from pathlib import Path

import pandas as pd
import scipy.io
import scipy.sparse as sp


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--in_dir", required=True, help="Directory written by the R step")
    p.add_argument("--out_pkl", required=True)
    p.add_argument("--project", default="scenicplus_run")
    args = p.parse_args()

    in_dir = Path(args.in_dir)

    # peaks x cells
    counts = scipy.io.mmread(str(in_dir / "atac_counts.mtx")).tocsr()
    barcodes = pd.read_csv(in_dir / "atac_barcodes.tsv", header=None)[0].tolist()
    regions = pd.read_csv(in_dir / "atac_regions.tsv", header=None)[0].tolist()
    meta = pd.read_csv(in_dir / "cell_metadata.tsv", sep="\t").set_index("barcode")
    meta = meta.loc[barcodes]

    counts_df = pd.DataFrame.sparse.from_spmatrix(
        counts, index=regions, columns=barcodes
    )

    from pycisTopic.cistopic_class import create_cistopic_object

    # tag_cells=False is NOT the default and it matters.
    #
    # create_cistopic_object defaults tag_cells=True, which appends
    # "___<project>" to every cell name -- measured:
    #
    #   default          ['AAACAGCCAAGGAATC-1___scenicplus_run', ...]
    #   tag_cells=False  ['AAACAGCCAAGGAATC-1', ...]
    #
    # That tag exists to keep barcodes distinct when SEVERAL samples are merged
    # into one cisTopic object. Here everything comes from a single Seurat
    # object whose barcodes are already unique, so it disambiguates nothing --
    # and it silently breaks step 6, which intersects these names with the RNA
    # AnnData's obs_names (untagged, from the same object) and dies with
    # "No cells found which are present in both assays". The two sides came
    # from one file and still did not match.
    cto = create_cistopic_object(
        fragment_matrix=counts_df,
        cell_names=barcodes,
        region_names=regions,
        project=args.project,
        tag_cells=False,
    )
    cto.add_cell_data(meta)

    out = Path(args.out_pkl)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(out.suffix + ".partial")
    with open(tmp, "wb") as fh:
        pickle.dump(cto, fh)
    tmp.replace(out)
    print(f"[create_cistopic] Wrote: {out}  "
          f"(n_regions={len(cto.region_names)}, n_cells={len(cto.cell_names)})")


if __name__ == "__main__":
    main()
