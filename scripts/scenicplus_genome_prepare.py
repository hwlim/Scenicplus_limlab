#!/usr/bin/env python
"""Validate a supplied genome annotation + chromsizes pair, then place them.

WHY THIS EXISTS SEPARATELY FROM STEP 7. `scenicplus prepare_data
download_genome_annotations` cannot produce chromsizes for anyone -- it queries
NCBI's retired `db=genome` -- and for mouse it silently reports whatever
assembly Ensembl serves today, which is GRCm39, against mm10 data. So supplying
both files is the normal path, and the only question left is whether the pair
you supplied is the RIGHT pair. That is what this checks, before the expensive
stages consume them rather than after.

Four checks, in the order that a failure is cheapest to understand:

  1. shape      -- chromsizes needs a tab-separated header row; the annotation
                   needs the seven columns get_search_space reads.
  2. assembly   -- chromosome 1's length IS the assembly. 195,471,971 is mm10
                   and 195,154,279 is GRCm39, so one number catches the failure
                   that otherwise reaches the results unannounced.
  3. naming     -- UCSC `chr1` vs Ensembl `1`, across all THREE files that
                   get_search_space joins, including the ATAC peaks. A
                   disagreement there surfaces today as a pandas KeyError two
                   stages later.
  4. overlap    -- the annotation's chromosomes against the peaks' actual
                   chromosomes. Agreeing on style is not the same as agreeing.

Then it copies, never before: a partial copy that passes for complete is how
the driver's version of this went wrong twice.
"""
import argparse
import hashlib
import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

# Reuse the driver's own checks rather than restating them. One source for what
# "correctly shaped" means, so the two drivers cannot drift apart on it.
sys.path.insert(0, str(Path(__file__).resolve().parent))
import scenicplus_06_grn_stage as gs  # noqa: E402

ANNOTATION_COLUMNS = ["Chromosome", "Start", "End", "Strand", "Gene",
                      "Transcription_Start_Site", "Transcript_type"]

# chr1's length, per assembly. hg38 measured from the reference pair that
# produced the validated PBMC run; mm10 from BSgenome.Mmusculus.UCSC.mm10.
EXPECTED_CHR1 = {"hg38": 248956422, "mm10": 195471971}

# Never used to gate -- only to make a mismatch legible. An observed length
# landing on one of these names the mistake instead of just reporting a number.
KNOWN_CHR1 = {248956422: "hg38/GRCh38", 249250621: "hg19/GRCh37",
              195471971: "mm10/GRCm38", 195154279: "GRCm39"}


def die(*lines):
    for ln in lines:
        print(f"[genome] {ln}", file=sys.stderr)
    sys.exit(1)


def sha256(path, limit=None):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_header(path):
    with open(path) as fh:
        return fh.readline().rstrip("\n").split("\t")


def chrom_lengths(chromsizes):
    out = {}
    with open(chromsizes) as fh:
        next(fh, None)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[2].strip().isdigit():
                out[f[0].strip()] = int(f[2])
    return out


def region_chroms(atac_regions):
    """Chromosome names from step 01's atac_regions.tsv ("chr1:100-200").

    No header: the file is one region per line, written by the R step.
    """
    chroms = set()
    with open(atac_regions) as fh:
        for line in fh:
            s = line.strip()
            if s:
                chroms.add(s.split(":")[0])
    return chroms


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--annotation", required=True)
    p.add_argument("--chromsizes", required=True)
    p.add_argument("--atac_regions", required=True)
    p.add_argument("--out_annotation", required=True)
    p.add_argument("--out_chromsizes", required=True)
    p.add_argument("--record", required=True)
    p.add_argument("--species", required=True)
    p.add_argument("--assembly", default="")
    args = p.parse_args()

    # --- 1. shape ------------------------------------------------------------
    gs._check_chromsizes_shape(args.chromsizes)      # exits with its own message

    head = read_header(args.annotation)
    missing = [c for c in ANNOTATION_COLUMNS if c not in head]
    if missing:
        die(f"{args.annotation} is missing column(s) {missing}.",
            f"  header:   {head}",
            f"  required: {ANNOTATION_COLUMNS}",
            "  get_search_space reads these by name, and an annotation left by a",
            "  failed step 7 is missing the conversion as well as the columns.")

    lengths = chrom_lengths(args.chromsizes)
    if not lengths:
        die(f"{args.chromsizes} has a header but no usable rows.")

    # --- 2. assembly ---------------------------------------------------------
    # chr1 under either naming, because the naming check has not run yet.
    chr1 = lengths.get("chr1", lengths.get("1"))
    assembly = args.assembly.strip()
    expected = EXPECTED_CHR1.get(assembly)
    if chr1 is None:
        print("[genome] WARNING: no chromosome 1 in the chromsizes, so the "
              "assembly cannot be confirmed.", file=sys.stderr)
    elif expected and chr1 != expected:
        looks = KNOWN_CHR1.get(chr1, "an assembly not in this table")
        die(f"chromosome 1 is {chr1:,} bp, but {assembly} is {expected:,} bp.",
            f"  {chr1:,} is {looks}.",
            "  Chromosome 1's length IS the assembly. This is the failure that",
            "  does not announce itself: names still convert, the search space",
            "  still builds, and the peak-gene links come out quietly wrong.",
            "  Rebuild the pair with scripts/scenicplus_make_genome_files.R, or",
            "  correct input.assembly if the data really is the other one.")
    elif expected:
        print(f"[genome] chromosome 1 = {chr1:,} bp, consistent with {assembly}")
    else:
        looks = KNOWN_CHR1.get(chr1)
        print(f"[genome] chromosome 1 = {chr1:,} bp"
              + (f" ({looks})" if looks else "")
              + f"; input.assembly is {assembly or 'unset'}, so nothing to check it "
                "against.", file=sys.stderr)

    # --- 3. naming, across all three files -----------------------------------
    ann_chroms = gs._tsv_chroms(args.annotation)
    peaks = region_chroms(args.atac_regions)
    styles = {
        "annotation": gs._chrom_style(ann_chroms),
        "chromsizes": gs._chrom_style(list(lengths)),
        "ATAC peaks": gs._chrom_style(sorted(peaks)),
    }
    if len(set(styles.values())) != 1:
        die("the three files get_search_space joins do not agree on chromosome "
            "naming:",
            *[f"  {k:<12} {v:<10} e.g. {next(iter(s), '?')}"
              for (k, v), s in zip(styles.items(),
                                   (ann_chroms, list(lengths), sorted(peaks)))],
            "  Nothing can overlap across a naming mismatch. An Ensembl-style",
            "  annotation is what a failed step 7 leaves behind -- its UCSC",
            "  conversion lives in the same branch that fetches chromsizes.")
    print(f"[genome] chromosome naming: {next(iter(styles.values()))} in all three")

    # --- 4. overlap with the data itself -------------------------------------
    ann_all = set(gs._tsv_chroms(args.annotation, limit=10**9))
    shared = peaks & ann_all
    if not shared:
        die("no chromosome appears in BOTH the annotation and the ATAC peaks.",
            f"  peaks:      {sorted(peaks)[:6]} ({len(peaks)} total)",
            f"  annotation: {sorted(ann_all)[:6]} ({len(ann_all)} total)",
            "  They agree on naming style but name different things, so the",
            "  search space would be empty.")
    only_peaks = sorted(peaks - ann_all)
    print(f"[genome] {len(shared)}/{len(peaks)} peak chromosomes are in the "
          f"annotation")
    if only_peaks:
        print(f"[genome] not in the annotation, so their peaks get no genes: "
              f"{only_peaks[:10]}{' ...' if len(only_peaks) > 10 else ''}")

    # --- place, and only now -------------------------------------------------
    for src, dst in ((args.annotation, args.out_annotation),
                     (args.chromsizes, args.out_chromsizes)):
        Path(dst).parent.mkdir(parents=True, exist_ok=True)
        tmp = Path(str(dst) + ".partial")
        shutil.copyfile(src, tmp)
        tmp.replace(dst)
        print(f"[genome] {src} -> {dst}")

    # --- the record ----------------------------------------------------------
    # Nothing today records which assembly a run used, which is exactly how the
    # GRCm39 problem stays invisible. A provenance bundle can answer it from
    # this file instead of from memory.
    record = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "species": args.species,
        "assembly": assembly or None,
        "chr1_bp": chr1,
        "chr1_matches": KNOWN_CHR1.get(chr1),
        "chromosome_naming": next(iter(styles.values())),
        "n_chromosomes": len(lengths),
        "n_annotation_rows": sum(1 for _ in open(args.annotation)) - 1,
        "n_annotation_chromosomes": len(ann_all),
        "peak_chromosomes": len(peaks),
        "peak_chromosomes_annotated": len(shared),
        "peak_chromosomes_unannotated": only_peaks,
        "source": {
            "annotation": {"path": str(Path(args.annotation).resolve()),
                           "sha256": sha256(args.annotation)},
            "chromsizes": {"path": str(Path(args.chromsizes).resolve()),
                           "sha256": sha256(args.chromsizes)},
        },
    }
    Path(args.record).parent.mkdir(parents=True, exist_ok=True)
    Path(args.record).write_text(json.dumps(record, indent=2) + "\n")
    print(f"[genome] recorded: {args.record}")


if __name__ == "__main__":
    main()
