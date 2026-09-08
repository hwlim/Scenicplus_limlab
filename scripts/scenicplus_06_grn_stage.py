#!/usr/bin/env python
"""
Run a single SCENIC+ GRN-inference stage as a direct `scenicplus` CLI call.

This replaces the inner snakemake (and the old init_inner.sh + patch_config.py
scaffolding). Each stage below is a verbatim transcription of the corresponding
rule's `shell:` block in the SCENIC+-generated Snakefile, with arguments sourced
from *our* top-level config.yaml instead of a patched snakemake config. The
master driver (scenicplus_run_pipeline.sh) invokes this once per stage so every
stage gets its own sentinel + .cfgsha + cascade treatment.

Intermediate files all live under --scplus_out (paths mirror the old
patch_config.py `output_data` map, so nothing downstream needs to change).

Usage:
    scenicplus_06_grn_stage.py --stage <name> --config <config.yaml> \
        --scplus_out <dir> [--cistopic_obj PKL] [--adata H5AD] \
        [--region_sets DIR]

Stages (in dependency order):
    prepare_gex_acc  genome_annot  search_space  cistarget  dem
    prepare_menr  tf_to_gene  region_to_gene  egrn_direct  egrn_extended
    aucell_direct  aucell_extended  scplus_mudata
"""
import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import yaml


def toks(val) -> list[str]:
    """Flatten a config value into CLI tokens.

    Lists become one token per element; whitespace-delimited strings (e.g.
    search-space "1000 150000", which the CLI reads as nargs=2) are split;
    scalars become a single token. Mirrors how snakemake expanded these
    params into the shell command.
    """
    if isinstance(val, (list, tuple)):
        return [str(x) for x in val]
    if isinstance(val, str) and (" " in val.strip()):
        return val.split()
    return [str(val)]


def out_paths(scplus_out: Path) -> dict:
    """Intermediate/output file paths — mirrors patch_config.py output_data."""
    o = scplus_out
    return {
        "combined_GEX_ACC_mudata": o / "ACC_GEX.h5mu",
        "dem_result_fname":        o / "dem_results.hdf5",
        "ctx_result_fname":        o / "ctx_results.hdf5",
        "output_fname_dem_html":   o / "dem_results.html",
        "output_fname_ctx_html":   o / "ctx_results.html",
        "cistromes_direct":        o / "cistromes_direct.h5ad",
        "cistromes_extended":      o / "cistromes_extended.h5ad",
        "tf_names":                o / "tf_names.txt",
        "genome_annotation":       o / "genome_annotation.tsv",
        "chromsizes":              o / "chromsizes.tsv",
        "search_space":            o / "search_space.tsv",
        "tf_to_gene_adjacencies":  o / "tf_to_gene_adj.tsv",
        "region_to_gene_adjacencies": o / "region_to_gene_adj.tsv",
        "eRegulons_direct":        o / "eRegulons_direct.tsv",
        "eRegulons_extended":      o / "eRegulons_extended.tsv",
        "AUCell_direct":           o / "AUCell_direct.h5mu",
        "AUCell_extended":         o / "AUCell_extended.h5mu",
        "scplus_mdata":            o / "scplusmdata.h5mu",
    }


def build_argv(stage: str, cfg: dict, out: dict, args) -> list[str]:
    """Return the full `scenicplus ...` argv for one stage."""
    sp = cfg["scenicplus"]
    g = cfg["grn"]
    inp = cfg["input"]
    n_cpu = str(cfg["resources"]["n_cpu"])
    seed = str(cfg["resources"]["seed"])
    temp_dir = str(Path(cfg["output"]["tmp"]).resolve())
    species = inp["species"]

    ctx_db = str(Path(inp["ctx_db"]).resolve())
    dem_db = str(Path(inp["dem_db"]).resolve())
    motif_ann = str(Path(inp["motif_annotations"]).resolve())

    def req(p, what):
        if not p:
            sys.exit(f"[grn_stage] stage '{stage}' requires --{what}")
        return str(Path(p).resolve())

    if stage == "prepare_gex_acc":
        av = [
            "scenicplus", "prepare_data", "prepare_GEX_ACC",
            "--cisTopic_obj_fname", req(args.cistopic_obj, "cistopic_obj"),
            "--GEX_anndata_fname", req(args.adata, "adata"),
            "--out_file", str(out["combined_GEX_ACC_mudata"]),
            "--bc_transform_func", str(sp["bc_transform_func"]),
        ]
        if not sp["is_multiome"]:
            av += [
                "--is_not_multiome",
                # patch_config.py hard-coded these for the non-multiome branch.
                "--key_to_group_by", "",
                "--nr_cells_per_metacells", "10",
            ]
        return av

    if stage == "genome_annot":
        return [
            "scenicplus", "prepare_data", "download_genome_annotations",
            "--species", species,
            "--biomart_host", inp["biomart_host"],
            "--genome_annotation_out_fname", str(out["genome_annotation"]),
            "--chromsizes_out_fname", str(out["chromsizes"]),
        ]

    if stage == "search_space":
        # NOTE: the generated Snakefile spells this subcommand "search_spance".
        # Transcribed verbatim; verify against the installed CLI and fix here
        # (single source) if it is actually "search_space".
        return [
            "scenicplus", "prepare_data", "search_spance",
            "--multiome_mudata_fname", str(out["combined_GEX_ACC_mudata"]),
            "--gene_annotation_fname", str(out["genome_annotation"]),
            "--chromsizes_fname", str(out["chromsizes"]),
            "--out_fname", str(out["search_space"]),
            "--upstream", *toks(sp["search_space_upstream"]),
            "--downstream", *toks(sp["search_space_downstream"]),
            "--extend_tss", *toks(sp["search_space_extend_tss"]),
        ]

    if stage == "cistarget":
        return [
            "scenicplus", "grn_inference", "motif_enrichment_cistarget",
            "--region_set_folder", req(args.region_sets, "region_sets"),
            "--cistarget_db_fname", ctx_db,
            "--output_fname_cistarget_result", str(out["ctx_result_fname"]),
            "--temp_dir", temp_dir,
            "--species", species,
            "--fr_overlap_w_ctx_db", str(sp["fraction_overlap_w_ctx_database"]),
            "--auc_threshold", str(sp["ctx_auc_threshold"]),
            "--nes_threshold", str(sp["ctx_nes_threshold"]),
            "--rank_threshold", str(sp["ctx_rank_threshold"]),
            "--path_to_motif_annotations", motif_ann,
            "--annotation_version", str(sp["annotation_version"]),
            "--motif_similarity_fdr", str(sp["motif_similarity_fdr"]),
            "--orthologous_identity_threshold", str(sp["orthologous_identity_threshold"]),
            "--annotations_to_use", *toks(sp["annotations_to_use"]),
            "--write_html",
            "--output_fname_cistarget_html", str(out["output_fname_ctx_html"]),
            "--n_cpu", n_cpu,
        ]

    if stage == "dem":
        av = [
            "scenicplus", "grn_inference", "motif_enrichment_dem",
            "--region_set_folder", req(args.region_sets, "region_sets"),
            "--dem_db_fname", dem_db,
            "--output_fname_dem_result", str(out["dem_result_fname"]),
            "--temp_dir", temp_dir,
            "--species", species,
            "--fraction_overlap_w_dem_database", str(sp["fraction_overlap_w_dem_database"]),
            "--max_bg_regions", str(sp["dem_max_bg_regions"]),
        ]
        if sp["dem_balance_number_of_promoters"]:
            # balanced branch also consumes the genome annotation.
            av += [
                "--genome_annotation", str(out["genome_annotation"]),
                "--balance_number_of_promoters",
                "--promoter_space", str(sp["dem_promoter_space"]),
            ]
        av += [
            "--adjpval_thr", str(sp["dem_adj_pval_thr"]),
            "--log2fc_thr", str(sp["dem_log2fc_thr"]),
            "--mean_fg_thr", str(sp["dem_mean_fg_thr"]),
            "--motif_hit_thr", str(sp["dem_motif_hit_thr"]),
            "--path_to_motif_annotations", motif_ann,
            "--annotation_version", str(sp["annotation_version"]),
            "--motif_similarity_fdr", str(sp["motif_similarity_fdr"]),
            "--orthologous_identity_threshold", str(sp["orthologous_identity_threshold"]),
            "--annotations_to_use", *toks(sp["annotations_to_use"]),
            "--write_html",
            "--output_fname_dem_html", str(out["output_fname_dem_html"]),
            "--seed", seed,
            "--n_cpu", n_cpu,
        ]
        return av

    if stage == "prepare_menr":
        return [
            "scenicplus", "prepare_data", "prepare_menr",
            "--paths_to_motif_enrichment_results",
            str(out["dem_result_fname"]), str(out["ctx_result_fname"]),
            "--multiome_mudata_fname", str(out["combined_GEX_ACC_mudata"]),
            "--out_file_tf_names", str(out["tf_names"]),
            "--out_file_direct_annotation", str(out["cistromes_direct"]),
            "--out_file_extended_annotation", str(out["cistromes_extended"]),
            "--direct_annotation", *toks(sp["direct_annotation"]),
            "--extended_annotation", *toks(sp["extended_annotation"]),
        ]

    if stage == "tf_to_gene":
        return [
            "scenicplus", "grn_inference", "TF_to_gene",
            "--multiome_mudata_fname", str(out["combined_GEX_ACC_mudata"]),
            "--tf_names", str(out["tf_names"]),
            "--temp_dir", temp_dir,
            "--out_tf_to_gene_adjacencies", str(out["tf_to_gene_adjacencies"]),
            "--method", str(g["tf_to_gene_importance_method"]),
            "--n_cpu", n_cpu,
            "--seed", seed,
        ]

    if stage == "region_to_gene":
        return [
            "scenicplus", "grn_inference", "region_to_gene",
            "--multiome_mudata_fname", str(out["combined_GEX_ACC_mudata"]),
            "--search_space_fname", str(out["search_space"]),
            "--temp_dir", temp_dir,
            "--out_region_to_gene_adjacencies", str(out["region_to_gene_adjacencies"]),
            "--importance_scoring_method", str(g["region_to_gene_importance_method"]),
            "--correlation_scoring_method", str(g["region_to_gene_correlation_method"]),
            "--n_cpu", n_cpu,
        ]

    if stage in ("egrn_direct", "egrn_extended"):
        extended = stage == "egrn_extended"
        cistromes = out["cistromes_extended"] if extended else out["cistromes_direct"]
        ereg_out = out["eRegulons_extended"] if extended else out["eRegulons_direct"]
        av = ["scenicplus", "grn_inference", "eGRN"]
        if extended:
            av += ["--is_extended"]
        av += [
            "--TF_to_gene_adj_fname", str(out["tf_to_gene_adjacencies"]),
            "--region_to_gene_adj_fname", str(out["region_to_gene_adjacencies"]),
            "--cistromes_fname", str(cistromes),
            "--ranking_db_fname", ctx_db,
            "--eRegulon_out_fname", str(ereg_out),
            "--temp_dir", temp_dir,
            "--order_regions_to_genes_by", str(g["order_regions_to_genes_by"]),
            "--order_TFs_to_genes_by", str(g["order_TFs_to_genes_by"]),
            "--gsea_n_perm", str(g["gsea_n_perm"]),
            "--quantiles", *toks(g["quantile_thresholds_region_to_gene"]),
            "--top_n_regionTogenes_per_gene", *toks(g["top_n_regionTogenes_per_gene"]),
            "--top_n_regionTogenes_per_region", *toks(g["top_n_regionTogenes_per_region"]),
            "--min_regions_per_gene", str(g["min_regions_per_gene"]),
            "--rho_threshold", str(g["rho_threshold"]),
            "--min_target_genes", str(g["min_target_genes"]),
            "--n_cpu", n_cpu,
        ]
        return av

    if stage in ("aucell_direct", "aucell_extended"):
        extended = stage == "aucell_extended"
        ereg = out["eRegulons_extended"] if extended else out["eRegulons_direct"]
        auc_out = out["AUCell_extended"] if extended else out["AUCell_direct"]
        return [
            "scenicplus", "grn_inference", "AUCell",
            "--eRegulon_fname", str(ereg),
            "--multiome_mudata_fname", str(out["combined_GEX_ACC_mudata"]),
            "--aucell_out_fname", str(auc_out),
            "--n_cpu", n_cpu,
        ]

    if stage == "scplus_mudata":
        return [
            "scenicplus", "grn_inference", "create_scplus_mudata",
            "--multiome_mudata_fname", str(out["combined_GEX_ACC_mudata"]),
            "--e_regulon_auc_direct_mudata_fname", str(out["AUCell_direct"]),
            "--e_regulon_auc_extended_mudata_fname", str(out["AUCell_extended"]),
            "--e_regulon_metadata_direct_fname", str(out["eRegulons_direct"]),
            "--e_regulon_metadata_extended_fname", str(out["eRegulons_extended"]),
            "--out_file", str(out["scplus_mdata"]),
        ]

    sys.exit(f"[grn_stage] unknown stage: {stage}")


STAGES = [
    "prepare_gex_acc", "genome_annot", "search_space", "cistarget", "dem",
    "prepare_menr", "tf_to_gene", "region_to_gene", "egrn_direct",
    "egrn_extended", "aucell_direct", "aucell_extended", "scplus_mudata",
]


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--stage", required=True, choices=STAGES)
    p.add_argument("--config", required=True)
    p.add_argument("--scplus_out", required=True)
    p.add_argument("--cistopic_obj")
    p.add_argument("--adata")
    p.add_argument("--region_sets")
    args = p.parse_args()

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)

    scplus_out = Path(args.scplus_out).resolve()
    scplus_out.mkdir(parents=True, exist_ok=True)
    out = out_paths(scplus_out)

    # genome_annot is the one stage that can be satisfied from disk, and the one
    # that lies about having succeeded. Handled before dispatch.
    if args.stage == "genome_annot" and _use_supplied_annotations(cfg, out):
        return 0

    argv = build_argv(args.stage, cfg, out, args)
    print("[grn_stage] " + " ".join(argv), flush=True)
    try:
        subprocess.run(argv, check=True)
        if args.stage == "genome_annot":
            _assert_genome_annot_complete(cfg, out)
    except subprocess.CalledProcessError:
        # prepare_GEX_ACC fails with "No cells found which are present in both
        # assays, check input and consider using `bc_transform_func`" and does
        # not say what the two sets of names LOOK like, which is the only thing
        # a reader needs. Print that here rather than make them write this
        # snippet themselves at the worst moment.
        if args.stage == "prepare_gex_acc":
            _diagnose_barcodes(args, cfg)
        raise
    return 0


def _use_supplied_annotations(cfg: dict, out: dict) -> bool:
    """Copy user-provided annotation/chromsizes into place; True if handled.

    `download_genome_annotations` needs Ensembl BioMart AND NCBI E-utilities.
    Behind a proxy that gates either, there is no way through it -- hence the
    escape hatch. Both files come out of one call, so both must be supplied
    together; taking one from disk and the other from the network would still
    need the network.
    """
    inp = cfg.get("input", {})
    ann, chrom = inp.get("genome_annotation") or "", inp.get("chromsizes") or ""
    if not ann and not chrom:
        return False
    if bool(ann) != bool(chrom):
        have, miss = ("genome_annotation", "chromsizes") if ann else ("chromsizes", "genome_annotation")
        sys.exit(f"[grn_stage] input.{have} is set but input.{miss} is not. Both "
                 f"come from one download, so supplying one still requires the "
                 f"network for the other -- set both, or neither.")
    pairs = (("genome_annotation", ann, out["genome_annotation"]),
             ("chromsizes", chrom, out["chromsizes"]))
    # Validate BOTH before copying EITHER: copying the first and then failing on
    # the second leaves half the stage's outputs in place, which the next run
    # has no way to tell from a good one.
    for key, src, _ in pairs:
        if not Path(src).expanduser().is_file():
            sys.exit(f"[grn_stage] input.{key} = {src!r} does not exist")
    for key, src, dst in pairs:
        p = Path(src).expanduser()
        shutil.copyfile(p, dst)
        print(f"[grn_stage] input.{key}: {p} -> {dst}", flush=True)
    _check_chromsizes_shape(out["chromsizes"])
    return True


def _check_chromsizes_shape(path) -> None:
    """A hand-made chromsizes file is easy to get subtly wrong; say so here.

    search_space reads it with pd.read_table, so it is TAB-separated WITH a
    header row and needs the columns SCENIC+ writes: Chromosome, Start, End.
    A UCSC .chrom.sizes file has neither the header nor the Start column, and
    pandas will happily read it as a two-column frame whose header is the first
    chromosome -- which fails much later, inside get_search_space.
    """
    try:
        first = Path(path).read_text().splitlines()[0].split("\t")
    except (OSError, IndexError) as e:
        sys.exit(f"[grn_stage] cannot read chromsizes {path}: {e}")
    need = ["Chromosome", "Start", "End"]
    if first[:3] != need:
        sys.exit(
            f"[grn_stage] {path} has header {first!r}, expected {need!r}.\n"
            f"  It is read with pandas.read_table, so it needs a TAB-separated\n"
            f"  header row. A raw UCSC .chrom.sizes converts with:\n"
            f"    awk 'BEGIN{{OFS=\"\\t\"; print \"Chromosome\",\"Start\",\"End\"}} "
            f"{{print $1,0,$2}}' hg38.chrom.sizes > chromsizes.tsv")


def _assert_genome_annot_complete(cfg: dict, out: dict) -> None:
    """download_genome_annotations exits 0 without writing chromsizes.

    Its helper wraps the NCBI assembly-report lookup in a bare `except
    Exception`, prints "Chromosome sizes will not be returned", and returns the
    annotation alone; the CLI then logs "Chrosomome sizes was not found, please
    provide this information manually" and exits 0. The stage's sentinel is
    genome_annotation.tsv, which WAS written -- so the driver marks step 7 done
    and step 8 is the first to notice, by which point re-running step 7 is both
    skipped and futile. Fail here instead.
    """
    if Path(out["chromsizes"]).is_file():
        _check_chromsizes_shape(out["chromsizes"])
        return
    sys.exit(
        f"[grn_stage] the download exited 0 but wrote no {out['chromsizes']}.\n"
        f"  That is its documented behaviour when the NCBI assembly-report\n"
        f"  lookup fails (it catches every exception and returns the gene\n"
        f"  annotation alone), so it will not succeed on a retry from behind\n"
        f"  the same proxy. Supply the file instead -- in config.yaml:\n"
        f"      input:\n"
        f"        genome_annotation: {out['genome_annotation']}\n"
        f"        chromsizes: /path/to/chromsizes.tsv\n"
        f"  The genome annotation above was written and is reusable. To build\n"
        f"  chromsizes for hg38:\n"
        f"    curl -O https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes\n"
        f"    awk 'BEGIN{{OFS=\"\\t\"; print \"Chromosome\",\"Start\",\"End\"}} "
        f"{{print $1,0,$2}}' \\\n"
        f"        hg38.chrom.sizes > chromsizes.tsv")


def _diagnose_barcodes(args, cfg) -> None:
    """Show why ACC and GEX barcodes did not intersect, and the fix."""
    try:
        import pickle
        import anndata
        with open(args.cistopic_obj, "rb") as fh:
            acc = list(pickle.load(fh).cell_names)
        gex = list(anndata.read_h5ad(args.adata, backed="r").obs_names)
    except Exception as e:                      # never mask the real failure
        print(f"[grn_stage] (barcode diagnosis unavailable: {e})", file=sys.stderr)
        return

    fn = cfg.get("scenicplus", {}).get("bc_transform_func", "lambda x: x")
    print("\n[grn_stage] barcode shapes:", file=sys.stderr)
    print(f"[grn_stage]   ACC (cisTopic): {len(acc)} cells, e.g. {acc[:2]}", file=sys.stderr)
    print(f"[grn_stage]   GEX (anndata) : {len(gex)} cells, e.g. {gex[:2]}", file=sys.stderr)
    print(f"[grn_stage]   bc_transform_func in config: {fn!r}", file=sys.stderr)
    print("[grn_stage]   (it is applied to the GEX names, to map them onto ACC)",
          file=sys.stderr)

    overlap = len(set(acc) & set(gex))
    if overlap:
        print(f"[grn_stage]   {overlap} name(s) already match -- the mismatch is "
              f"elsewhere.", file=sys.stderr)
        return
    # The usual cause: pycisTopic's create_cistopic_object defaults to
    # tag_cells=True and appends "___<project>" to every ACC name.
    a0, g0 = acc[0], gex[0]
    if "___" in a0 and a0.split("___")[0] == g0:
        suffix = "___" + a0.split("___", 1)[1]
        print(f"[grn_stage]   ACC names carry the pycisTopic sample tag {suffix!r};",
              file=sys.stderr)
        print("[grn_stage]   GEX names do not. Either rebuild step 3 without the",
              file=sys.stderr)
        print("[grn_stage]   tag (--from 3, redoes topic modeling), or set in",
              file=sys.stderr)
        print("[grn_stage]   config.yaml under scenicplus::", file=sys.stderr)
        print(f'[grn_stage]     bc_transform_func: \'lambda x: x + "{suffix}"\'',
              file=sys.stderr)
    else:
        print("[grn_stage]   No shared names and no recognised tag pattern. The two",
              file=sys.stderr)
        print("[grn_stage]   sides came from different cell sets, or one was",
              file=sys.stderr)
        print("[grn_stage]   renamed. Compare the samples they were built from.",
              file=sys.stderr)


if __name__ == "__main__":
    sys.exit(main())
