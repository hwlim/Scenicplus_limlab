#!/usr/bin/env python
"""
Patch the SCENIC+ snakemake config.yaml (produced by
`scenicplus init_snakemake`) in place, filling input/output paths and
parameters from our top-level config.yaml.
"""
import argparse
from pathlib import Path

import yaml


def patch_config(scplus_cfg_path: Path, cfg: dict, paths: dict) -> None:
    with open(scplus_cfg_path) as fh:
        sc_cfg = yaml.safe_load(fh)

    sc_cfg["input_data"] = {
        "cisTopic_obj_fname": str(Path(paths["cistopic_obj"]).resolve()),
        "GEX_anndata_fname":  str(Path(paths["adata"]).resolve()),
        "region_set_folder":  str(Path(paths["region_sets"]).resolve()),
        "ctx_db_fname":       str(Path(cfg["input"]["ctx_db"]).resolve()),
        "dem_db_fname":       str(Path(cfg["input"]["dem_db"]).resolve()),
        "path_to_motif_annotations": str(Path(cfg["input"]["motif_annotations"]).resolve()),
    }

    out_root = Path(paths["scplus_out"]).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    sc_cfg["output_data"] = {
        "combined_GEX_ACC_mudata":     str(out_root / "ACC_GEX.h5mu"),
        "dem_result_fname":            str(out_root / "dem_results.hdf5"),
        "ctx_result_fname":            str(out_root / "ctx_results.hdf5"),
        "output_fname_dem_html":       str(out_root / "dem_results.html"),
        "output_fname_ctx_html":       str(out_root / "ctx_results.html"),
        "cistromes_direct":            str(out_root / "cistromes_direct.h5ad"),
        "cistromes_extended":          str(out_root / "cistromes_extended.h5ad"),
        "tf_names":                    str(out_root / "tf_names.txt"),
        "genome_annotation":           str(out_root / "genome_annotation.tsv"),
        "chromsizes":                  str(out_root / "chromsizes.tsv"),
        "search_space":                str(out_root / "search_space.tsv"),
        "tf_to_gene_adjacencies":      str(out_root / "tf_to_gene_adj.tsv"),
        "region_to_gene_adjacencies":  str(out_root / "region_to_gene_adj.tsv"),
        "eRegulons_direct":            str(out_root / "eRegulons_direct.tsv"),
        "eRegulons_extended":          str(out_root / "eRegulons_extended.tsv"),
        "AUCell_direct":               str(out_root / "AUCell_direct.h5mu"),
        "AUCell_extended":             str(out_root / "AUCell_extended.h5mu"),
        "scplus_mdata":                str(out_root / "scplusmdata.h5mu"),
    }

    sc_cfg["params_general"] = {
        "temp_dir": str(Path(cfg["output"]["tmp"]).resolve()),
        "n_cpu": cfg["resources"]["n_cpu"],
        "seed":  cfg["resources"]["seed"],
    }

    sp = cfg["scenicplus"]
    sc_cfg["params_data_preparation"] = {
        "bc_transform_func":      sp["bc_transform_func"],
        "is_multiome":            sp["is_multiome"],
        "key_to_group_by":        "",
        "nr_cells_per_metacells": 10,
        "direct_annotation":      sp["direct_annotation"],
        "extended_annotation":    sp["extended_annotation"],
        "species":                cfg["input"]["species"],
        "biomart_host":           cfg["input"]["biomart_host"],
        "search_space_upstream":   sp["search_space_upstream"],
        "search_space_downstream": sp["search_space_downstream"],
        "search_space_extend_tss": sp["search_space_extend_tss"],
    }

    sc_cfg["params_motif_enrichment"] = {
        "species":             cfg["input"]["species"],
        "annotation_version":  sp["annotation_version"],
        "motif_similarity_fdr": sp["motif_similarity_fdr"],
        "orthologous_identity_threshold": sp["orthologous_identity_threshold"],
        "annotations_to_use":  sp["annotations_to_use"],
        "fraction_overlap_w_dem_database": sp["fraction_overlap_w_dem_database"],
        "dem_max_bg_regions":  sp["dem_max_bg_regions"],
        "dem_balance_number_of_promoters": sp["dem_balance_number_of_promoters"],
        "dem_promoter_space":  sp["dem_promoter_space"],
        "dem_adj_pval_thr":    sp["dem_adj_pval_thr"],
        "dem_log2fc_thr":      sp["dem_log2fc_thr"],
        "dem_mean_fg_thr":     sp["dem_mean_fg_thr"],
        "dem_motif_hit_thr":   sp["dem_motif_hit_thr"],
        "fraction_overlap_w_ctx_database": sp["fraction_overlap_w_ctx_database"],
        "ctx_auc_threshold":   sp["ctx_auc_threshold"],
        "ctx_nes_threshold":   sp["ctx_nes_threshold"],
        "ctx_rank_threshold":  sp["ctx_rank_threshold"],
    }

    g = cfg["grn"]
    sc_cfg["params_inference"] = {
        "tf_to_gene_importance_method":     g["tf_to_gene_importance_method"],
        "region_to_gene_importance_method": g["region_to_gene_importance_method"],
        "region_to_gene_correlation_method": g["region_to_gene_correlation_method"],
        "order_regions_to_genes_by":        g["order_regions_to_genes_by"],
        "order_TFs_to_genes_by":            g["order_TFs_to_genes_by"],
        "gsea_n_perm":                      g["gsea_n_perm"],
        "quantile_thresholds_region_to_gene": g["quantile_thresholds_region_to_gene"],
        "top_n_regionTogenes_per_gene":     g["top_n_regionTogenes_per_gene"],
        "top_n_regionTogenes_per_region":   g["top_n_regionTogenes_per_region"],
        "min_regions_per_gene":             g["min_regions_per_gene"],
        "rho_threshold":                    g["rho_threshold"],
        "min_target_genes":                 g["min_target_genes"],
    }

    with open(scplus_cfg_path, "w") as fh:
        yaml.safe_dump(sc_cfg, fh, sort_keys=False)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config", required=True, help="Top-level config.yaml")
    p.add_argument("--scplus_cfg", required=True,
                   help="SCENIC+ snakemake config.yaml to patch in place")
    p.add_argument("--cistopic_obj", required=True)
    p.add_argument("--adata", required=True)
    p.add_argument("--region_sets", required=True)
    p.add_argument("--scplus_out", required=True)
    args = p.parse_args()

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)

    patch_config(
        Path(args.scplus_cfg),
        cfg,
        paths={
            "cistopic_obj": args.cistopic_obj,
            "adata":        args.adata,
            "region_sets":  args.region_sets,
            "scplus_out":   args.scplus_out,
        },
    )


if __name__ == "__main__":
    main()
