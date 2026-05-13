#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Initialize the SCENIC+ snakemake folder via `scenicplus init_snakemake`,
# then patch its config.yaml with paths and parameters from our top-level
# config.yaml. The patched directory is what `snakemake --cores N` runs.
# -----------------------------------------------------------------------------
set -euo pipefail

usage() {
    cat >&2 <<EOF
Usage: $(basename "$0") \\
    --config <top_config.yaml> \\
    --out_dir <scplus_pipeline_dir> \\
    --cistopic_obj <pkl> \\
    --adata <h5ad> \\
    --region_sets <dir> \\
    --scplus_out <dir>
EOF
    exit 2
}

CONFIG=""; OUT_DIR=""; CISTOPIC_OBJ=""; ADATA=""; REGION_SETS=""; SCPLUS_OUT=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)        CONFIG="$2";       shift 2;;
        --out_dir)       OUT_DIR="$2";      shift 2;;
        --cistopic_obj)  CISTOPIC_OBJ="$2"; shift 2;;
        --adata)         ADATA="$2";        shift 2;;
        --region_sets)   REGION_SETS="$2";  shift 2;;
        --scplus_out)    SCPLUS_OUT="$2";   shift 2;;
        -h|--help)       usage;;
        *) echo "[scenicplus_06_init_inner] unknown arg: $1" >&2; usage;;
    esac
done

if [[ -z "$CONFIG" || -z "$OUT_DIR" || -z "$CISTOPIC_OBJ" \
   || -z "$ADATA" || -z "$REGION_SETS" || -z "$SCPLUS_OUT" ]]; then
    echo "[scenicplus_06_init_inner] ERROR: missing required argument(s)" >&2
    usage
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "$OUT_DIR"
rm -rf "$OUT_DIR/Snakemake"

echo "[scenicplus_06_init_inner] scenicplus init_snakemake --out_dir $OUT_DIR"
scenicplus init_snakemake --out_dir "$OUT_DIR"

CFG_PATH="$OUT_DIR/Snakemake/config/config.yaml"
if [[ ! -f "$CFG_PATH" ]]; then
    echo "[scenicplus_06_init_inner] ERROR: $CFG_PATH not created by init_snakemake." >&2
    exit 1
fi

python "$SCRIPT_DIR/scenicplus_06_patch_config.py" \
    --config       "$CONFIG" \
    --scplus_cfg   "$CFG_PATH" \
    --cistopic_obj "$CISTOPIC_OBJ" \
    --adata        "$ADATA" \
    --region_sets  "$REGION_SETS" \
    --scplus_out   "$SCPLUS_OUT"

echo "[scenicplus_06_init_inner] Patched: $CFG_PATH"
