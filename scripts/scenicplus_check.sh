#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Preflight check for the SCENIC+ pipeline.
#
# Verifies that the active Python and R environments have the packages the
# pipeline will need at runtime. Exits non-zero with a clear list of missing
# packages, so users fail fast instead of mid-pipeline.
#
# Usage:
#   scenicplus_check.sh                  # check whatever shell is active
#   SCENICPLUS_ENV=scenicplus scenicplus_check.sh
# -----------------------------------------------------------------------------
set -euo pipefail

# Optional: activate a conda env first, mirroring the launcher behavior.
#if [[ -n "${SCENICPLUS_ENV:-}" ]]; then
#    # shellcheck disable=SC1091
#    #source "$(conda info --base)/etc/profile.d/conda.sh"
#    #source activate "$SCENICPLUS_ENV"
#fi
#
FAIL=0

echo "[preflight] python: $(command -v python || echo MISSING)"
if ! command -v python >/dev/null 2>&1; then
    echo "[preflight] ERROR: python not on PATH" >&2
    FAIL=1
else
    python - <<'PY' || FAIL=1
import importlib, sys
required = [
    "yaml", "numpy", "pandas", "scipy", "matplotlib", "networkx",
    "pyranges", "scanpy", "anndata", "mudata",
    "scenicplus", "pycisTopic", "pycistarget",
]
missing = []
for m in required:
    try:
        print("Checking: " + m, file=sys.stderr)
        importlib.import_module(m)
    except Exception as e:
        missing.append(f"{m} ({type(e).__name__}: {e})")
if missing:
    print("[preflight] MISSING Python packages:", file=sys.stderr)
    for m in missing:
        print(f"  - {m}", file=sys.stderr)
    sys.exit(1)
print("[preflight] Python packages OK")
PY
fi

echo "[preflight] Rscript: $(command -v Rscript || echo MISSING)"
if ! command -v Rscript >/dev/null 2>&1; then
    echo "[preflight] ERROR: Rscript not on PATH" >&2
    FAIL=1
else
    Rscript - <<'R' || FAIL=1
required <- c("Seurat", "Matrix", "Signac", "optparse")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
    cat("[preflight] MISSING R packages:\n", file = stderr())
    for (m in missing) cat(sprintf("  - %s\n", m), file = stderr())
    quit(status = 1)
}
cat("[preflight] R packages OK\n")
R
fi

if [[ "$FAIL" -ne 0 ]]; then
    echo "[preflight] FAILED. Install the missing packages (see environment.yml) and retry." >&2
    exit 1
fi
echo "[preflight] All required packages are available."
