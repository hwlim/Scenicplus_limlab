#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Preflight check for the SCENIC+ pipeline.
#
# Verifies that the ACTIVE python and R can do what the pipeline will ask of
# them, and exits non-zero naming what cannot, so a run fails in seconds rather
# than at step 6 after topic modeling has already burned an hour.
#
#   scenicplus_check.sh
#
# WHY THIS CHECKS SUBMODULES AND NOT JUST PACKAGES
# ------------------------------------------------
# It used to `import pycistarget` and call that good enough. Measured, that
# imports NOTHING of what fails:
#
#     import pycistarget                                     -> OK
#       pycistarget.motif_enrichment_cistarget loaded?  False
#       IPython loaded?                                 False
#       sqlite3 loaded?                                 False
#
# The package __init__ is light. The chain that actually breaks only runs when
# the CLI imports a submodule:
#
#     scenicplus.cli.commands
#       -> pycistarget.motif_enrichment_cistarget
#         -> pycistarget.motif_enrichment_result
#           -> IPython.display
#             -> IPython.core.history
#               -> sqlite3
#                 -> _sqlite3 -> libicui18n.so.78 -> libstdc++.so.6
#
# and the last link is where a conda env on an old host comes apart:
#
#     ImportError: /lib64/libstdc++.so.6: version `CXXABI_1.3.15' not found
#       (required by <prefix>/lib/python3.11/lib-dynload/../.././libicui18n.so.78)
#
# The env's ICU has DT_RPATH `$ORIGIN/`, which the loader searches BEFORE
# LD_LIBRARY_PATH -- so the only way it reaches /lib64 is if <prefix>/lib has no
# libstdc++.so.6 at all. On this pipeline's history that means a package that
# extracted only partially on the network filesystem: conda-meta says installed,
# the files are not there. Which is why the ABI check below looks at the FILE and
# not just at the import.
#
# This check PASSED on the env that then failed at step 6, which is the whole
# reason the header above is this long: it was testing a weaker property than
# the pipeline needs, and a green preflight is worse than none when it is wrong.
#
# Every check here is one the pipeline will make anyway. The point is to make
# them all at the start, cheaply, on the machine that will run the work -- the
# launchers run this from the driver, inside the job, not on the submit host.
# -----------------------------------------------------------------------------
set -uo pipefail          # NOT -e: every check must run, so the report is whole

FAIL=0

echo "[preflight] python:  $(command -v python || echo MISSING)"
echo "[preflight] Rscript: $(command -v Rscript || echo MISSING)"

# --- the C++ ABI, before anything imports ------------------------------------
# Reported even when the imports pass: a run that works today on one node class
# and not on another is this, and the two numbers below are the evidence.
if command -v python >/dev/null 2>&1; then
    python - <<'PY' || FAIL=1
import os, subprocess, sys, sysconfig
prefix = sys.prefix
libdir = os.path.join(prefix, "lib")
so = os.path.join(libdir, "libstdc++.so.6")
if not os.path.exists(so):
    print(f"[preflight] ERROR: {so} is MISSING.", file=sys.stderr)
    print("[preflight]   The env's own C++ runtime is not there, so anything in", file=sys.stderr)
    print("[preflight]   it that needs libstdc++ falls through to the system one,", file=sys.stderr)
    print("[preflight]   which on an older host is too old:", file=sys.stderr)
    print("[preflight]     ImportError: /lib64/libstdc++.so.6: version `CXXABI_1.3.15' not found", file=sys.stderr)
    print("[preflight]   Repair (a partial extraction leaves conda-meta claiming", file=sys.stderr)
    print("[preflight]   it is installed, so a plain install is a no-op):", file=sys.stderr)
    print(f"[preflight]     conda install -p {prefix} -c conda-forge --force-reinstall libstdcxx-ng", file=sys.stderr)
    sys.exit(1)
import re
# Only the NUMBERED versions sort. libstdc++ also exports CXXABI_FLOAT128,
# CXXABI_TM_1 and friends, and feeding those to int() crashes the check on a
# perfectly healthy env -- which is exactly what the first version of this
# block did. Reporting is all this is for; a missing FILE is the only thing
# here worth failing on, and that was handled above.
vers = []
try:
    out = subprocess.run(["strings", so], capture_output=True, text=True,
                         timeout=120).stdout
    for line in out.splitlines():
        m = re.fullmatch(r"CXXABI_(\d+(?:\.\d+)*)", line)
        if m:
            vers.append(tuple(int(x) for x in m.group(1).split(".")))
except Exception as e:                      # strings absent on a bare image
    print(f"[preflight] (could not read {so}: {e})", file=sys.stderr)
top = ".".join(str(x) for x in max(vers)) if vers else "unknown"
print(f"[preflight] libstdc++: {os.path.realpath(so)}  (max CXXABI_{top})")
PY
fi

# --- python: the submodules the CLI actually imports -------------------------
if ! command -v python >/dev/null 2>&1; then
    echo "[preflight] ERROR: python not on PATH" >&2
    FAIL=1
else
    python - <<'PY' || FAIL=1
import importlib, sys

# Plain packages: cheap, and a missing one is unambiguous.
packages = [
    "yaml", "numpy", "pandas", "scipy", "matplotlib", "networkx",
    "pyranges", "scanpy", "anndata", "mudata",
    "scenicplus", "pycisTopic", "pycistarget",
]
# The submodules that carry the real import chains. `import pycistarget` does
# not reach any of these -- measured, see this file's header.
submodules = [
    "sqlite3",                                  # the libstdc++ canary
    "pycistarget.motif_enrichment_cistarget",   # -> IPython -> sqlite3
    "pycistarget.motif_enrichment_dem",
    "scenicplus.cli.commands",                  # what every GRN stage imports
    "pycisTopic.lda_models",                    # step 4
    "ray",                                      # step 4's parallelism
]

bad = []
for m in packages + submodules:
    try:
        importlib.import_module(m)
    except Exception as e:
        bad.append((m, f"{type(e).__name__}: {e}"))

if bad:
    print("[preflight] python imports that FAILED:", file=sys.stderr)
    for m, why in bad:
        print(f"  - {m}\n      {why}", file=sys.stderr)
    sys.exit(1)
print(f"[preflight] python: {len(packages)} package(s) + {len(submodules)} "
      f"submodule chain(s) OK")
PY
fi

# --- the console script, run as the pipeline runs it -------------------------
# scripts/scenicplus_06_grn_stage.py execs `scenicplus` as a COMMAND, so PATH
# resolution and the entry point matter, not just importability.
if ! command -v scenicplus >/dev/null 2>&1; then
    echo "[preflight] ERROR: 'scenicplus' is not on PATH." >&2
    echo "[preflight]   Every GRN stage execs it as a command. If the package" >&2
    echo "[preflight]   imports but this is missing, pip installed it outside the" >&2
    echo "[preflight]   env -- see install_local.sh." >&2
    FAIL=1
elif ! scenicplus --help >/dev/null 2>&1; then
    echo "[preflight] ERROR: '$(command -v scenicplus) --help' exits non-zero:" >&2
    scenicplus --help 2>&1 | tail -20 | sed 's/^/    /' >&2
    FAIL=1
else
    echo "[preflight] scenicplus CLI: $(command -v scenicplus) OK"
fi

# --- R -----------------------------------------------------------------------
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
    echo "[preflight] FAILED -- fix the above before running. Set" >&2
    echo "[preflight] SCENICPLUS_SKIP_CHECK=1 to run anyway, at the cost of" >&2
    echo "[preflight] finding out at whichever step first needs the broken thing." >&2
    exit 1
fi
echo "[preflight] all checks passed"
