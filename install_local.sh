#!/usr/bin/env bash
# Build the SCENIC+ environment: conda layer from environment.local.yml, then
# the pip layer in two phases.
#
# Verified 2026-09-07 on linux-64 (WSL2, micromamba). Result: scenicplus 1.0a2,
# pycistarget, pycisTopic, snakemake 8.5.5, Seurat 5.5.1, Signac 1.17.1.
#
#   ./install_local.sh [PREFIX]      # default: ./scenicplus_env
#
#   SCP_COPY=1  ./install_local.sh ...   # copy instead of hardlink (SMB/CIFS,
#                                        # lustre, or anywhere hardlinks fail)
#   SCP_FRESH=1 ./install_local.sh ...   # delete an existing prefix first
#   SCP_CLEAN_PKGS=1 ./install_local.sh ...  # drop partially-extracted packages
#                                        # from the cache (CondaVerificationError
#                                        # "... appears to be corrupted")
#
# Re-running is safe: an existing env is updated in place, and phase 1 is
# skipped when pybedtools 0.9.1 is already installed. An interrupted install
# can simply be re-run.
#
# Works with micromamba, mamba or conda -- whichever is found.
#
# WHY TWO PIP PHASES. scenicplus ships a fully-frozen pin set, and two of the
# pinned versions predate python 3.11 and have only an sdist on PyPI:
#
#   pybedtools==0.9.1 -- must compile (no py311 conda build exists either), and
#     its sdist does not declare setuptools, so pip's ISOLATED build env has
#     none and the build dies with "setuptools was not found".
#     -> phase 1, --no-build-isolation, using this env's setuptools + cython.
#
#   loomxpy (git dep) -- needs the `poetry.masonry.api` build backend, which is
#     NOT in this env, so pip must fetch it into an isolated build env.
#     -> phase 2, isolation left ON (the default).
#
# The two requirements are contradictory in a single `pip install`, which is why
# conda's `pip:` section cannot express this and this script exists.
#
# CC/CXX are exported explicitly because conda ships compilers under prefixed
# names and sets CC/CXX only on `conda activate` -- and this pipeline
# deliberately does not activate (commit ae98960). Without them the pybedtools
# build looks for a bare `g++` and fails. Harmless if the system has its own.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_PREFIX="${1:-$HERE/scenicplus_env}"
YML="$HERE/environment.local.yml"

CONDA=""
for c in micromamba mamba conda; do
    if command -v "$c" >/dev/null 2>&1; then CONDA="$c"; break; fi
done
[[ -n "$CONDA" ]] || { echo "ERROR: need micromamba, mamba or conda on PATH" >&2; exit 1; }
echo "### using $CONDA -> $ENV_PREFIX"

# A shared/system conda install often has a package cache the user cannot write
# ("Could not open lockfile .../pkgs/cache/cache.lock"). Point the cache at a
# writable location rather than fighting it.
if [[ -z "${CONDA_PKGS_DIRS:-}" ]]; then
    export CONDA_PKGS_DIRS="${TMPDIR:-$HOME/tmp}/conda-pkgs-$USER"
    mkdir -p "$CONDA_PKGS_DIRS"
    echo "### CONDA_PKGS_DIRS=$CONDA_PKGS_DIRS (override by exporting it yourself)"
fi

# SMB/CIFS and other network filesystems break conda's default install method:
# it HARDLINKS from the package cache into the env, and those filesystems either
# refuse hardlinks or apply surprising ownership. Set SCP_COPY=1 to copy files
# instead -- slower and larger, but it survives. (Symptom without it: permission
# or "Operation not permitted" errors partway through the link step.)
COPY_ARGS=()
if [[ "${SCP_COPY:-0}" == "1" ]]; then
    export CONDA_ALWAYS_COPY=true          # conda / mamba
    COPY_ARGS=(--always-copy)              # micromamba
    echo "### SCP_COPY=1 -> copying instead of hardlinking (network filesystem mode)"
fi

# An interrupted install (killed job, filesystem hiccup) can leave PARTIALLY
# EXTRACTED package directories in the cache. conda then verifies each against
# its manifest and aborts with, e.g.
#
#   CondaVerificationError: The package for r-base located at
#   <pkgs>/r-base-4.5.3-h502d0c9_3 appears to be corrupted. The path
#   'share/man/man1/Rscript.1' specified in the package manifest cannot be found.
#
# The downloaded ARCHIVES are usually intact -- only the unpacked directories
# are truncated. SCP_CLEAN_PKGS=1 removes the unpacked directories and keeps the
# archives, so packages re-extract without re-downloading.
if [[ "${SCP_CLEAN_PKGS:-0}" == "1" ]]; then
    for d in ${CONDA_PKGS_DIRS//:/ }; do
        [[ -d "$d" ]] || continue
        n=$(find "$d" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | wc -l)
        echo "### SCP_CLEAN_PKGS=1 -> removing $n extracted package dir(s) from $d (archives kept)"
        find "$d" -mindepth 1 -maxdepth 1 -type d -exec rm -rf {} + 2>/dev/null || true
    done
fi

echo "### phase 0: conda layer"
# Creating an env FROM A YAML differs by tool, and getting it wrong is silent-ish:
#   micromamba create -f env.yml      -- accepts a YAML directly
#   conda/mamba env create -f env.yml -- NOTE the `env` subcommand
# `conda create -f` means a SPEC FILE (one package per line), so passing a YAML
# there makes conda read the PATH as a package name and report
# "<path> does not exist (perhaps a typo or a missing channel)".
#
# RESUME: an interrupted install leaves a partial prefix, and `create` then
# refuses it ("prefix already exists"). If one is there, install/update INTO it
# instead, which is idempotent -- already-satisfied packages are no-ops. Use
# --fresh to start over instead.
if [[ "${SCP_FRESH:-0}" == "1" && -d "$ENV_PREFIX" ]]; then
    echo "### SCP_FRESH=1 -> removing existing $ENV_PREFIX"
    rm -rf "$ENV_PREFIX"
fi

if [[ -d "$ENV_PREFIX/conda-meta" ]]; then
    echo "### existing env found at $ENV_PREFIX -- updating in place (resume)"
    case "$(basename "$CONDA")" in
        micromamba) "$CONDA" install -y -p "$ENV_PREFIX" -f "$YML" "${COPY_ARGS[@]}" ;;
        *)          "$CONDA" env update -p "$ENV_PREFIX" -f "$YML" ;;
    esac
elif [[ -d "$ENV_PREFIX" ]]; then
    # A directory with no conda-meta is a prefix that never got far enough to be
    # an env -- an interrupted first attempt. `create` would refuse it.
    echo "### $ENV_PREFIX exists but is not a conda env (interrupted?) -- creating into it"
    case "$(basename "$CONDA")" in
        micromamba) "$CONDA" create -y -p "$ENV_PREFIX" -f "$YML" "${COPY_ARGS[@]}" ;;
        *)          "$CONDA" env create --force -p "$ENV_PREFIX" -f "$YML" ;;
    esac
else
    case "$(basename "$CONDA")" in
        micromamba) "$CONDA" create -y -p "$ENV_PREFIX" -f "$YML" "${COPY_ARGS[@]}" ;;
        *)          "$CONDA" env create -p "$ENV_PREFIX" -f "$YML" ;;
    esac
fi

PY="$ENV_PREFIX/bin/python"
[[ -x "$PY" ]] || { echo "ERROR: $PY missing after create" >&2; exit 1; }

# Conda's compilers are prefixed; fall back to the system ones if absent.
for cc in "$ENV_PREFIX"/bin/*-cc; do [[ -x "$cc" ]] && export CC="$cc"; done
for cxx in "$ENV_PREFIX"/bin/*-c++; do [[ -x "$cxx" ]] && export CXX="$cxx"; done
export CPATH="${CPATH:+$CPATH:}$ENV_PREFIX/include"
echo "### CC=${CC:-<system>}  CXX=${CXX:-<system>}"

echo "### phase 1: pybedtools==0.9.1 (compiles; no isolation)"
if "$PY" -c 'import pybedtools,sys; sys.exit(0 if pybedtools.__version__=="0.9.1" else 1)' 2>/dev/null; then
    echo "  already at 0.9.1 -- skipping the compile"
else
    "$PY" -m pip install --no-build-isolation "pybedtools==0.9.1"
fi

echo "### phase 2: scenicplus (isolation ON)"
"$PY" -m pip install "scenicplus @ git+https://github.com/aertslab/scenicplus.git"

echo "### verify"
# Run from / so the env PREFIX DIRECTORY cannot be picked up as a namespace
# package -- `import scenicplus` from the parent dir "succeeds" with __file__
# None and proves nothing.
cd /
"$PY" - <<'PY'
import importlib.util as u
bad = 0
for m in ("scenicplus", "pycistarget", "pycisTopic", "pybedtools",
          "anndata", "mudata", "snakemake", "yaml"):
    s = u.find_spec(m)
    if s is None or s.origin is None:
        bad += 1
        print(f"  {m}: {'MISSING' if s is None else 'namespace-pkg (NOT a real install)'}")
    else:
        print(f"  {m}: ok")
raise SystemExit(1 if bad else 0)
PY
"$ENV_PREFIX/bin/Rscript" -e 'for (p in c("Seurat","Signac","Matrix","optparse")) cat(sprintf("  %-9s %s\n", p, as.character(packageVersion(p))))' 2>&1 | grep -v '^Loading'
"$ENV_PREFIX/bin/scenicplus" --help | head -2

cat <<EOF

### done. To use:
    export SCENICPLUS_PATH=$HERE
    export PATH=\$SCENICPLUS_PATH/scripts:$ENV_PREFIX/bin:\$PATH
    export SCRNA_CONDA_ENV=$ENV_PREFIX     # scrna.tool.sh-style prefix
EOF
