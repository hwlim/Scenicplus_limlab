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
#   SCP_FROM=1  ./install_local.sh ...   # skip phase 0; the conda layer is done
#   SCP_FROM=2  ./install_local.sh ...   # skip phases 0-1; only pip scenicplus
#
# SCP_FROM is for when phase 0 has already SUCCEEDED and re-solving it is the
# expensive part -- an old conda over a network filesystem can take longer to
# re-verify a finished env than the pip layer takes to run. It refuses to skip
# phase 0 if <prefix>/bin/python is absent, so it cannot silently proceed
# against an env that was never built.
#
# Re-running is safe: an existing env is updated in place, and phase 1 is
# skipped when pybedtools 0.9.1 is already installed. An interrupted install
# can simply be re-run.
#
# IF $ENV_PREFIX/bin/scenicplus IS MISSING after an otherwise clean run, pip
# installed the package outside this prefix -- almost always the user site,
# because ~/.config/pip/pip.conf says `user = true` or PIP_USER is set. Repair:
#
#   SCP_FROM=2 ./install_local.sh <prefix>      # pip layer only, ~2 min
#
# See the PIP_USER / PYTHONNOUSERSITE block below for why the import checks used
# to pass anyway.
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

SCP_FROM="${SCP_FROM:-0}"

if [[ "$SCP_FROM" -gt 0 ]]; then
    echo "### SCP_FROM=$SCP_FROM -> skipping phase 0 (conda layer)"
    # Skipping the conda layer means trusting an env that is already there. Say
    # so loudly if it is not, rather than failing later inside pip with a
    # confusing message about a missing interpreter.
    if [[ ! -x "$ENV_PREFIX/bin/python" ]]; then
        echo "ERROR: SCP_FROM=$SCP_FROM skips phase 0, but $ENV_PREFIX/bin/python" >&2
        echo "       does not exist. Run without SCP_FROM to build the conda layer first." >&2
        exit 1
    fi
else
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

fi

PY="$ENV_PREFIX/bin/python"
[[ -x "$PY" ]] || { echo "ERROR: $PY missing after create" >&2; exit 1; }

# pip must install INTO $ENV_PREFIX, never into the USER SITE
# (~/.local/lib/pythonX.Y/site-packages). Two things send it there on a shared
# cluster, and neither is loud:
#
#   * `user = true` in ~/.config/pip/pip.conf (or PIP_USER=1 in the
#     environment) -- pip installs the package under ~/.local and puts the
#     console scripts in ~/.local/bin. `import scenicplus` then works while
#     $ENV_PREFIX/bin/scenicplus never appears;
#   * a copy ALREADY in the user site -- pip reports "Requirement already
#     satisfied" and installs nothing at all, with the same result.
#
# The user site sits BEFORE the env's site-packages on sys.path (measured, not
# assumed), so imports resolve to it either way and every import-based check
# passes on an env that does not actually contain the package. PIP_USER=0
# overrides a pip.conf that says otherwise (measured: with `user = true` in the
# config, `pip install six` landed in the user site; with PIP_USER=0 layered on
# top, the same command landed in the env prefix). PYTHONNOUSERSITE=1 hides the
# user site from both pip's already-satisfied check and our own verify.
export PIP_USER=0
export PYTHONNOUSERSITE=1

# Conda's compilers are prefixed; fall back to the system ones if absent.
for cc in "$ENV_PREFIX"/bin/*-cc; do [[ -x "$cc" ]] && export CC="$cc"; done
for cxx in "$ENV_PREFIX"/bin/*-c++; do [[ -x "$cxx" ]] && export CXX="$cxx"; done
export CPATH="${CPATH:+$CPATH:}$ENV_PREFIX/include"
echo "### CC=${CC:-<system>}  CXX=${CXX:-<system>}"

if [[ "$SCP_FROM" -gt 1 ]]; then
    echo "### SCP_FROM=$SCP_FROM -> skipping phase 1 (pybedtools)"
else
echo "### phase 1: pybedtools==0.9.1 (compiles; no isolation)"
if "$PY" -c 'import pybedtools,sys; sys.exit(0 if pybedtools.__version__=="0.9.1" else 1)' 2>/dev/null; then
    echo "  already at 0.9.1 -- skipping the compile"
else
    "$PY" -m pip install --no-build-isolation "pybedtools==0.9.1"
fi

fi

echo "### phase 2: scenicplus (isolation ON)"
"$PY" -m pip install "scenicplus @ git+https://github.com/aertslab/scenicplus.git"

echo "### verify"
# Run from / so the env PREFIX DIRECTORY cannot be picked up as a namespace
# package -- `import scenicplus` from the parent dir "succeeds" with __file__
# None and proves nothing.
cd /
# "importable" is NOT the same as "installed here": the module can resolve from
# the user site or from PYTHONPATH, so each origin is checked against the
# prefix. Without that, an install that went to ~/.local passes every line.
SCP_PREFIX="$ENV_PREFIX" "$PY" - <<'PY'
import importlib.util as u, os
prefix = os.path.realpath(os.environ["SCP_PREFIX"]) + os.sep
bad = 0
for m in ("scenicplus", "pycistarget", "pycisTopic", "pybedtools",
          "anndata", "mudata", "snakemake", "yaml"):
    s = u.find_spec(m)
    if s is None or s.origin is None:
        bad += 1
        print(f"  {m}: {'MISSING' if s is None else 'namespace-pkg (NOT a real install)'}")
        continue
    origin = os.path.realpath(s.origin)
    if not origin.startswith(prefix):
        bad += 1
        print(f"  {m}: OUTSIDE THE ENV -> {origin}")
    else:
        print(f"  {m}: ok")
raise SystemExit(1 if bad else 0)
PY
"$ENV_PREFIX/bin/Rscript" -e 'for (p in c("Seurat","Signac","Matrix","optparse")) cat(sprintf("  %-9s %s\n", p, as.character(packageVersion(p))))' 2>&1 | grep -v '^Loading'

# The console script, and a repair when it is absent. scenicplus 1.0a2 declares
#   [project.scripts] scenicplus = "scenicplus.cli.scenicplus:main"
# so pip writes $ENV_PREFIX/bin/scenicplus whenever it installs the package into
# this prefix. Absent means the install went somewhere else, or pip decided it
# was already satisfied and wrote nothing. This is not cosmetic: every GRN stage
# in scripts/scenicplus_06_grn_stage.py execs `scenicplus` as a command.
if [[ ! -x "$ENV_PREFIX/bin/scenicplus" ]]; then
    echo "### $ENV_PREFIX/bin/scenicplus missing -- reinstalling scenicplus alone"
    echo "###   (--no-deps: the pinned dependency layer above must NOT be redone;"
    echo "###    a plain --force-reinstall would try to rebuild pybedtools 0.9.1"
    echo "###    WITH build isolation and fail)"
    "$PY" -m pip install --force-reinstall --no-deps \
        "scenicplus @ git+https://github.com/aertslab/scenicplus.git"
fi
if [[ ! -x "$ENV_PREFIX/bin/scenicplus" ]]; then
    echo "ERROR: $ENV_PREFIX/bin/scenicplus is still missing after a forced" >&2
    echo "       reinstall. pip put the package somewhere this prefix cannot see," >&2
    echo "       or could not write to $ENV_PREFIX/bin. Diagnosis:" >&2
    env -u PYTHONNOUSERSITE SCP_PREFIX="$ENV_PREFIX" "$PY" - >&2 <<'PY'
import importlib.util as u, os, site, sys
s = u.find_spec("scenicplus")
print("       with the user site ENABLED, scenicplus resolves to:",
      s.origin if s else "NOT FOUND")
print("       user site:", site.getusersitepackages())
print("       env  site:", os.path.join(os.environ["SCP_PREFIX"], "lib",
                                        "python%d.%d" % sys.version_info[:2],
                                        "site-packages"))
PY
    ls -d "$HOME"/.local/bin/scenicplus 2>/dev/null \
        | sed 's/^/       console script landed here instead: /' >&2
    echo "       Check for 'user = true' in ~/.config/pip/pip.conf and for" >&2
    echo "       PIP_USER in the environment, then re-run with SCP_FROM=2." >&2
    exit 1
fi
"$ENV_PREFIX/bin/scenicplus" --help | head -2

cat <<EOF

### done. To use:
    export SCENICPLUS_PATH=$HERE
    export PATH=\$SCENICPLUS_PATH/scripts:$ENV_PREFIX/bin:\$PATH
    export SCRNA_CONDA_ENV=$ENV_PREFIX     # scrna.tool.sh-style prefix
    export PYTHONNOUSERSITE=1              # see below
EOF

cat >&2 <<'EOF'

PYTHONNOUSERSITE=1 is not optional hygiene. ~/.local/lib/pythonX.Y/site-packages
comes BEFORE the env's site-packages on sys.path, so any copy of scenicplus,
pycisTopic or pycistarget left there by an earlier `pip install --user` silently
WINS over the one in this env -- at run time, in every step, with no message.
EOF
