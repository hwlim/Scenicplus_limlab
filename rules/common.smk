# =============================================================================
# Helpers for the SCENIC+ workflow. NO RULES LIVE HERE, deliberately.
#
# Everything in this file answers "where does X go" or "what is this species
# called in that tool's vocabulary". A rule belongs in the modality-ish file
# that owns it (prepare / genome / grn / report), so that reading one of those
# files shows the whole of a stage rather than half of it.
#
# The one invariant this file exists to enforce: NO .smk FILE HARDCODES A
# DIRECTORY. Paths come from stage_path() and logs from log_path(). The same
# rule in scRNA_LimLab_Snake was learned the expensive way, when 71 hardcoded
# interpolations had to be found and changed at once.
# =============================================================================

import shlex

# --- Stage taxonomy ----------------------------------------------------------
# Numbered by PIPELINE STAGE, never by tool. SCENIC+ is linear and single
# sample, so it needs fewer stages than the multi-modality pipeline this is
# modelled on. Steps 3 to 5 collapse into `3.cistopic` because they are one
# object being refined, the same reason `3.samples` collapses over there.
STAGES = {
    "input":    "0.input",       # the Seurat .rds, or a symlink to it
    "export":   "1.export",      # seurat_export/: mtx, barcodes, regions, metadata,
                                 #   and one embedding_<name>.tsv per reduction
    "anndata":  "2.anndata",     # rna.h5ad
    "cistopic": "3.cistopic",    # cistopic_obj.pkl, *_with_topics.pkl, region_sets/
    "grn":      "4.grn",         # ACC_GEX.h5mu ... scplusmdata.h5mu
    "analysis": "5.analysis",    # tsv/ plots/
    "qc":       "QC",            # model selection, cell and region counts per stage
    "logs":     "logs",
    "provenance": "provenance",
}


def stage_dir(stage):
    """The directory for a named stage. Fails loudly on a typo'd stage name."""
    if stage not in STAGES:
        raise KeyError(
            f"unknown stage {stage!r}. Known: {', '.join(sorted(STAGES))}. "
            f"Add it to STAGES rather than writing the directory into a rule.")
    return STAGES[stage]


def stage_path(stage, *parts):
    """A path inside a stage: stage_path('grn', 'search_space.tsv')."""
    return os.path.join(stage_dir(stage), *parts)


def log_path(rule_name):
    """One log per rule, all in one place.

    `bsub -o` APPENDS and snakemake's {jobid} restarts at 0 each run, so a log
    file can hold several runs' worth of epilogue. Read the LAST block. This is
    why rule bodies use `2>&1 | tee {log}` with pipefail rather than `> {log}`:
    truncation would race the LSF epilogue that keeps being appended.
    """
    return os.path.join(stage_dir("logs"), f"{rule_name}.log")


def opt_arg(flag, value):
    """`--flag value`, or nothing at all when the value is empty.

    Do NOT write `--flag {params.x:q}` for a value that can be empty.
    Snakemake's `:q` renders an empty string as NOTHING, not as `''` --
    measured:

        params: empty="", full="x"
        shell:  "echo A {params.empty:q} B {params.full:q} C"
        ->      echo A  B x C

    so the flag loses its argument and swallows whatever came next. R01 died
    exactly that way: `--celltype_scope '' --reduction ''` reached optparse as
    `--celltype_scope --reduction`, and optparse reported that
    `celltype_scope` requires an argument -- an error about the flag AFTER the
    one that was actually wrong.

    Every option this is used for defaults to empty in the script itself, so
    omitting the flag is exactly equivalent to passing an empty one.
    """
    v = "" if value is None else str(value)
    return f"{flag} {shlex.quote(v)}" if v else ""


def script_path(name):
    """A shipped script, addressed through the workflow's own directory.

    DECLARE IT AS AN `input:`, NOT A `params:`. Snakemake's `code` rerun
    trigger covers a rule's own text, NOT the content of a script the rule
    shells out to -- measured: edit the script a rule calls and snakemake
    reports "Nothing to be done". As an input, the same edit reports "updated
    input files" and the rule re-runs.

    That distinction is the whole reason this workflow exists. The bash driver
    hashes CONFIG and not CODE, which is how a fixed step 3 never re-ran and
    tagged cell names survived the fix. Putting the script in `params` would
    have reproduced that bug faithfully.

    Not $SCENICPLUS_PATH: the workflow already knows where it lives, and a rule
    that depends on an environment variable is a rule that behaves differently
    depending on who launched it.
    """
    p = os.path.join(workflow.basedir, "scripts", name)
    if not os.path.exists(p):
        raise FileNotFoundError(f"no such script: {p}")
    return p


def cfg_params(*keys):
    """The config values a rule's script reads, as a params dict.

    THIS IS THE RERUN TRIGGER, and it is not decoration. The step scripts take
    `--config` and read the file themselves, so from Snakemake's point of view
    the config is neither an input nor a param, and editing it changes NOTHING.
    That exact shape is a live bug in the sibling pipeline, where one modality
    passes config values inline and the other passes only the filename: a config
    edit re-runs the first and is a silent no-op for the second.

    Declaring the values here puts them in the `params` trigger, so changing
    `cistopic.n_topics` re-runs topic modeling and nothing else. The key lists
    come from the bash driver's own `.cfgsha` slices, which already had to work
    this out per step -- reusing them means the two drivers agree about what
    each step depends on.

    Dotted keys, missing ones omitted, so a config without an optional key
    hashes the same as it did before that key existed.
    """
    out = {}
    for k in keys:
        cur, ok = config, True
        for part in k.split("."):
            if isinstance(cur, dict) and part in cur:
                cur = cur[part]
            else:
                ok = False
                break
        if ok:
            out[k] = cur
    return out


def n_cpu():
    """Threads for the steps that take --n_cpu. One place, so a rule and the
    value it passes to the script cannot disagree."""
    return int(config.get("resources", {}).get("n_cpu", 1))


# --- Species vocabularies ----------------------------------------------------
# `input.species` feeds calls that want DIFFERENT spellings of the same
# organism: the genome-annotation step wants `mmusculus`, the motif-enrichment
# steps want `mus_musculus`. Today that only works because the motif annotation
# file is always supplied, so the species argument is never consulted there
# (pycistarget/utils.py:98). One table, one lookup, so the two cannot disagree.
#
# Keyed by the spelling the config already uses, because this workflow reads the
# SAME config as the bash driver. SnakemakePlan.md's "one spelling, derived not
# repeated" is a config migration and belongs with the genome work in I2, not
# here, where it would break every existing workspace for no gain yet.
SPECIES = {
    "hsapiens": {
        "pycistarget": "homo_sapiens",
        "assembly":    "hg38",
        "ensdb":       "EnsDb.Hsapiens.v86",
        "bsgenome":    "BSgenome.Hsapiens.UCSC.hg38",
        "motif_tbl":   "hgnc",
    },
    "mmusculus": {
        "pycistarget": "mus_musculus",
        "assembly":    "mm10",
        "ensdb":       "EnsDb.Mmusculus.v79",
        "bsgenome":    "BSgenome.Mmusculus.UCSC.mm10",
        "motif_tbl":   "mgi",
    },
}


def species_info(cfg):
    """Look up the current species, or refuse with the reason.

    The schema allows the four species the SCENIC+ CLI names, because the bash
    driver passes the value straight through and one config serves both
    drivers. This table carries only the two the lab runs. Guessing the EnsDb
    and BSgenome package names for the others would be inventing facts, and a
    wrong genome reference is the failure that does not announce itself, so an
    unsupported species stops here instead.
    """
    s = cfg["input"]["species"]
    if s not in SPECIES:
        raise ValueError(
            # NOT __file__: inside an included .smk that resolves to snakemake's
            # own module path, which sends the reader to the wrong file.
            f"input.species = {s!r} has no entry in SPECIES "
            f"(rules/common.smk).\n"
            f"  Supported here: {', '.join(sorted(SPECIES))}.\n"
            f"  The config schema is deliberately wider, because the bash "
            f"driver passes species straight to the CLI. To add one, fill in "
            f"its row -- do not guess the EnsDb or BSgenome package name.")
    return SPECIES[s]


# --- Resource tiers ----------------------------------------------------------
# Memory, runtime and THREADS travel together, on purpose. Setting them apart is
# how scRNA_LimLab_Snake ended up submitting `-n {threads}` with threads
# defaulting to 1 while the script forked 16 workers: a 16x oversubscription
# that no single file made visible.
#
# PROVISIONAL. There is exactly one tier today, and it holds what the single
# bsub'd driver job currently reserves for all twenty steps at once. That is not
# a measurement of any step; it is the status quo, expressed so that rules can
# start asking for resources by name. I5 is where per-rule tiers get real
# numbers, from `logs/lsf/driver_*.out` and a run that records peak RSS per
# step. Adding a tier with an invented number before then would look like
# evidence.
TIERS = {
    "default": {"mem_mb": 128000, "runtime": 72 * 60, "threads": 16},
}


def tier(name, cfg=None):
    if name not in TIERS:
        raise KeyError(f"unknown resource tier {name!r}. Known: {', '.join(TIERS)}")
    t = dict(TIERS[name])
    # resources.n_cpu is the number the STEP SCRIPTS use internally, so a rule's
    # thread count must not exceed it or the reservation lies about the shape of
    # the work.
    if cfg is not None:
        t["threads"] = min(t["threads"], cfg.get("resources", {}).get("n_cpu", t["threads"]))
    return t


def mem(name, cfg=None):
    return tier(name, cfg)["mem_mb"]


def rt(name, cfg=None):
    return tier(name, cfg)["runtime"]


def cpus(name, cfg=None):
    return tier(name, cfg)["threads"]
