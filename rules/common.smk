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


def pin_env():
    """The environment that makes a step's output reproducible. MEASURED.

    Two independent sources of run-to-run variation, both found by comparing a
    driver run against a workflow run of the same data:

    **PYTHONHASHSEED.** Python randomises string hashing per PROCESS, and
    SCENIC+ does `list(set(names))` in `utils.py` (lines 394, 404, 405), so that
    list comes out in a different order every invocation -- on the same host,
    with the same input. Pandas sorts are stable, so a permuted input permutes
    the ties, and a top-N cut over tied values then keeps different rows. That
    is the difference seen in the cistromes: same entries, permuted order.
    Freezing the seed makes runs comparable; it does NOT make the ordering
    meaningful, so an eRegulon set still deserves a stability caveat.

    **The BLAS thread variables.** Reduction order follows thread count and
    thread count defaults to the host's core count, so a 48-core and a 64-core
    node give correlations that differ in the last bit. Measured: `rho` differed
    by 1e-16 on 9 rows of 1,428,119, while the seeded gradient boosting was
    bit-identical across four hosts. Pinning made two different hosts agree
    exactly.

    Pinned to `resources.n_cpu`, NOT to the rule's `threads`, on purpose. The
    bash driver runs all twenty steps in one job under one thread count, so
    matching that is what makes the two drivers comparable at all. The cost is
    that a one-slot rule may run BLAS with more threads than it reserved; every
    rule where that applies finishes in under three minutes.
    """
    n = n_cpu()
    return (
        "export PYTHONHASHSEED=0"
        f" OMP_NUM_THREADS={n}"
        f" OPENBLAS_NUM_THREADS={n}"
        f" MKL_NUM_THREADS={n}"
        f" NUMEXPR_NUM_THREADS={n}; "
    )


def shell_prefix():
    """What every rule body runs before its own command.

    `set -o pipefail` and nothing else from the shell side: rule bodies end
    `2>&1 | tee {log}`, and WITHOUT pipefail that pipeline's status is tee's,
    which is 0 -- a failing script reporting success. Not `-e`, because bodies
    chain with `&&`; not `-u`, because site `module` functions dereference unset
    variables.

    One prefix rather than eighteen rule edits, so the pinning cannot be applied
    to some rules and not others -- which would be worse than not applying it,
    since the run would look reproducible and not be.
    """
    return "set -o pipefail; " + pin_env()


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
# MEASURED, 2026-09-09/10, from one full Snakemake run of all 18 step rules on
# the PBMC multiome fixture. Every figure behind these tiers is in
# `tests/measured_resources.tsv`, with the run's provenance and its two caveats
# at the top of that file; `tests/test_resources.py` re-derives every assignment
# below from it and fails if a tier stops covering its rule.
#
# TWO AXES, NOT ONE, because the measurement says they do not correlate. R04
# runs 78 minutes in 21 GB; R09 finishes in 11 minutes and wants 149 GB. Folding
# those into one ladder forces every long rule to buy memory it never touches,
# or every large one to buy hours. So a rule names a memory tier and a time tier
# separately, and `tests/test_resources.py` checks each against its own figure.
#
# THREADS ARE NOT A TIER. They come from `n_cpu()`, which reads
# `resources.n_cpu` -- the same number the step scripts pass to the tools that
# fork. A tier carrying a thread count could disagree with it, and a reservation
# that disagrees with the work is the 16x oversubscription scRNA_LimLab_Snake
# spent a release chasing. One source, config, for both sides.

# Memory ladder, MB, named for its own size in decimal GB. The names are not
# t-shirt sizes on purpose: a rule asks for its tier BY RULE NAME (see `mem()`),
# so these labels appear nowhere but this file, and a label that states its own
# number cannot drift from it the way "large" can.
MEM_TIERS = {
    "4g":     4000,
    "16g":   16000,
    "32g":   32000,
    "48g":   48000,
    "96g":   96000,
    "128g": 128000,
    "192g": 192000,
    "256g": 256000,
}

# Wall-clock ladder, MINUTES -- the unit snakemake's `runtime` resource uses and
# the unit `-W` takes in the profile. Only two, because only one rule in the
# whole workflow is slow: R04 at 78 minutes against a worst case of 12 elsewhere.
# A third tier would be an invented number.
TIME_TIERS = {
    "quick":   60,
    "normal": 240,
}

# The headroom each tier assignment must clear, and why it differs by rule.
#
#   SCALES (3x)  memory grows with cells, regions or genes, and the fixture is
#                ~11k cells while a real cohort is several times that. Three is
#                the factor the assignments below are checked against.
#   FIXED (1.5x) memory is set by the cisTarget databases being read (32.8 GB
#                of ctx feathers, 12.9 GB of dem), not by the experiment, so a
#                bigger dataset does not move it. R09 and R10 only.
#
# R09 at 1.5x already asks for 256 GB. Applying the scaling factor there would
# ask for half a terabyte to guard against growth that cannot happen.
MEM_HEADROOM = {"scales": 3.0, "fixed": 1.5}

# How each rule's memory grows, and the tiers it gets. UNMEASURED rules are
# marked; R19 and R20 did not exist when the run above was made.
RULE_TIERS = {
    #                       mem      time      growth
    "R01_seurat_export":   ("16g",  "quick",  "scales"),
    "R02_build_anndata":   ("4g",   "quick",  "scales"),
    "R03_create_cistopic": ("32g",  "quick",  "scales"),
    "R04_topic_modeling":  ("96g",  "normal", "scales"),
    "R05_region_sets":     ("128g", "quick",  "scales"),
    "R06_prepare_gex_acc": ("32g",  "quick",  "scales"),
    "R07_genome_annot":    ("4g",   "quick",  "scales"),
    "R08_search_space":    ("16g",  "quick",  "scales"),
    "R09_cistarget":       ("256g", "quick",  "fixed"),
    "R10_dem":             ("192g", "quick",  "fixed"),
    "R11_prepare_menr":    ("32g",  "quick",  "scales"),
    "R12_tf_to_gene":      ("32g",  "quick",  "scales"),
    "R13_region_to_gene":  ("32g",  "quick",  "scales"),
    "R14_egrn_direct":     ("128g", "quick",  "scales"),
    "R15_egrn_extended":   ("128g", "quick",  "scales"),
    "R16_aucell_direct":   ("48g",  "quick",  "scales"),
    "R17_aucell_extended": ("48g",  "quick",  "scales"),
    "R18_scplus_mudata":   ("16g",  "quick",  "scales"),
    # MEASURED 2026-09-11, the first run that included them. Both were guessed
    # by analogy to R18 in I5; the guess held for R19 (3978 MB, 4.0x of 16g) and
    # was two rungs and a whole time tier too generous for R20, which reserved
    # 240 minutes for a 95-second rule. Retiered to what the data says.
    #
    # They are close enough to each other to share tiers for a real reason: both
    # read the same `scplusmdata.h5mu` and neither holds much beyond it -- R20's
    # figures are drawn and freed one at a time.
    "R19_postprocess_tsv": ("16g",  "quick",  "scales"),
    "R20_visualize":       ("16g",  "quick",  "scales"),
    # R21 reads TSVs with csv.reader and base64s a handful of PNGs. It
    # holds one figure at a time, so its peak follows the LARGEST FIGURE
    # rather than the dataset -- and the largest here is ~12 Mpx. 4g is
    # generous for that and does not scale with cells.
    "R21_report":          ("4g",   "quick",  "fixed"),
}


def mem(rule_name):
    """Memory in MB for a rule, by name. Unknown rule = hard error, not a default.

    A rule that falls through to the profile's `default-resources` gets 8000 MB,
    which is under the measured peak of twelve of the eighteen steps. Silence
    there would read as "sized" and mean "8 GB".
    """
    return MEM_TIERS[_rule_tiers(rule_name)[0]]


def rt(rule_name):
    """Wall-clock minutes for a rule, by name."""
    return TIME_TIERS[_rule_tiers(rule_name)[1]]


def _rule_tiers(rule_name):
    if rule_name not in RULE_TIERS:
        raise KeyError(
            f"rule {rule_name!r} has no resource tier. Add it to RULE_TIERS in "
            f"rules/common.smk with a measurement in tests/measured_resources.tsv, "
            f"or tests/test_resources.py will fail. Known: {', '.join(sorted(RULE_TIERS))}"
        )
    m, t, _ = RULE_TIERS[rule_name]
    for table, key in ((MEM_TIERS, m), (TIME_TIERS, t)):
        if key not in table:
            raise KeyError(f"rule {rule_name!r} names unknown tier {key!r}")
    return RULE_TIERS[rule_name]
