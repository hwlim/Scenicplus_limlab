# Plan: a Snakemake workflow for SCENIC+

Written 2026-09-09, after the flattened bash driver completed all 20 steps twice
on CCHMC (§7 of `RUNBOOK.md`). This is a plan, not a change. Nothing here is
built.

Modelled on `scRNA_LimLab_Snake`, which is the same lab's answer to the same
shape of problem: `Snakefile` + `rules/*.smk` + `schemas/config.schema.yaml` +
`Template/config.yml` + `profiles/lsf/` + a `*.run.sh` runner.

---

## What changed after this was written

Recorded here rather than edited into the body, so the plan still reads as it
was reasoned and the deltas are visible.

**Mouse has run** (2026-09-09). The gating baseline is no longer human-only, so
I2's mm10 check can be diffed against a real run rather than against the
generator's own output. What it also produced is the defect below.

**`input.reduction` exists, and three steps read it** (`57ecceb`). It names the
Seurat reduction that becomes `X_umap`, the layout every figure is drawn on.
Consequences for this plan, all small but none automatic:

- The schema must carry it, as a plain string with no `default:`, like every
  other key. An unknown reduction is caught by R01 against the object, not by
  the schema, which cannot know what an .rds contains.
- R01 and R02 take it as a param. Changing it must re-run them, which
  Snakemake's `params` trigger gives for free and the `.cfgsha` list gave only
  because it was edited by hand in three places.
- **R20 gains a second input**: the embedding TSV under `1.export/`, not just
  the MuData from `4.grn/`. The eRegulon object is concatenated from the AUC
  modalities and inherits no embedding, so the layout cannot arrive through the
  DAG's main spine. This is the one place where a rule reaches back across
  stages, and it is deliberate: it is also what makes redrawing a finished run
  cheap.

**A per-workspace `run.sh` is how runs are actually launched**, and the
preflight runs on the submitting host before anything is queued. That is
evidence for Decision 5 rather than a change to it; see the note there.

**`quickstart.md` now exists** as the followed path, with RUNBOOK as the
reference behind it. I8 retires the bash driver, so it must rewrite both.

**`input.biomart_host` goes with the download.** It is read in exactly one
place, the step 7 argv, and that call is skipped whenever `genome_annotation`
and `chromsizes` are supplied. Since supplying them is the normal path, and the
download it configures cannot produce chromsizes for anyone and reports the
wrong assembly for mouse, the key only matters on a path that cannot finish.
I2 replaces that path with `rule genome_files`, and `biomart_host` should be
deleted in the same increment rather than left as a second `input.assembly`.

**One template, not two.** The layout below borrowed `Template/config.yml` from
`scRNA_LimLab_Snake`. Taking it would give this repo two templates during the
whole transition, one per driver, with nothing holding them in step. The
workflow instead reads the same `config/config.yaml` the bash driver reads, so
one workspace runs under either and the schema validates both. That is also what
makes I1's gate meaningful: the same config, two drivers, compare the outputs.

---

## What this is NOT

**Not a rewrite of the science.** The eight step scripts stay as they are, and
`scenicplus_06_grn_stage.py` in particular stays: it is the single place that
knows each of the 13 flattened CLI stages' flags, and it took a run through all
20 steps to get right. Rules will call it with `--stage X`, exactly as the bash
driver does.

**Not a re-flattening.** The inner-snakemake removal (`f71405a`) stands. What
gets replaced is the *scheduling* layer — `run_step`'s sentinel + `.cfgsha` +
cascade logic in `scripts/scenicplus_run_pipeline.sh` — with Snakemake's DAG.

---

## Why bother: three things the bash driver cannot do

### 1. Every step gets one job's worth of resources

`scripts/scenicplus_run_lsf.sh:85-95` bsubs the **whole driver as a single
job** — one `bsub` wrapping `/bin/bash -c "$DRIVER"`. `CLAUDE.md:60` says that
job is "sized for the heaviest GRN stages", which is the only thing it can be.
So step 20 (a few plots) holds the same 16 cores and 128 GB as `cistarget`, for
hours.

The steps are not remotely alike:

| step | shape |
|---|---|
| 4 topic_modeling | memory ∝ nnz × 120 B × *concurrent models*; wants few cores, lots of RAM |
| 9/10 cistarget, dem | I/O bound reading 32.8 GB + 12.9 GB of feathers |
| 12/13 tf_to_gene, region_to_gene | GBM, genuinely CPU-parallel |
| 1–3, 19, 20 | minutes, one core |

Per-rule `resources` plus the LSF profile right-sizes each. On a shared cluster
that is the difference between one 72-hour reservation and a dozen small ones.

### 2. Rerun triggers that are not hand-rolled

`run_step` compares a sentinel's `.cfgsha` sidecar against a hash of the config
keys that step reads, and cascades forward once anything runs. That is a
reimplementation of Snakemake's `params` trigger, and it has a gap the real one
does not: **it hashes config, not code**. A step whose *script* changed is
still "fresh". That bit us this session — `scenicplus_03_create_cistopic.py`
gained `tag_cells=False` and no workspace re-ran it.

Snakemake's `code` trigger covers that, and its `input`/`software-env` triggers
cover two more.

> **Carry over the hard-won warning.** `scRNA_LimLab_Snake/CLAUDE.md` documents
> it and it will apply here identically: **never delete `.snakemake/`** from a
> workspace you intend to rerun. Four of the five triggers are *differential* —
> they compare against `.snakemake/metadata/`. Delete it and only mtime
> survives, and a config edit silently does nothing. Recovery is `--forceall`,
> never `--touch`.

### 3. A config that rejects what it does not understand

`input.assembly` is dead config today: nothing reads it, and setting `mm10`
there does nothing while looking like it does something. A JSON-Schema with
`additionalProperties: false` — the `schemas/config.schema.yaml` +
`snakemake.utils.validate()` pattern — makes an unknown or misspelled key fail
at DAG build, before a job is submitted.

---

## Layout

Marked against what exists as of 2026-09-09: [x] built, [ ] pending, [-]
superseded.

```
[x] Snakefile                  orchestration only: configfile, validate(), the
                               species resolution, includes, rule all,
                               onstart/onsuccess/onerror
[x] rules/common.smk           helpers, NO rules: STAGES + stage_path(),
                               resource tiers, the species table, log_path()
[x] rules/prepare.smk          R01-R05
[x] rules/genome.smk           R07: CHECKS a supplied pair, does not generate
[x] rules/grn.smk              R06, R08-R18 (via scenicplus_06_grn_stage.py)
[x] rules/report.smk           R19-R20. The report.html rule is I6.
[x] schemas/config.schema.yaml the contract
[-] Template/config.yml        SUPERSEDED: one template, not two. The workflow
                               reads the same config/config.yaml the bash
                               driver reads. See "What changed", above.
[x] profiles/lsf/              config.yaml + lsf-status.sh, adapted from
                               scRNA_LimLab_Snake's. Early, so the executor
                               could be tested; the RESOURCE numbers are I5.
[x] scripts/                   unchanged: the 8 step scripts, the helper, the
                               genome-file generator, scenicplus_check.sh
[x] scripts/scenicplus.run.sh  the runner (-n, -j, -p, -f N, --lsf)
[ ] scripts/scenicplus.init.sh scaffold a workspace. NOT renamed: today's
                               scenicplus_init.sh serves both drivers and
                               writes the one config they share. Renaming it
                               is I8's job, with the driver it belongs to.
[x] tests/                     test_config_schema.py, dryrun.sh,
                               cluster_smoke.sh + cluster_smoke/Snakefile
[ ] docs/TODO.md, REFACTOR.md  deferred. Development.md and this file carry
                               the record while the workflow is one increment
                               old; splitting them now would be filing
                               cabinets for four documents.
```

### Stage taxonomy

`scRNA_LimLab_Snake` numbers directories by *pipeline stage*, never by
modality. SCENIC+ is linear and single-sample, so it needs fewer:

```
0.input/       the Seurat .rds (or a symlink to it)
1.export/      seurat_export/  (mtx, barcodes, regions, cell_metadata)
2.anndata/     rna.h5ad
3.cistopic/    cistopic_obj.pkl, cistopic_obj_with_topics.pkl, region_sets/
4.grn/         ACC_GEX.h5mu, search_space.tsv, ctx/dem results, adjacencies,
               eRegulons, AUCell, scplusmdata.h5mu
5.analysis/    tsv/  plots/
QC/            model-selection plots, cell/region counts per stage
logs/          one per rule
provenance/<ts>_<status>/
```

Steps 3–5 collapse into one stage because they are one object being refined,
which is how `3.samples/` works over there.

---

## The genome files, per species and assembly

This is the part that most needs to stop being manual, and it is where the plan
differs most from the current driver.

### What is wrong today

`scenicplus prepare_data download_genome_annotations` cannot produce chromsizes
at all — it queries NCBI's retired `db=genome` (upstream #640). Worse, for
mouse it reports whatever assembly Ensembl serves *today*, which is **GRCm39**,
against mm10/GRCm38 data — and that failure is silent. `RUNBOOK.md` §2b has the
detail. The current answer is a hard stop at step 7 and two files built by hand.

### What the rule should do

**Genome files are reference data, not run outputs.** They depend only on
species + assembly, so they belong beside the cisTarget databases and are shared
across analyses:

```yaml
reference:
  dir: /data/limlab/Resource/scenicplus        # shared, built once per assembly
input:
  species: mouse                                # ONE spelling; see below
  assembly: mm10                                # no longer dead config
```

```
rule genome_files:
    output:
        annotation = "{ref}/{assembly}/genome_annotation.tsv",
        chromsizes = "{ref}/{assembly}/chromsizes.tsv",
        record     = "{ref}/{assembly}/assembly.json",
    params:   species, assembly, ensdb package + version
    script:   scripts/scenicplus_make_genome_files.R
```

`scenicplus_make_genome_files.R` already exists and is verified: it derives both
files from a **pinned EnsDb**, so the assembly is the one chosen rather than the
one Ensembl ships. Its human output matched the chromsizes from the working PBMC
run on all 25 shared chromosomes.

Three additions the rule needs:

1. **`assembly.json`, a record.** Species, assembly, EnsDb package **and its
   version**, BSgenome package, chr1 length, and the generation date. Nothing
   today records which assembly an annotation came from — which is exactly how
   the GRCm39 problem stays invisible. With the record, a provenance bundle can
   answer "which assembly was this run on" from a file rather than from memory.
   `chr1 = 195,471,971` is mm10; `195,154,279` is GRCm39. One number settles it.

2. **An override that is still checked.** If `input.genome_annotation` /
   `input.chromsizes` are set, the rule is bypassed and the files are plain
   inputs — but the naming and shape checks from
   `scenicplus_06_grn_stage.py` (`_chrom_style`, `_check_chromsizes_shape`) move
   into a small validation rule so a supplied pair is checked the same way a
   generated one is.

3. **A cross-check against the data.** A rule that compares the annotation's
   chromosome names and lengths against the ATAC peaks' seqnames from step 1,
   and fails when they disagree. Today that disagreement surfaces at step 8 as
   a pandas `KeyError` about missing columns. The check is cheap and belongs
   before the expensive stages, not after them.

### One species spelling, derived not repeated

`input.species` today feeds three CLI calls that want **two different
vocabularies**: step 7 wants `mmusculus`, steps 9–10 want `mus_musculus`. It
works only because `load_motif_annotations` consults the species solely when no
annotation file is given (`pycistarget/utils.py:98`), and one is always set.

A table in `rules/common.smk` should own this, keyed by one validated enum:

| `species` | ensembl | pycistarget | assembly | EnsDb | BSgenome | motif .tbl |
|---|---|---|---|---|---|---|
| `human` | hsapiens | homo_sapiens | hg38 | EnsDb.Hsapiens.v86 | ...UCSC.hg38 | hgnc |
| `mouse` | mmusculus | mus_musculus | mm10 | EnsDb.Mmusculus.v79 | ...UCSC.mm10 | mgi |

The user writes `species: mouse`. Everything else is looked up. That removes a
whole class of "the two spellings disagreed" — the failure mode this project has
hit eight times in other guises.

---

## Increments

Each leaves the pipeline runnable, and each is gated against the **completed
PBMC-400 workspace**, which is the great advantage of doing this now rather than
first: there is a trusted end-to-end result to diff against.

| | scope | gate |
|---|---|---|
| **I0** | skeleton: Snakefile, common.smk, schema, Template config, runner. No rules yet beyond `all`. | schema accepts the template, rejects a typo'd key; `-n` parses |
| **I1** | `rules/prepare.smk` — steps 1–5 | rerun in a copy of the PBMC workspace: outputs byte-identical to the bash driver's |
| **I2** | `rules/genome.smk` + `assembly.json` + the peak/annotation cross-check | build hg38, diff against the files the working run used; build mm10, confirm chr1 = 195,471,971 |
| **I3** | `rules/grn.smk` — steps 6, 8–18 | full run from `3.cistopic/`; `search_space.tsv` and `scplusmdata.h5mu` match |
| **I4** | `rules/report.smk` — steps 19–20 | figures render; the RSS PNG is ~12 Mpx, not 301 |
| **I5** | per-rule resources + LSF profile | a real cluster run; compare wall-clock and peak RSS per step against the single-job baseline |
| **I6** | `report.html`, and `rule all` switched to it | renders for a complete run; a run that fails mid-DAG still leaves usable logs + bundle |
| **I7** | provenance bundle (config.used, git, logs, `lsf_jobs.tsv`, `assembly.json`) | bundle from a failed run contains the failing log |
| **I8** | retire `scenicplus_run_pipeline.sh` | the runner is the only entry point; `quickstart.md` and RUNBOOK both rewritten |

### I0 is built and gated, 2026-09-09

`Snakefile`, `rules/common.smk`, `schemas/config.schema.yaml`,
`scripts/scenicplus.run.sh`, `tests/test_config_schema.py`, `tests/dryrun.sh`.
No step rules yet, so `rule all` targets nothing and says so out loud at
`onstart` rather than reporting a silent success. The bash driver is untouched
and remains the way to run anything real.

    python tests/test_config_schema.py        # 15 refusals, 6 shapes that stay legal
    tests/dryrun.sh                           # parses, and refuses 5 bad configs

Both gates were themselves made to fail before being believed. Weakening the
schema in two separate places turned exactly the matching case red; commenting
out the Snakefile's `validate()` call turned the three schema-dependent cases
red and left the species and runner cases green, which is the discrimination
that makes a green run mean something.

Two things the build settled that the plan had not:

- **The species table refuses what the schema allows.** The schema permits the
  four species the SCENIC+ CLI names, because one config serves both drivers
  and the bash driver passes the value straight through. `rules/common.smk`
  carries reference rows for the two the lab runs, and stops for the others
  rather than guessing an EnsDb or BSgenome package name. A wrong genome
  reference is the failure that does not announce itself.
- **The runner passes no `--configfile`.** Snakemake EXTENDS the `configfile:`
  directive with a command-line one rather than replacing it, so the two get
  merged and a key missing from the second silently keeps the first one's
  value. The workspace's config is the config.

### I1 is built and gated, 2026-09-09

`rules/prepare.smk`, R01-R05, calling the same scripts with the same flags.
`rule all` now targets the region sets.

**The gate the plan asked for, run on the real PBMC-400 fixture.** Steps 1-3
rebuilt from the object into the new stage layout and compared against what the
bash driver produced in `scenicplus-pbmc400/`:

| | |
|---|---|
| 16 export files, all 7 embeddings included | byte-identical |
| `rna.h5ad` | byte-identical |
| `summary.txt` | differs: the two provenance lines added earlier today |
| `cistopic_obj.pkl` | differs by 19,533 bytes |

The pickle difference is explained rather than tolerated: 19,533 = 17 bytes x
1149 cells, and 17 is the length of `___scenicplus_run`. The driver's file has
TAGGED cell names because it was written before the `tag_cells=False` fix. The
workflow's file is the corrected one. Steps 4 and 5 could not run here (Ray's
plasma socket does not work in this sandbox) and need the cluster.

**The rerun triggers work, which is what the increment was for.** Changing
`input.reduction` re-ran R01 with reason *params have changed since last
execution*, and R02 followed as *input files updated by another job*. Nothing
else moved.

### I2 is built and gated, 2026-09-09 -- but R07 does NOT generate

`rules/genome.smk` + `scripts/scenicplus_genome_prepare.py`. The rule block
sketched above, with `script: scenicplus_make_genome_files.R`, is SUPERSEDED,
and by a measurement rather than a preference: the SCENIC+ environment has none
of the R stack that script needs. `EnsDb.Hsapiens.v86`, `EnsDb.Mmusculus.v79`,
both BSgenome packages, `ensembldb` and `AnnotationFilter` are all absent. A
generating rule could not run there, and adding two genomes to an environment
whose pin set took a week to settle is not a trade this increment should make.

So generation stays what it already was: a tool, run once per assembly, in the
`scRNA_LimLab_Snake` environment. That is also the right shape -- the files
depend only on species and assembly, exactly like the 45.7 GB of cisTarget
databases nobody expects the pipeline to build. Decision 1's shared reference
directory is satisfied by pointing `input.genome_annotation` /
`input.chromsizes` into one, with no new config key.

**What R07 does instead is the part that was never done at all: four checks, in
the order that a failure is cheapest to understand.**

1. **Shape.** Chromsizes needs its tab-separated header; the annotation needs
   the seven columns `get_search_space` reads by name.
2. **Assembly.** Chromosome 1's length IS the assembly, so one number catches
   the failure that otherwise reaches the results in silence. A GRCm39
   chromsizes under `assembly: mm10` is refused and told it is GRCm39.
3. **Naming.** UCSC against Ensembl, across all THREE files the search space
   joins -- including this run's actual peaks.
4. **Overlap.** The annotation's chromosomes against the peaks' chromosomes.
   Agreeing on style is not the same as agreeing, and style agreement is all
   the existing diagnostic could check.

Then it copies, never before, and writes `QC/assembly.json`: species, assembly,
chr1 length and what that length matches, naming style, chromosome and
transcript counts, which peak chromosomes have no annotation, and the sha256 of
both sources. Nothing recorded that before, which is how the GRCm39 problem
stayed invisible for four months.

**The consequence to accept:** R07 now depends on R01, because the peak
chromosomes come from the export. The driver's step 7 had no such dependency --
and no such check.

**Gated.** The real hg38 pair from the validated PBMC run passes and produces
the record (chr1 = 248,956,422, UCSC throughout, 23 of 23 peak chromosomes
annotated). `tests/genome_checks.sh` then makes every check fail on purpose,
including that a refused pair leaves nothing behind. Two of those cases
initially "passed" because the HARNESS was wrong -- `${2:-chr}` substitutes on
an empty argument as well as a missing one, so the Ensembl fixture came out
UCSC. Reading a gate's failures before its successes is what caught it.

### I3 is built, 2026-09-09 -- and it recovers the parallelism flattening cost

`rules/grn.smk`: R06 and R08-R18, each one `scenicplus_06_grn_stage.py --stage X`
exactly as the driver calls it. R07 is I2's checking rule, not this file's.
Eighteen rules now exist, `--list` prints them in step order, and the full DAG
builds.

**The dependencies come from each stage's ARGUMENTS, not from the driver's
ordering**, and that is the whole gain. Verified by reading the graph itself
rather than by assertion:

```
R06 needs R02, R04          R12 needs R06, R11
R07 needs R01               R13 needs R06, R08
R08 needs R06, R07          R14/R15 need R11, R12, R13
R09 needs R05               R16 needs R06, R14
R10 needs R05               R17 needs R06, R15
R11 needs R06, R09, R10     R18 needs R06, R14, R15, R16, R17
```

Four pairs are independent and will run together: cistarget with dem,
tf_to_gene with region_to_gene, the two eGRN rules, the two AUCell rules. The
driver runs all thirteen in a line because flattening the inner snakemake gave
that up; declaring real inputs gets it back without re-introducing anything.

**One edge was not obvious and the driver never declared it.** R12 reads
`tf_names.txt`, which R11 writes. Sequential execution satisfied that by
accident. `tests/dryrun.sh` now checks the graph's SHAPE -- the 18 rules, that
edge, R08-after-R07, R13-after-R08, and the four independent pairs -- as
properties rather than as a snapshot, because a snapshot breaks on every
legitimate change and teaches people to re-bless it. Removing R12's `tf_names`
input turns exactly that check red.

**I3 IS VALIDATED ON THE CLUSTER, 2026-09-11.** Both drivers were run on the
same data with the reproducibility pinning in place, and they agree: every text
output matches by md5, and the remaining numeric differences are under 1e-9.

That closes the question the whole comparison existed to answer. The eRegulon
gap of 86 against 74 was never the scheduling layer. It was the environment:
`PYTHONHASHSEED` randomising the order of `list(set(...))` over names, and BLAS
thread count following the host's core count. Fix both and eighteen rules across
seven hosts reproduce one bsub'd job running twenty steps in a line.

Three things that are now established rather than argued:

- **The flattened stage scripts are driver-agnostic.** The same commands, given
  the same environment, give the same answers whoever schedules them.
- **The DAG is right**, including the edge the driver never declared (R12 needs
  R11's `tf_names.txt`). A wrong dependency would have shown up as a different
  answer, not just a different order.
- **Reproducibility was the precondition, not a nicety.** Without the pinning
  this comparison could not have been made at all, and the 1e-16 arithmetic
  difference that started the investigation turned out to be the smaller half.

Residual under 1e-9 on the binary outputs, with the text outputs byte-identical.
Not zero, so it is worth knowing it exists; far below anything a threshold in
this pipeline acts on.

**The genome pair is now mandatory**, a parse-time refusal rather than I2's
warning, because R08 cannot build a search space without it.

### Two things this plan got wrong, found by building it

**1. Snakemake's `code` trigger does NOT cover an external script.** The plan
says it does, and that was the second of the three reasons for the whole
exercise. Measured: a rule whose `shell:` calls a script, then the script is
edited, and snakemake reports *Nothing to be done*. The trigger covers the
rule's own text, not what the rule shells out to.

The fix is one line per rule: declare the script as an `input:` rather than a
`params:`. Then the same edit reports *updated input files* and the rule
re-runs. Every rule in `prepare.smk` does this, and `script_path()` says why.
Had it stayed in `params`, this workflow would have reproduced the driver's
exact blind spot -- the one that let a fixed step 3 never re-run -- while
claiming to have removed it.

**2. `:q` on an empty param renders as NOTHING, not as `''`.** So
`--celltype_scope {params.scope:q} --reduction ...` becomes
`--celltype_scope --reduction`, and the flag swallows the next one. R01 died
that way on its first run, with optparse blaming `celltype_scope` for a problem
in the value after it. `opt_arg()` omits the flag entirely instead; every option
it is used for defaults to empty in the script, so the two are equivalent.

A third, smaller: `snakemake --list` prints each rule's docstring after its
name, so the runner's step-number lookup has to take the first field. It did
not, and handed a whole docstring to `--forcerun`. Caught by the check that
asserts `-f 1` RESOLVES, which exists because asserting only that `-f 7`
refuses would have passed while the lookup was broken.

### I4 is built, 2026-09-11 -- all twenty steps are now rules

`rules/report.smk`: R19 and R20. Both read only `scplusmdata.h5mu`, so they are
independent and run together; R20 additionally reads step 01's embedding, which
is what lets a finished run be redrawn without recomputing anything.

**Which outputs are declared, and why not all of them.** The plan's rule is to
declare everything including side effects, because a stage that exits 0 without
its outputs is what Snakemake catches for free. Two kinds resist it. Names that
are data-dependent cannot be declared at DAG-build time: `02_umap_eRegulon_<name>`
is one figure per top eRegulon. And `03_rss_per_celltype` is computed inside a
try/except that prints and continues, so a sparse cell type can legitimately
leave it out; declaring it would turn a warning into a failed run. Everything
else is declared, which is what makes the plotnine regression -- the two
heatmap-dotplots that crashed on `.savefig` -- a MissingOutput failure rather
than a run that finishes without them. Declaring those two is evidence-based:
the validated run produced both.

**`tests/output_names.py` guards the rule-versus-script name contract**, which
nothing else does and a dry run cannot: `-n` never executes a script, so it
cannot know what the script would have written. A declared name the script never
writes fails at the END of a run, after every expensive stage has succeeded. The
gate compares two independent sources, the declarations and the write calls, and
was made to fail in both directions -- a typo'd declaration, and an
unconditional figure left undeclared.

Writing it caught me mis-reading the code: a grep suggested the extended
dotplot was written as `04_heatmap_dotplot_extended`, and I nearly reported that
as a defect. The number is chosen by a conditional on the following line, so
`05_` is correct. That is the argument for the gate rather than for reading
carefully.

**Not gated here:** R19 and R20 need a real `scplusmdata.h5mu`, so the figures
have not been rendered by this workflow. What is established is the DAG, the
parameter slices, the commands, and that every declared name is one the scripts
write.

### I5 has a prerequisite nobody has bought yet: an executor plugin

Measured 2026-09-09 against the environment `install_cchmc.sh` builds.

Snakemake's version here is not a choice. `scenicplus` depends on
`snakemake==8.5.5`, so it arrives transitively through pip and cannot be moved
without breaking that pin (`environment.cchmc.yml:49-50`). Snakemake 8 removed
`--cluster` and moved cluster submission into executor PLUGINS, and the
environment ships none:

    $ snakemake --executor cluster-generic --help
    invalid choice: 'cluster-generic' (choose from 'local', 'dryrun', 'touch')

    installed: snakemake, snakemake-interface-{common,executor-plugins,
               report-plugins,storage-plugins}     <- interfaces only

So **as shipped, this environment can run the workflow locally and nowhere
else**, and I5 -- the increment the whole plan exists for, per-rule right-sized
LSF jobs -- cannot start until an executor plugin is added.

**Both candidates have a compatible release, checked against PyPI 2026-09-09.**
The binding constraint is not snakemake's own range but scenicplus's, which
pins the interface packages EXACTLY:

    scenicplus 1.0a2 requires:  snakemake==8.5.5
                                snakemake-interface-common==1.17.1
                                snakemake-interface-executor-plugins==8.2.0
                                snakemake-interface-report-plugins==1.0.0
                                snakemake-interface-storage-plugins==3.1.1

So any plugin wanting `snakemake-interface-executor-plugins >=9.0.0` would
force pip to move a package scenicplus pins with `==`, which breaks the
scenicplus install rather than merely upgrading something.

| plugin | newest usable | why the newer ones are not |
|---|---|---|
| `snakemake-executor-plugin-cluster-generic` | **1.0.8** | 1.0.9 requires interface-executor-plugins >=9.0.0 |
| `snakemake-executor-plugin-lsf` | **0.2.0** | 0.2.1 onward require >=9.0.0; 0.3.x also require snakemake >=9 |

Neither of those two moves a pinned package: both accept
interface-executor-plugins 8.2.0 and interface-common 1.17.1 as installed.

**Prefer cluster-generic 1.0.8 if this is done at all.** It wraps a submit
command, which is exactly what `scenicplus_run_lsf.sh` already does with
`bsub`, so the LSF specifics stay in this repo -- including the `-M` versus
`rusage[mem]` lesson, which a third-party plugin would own instead. The
dedicated LSF plugin's usable version, 0.2.0, is the oldest release of its
line and would be the least-exercised code in the stack.

**Tried, 2026-09-09: cluster-generic 1.0.8 installs cleanly and works.** On the
WSL build of this environment, `pip install --dry-run` reported every
dependency already satisfied and "Would install" exactly one package; the real
install moved nothing, all five pinned versions and scenicplus 1.0a2 were
unchanged afterwards, and `scenicplus_check.sh` still passed. `--executor` then
offers `cluster-generic` alongside local, dryrun and touch.

`profiles/lsf/` and `tests/cluster_smoke.sh` came out of that, ahead of I5's
proper resource work, because a profile nobody can test is not worth writing.
The profile is adapted from `scRNA_LimLab_Snake`'s, which has run on this
cluster, and keeps its two expensive comments: `-M` alongside `rusage[mem]`,
and why a status command is not optional.

**PROVEN ON THE CLUSTER, 2026-09-09.** `tests/cluster_smoke.sh --lsf` passed
every check, including the two local mode cannot make:

    ok   b_bigger was allocated 4 slots
    ok   a failing job was detected and reported (exit 1)

The first says `threads` reaches `bsub -n` and LSF honours it, so the
sixteen-fold oversubscription the sibling repo hit cannot happen here by
construction. The second says `lsf-status.sh` works against real `bjobs`
output: bsub returns immediately, so nothing but the status command could have
noticed that job die, and the alternative to noticing is a workflow that waits
forever.

**So I5 is unblocked.** The executor, the profile, the submit template and the
status probe are all exercised end to end on the target cluster, with a
four-job DAG that needs no data. What remains for I5 is the part this never
touched: what each of the twenty steps should actually ASK for.

Two branches of `lsf-status.sh` remain unexercised, both narrow: the `bhist`
fallback for a job that leaves `bjobs` before the next poll, and the `UNKWN`
transient. Neither can be triggered on demand. They matter only for jobs that
finish inside one poll interval or lose contact with their host, and the code
treats both conservatively.

I0–I4 are mechanical; I5 is where the payoff lands. Stopping after I4 is
coherent (correct, still one big job); after I5 (right-sized jobs, no report);
after I6 (a run someone can read). I7–I8 are the tidy-up.

---

## Traps to carry across, not rediscover

Each of these cost a real run in one repo or the other:

- **`.snakemake/` is not a cache.** Never delete it. `--forceall` to recover,
  never `--touch`.
- **`2>&1 | tee {log}`, with `shell.prefix("set -o pipefail; ")`.** `bsub -o`
  *appends*, so `> {log}` gets truncated by a rerun while the LSF epilogue keeps
  growing. Without `pipefail`, a failing script through `tee` exits 0.
- **Declare `threads:` AND pass it.** Setting `resources.cpus_per_task` while
  the submit command interpolates `-n {threads}` gave 16× oversubscription over
  there. The same class as `rusage[mem]` vs `-M` here.
- **`additionalProperties: false`, and no `default:` in the schema.** Defaults
  in two places drift; the code's `.get()` defaults are the ones that run.
- **A step that exits 0 without its outputs.** Snakemake catches this for free
  (`MissingOutputException`) — which is precisely how upstream #611 saw the
  chromsizes bug that our exit-0 driver could not. Declare every output,
  including the ones currently written as side effects.

---

## Decisions

Settled 2026-09-09. Each was open in the first draft; none is a guess.

### 1. Genome files live in a shared reference directory, named in the config

```yaml
reference:
  dir: /data/limlab/Resource/scenicplus     # built once per assembly
```

They depend only on species + assembly, so they belong beside the cisTarget
databases rather than copied into every workspace — those are 45.7 GB and the
genome files should follow the same rule.

### 2. There will be a `report.html`, after the resource work

`scRNA_LimLab_Snake` has one and it is how anyone actually looks at a run;
SCENIC+ currently emits loose PNGs and TSVs. It lands as **I6**, after per-rule
resources (I5), because a report of a run that took one oversized job says less
than a report of one that was scheduled properly.

### 3. Rule names are `R01_*` … `R20_*` — the step number survives

This keeps the number as the interface while giving Snakemake real rule names.
The `R` prefix is load-bearing, not decoration:

```
rule 07_genome_annot:      SyntaxError: invalid decimal literal
rule R-07_genome_annot:    SyntaxError   (a hyphen is a minus sign)
rule R07_genome_annot:     accepted
```

Snakemake compiles rule names as **Python identifiers**, so a name cannot begin
with a digit and cannot contain a hyphen. Tested, not assumed.

Zero-padding matters too: `R01`…`R20` sort in step order as strings, so
`--list` prints the pipeline in the order it runs. Verified that all three
targeting flags take the name and behave:

```
snakemake --list                              # R01_, R02_, R07_ in order
snakemake --forcerun R02_build_anndata        # cascades to everything downstream
snakemake --until R02_build_anndata           # stops there
```

`--forcerun R07_genome_annot` is therefore the direct replacement for
`--from 7`, and the runner can keep a `-f 7` shorthand that maps a number onto
the matching rule name.

### 4. `rule all` targets the report

`report.html` is the default target, so an ordinary `snakemake` run is not done
until the run is readable. That makes the report non-optional by construction
rather than by discipline — the failure mode it prevents is a green run whose
outputs nobody opened, which this pipeline has already produced once (the
301-megapixel RSS figure was "successful" for a day).

Consequence to accept: **until I6 lands, `rule all` targets
`5.analysis/` — the TSVs and plots from R19-R20.** The switch happens with the
report rule, not before, or every increment up to it fails its own target.

Consequence to design around: a *partial* run must still produce a report.
Snakemake will not build `report.html` if an upstream rule failed, so the
report rule cannot be the only way to see what happened — that is what the
provenance bundle (I7) and the per-rule logs are for. The report is the default
target, not the only artifact.

### 5. The preflight runs in the runner, not as a rule

`scenicplus_check.sh` has no outputs, so as a rule it would either run on every
invocation or need a sentinel that lies about when it last passed. The runner
executes it before snakemake starts — which is also where it belongs, since half
of what it checks (which libstdc++ actually loaded, what pip's config says) is a
property of the *environment the job is running in*, and is worth knowing before
a DAG is built rather than as one node inside it.

Note that this is the one place the current bash driver already got right, for
the wrong reason: the check is in the driver because there was nowhere else to
put it. Here it stays there on purpose.

**Strengthened 2026-09-09 by how the cluster runs actually work.**
`scenicplus_run_lsf.sh` runs the preflight on the SUBMITTING host, before
anything is queued, and LSF then carries that same environment into the job. So
the check is not merely convenient in the runner: the runner is the only place
that can check the environment the work will run in *before* committing a
reservation to the queue. Under Snakemake this matters more, not less, because a
bad environment would otherwise be discovered once per rule, in twenty separate
jobs, each after its own queue wait. Keep it in `scenicplus.run.sh`, ahead of
the snakemake invocation, and keep `SCENICPLUS_SKIP_CHECK` as the deliberate
opt-out.

---

## Future work: make the result ORDER-INVARIANT, not merely reproducible

Parked 2026-09-11. The pinning makes two runs agree; it does not make the answer
independent of an arbitrary order. The difference matters because a pinned run is
reproducible by construction while still resting on whichever order a set
happened to iterate in.

The tie-sensitive line is in the eGRN builder, and it is NOT a top-N cut:

    TF2G_adj_relevant_pos.loc[TF].set_index('target')[order_TFs_to_genes_by]
        .sort_values(ascending=False)

That ranking goes to GSEA as `rnk`, and the enrichment score depends on where
set members sit in it. `sort_values` is stable, so tied importances keep input
order. Two candidate interventions, which are not interchangeable:

1. **A total order.** Break ties by a deterministic key -- the gene or region
   name -- wherever a sort feeds a ranking or a cut. Gives INVARIANCE: the
   answer stops depending on input order at all, and the hash-seed pin stops
   being load-bearing for this class. Does not change which items are best.
2. **Include every tied member at a top-N cut**, even past N. Defensible, and
   the idea that prompted this, but it applies only to the cuts and it changes
   what a module CONTAINS rather than just stabilising it. A science change, to
   be argued on its own merits.

Both live in the SCENIC+ package, so either means carrying a patch or
upstreaming one. **Measure before engineering:** gradient-boosting importances
are continuous, so exact ties should be rare unless there is a pile at zero.

    awk -F'\t' 'NR>1 {c[$3]++} END {
      t=0; for (v in c) if (c[v] > 1) t += c[v]
      printf "rows sharing an importance: %d of %d (%.2f%%)\n", t, NR-1, 100*t/(NR-1) }' \
      tf_to_gene_adj.tsv

Under a percent and ties were never the mechanism, which would leave the
unexplained cistrome motif CATEGORY counts (363 vs 364 direct, 425 vs 423
extended) as the remaining suspect -- no reordering can change a count, and
neither intervention above would address it.

---

## Still genuinely open

- **`QC/` contents.** The stage taxonomy reserves it; what belongs there beyond
  the topic-model selection plots is not decided.
