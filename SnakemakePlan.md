# Plan: a Snakemake workflow for SCENIC+

Written 2026-09-09, after the flattened bash driver completed all 20 steps twice
on CCHMC (§7 of `RUNBOOK.md`). This is a plan, not a change. Nothing here is
built.

Modelled on `scRNA_LimLab_Snake`, which is the same lab's answer to the same
shape of problem: `Snakefile` + `rules/*.smk` + `schemas/config.schema.yaml` +
`Template/config.yml` + `profiles/lsf/` + a `*.run.sh` runner.

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

```
Snakefile                  orchestration only: configfile, validate(), the
                           species/assembly resolution, includes, rule all,
                           onsuccess/onerror -> provenance
rules/common.smk           helpers, NO rules: STAGES + stage_path(), resource
                           tiers, the species table, log_path()
rules/prepare.smk          steps 1-5
rules/genome.smk           the genome files + the assembly record  <-- see below
rules/grn.smk              steps 6, 8-18 (all via scenicplus_06_grn_stage.py)
rules/report.smk           steps 19-20
schemas/config.schema.yaml the contract
Template/config.yml        the shipped template, matching the schema key for key
profiles/lsf/config.yaml   + lsf-status.sh, as scRNA_LimLab_Snake has
scripts/                   unchanged: the 8 step scripts, the helper, the
                           genome-file generator, scenicplus_check.sh
scripts/scenicplus.run.sh  the runner (LSF default, -l local, -n dry run, -j, -p)
scripts/scenicplus.init.sh scaffold a workspace (today's scenicplus_init.sh)
tests/                     schema tests + a dry-run matrix, as over there
docs/TODO.md, REFACTOR.md  open work and the record
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
| **I6** | provenance bundle (config.used, git, logs, `lsf_jobs.tsv`, `assembly.json`) | bundle from a failed run contains the failing log |
| **I7** | retire `scenicplus_run_pipeline.sh` | the runner is the only entry point; RUNBOOK rewritten |

I0–I4 are mechanical; I5 is where the payoff lands. Stopping after I4 is
coherent (correct, still one big job); stopping after I5 is coherent (right-sized
jobs, no provenance).

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

## Open questions

1. **Shared reference dir or per-workspace?** Recommended shared
   (`reference.dir`), since genome files depend only on species+assembly and the
   cisTarget databases already live that way. Per-workspace copies would be 46 GB
   each.
2. **A `report.html`?** `scRNA_LimLab_Snake` has one and it is the main way
   anyone looks at a run. SCENIC+ currently emits loose PNGs and TSVs. Probably
   worth it, but after I5.
3. **Does the workspace keep `--from N` ergonomics?** Snakemake's equivalent is
   `--forcerun <rule>`, which cascades the same way. The runner can map a
   `-f <rule>` flag onto it, but the step *numbers* stop being the interface —
   rule names take over. Worth deciding before the RUNBOOK is rewritten.
4. **Where does `scenicplus_check.sh` run?** As a rule with no outputs it would
   run every invocation. Probably stays in the runner, before snakemake starts.
