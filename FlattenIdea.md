# FlattenIdea — removing the nested SCENIC+ Snakemake

Notes for revisiting in another Claude Code session / device.

## The current reality (important reframing)

Despite the "avoid nested Snakemake" framing, **the outer pipeline is not
Snakemake anymore.** `scripts/scenicplus_run_pipeline.sh` is a hand-rolled bash
driver with sentinel-based skip logic (`step_is_fresh`, `.cfgsha` sidecars,
cascade re-runs). CLAUDE.md says it "replaces snakemake."

The *only* Snakemake in the whole pipeline is the one SCENIC+ spawns internally
at step 07:

```
run_step 7 run_scenicplus ...
    cd '$SCPLUS_PIPELINE/Snakemake' && snakemake --cores N --rerun-incomplete
```

So there is no nesting today — but there is one opaque `snakemake` child
process the driver can't see into (no per-stage sentinels, no per-stage config
hashing, no partial-resume through *our* logic). That black box is the thing
worth removing. Right now the driver hashes only `resources.n_cpu` for step 07,
so changing any GRN parameter either re-runs the entire inner snakemake or
nothing.

## What the inner Snakemake actually is

`scenicplus init_snakemake` generates a Snakefile that is just a DAG of shell
rules calling two CLI entrypoints — `scenicplus prepare_data …` and
`scenicplus grn_inference …`. Every node in that DAG maps one-to-one to an
`output_data` file already enumerated in `scripts/scenicplus_06_patch_config.py`.
The generated Snakefile is present in the repo at `Snakemake/Snakefile` — the
DAG below was traced from its actual `input`/`output` wiring (13 rules, not the
earlier ~9-rule approximation; earlier notes omitted `prepare_menr` and
mis-attributed the genome/chromsizes edges).

**True dependency structure, by level (what runs in parallel):**

```
LEVEL 0  (roots — fire simultaneously)
├─ prepare_GEX_ACC ──────► combined_GEX_ACC_mudata (ACC_GEX.h5mu)
├─ download_genome_annot ─► genome_annotation.tsv, chromsizes.tsv
├─ motif_cistarget ──────► ctx_result.hdf5     ★ usually the slowest single step
└─ motif_dem ────────────► dem_result.hdf5     (waits on genome_annot ONLY in
                                                 "balance_number_of_promoters" mode)
LEVEL 1
├─ get_search_space   ◄─ GEX_ACC + genome_annot + chromsizes    → search_space.tsv
└─ prepare_menr       ◄─ cistarget + dem + GEX_ACC              → tf_names, cistromes_{direct,extended}
LEVEL 2  (parallel branches)
├─ tf_to_gene         ◄─ GEX_ACC + prepare_menr(tf_names)       → tf_to_gene_adj.tsv
└─ region_to_gene     ◄─ GEX_ACC + search_space                → region_to_gene_adj.tsv
LEVEL 3  (parallel)
├─ eGRN_direct        ◄─ tf_to_gene + region_to_gene + cistromes_direct   → eRegulons_direct.tsv
└─ eGRN_extended      ◄─ tf_to_gene + region_to_gene + cistromes_extended → eRegulons_extended.tsv
LEVEL 4  (parallel)
├─ AUCell_direct      ◄─ eGRN_direct + GEX_ACC                  → AUCell_direct.h5mu
└─ AUCell_extended    ◄─ eGRN_extended + GEX_ACC               → AUCell_extended.h5mu
LEVEL 5
└─ scplus_mudata      ◄─ AUCell_{direct,extended} + eGRN_{direct,extended} + GEX_ACC
                                                                → scplusmdata.h5mu (rule all)
```

Parallelism the inner snakemake exploits:

- **Level 0 is the big win**: cistarget, dem, prepare_GEX_ACC, and
  genome_annotations have no inter-dependencies and run at once. cistarget (★)
  is typically the longest job, so dem + prepare + genome-download hide behind
  it entirely.
- Three more fork points: tf_to_gene ∥ region_to_gene (L2), and the
  direct/extended split at eGRN (L3) and AUCell (L4).
- **Critical path** (bounds wall time regardless of cores):
  `cistarget → prepare_menr → tf_to_gene → eGRN_direct → AUCell_direct →
  scplus_mudata` (6 nodes). Everything else can overlap beside it.

**Gotchas before transcribing shell blocks verbatim:**

1. Likely typo in the generated Snakefile — `get_search_space` calls
   `scenicplus prepare_data search_spance` ("spance", line ~250). Verify whether
   the real CLI subcommand is misspelled too or this is a generator bug.
2. The `dem → genome_annotation` edge is **conditional** (only under
   `balance_number_of_promoters`). Honor it in balanced mode when grouping.

## Recommendation: flatten, don't nest deeper (Option A)

Fits this repo's philosophy (bash driver + sentinels, no Snakemake). Replace
step 07's `snakemake --cores` call with direct `scenicplus` CLI stage calls,
each promoted to its own `run_step`. Delete the `init_snakemake` +
`patch_config` machinery (steps 06–07) and inline the DAG as native steps:

```
06 prepare_gex_acc        S=ACC_GEX.h5mu
07 genome_annotation      S=genome_annotation.tsv   (parallel-eligible)
08 search_space           S=search_space.tsv
09 motif_cistarget        S=ctx_results.hdf5
10 motif_dem              S=dem_results.hdf5
11 tf_to_gene             S=tf_to_gene_adj.tsv
12 region_to_gene         S=region_to_gene_adj.tsv
13 egrn                   S=eRegulons_direct.tsv
14 aucell                 S=AUCell_direct.h5mu
15 assemble_mudata        S=scplusmdata.h5mu
16 postprocess_tsv / 17 visualize  (current 08/09)
```

Why this is the right call:

- **One process model.** Every SCENIC+ stage gets the same sentinel + `.cfgsha`
  + cascade treatment as steps 01–05. Real per-stage incremental resume driven
  by *our* config keys (change one GRN param → only downstream stages re-run).
- **No second scheduler.** Drop `snakemake` from `environment.yml`, plus the
  `init_snakemake` / `scenicplus_06_patch_config.py` / `scenicplus_06_init_inner.sh`
  layer (~200 lines that exist only to generate/rewrite a config for a
  scheduler we're eliminating).
- **Cluster-friendly.** LSF today `bsub`s one fat job sized for the inner
  snakemake's `--cores N`. Flattened, each stage is a discrete step that could
  later be submitted as its own right-sized LSF job — impossible while buried
  inside a child snakemake.

Cost / risk:

- The generated Snakefile is the authoritative source for exact subcommand
  names, flags, and DAG edges (some stages read multiple upstreams). Lift them
  verbatim, don't guess. **Real implementation step:** run
  `scenicplus init_snakemake` once, open `Snakemake/workflow/Snakefile`, and
  transcribe each rule's `shell:` block into a `run_step`.
  (Couldn't read it in the originating session — `scenicplus` wasn't importable
  in that env.)
- Whatever intra-stage parallelism the inner snakemake exploited (e.g.
  cistarget vs dem concurrently) becomes sequential unless added explicitly —
  usually fine since each `scenicplus` stage is already multi-threaded via
  `n_cpu`.

## Lighter alternative (Option B)

Keep step 07 but split into two driver steps instead of one snakemake call:

```
06 scenicplus prepare_data   --config <patched>   → intermediate files
07 scenicplus grn_inference  --config <patched>   → scplusmdata.h5mu
```

Still removes the `snakemake` invocation (both CLIs run a full stage without the
scheduler) with far less churn, but only 2 sentinels instead of ~10, so
incremental resume stays coarse. Pragmatic middle ground.

## Middle ground: strategic grouping (Option A′) — recommended sweet spot

"2 vs ~10 stages" is a false binary. The inner Snakefile's rules are **already
CLI subcommands** (each rule's `shell:` is a `scenicplus prepare_data …` /
`scenicplus grn_inference …` call). So "A's transcription cost" is just *copy
the subcommand + flags out of the Snakefile* — and once you've paid that, **how
many `run_step`s you wrap them into is a free choice.** The real effort cliff is
only B → anything-finer (B needs no Snakefile reading). Beyond that, landing on
a good stage count K is nearly free.

Group along the DAG's level boundaries so no group ever depends on a *later*
group (sentinels stay monotonic, exactly as the driver expects). A clean
5-stage cut:

```
06 prepare     = prepare_GEX_ACC + genome_annot + search_space    S=search_space.tsv
07 motif       = cistarget + dem + prepare_menr                   S=cistromes_direct   ★ protect this
08 adjacency   = tf_to_gene + region_to_gene                      S=region_to_gene_adj.tsv
09 egrn        = eGRN_direct + eGRN_extended                      S=eRegulons_direct.tsv
10 aucell+asm  = AUCell_direct + AUCell_extended + scplus_mudata  S=scplusmdata.h5mu
```

- ~90% of A's resume benefit for ~40% of the wrapper count. The expensive motif
  enrichment (★) sits behind its own `.cfgsha`: change an eGRN/AUCell knob →
  only 09/10 re-run; motif work untouched.
- Want the Level-0 parallelism back? Split 07 into `07a cistarget` / `07b dem`
  and background them; likewise `08a tf_to_gene` / `08b region_to_gene`
  (~7 stages). Recovers the cross-branch overlap B/monolithic-A lose.

### Why B is actually the *weakest* option for re-run cost

Worth spelling out (the doc previously just said "coarse resume"): the nested
snakemake gives free `--rerun-incomplete` rule-level resume today. Option B
replaces it with a single monolithic `grn_inference` CLI that does **not** skip
existing intermediates — so a crash 90% through re-does all of motif enrichment,
and changing any downstream knob recomputes cistarget/dem too. B therefore loses
**both** the cross-branch parallelism **and** the internal resume; only per-stage
`n_cpu` threading survives. A/A′ restore stage-level resume through our driver.

## The other direction (Option C)

The phrase "the main **snakemake** process" could instead mean going the
*opposite* way: turn the outer bash driver back into a real Snakefile and
`include:` SCENIC+'s rules so it's one flat DAG under one `snakemake`
invocation. Viable, but contradicts CLAUDE.md's current commitment (bash driver
*replacing* snakemake).

## Open decision (pick before writing code)

- **A** — Bash driver, inline ~13 native stages (finest resume; matches CLAUDE.md).
- **A′** — Bash driver, ~5 stages grouped along DAG levels (~90% of A's resume
  benefit, ~40% of the work). *Recommended sweet spot.*
- **B** — Bash driver, 2 coarse stages (minimal churn, but forfeits the inner
  snakemake's free rule-level resume — weakest for re-run cost).
- **C** — Convert outer to Snakemake, include SCENIC+ rules (contradicts current design).

## Files touched when implementing (reference)

- `scripts/scenicplus_run_pipeline.sh` — step definitions (lines ~216–235 are
  the current steps 06/07 to replace).
- `scripts/scenicplus_06_init_inner.sh` — delete (Option A).
- `scripts/scenicplus_06_patch_config.py` — delete or repurpose its
  input/output path map as the CLI arg source (Option A).
- `environment.yml` / `environment_macos.yml` — drop `snakemake` (Option A).
- `CLAUDE.md` / `README.md` — update the "inner snakemake is unavoidable" note.
