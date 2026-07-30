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
Approximate DAG:

```
prepare_GEX_ACC            → ACC_GEX.h5mu
download_genome_annot      → genome_annotation.tsv, chromsizes.tsv
get_search_space           → search_space.tsv
motif_enrichment_cistarget → ctx_results.hdf5, cistromes_*
motif_enrichment_dem       → dem_results.hdf5
TF_to_gene                 → tf_to_gene_adj.tsv
region_to_gene             → region_to_gene_adj.tsv
eGRN (direct/extended)     → eRegulons_direct.tsv, eRegulons_extended.tsv
AUCell (direct/extended)   → AUCell_direct.h5mu, AUCell_extended.h5mu
→ assemble                 → scplusmdata.h5mu
```

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

## The other direction (Option C)

The phrase "the main **snakemake** process" could instead mean going the
*opposite* way: turn the outer bash driver back into a real Snakefile and
`include:` SCENIC+'s rules so it's one flat DAG under one `snakemake`
invocation. Viable, but contradicts CLAUDE.md's current commitment (bash driver
*replacing* snakemake).

## Open decision (pick before writing code)

- **A** — Bash driver, inline ~10 native stages (best resume; matches CLAUDE.md). *Recommended.*
- **B** — Bash driver, 2 coarse stages (minimal churn, coarse resume).
- **C** — Convert outer to Snakemake, include SCENIC+ rules (contradicts current design).

## Files touched when implementing (reference)

- `scripts/scenicplus_run_pipeline.sh` — step definitions (lines ~216–235 are
  the current steps 06/07 to replace).
- `scripts/scenicplus_06_init_inner.sh` — delete (Option A).
- `scripts/scenicplus_06_patch_config.py` — delete or repurpose its
  input/output path map as the CLI arg source (Option A).
- `environment.yml` / `environment_macos.yml` — drop `snakemake` (Option A).
- `CLAUDE.md` / `README.md` — update the "inner snakemake is unavoidable" note.
