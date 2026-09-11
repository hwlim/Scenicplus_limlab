# Development Log

20260512: Initial test
  - /Volumes/parklabngs/Lim.0.Analysis/20240521-multiome-kidney/Seurat/merged_macs/ENP/ScenicPlus
  - Failed in 07_run_scenicplus.log
    ```bash
    26 [Tue May 12 22:19:35 2026]
    27 localrule download_genome_annotations:
    28     output: /data/parklabngs/Lim.0.Analysis/20240521-multiome-kidney/Seurat/merged_macs/ENP/ScenicPlus/results/scplus_out/genome_annotation.tsv, /data/parklabngs/Lim.0.Analysis/20240521-multiome-kidney/Seurat/merged_macs/ENP/ScenicPlus/results/scplus_out/chromsizes.tsv
    29     jobid: 8
    30     reason: Missing output files: /data/parklabngs/Lim.0.Analysis/20240521-multiome-kidney/Seurat/merged_macs/ENP/ScenicPlus/results/scplus_out/genome_annotation.tsv, /data/parklabngs/Lim.0.Analysis/20240521-multiome-kidney/Seurat/merged_macs/ENP/ScenicPlus/results/scplus_out/chromsizes.tsv
    31     resources: tmpdir=/scratch/limc8h
    32 
    33 /data/limlab/Resource/conda_env/scenicplus_limlab/lib/python3.11/site-packages/pybiomart/dataset.py:269: DtypeWarning: Columns (0) have mixed types. Specify dtype option on import or set low_memory=False.
    34   result = pd.read_csv(StringIO(response.text), sep='\t')
    35 2026-05-12 22:20:09,890 Download gene annotation INFO     Using genome: GRCm39
    36 Could not find Id on https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=genome&term=GRCm39
    37 Returning gene annotation without subestting for assembled chromosomesand converting to UCSC style. Please make sure that the chromosome namesin the returned object match with the chromosome names in the scplus_obj.Chromosome sizes will not be returned
    38 2026-05-12 22:20:10,163 SCENIC+      INFO     Chrosomome sizes was not found, please provide this information manually.
    39 2026-05-12 22:20:10,164 SCENIC+      INFO     Saving genome annotation to: /data/parklabngs/Lim.0.Analysis/20240521-multiome-kidney/Seurat/merged_macs/ENP/ScenicPlus/results/scplus_out/genome_annotation.tsv
    40 Waiting at most 5 seconds for missing files.
    41 MissingOutputException in rule download_genome_annotations in file /data/parklabngs/Lim.0.Analysis/20240521-multiome-kidney/Seurat/merged_macs/ENP/ScenicPlus/results/scplus_pipeline/Snakemake/workflow/Snakefile, line 221:
    42 Job 8  completed successfully, but some output files are missing. Missing files after 5 seconds. This might be due to filesystem latency. If that is the case, consider to increase the wait time with --latency-wait:
    43 /data/parklabngs/Lim.0.Analysis/20240521-multiome-kidney/Seurat/merged_macs/ENP/ScenicPlus/results/scplus_out/chromsizes.tsv
    44 Removing output files of failed job download_genome_annotations since they might be corrupted:
    45 /data/parklabngs/Lim.0.Analysis/20240521-multiome-kidney/Seurat/merged_macs/ENP/ScenicPlus/results/scplus_out/genome_annotation.tsv
    46 Shutting down, this might take some time.
    47 Exiting because a job execution failed. Look above for error message
    48 Complete log: .snakemake/log/2026-05-12T221935.788913.snakemake.log
    49 WorkflowError:
    50 At least one job did not complete successfully.
    ```

20260730: Flattened the inner snakemake (`f71405a`)
  - SCENIC+'s `grn_inference`/`prepare_data` steps ran through an opaque child
    snakemake spawned inside step 07. Replaced with 13 native driver stages
    dispatched by `scripts/scenicplus_06_grn_stage.py`, one `scenicplus` CLI
    call each. Trades that snakemake's intra-DAG parallelism for real per-stage
    sentinel / `.cfgsha` / cascade resume. Branch `flatten-inner-snakemake`.
  - Not run at this point. Only the flags were checked, against the generated
    Snakefile and the installed CLI.

20260907: A working environment, and steps 1-3 on real data
  - `environment.yml` never produced a usable env (its own header says it was
    never run-tested). `environment.cchmc.yml` + `install_cchmc.sh` do, in two
    pip phases: `pybedtools==0.9.1` needs `--no-build-isolation` (its sdist
    declares no setuptools) while `loomxpy` needs isolation ON (poetry backend).
    One `pip:` section cannot express both, which is why the installer exists.
  - Steps 1-3 run locally (WSL) on a PBMC-400 multiome fixture: 1149 cells,
    12210 genes, 61490 peaks, 25 cell types. Step 4 blocked by the sandbox --
    Ray cannot create its plasma socket. Not a code problem.
  - Also found by running step 3: it reported `n_regions` as the CELL count.

20260908: End to end on CCHMC, twice
  - **All 20 steps completed**, human/hg38, PBMC-400 fixture. Then a SECOND run
    from scratch on the fixed code, stopping once at step 7 to build the genome
    files, supplied through `input.genome_annotation` / `input.chromsizes`.
    That second run is the one that matters: the first still used a pre-fix
    step 3 plus a `bc_transform_func` workaround.
  - Getting there took EIGHT failures that stopped a run, plus a ninth found
    only by opening an output. Every one was two things that had to agree while
    the error named neither:

    | | stopped at | the disagreement | fix |
    |---|---|---|---|
    | 1 | install "succeeded", no `bin/scenicplus` | user site vs env site-packages | `c0c672e` |
    | 2 | pip took an sdist needing rust | login-node glibc vs compute-node glibc | `b6a6e05` |
    | 3 | OOM with 128 GB requested | `rusage[mem]` reserves, `-M` enforces | `61e184b` |
    | 4 | step 6 `CXXABI_1.3.15` | a wheel's libstdc++ vs conda's | `caf82cd` |
    | 5 | step 6 "no cells in both assays" | cisTopic's tagged barcodes vs the RNA AnnData's | `087ee71` |
    | 6 | step 7 exited 0 with no chromsizes | what it wrote vs what step 8 needed | `3f9b6c5` |
    | 7 | step 8 pandas `KeyError` | Ensembl vs UCSC chromosome names | `3e07d07` |
    | 8 | step 20 `'ggplot' has no attribute 'savefig'` | matplotlib's API vs plotnine's | `187bdf3` |
    | 9 | nothing -- a 301-megapixel PNG | `plot_rss` multiplies the figsize it is given | `9511ac7` |

  - #6 and #7 are one cause. `download_genome_annotations` derives chromsizes
    from NCBI E-utilities `db=genome`, a RETIRED database that answers HTTP 200
    with `<Count>0</Count>` for every term -- reproduced from two unrelated
    networks, so not a firewall. The same swallowed branch also does the UCSC
    chromosome-name conversion and the assembled-molecule subsetting. Reported
    upstream as aertslab/scenicplus#640 with a verified fix (`db=assembly`
    returns the same UID the code already needs); commented on the open #611.
    #476 is the identical failure, closed in 2024 on a per-user workaround with
    no diagnosis. No PR has ever touched the function.
  - THIS IS THE SAME FAILURE AS 20260512 ABOVE. The mouse run hit `db=genome`
    too -- and reported `Using genome: GRCm39` while the data is mm10/GRCm38,
    which is the half that would NOT have announced itself.
  - Outputs looked at, not merely produced: `search_space.tsv` (185,070 links,
    56,032 of 61,490 peaks, median TSS distance 53 kb, zero non-standard
    contigs) and the RSS plot, which recovers SPIB/BCL11A in naive B, LEF1 in
    naive CD4 T, KLF4 in classical monocytes, CEBPA/MAFB in intermediate
    monocytes, TBX21 in effector/MAIT -- consistent with the FigR result on the
    same data. 25 cell types over 1149 cells is too fine for RSS; collapse to
    lineage before trusting per-type values.
  - `scenicplus_make_genome_files.R` (`91c3df0`): both files from a pinned
    EnsDb, so the assembly is the one chosen rather than the one Ensembl serves
    today. Prints `chr1 = 195,471,971` for mm10; GRCm39's is 195,154,279. Its
    hg38 output matched the chromsizes from the working run on all 25 shared
    chromosomes. RUNBOOK section 2b covers the mouse setup.

20260909: Plan for a Snakemake workflow (`SnakemakePlan.md`, branch `test`)
  - Written now because there is finally a trusted end-to-end result to gate
    each increment against. Not a rewrite: the step scripts and
    `scenicplus_06_grn_stage.py` stay; only `run_step`'s sentinel/`.cfgsha`/
    cascade is replaced by Snakemake's DAG.
  - Three reasons. The LSF launcher bsubs the whole driver as ONE job sized for
    the heaviest stage, so step 20's plots hold 16 cores and 128 GB for hours.
    `run_step` hashes config but NOT code, so a step whose script changed still
    reads fresh -- which is why no workspace re-ran step 3 after `087ee71`. And
    `input.assembly` is dead config that looks live.
  - Open, not decided: whether step NUMBERS survive as the interface
    (`--from 7`) once rule names exist.

20260909: `input.reduction` -- every figure so far was drawn on a layout
          nobody chose
  - Found from the mouse run: a multi-sample INTEGRATED object, and the figures
    showed an unintegrated layout. The cause was not a wrong reduction being
    picked. It was that none was.
  - `scenicplus_02_build_anndata.py` decided which reduction became `X_umap` --
    the key every plotting call reads -- by name prefix: anything starting with
    "umap". Seurat produces `wnn.umap`, `rna.umap`, `atac.umap`,
    `umap.harmony`. None of them start with "umap". So the test matched nothing
    on any object this pipeline has ever been given, the fallback fired, and a
    fresh PCA-UMAP of the RNA matrix with no batch correction took the key. The
    human PBMC run's own log says so in one line that reads like an
    observation rather than a warning:

    ```
    [build_anndata] Imported reduction 'wnn.umap' -> obsm['X_wnn.umap']
    [build_anndata] No UMAP in Seurat object - computing one for plots.
    ```

  - **Scope: figures only, and that was checked rather than assumed.** No stage
    between 3 and 18 reads an embedding (grep for obsm/umap across the cisTopic,
    region-set and postprocess scripts returns nothing), and
    `scenicplus.is_multiome: true` sends `prepare_GEX_ACC` down the
    barcode-pairing branch, so the metacell path -- the one place SCENIC+ groups
    cells by anything but a label -- never runs. Topics come from fragment
    counts; DARs and RSS from `celltype_column`. No eRegulon, importance or AUC
    value changes.
  - Fixed as one named config key, `input.reduction`, read by three steps:
    step 01 validates it against `Reductions(obj)` and fails in the first
    minute if the object has no such reduction; step 02 maps every reduction to
    `X_<name>` and assigns the chosen one to `X_umap` LAST, so a reduction
    literally named "umap" cannot outrank it; step 20 reads the layout from
    step 01's `embedding_<name>.tsv` rather than hoping it survived nine
    intermediate files, and puts its name in every figure title.
  - Step 20 also stops trusting propagation for a second reason: the eRegulon
    object is concatenated from the AUC modalities alone and inherits no
    embedding, so the curated layout could never have reached those panels
    however step 02 chose it. Reading the TSV means a finished run can be
    redrawn with `--only 20`, without recomputing any of the GRN.
  - Gated on real data, both directions: on the PBMC export, `wnn.umap` gives an
    `X_umap` byte-equal to the source TSV, a name the object lacks exits 1 and
    writes nothing, and the unset case still works and now says what it did.
    Step 20's two failure paths were made to fail on purpose -- a missing file,
    and barcodes tagged `___pbmc` on one side only, which is the step-6 mismatch
    of 2026-09-08 in a new place. Confirmed the added `.cfgsha` key does NOT
    disturb an existing workspace: for a config without `input.reduction` the
    hash is identical before and after, and changes the moment it is set.
  - Confirmed on the cluster the same day: a named reduction produces the
    figure it names.
  - Branches consolidated: `test` merged into `main` and
    `flatten-inner-snakemake`, then deleted. All three had identical content;
    the only obstacle was one commit duplicated by an earlier cross-branch
    cherry-pick. Work continues on `flatten-inner-snakemake`. Entries above
    that name `test` describe where things happened at the time.
  - `quickstart.md`: the path a lab member follows, install through submission,
    with the CCHMC HPC and a personal machine as two named routes. The lab's
    pre-downloaded cisTarget databases and genome files under
    `/data/limlab/Resource/Scenicplus_db/` are recorded there, so nobody
    re-downloads 45.7 GB. Written in plain ASCII deliberately: the first draft
    used typographic dashes and section signs, and they did not survive being
    copied between editors. RUNBOOK keeps the reasoning and stops being the
    tutorial.
  - Both documents record how runs are ACTUALLY submitted rather than how the
    scripts could be used: `scenicplus_run_lsf.sh` called from a
    per-workspace `run.sh`, from
    a node of the class the env was built for. The node class is not a
    preference -- `run_lsf.sh` runs the preflight on the SUBMITTING host before
    anything is queued, so the shell you type in is the one that gets checked,
    and it is also the environment LSF carries into the job. The `_cchmc`
    launcher, which loads modules inside the job, has never been the path in
    use; the runbook now says so rather than implying either is equally
    travelled.

20260909: Snakemake workflow, increment I0 (`SnakemakePlan.md`)
  - Skeleton and config contract: `Snakefile`, `rules/common.smk`,
    `schemas/config.schema.yaml`, `scripts/scenicplus.run.sh`, plus
    `tests/test_config_schema.py` and `tests/dryrun.sh`. No step rules, so
    `rule all` targets nothing and SAYS SO at onstart rather than reporting a
    silent success. The bash driver is untouched and remains the way to run
    anything real.
  - Both gates were made to fail before being believed. Weakening the schema in
    two separate places turned exactly the matching case red; commenting out
    the Snakefile's `validate()` turned the three schema-dependent checks red
    and left the species and runner checks green.
  - Deviation from the plan's layout, recorded there: NO `Template/config.yml`.
    The workflow reads the same `config/config.yaml` the bash driver reads, so
    one workspace runs under either and nothing drifts between two templates.
  - The schema's key list was derived by walking the config, not typed out, and
    cross-checked against what the step scripts actually read.
  - Two things the build settled that the plan had not: the species table
    refuses what the schema allows, because guessing an EnsDb or BSgenome name
    for an unsupported species is inventing the one fact that fails silently;
    and the runner passes no `--configfile`, because snakemake EXTENDS the
    `configfile:` directive rather than replacing it, so a key missing from the
    second file silently keeps the first one's value.

20260909: the cluster executor, which I5 turns out to depend on
  - snakemake 8 removed `--cluster` and moved submission into executor plugins;
    this environment shipped none, so it could run the workflow locally and
    nowhere else. The version is not ours to choose: `scenicplus 1.0a2` pins
    `snakemake==8.5.5` AND all four interface packages with `==`, so anything
    wanting the executor interface at 9.0+ breaks the scenicplus install rather
    than upgrading it. That is why `snakemake` cannot be dropped from the
    environment even though the flattening removed every call to it -- the
    dependency is declarative, and the package still ships
    `scenicplus/snakemake/Snakefile` for a CLI command we no longer use.
  - Newest usable: `cluster-generic` **1.0.8** (1.0.9 wants interface >=9.0.0)
    and `lsf` **0.2.0** (0.2.1+ want >=9.0.0; 0.3.x also want snakemake >=9).
  - Installed cluster-generic 1.0.8 locally: `--dry-run` reported every
    dependency satisfied and one package to add; the real install moved
    nothing, all five pins and scenicplus unchanged, `scenicplus_check.sh`
    still green.
  - `profiles/lsf/` and `tests/cluster_smoke.sh` landed early, ahead of I5's
    resource work, because a profile nobody can test is not worth writing. The
    profile is adapted from scRNA_LimLab_Snake's, which has run on this
    cluster, keeping its two expensive comments: `-M` beside `rusage[mem]`, and
    why a status command is not optional.
  - **Local mode cannot check the thing that matters.** It substitutes a shell
    for bsub, so a failing job returns non-zero synchronously and is caught
    without `lsf-status.sh` ever being consulted. Whether a job LSF KILLS is
    detected -- as opposed to waited on forever -- is a cluster-only result, and
    so is whether bsub receives the threads a rule asked for.
  - **Both came back green the same day**, from `cluster_smoke.sh --lsf` on
    CCHMC: `b_bigger was allocated 4 slots`, and `a failing job was detected
    and reported`. Threads reach `bsub -n` and LSF honours them, so the
    sixteen-fold oversubscription from the sibling repo cannot happen here by
    construction; and `lsf-status.sh` reads real `bjobs` output correctly,
    since bsub returns immediately and nothing else could have noticed that job
    die. I5 is unblocked, and what is left of it is the part none of this
    touched: what each of the twenty steps should actually ask for.

20260909: Snakemake workflow, increment I1 -- steps 1-5 as rules
  - `rules/prepare.smk`: R01-R05, calling the same scripts with the same flags.
    `rule all` targets the region sets.
  - **Gated against the bash driver's own outputs**, PBMC-400 fixture. All 16
    export files byte-identical, every embedding included, and `rna.h5ad`
    byte-identical. `summary.txt` differs by the two provenance lines added
    earlier today. `cistopic_obj.pkl` differs by 19,533 bytes = 17 x 1149
    cells, and 17 is the length of `___scenicplus_run`: the driver's file
    predates the `tag_cells=False` fix, so the WORKFLOW's file is the correct
    one. Steps 4-5 need the cluster; Ray does not run in this sandbox.
  - Rerun triggers verified on real data: changing `input.reduction` re-ran R01
    (*params have changed*), R02 followed (*input files updated by another
    job*), nothing else moved.
  - **The plan was wrong about the `code` trigger, and it mattered.** Snakemake
    covers a rule's own text, NOT a script the rule shells out to -- measured:
    edit the script, get "Nothing to be done". Declaring the script as an
    `input:` closes it. Left in `params:`, this workflow would have reproduced
    the driver's exact blind spot, the one that let a fixed step 3 never re-run,
    while claiming to have removed it. That was one of the three reasons for
    the whole exercise.
  - **`:q` on an empty param renders as NOTHING, not `''`**, so
    `--celltype_scope {params.scope:q} --reduction ...` reached optparse as
    `--celltype_scope --reduction` and the flag swallowed the next one. R01 died
    on its first run with optparse blaming the wrong option. `opt_arg()` omits
    the flag instead.
  - `snakemake --list` prints each rule's docstring after its name, so the
    runner's step-number lookup handed a whole docstring to `--forcerun`. Found
    by the check asserting `-f 1` RESOLVES; asserting only that `-f 7` refuses
    would have stayed green with the lookup broken.

20260909: Snakemake workflow, increment I2 -- R07 checks, it does not generate
  - `rules/genome.smk` + `scripts/scenicplus_genome_prepare.py`. The plan's
    generating rule is SUPERSEDED by a measurement: the SCENIC+ env has NONE of
    the R stack `scenicplus_make_genome_files.R` needs -- both EnsDb packages,
    both BSgenome packages, ensembldb and AnnotationFilter all absent. A
    generating rule could not run there, and adding two genomes to a pin set
    that took a week to settle is not this increment's trade. Generation stays
    a tool, run once per assembly in the scRNA_LimLab_Snake env, which is also
    the right shape: the files depend only on species and assembly, exactly
    like the 45.7 GB of cisTarget databases nobody expects the pipeline to
    build. Decision 1's shared reference directory needs no new config key --
    `input.genome_annotation` / `input.chromsizes` already point into one.
  - **Four checks, which is the part nothing ever did.** Shape (the chromsizes
    header; the seven annotation columns get_search_space reads by name);
    ASSEMBLY (chromosome 1's length IS the assembly, so a GRCm39 chromsizes
    under `assembly: mm10` is refused and told it is GRCm39 -- the failure that
    reached the 2026-05 kidney run and announced nothing); naming, across all
    three files the search space joins INCLUDING this run's peaks; and overlap,
    because agreeing on style is not agreeing, which is all the old diagnostic
    could check. Copies only after all four pass.
  - `QC/assembly.json`: species, assembly, chr1 length and what it matches,
    naming style, chromosome and transcript counts, which peak chromosomes have
    no annotation, sha256 of both sources. Nothing recorded any of that, which
    is how the wrong assembly stayed invisible for four months.
  - Consequence accepted: R07 depends on R01, because the peak chromosomes come
    from the export. The driver's step 7 had neither the dependency nor a check.
  - Gated on the real hg38 pair (chr1 = 248,956,422, UCSC throughout, 23/23
    peak chromosomes annotated) plus `tests/genome_checks.sh`, which makes every
    check fire on purpose, including that a refused pair leaves nothing behind.
    **Two cases first "passed" because the HARNESS was wrong**: `${2:-chr}`
    substitutes on an empty argument as well as a missing one, so the Ensembl
    fixture came out UCSC and two refusals tested nothing. Use `${2-chr}`.

20260909: Snakemake workflow, increment I3 -- the GRN DAG, 18 rules
  - `rules/grn.smk`: R06 and R08-R18, each one
    `scenicplus_06_grn_stage.py --stage X` exactly as the driver calls it. R07
    is I2's checking rule. `--list` prints all eighteen in step order and the
    full DAG builds.
  - **Dependencies taken from each stage's ARGUMENTS, not the driver's
    ordering**, which is the entire gain. Four pairs are independent and now run
    together: cistarget with dem, tf_to_gene with region_to_gene, the two eGRN
    rules, the two AUCell rules. Flattening the inner snakemake gave that
    parallelism up (`CLAUDE.md` says so); declaring real inputs recovers it
    without bringing the inner snakemake back.
  - **One edge the driver never declared.** R12 reads `tf_names.txt`, which R11
    writes; sequential execution satisfied it by accident. `tests/dryrun.sh` now
    checks the graph's SHAPE -- 18 rules, that edge, R08-after-R07,
    R13-after-R08, and the four independent pairs -- as PROPERTIES rather than a
    snapshot, since a snapshot breaks on every legitimate change and teaches
    people to re-bless it.
  - The genome pair is now a parse-time requirement rather than a warning: R08
    cannot build a search space without it.
  - **What cannot be gated locally:** R09 and R10 read the 45.7 GB of cisTarget
    databases, so I3's real test is a cluster run. Established here: the DAG, the
    parameter slices, and that each command is the driver's.
  - Worth remembering from proving the DAG check could fail: the FIRST attempt
    tested nothing. Removing a line by its first match hit the PRODUCING rule's
    output, not the consuming rule's input, so the graph failed to build for an
    unrelated reason and the check never ran. The second attempt asserted the
    anchor was unique first. That is the replace-first-match trap, and it cost a
    minute here rather than a cluster run.

20260910: why two runs of the same data disagreed -- TWO causes, both measured
  - Comparing a bash-driver run against a Snakemake run of the SAME data (the
    full PBMC fixture) showed outputs matching through R11 and R13, then
    diverging: 86 -> 74 direct eRegulons, 104 -> 96 extended.
  - **Cause 1: BLAS thread count.** `tf_to_gene_adj.tsv` differed only in `rho`,
    by 1e-16, on 3-9 rows of 1,428,119; `importance` and `regulation` were
    bit-identical across four hosts, so the SEEDED gradient boosting is
    deterministic. `rho` comes from a correlation, i.e. a matrix product, whose
    reduction order follows thread count, which defaults to the host's core
    count (48 vs 64 here). PROVEN both ways: two runs on the SAME host are
    byte-identical, and pinning the thread variables made two DIFFERENT hosts
    byte-identical.
  - **Cause 2, and the bigger one: PYTHONHASHSEED.** Pinning threads did NOT
    make the eGRN step host-independent. Python randomises string hashing per
    PROCESS, and SCENIC+ does `list(set(...))` over names
    (`utils.py:394,404,405`), so that order changes on every invocation -- same
    host, same input. Pandas sorts are stable, so a permuted input permutes the
    ties and a top-N cut over ties then keeps different rows. This is what the
    cistrome comparison saw as "same entries, PERMUTED order".
  - **A 1e-16 difference on 9 rows cannot move 12 eRegulons.** An earlier
    reading of this blamed the floating point and was wrong; the ordering is the
    plausible mechanism and the arithmetic is a side issue.
  - Fix: `shell_prefix()` in `rules/common.smk` puts `PYTHONHASHSEED=0` and the
    four `*_NUM_THREADS` variables into the GLOBAL shell prefix, so every rule
    inherits them. Not per-rule: pinning applied to some rules and not others
    yields a run that LOOKS reproducible and is not. Pinned to
    `resources.n_cpu`, not the rule's `threads`, because the bash driver runs
    every step in one job under one thread count and the two are only comparable
    if they agree.
  - For the bash driver, exporting the same variables in the workspace `run.sh`
    is sufficient: one bsub job, and LSF carries the submission environment in.
    `quickstart.md`'s runner skeleton now does it.
  - Gated behaviourally, not by grep: `tests/dryrun.sh` includes the real
    `common.smk`, calls the real `shell_prefix()`, and reads the environment
    from inside a running rule -- checking the variables are present, that the
    thread value FOLLOWS n_cpu (config says 7, env must say 7), and that two
    interpreters agree on a string's hash. Removing the seed turns it red.
  - **Freezing the seed buys reproducibility, not correctness.** Which
    eRegulons survive depends on the order an unordered set happened to iterate
    in. Report an eRegulon list with a stability caveat, and prefer several runs
    under different seeds to see which members are always present.
  - Still open: two eGRN runs on ONE host, which is the only comparison that
    separates the two causes for that step. Also unexplained, and NOT ordering:
    the cistrome motif category counts differ (363 vs 364 direct, 425 vs 423
    extended). A permutation cannot change a count.

20260911: the two drivers agree -- the comparison is closed
  - Both drivers re-run on the same data WITH the pinning
    (`PYTHONHASHSEED=0` + the `*_NUM_THREADS` set). Result: every text output
    matches by md5, remaining numeric differences under 1e-9.
  - **Figures closed the same day**, the last artifact class outstanding: the
    PNGs match by md5 AND were opened and read side by side. The hash says the
    same data was drawn; the reading says it is the same picture, not the same
    bytes by coincidence.
  - **The PDFs do NOT match, and that is not a defect.** Each figure is written
    twice (`scenicplus_08_visualize.py:53-62`) and matplotlib stamps
    `/CreationDate` into a PDF, so two runs differ in those bytes with identical
    content. Measured on matplotlib 3.6.3: same code, two runs, PNG md5
    identical, PDF md5 different, `/CreationDate` present. `pin_env()` does not
    address it because a timestamp is not variation in the DATA. Compare PNGs.
  - So the 86 -> 74 eRegulon gap was never the scheduling layer. It was the
    environment, and of the two causes the ORDERING one mattered more than the
    arithmetic: `list(set(...))` over names (`utils.py:394,404,405`) reordered
    per process, stable pandas sorts turned that into different tie-breaks, and
    a top-N cut kept different rows. The 1e-16 `rho` difference that started the
    investigation was the smaller half.
  - **I3 is validated on the cluster.** Eighteen rules across seven hosts
    reproduce one bsub'd job walking twenty steps in a line. That also confirms
    the DAG, including the edge the driver never declared (R12 needs R11's
    `tf_names.txt`) -- a wrong dependency would have shown as a different
    answer, not merely a different order.
  - Worth keeping: reproducibility was the PRECONDITION for this comparison, not
    a nicety. Without the pinning there was no way to tell a scheduling defect
    from environment noise, and the first attempt to explain the gap blamed the
    arithmetic and was wrong.
  - Residual under 1e-9 on the binary outputs is not zero. It is far below any
    threshold this pipeline acts on, but it is recorded rather than rounded away.

20260911: Snakemake workflow, increment I5 -- per-rule resources, from measurement
  - All twenty rules now name their own memory and wall-clock tier. The evidence
    is `tests/measured_resources.tsv`: LSF accounting for all 18 step rules from
    the 2026-09-09/10 run, derived out of `logs/lsf/*.out` rather than typed in.
    `tests/test_resources.py` re-derives every assignment from it.
  - **The single reservation was wrong in BOTH directions at once.** `cistarget`
    peaked at 152566 MB against the 128000 MB the job reserved -- 1.19x over --
    while nineteen other steps sat well inside it. 15 of 20 rules now reserve
    less than that uniform figure, 2 reserve more, and 12 ask for 32 GB or under.
  - **Why nobody noticed: `-M` is enforced PER PROCESS.** LSF's `Max Memory` is
    the whole process tree, so a 28-process job totalling 149 GB never trips a
    128 GB per-process ceiling. Exit 0, nothing logged, number still wrong. Same
    failure class as a silently-defaulted config key.
  - **Two axes, because the data says they do not correlate.** R04 runs 78
    minutes in 21.7 GB; R09 finishes in 11 minutes and wants 149 GB. One ladder
    makes every long rule buy memory or every large one buy hours. Only two time
    tiers exist, because only one rule is slow and a third would be invented.
  - **Headroom 3x where memory scales with the experiment, 1.5x for R09 and R10**,
    whose memory is set by the 32.8 GB + 12.9 GB of cisTarget feathers they read
    and does not move with cohort size. At 3x, R09 would ask for half a terabyte
    to guard against growth that cannot happen.
  - **The plan's own prediction about step 4 was wrong**, and the table is left
    standing so that stays visible: it forecast "wants few cores, lots of RAM",
    and topic modeling wants 21.7 GB -- five rules want more -- while being the
    only slow step in the workflow.
  - **The arithmetic does not support a savings claim.** Reserved MB-hours fall
    only 20% (289956 -> 232909, against 78358 used) because R04 holds two thirds
    of the runtime and its tier barely moved. The case is correctness and
    schedulability. Saying otherwise would not survive the numbers.
  - Gated locally and NOT yet run on a cluster: six source checks plus a new
    `dryrun.sh` 7b that opens each job snakemake actually RESOLVED -- a
    declaration can read correctly and still not take effect. All four
    mutation-tested (delete a rule's `resources:`; point R15 at R14's tier; drop
    R09 one rung; each turned exactly the intended check red).
  - **Two things to settle before the first cluster run.** R09 asks for 256000
    MB -- confirm a node in the queue has it, or it pends instead of failing.
    And R19/R20 carry UNMEASURED tiers: they did not exist when the run above was
    made, so they are sized by analogy to R18 and must be replaced with figures.
  - One check found its own bug while being mutation-tested: the first 7b regex
    spanned from a rule header to the NEXT job's `resources:` line, so a rule
    with none reported its neighbour's numbers under its own name -- a true
    failure with a false explanation. Split into blocks first.

20260911: Snakemake workflow, increment I6 -- report.html, and rule all targets it
  - One self-contained HTML page per run: run + genome + eRegulon tables +
    figures + LSF accounting + logs + the config as read. `rule all` targets it,
    so an ordinary run is not finished until the run is READABLE. That closes
    the failure this pipeline already had once -- a green run whose 301-megapixel
    RSS figure nobody could open, unnoticed for a day.
  - **R21_report is a rule, not a step.** 20 steps plus one; the bash driver has
    no equivalent. `-f 21` forces a redraw, which is what anyone actually wants
    to force.
  - **Standard library only**, and deliberately so rather than by necessity:
    jinja2 is guaranteed (snakemake requires it) and pandas is a scenicplus
    dep. The reason to use neither is the sibling repo's lesson -- an undeclared
    `.Rmd` param is a hard render failure no dry run catches, surfacing at the
    END of a cluster run. A function returning a string cannot fail that way.
  - **Figures embed as data URIs under a budget, and LINK above it**, saying
    which and why. One figure here is 12 Mpx; an unbudgeted embed makes a page
    no browser opens, which is the same failure class it exists to surface.
  - **The gate aims at the MISSING half.** `tests/report_render.sh` renders a
    deliberately PARTIAL workspace -- two of four tables absent, a figure family
    missing, an empty log, an assembly MISMATCH -- and asserts what the page
    says is wrong. A report that silently omits the section whose data never
    arrived looks finished and is not.
  - **Two defects found by writing the gate, both mine.** The LSF field regex
    was anchored `^` but compiled WITHOUT `re.M`, so every accounting number
    read 0 while the unanchored host regex kept working: a compute table of
    zeroes beside correct node names. And the gate's own missing-table assertion
    passed while the report skipped the table, because the summary row still
    carried the filename -- a check that could not report the problem it was
    written for. Both fixed, three mutations now turn the gate red.
  - **An existing test's assumption was falsified.** `dryrun.sh` asserted `-f 21`
    refuses, commenting that 21 "stays out of range however far the increments
    get". I6 added R21. The out-of-range number is now DERIVED from `--list`.
  - Accepted and now said out loud in `onerror`: snakemake does not build a
    target whose inputs failed, so a partial run leaves NO report. Its absence
    after a failure would otherwise read as a second problem.
  - **Not cluster-run.** No report has been produced from a real `5.analysis/`.

20260911: the rerun finished cleanly -- and it did NOT validate I5
  - Reported by the operator as a clean run, and I recorded it as I5's cluster
    validation. **That was wrong, and the epilogues say so.** R19 and R20 both
    report `Total Requested Memory: 128000.00 MB`, which is the PRE-I5 uniform
    tier verbatim; under I5 they request 16000 and 32000. So this run was
    launched from a clone that did not have I5. The per-rule tiers remain
    BUILT, NOT VALIDATED, and RUNBOOK's status row is corrected to say so.
  - Worth keeping as a shape: "the run was clean" and "the run exercised the
    change" are different claims, and only the second one needs evidence FROM
    INSIDE the run. The reservation an epilogue reports is that evidence, and
    it is the first thing to read after any resource change.
  - **It DID supply the two measurements I5 shipped without**, and those are
    valid regardless: what a rule uses is independent of what it reserved. R19
    3978 MB / 46 s, R20 4132 MB / 95 s. Both now in
    `tests/measured_resources.tsv`, so `test_resources.py` checks them instead
    of naming them unmeasured.
  - **One guess held, one did not.** R19's analogy to R18 was right (4.0x of
    16g). R20's was two memory rungs and a whole time tier too generous -- it
    reserved 240 minutes for a rule that runs in 95 seconds. Retiered
    32g/normal -> 16g/quick. That is the argument for measuring rather than
    reasoning about shape: the wrong guess was the one that FELT better
    justified, because R20 draws figures and figures sound expensive.
  - Node and slots are recorded as `unknown` for both: the operator supplied the
    resource block without the host line. Nothing reads those columns, and
    inventing a node name would put an unmeasured fact in a file whose whole
    purpose is measurement.
  - R21_report still needs the same treatment after the first run including it.

Status: end-to-end on human/hg38 small-scale PBMC, and on mouse (reported
2026-09-09; artifacts not inspected here). The PARAMETERS have never been
examined: the topic-count sweep, the DAR thresholds and the search-space width
are as shipped. A green run says the plumbing holds, not that the numbers mean
anything -- and every figure produced before today was drawn on a layout nobody
chose, which is the same lesson arriving through the output rather than through
a crash.
