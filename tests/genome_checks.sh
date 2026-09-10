#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Gate for scenicplus_genome_prepare.py, the R07 checks.
#
#   tests/genome_checks.sh
#
# Needs SCENICPLUS_PATH and a python with the pipeline's scripts importable.
# Builds its own tiny annotation / chromsizes / peaks files, so it needs no
# reference data and no network.
#
# Every case here is a REFUSAL except the first. That is the point: these checks
# exist to stop a wrong genome reaching the expensive stages, and the only
# evidence that they can is watching each one fire. The assembly case matters
# most -- it is the mm10-versus-GRCm39 failure that reached a real kidney run
# and announced nothing.
# -----------------------------------------------------------------------------
set -uo pipefail

: "${SCENICPLUS_PATH:?set SCENICPLUS_PATH to the repository root}"
PY="${PYTHON:-python}"
PREP="$SCENICPLUS_PATH/scripts/scenicplus_genome_prepare.py"

W="$(mktemp -d)"; trap 'rm -rf "$W"' EXIT
FAIL=0
say() { printf '  %-4s %s\n' "$1" "$2"; [[ "$1" == FAIL ]] && FAIL=1; return 0; }

# --- fixtures ---------------------------------------------------------------
# hg38's chr1 is 248,956,422. Anything else under assembly hg38 is a different
# genome wearing its name.
#
# `${2-chr}`, NOT `${2:-chr}`. The colon form substitutes on an EMPTY argument as
# well as a missing one, so passing "" to get Ensembl-style names silently
# produced UCSC ones instead, and two refusal cases passed while testing
# nothing. The harness was wrong, not the script -- which is why a gate's
# failures are worth reading before its successes.
mk_chromsizes() {  # mk_chromsizes <file> <chr1 length> [prefix]
    local pre="${3-chr}"
    { printf 'Chromosome\tStart\tEnd\n'
      printf '%s1\t0\t%s\n' "$pre" "$2"
      printf '%s2\t0\t242193529\n' "$pre"; } > "$1"
}
mk_annotation() {  # mk_annotation <file> [prefix] [drop-column]
    local pre="${2-chr}" drop="${3:-}"
    if [[ "$drop" == "drop" ]]; then
        printf 'Chromosome\tStart\tEnd\tStrand\tGene\n' > "$1"
        printf '%s1\t1000\t2000\t+\tAAA\n' "$pre" >> "$1"
    else
        printf 'Chromosome\tStart\tEnd\tStrand\tGene\tTranscription_Start_Site\tTranscript_type\n' > "$1"
        printf '%s1\t1000\t2000\t+\tAAA\t1000\tprotein_coding\n' "$pre" >> "$1"
        printf '%s2\t3000\t4000\t-\tBBB\t4000\tprotein_coding\n' "$pre" >> "$1"
    fi
}
mk_peaks() {       # mk_peaks <file> [prefix]
    local pre="${2-chr}"
    printf '%s1:100-200\n%s1:300-400\n%s2:500-600\n' "$pre" "$pre" "$pre" > "$1"
}

run() {  # run <annotation> <chromsizes> <peaks> <assembly> -> exit code, log in $W/out
    "$PY" "$PREP" --annotation "$1" --chromsizes "$2" --atac_regions "$3" \
        --out_annotation "$W/out_ann.tsv" --out_chromsizes "$W/out_chrom.tsv" \
        --record "$W/record.json" --species hsapiens --assembly "$4" \
        >"$W/out" 2>&1
}

echo "checking: $PREP"
echo

# --- 1. a correct pair must pass --------------------------------------------
mk_chromsizes "$W/cs.tsv" 248956422
mk_annotation "$W/ann.tsv"
mk_peaks      "$W/peaks.tsv"
if run "$W/ann.tsv" "$W/cs.tsv" "$W/peaks.tsv" hg38; then
    [[ -s "$W/out_ann.tsv" && -s "$W/record.json" ]] \
        && say ok "a correct pair passes, and both files plus the record are written" \
        || say FAIL "passed but did not write its outputs"
else
    say FAIL "a correct pair was refused"; sed -n '1,8p' "$W/out"
fi

echo
echo "  each of these MUST be refused:"

# --- 2. the assembly fingerprint -------------------------------------------
# GRCm39's chr1 under the name mm10: the exact shape of the silent failure.
mk_chromsizes "$W/cs_39.tsv" 195154279
mk_annotation "$W/ann_m.tsv"; mk_peaks "$W/peaks_m.tsv"
if run "$W/ann_m.tsv" "$W/cs_39.tsv" "$W/peaks_m.tsv" mm10; then
    say FAIL "GRCm39 chromsizes accepted while the config said mm10"
elif grep -q "GRCm39" "$W/out"; then
    say ok "a GRCm39 chr1 under assembly mm10 is refused, and NAMED as GRCm39"
else
    say FAIL "refused, but without naming what the length actually is"
fi

# --- 3. a raw UCSC .chrom.sizes, with no header -----------------------------
printf 'chr1\t248956422\nchr2\t242193529\n' > "$W/cs_raw.tsv"
run "$W/ann.tsv" "$W/cs_raw.tsv" "$W/peaks.tsv" hg38 \
    && say FAIL "a headerless chromsizes was accepted" \
    || say ok "a raw UCSC .chrom.sizes (no header) is refused"

# --- 4. an annotation missing a column get_search_space reads ---------------
mk_annotation "$W/ann_short.tsv" chr drop
if run "$W/ann_short.tsv" "$W/cs.tsv" "$W/peaks.tsv" hg38; then
    say FAIL "an annotation missing two columns was accepted"
elif grep -q "Transcription_Start_Site" "$W/out"; then
    say ok "a missing annotation column is refused, and the column is named"
else
    say FAIL "refused, but without naming the missing column"
fi

# --- 5. naming disagreement across the three files -------------------------
# Ensembl-style annotation against UCSC peaks: what a failed step 7 leaves.
mk_annotation "$W/ann_ens.tsv" ""
if run "$W/ann_ens.tsv" "$W/cs.tsv" "$W/peaks.tsv" hg38; then
    say FAIL "an Ensembl-style annotation was accepted against UCSC peaks"
elif grep -q "naming" "$W/out"; then
    say ok "a chromosome-naming disagreement is refused"
else
    say FAIL "refused, but not for the naming reason"
fi

# --- 6. same style, different contigs --------------------------------------
# The case style agreement cannot catch: both UCSC, nothing in common.
printf 'Chromosome\tStart\tEnd\tStrand\tGene\tTranscription_Start_Site\tTranscript_type\nchrZZ\t1\t2\t+\tX\t1\tprotein_coding\n' > "$W/ann_zz.tsv"
if run "$W/ann_zz.tsv" "$W/cs.tsv" "$W/peaks.tsv" hg38; then
    say FAIL "an annotation sharing no chromosome with the peaks was accepted"
elif grep -q "no chromosome appears in BOTH" "$W/out"; then
    say ok "an annotation sharing no chromosome with the peaks is refused"
else
    say FAIL "refused, but not for the overlap reason"
fi

# --- 7. nothing is placed by a failing run ---------------------------------
# Validate-then-copy, in that order. A half-written pair that looks complete is
# how the driver's version of this went wrong twice.
rm -f "$W/out_ann.tsv" "$W/out_chrom.tsv"
run "$W/ann_ens.tsv" "$W/cs.tsv" "$W/peaks.tsv" hg38 >/dev/null 2>&1
[[ ! -e "$W/out_ann.tsv" && ! -e "$W/out_chrom.tsv" ]] \
    && say ok "a refused pair leaves NOTHING behind" \
    || say FAIL "a refused pair still wrote an output"

echo
[[ "$FAIL" -ne 0 ]] && { echo "FAILED"; exit 1; }
echo "all checks passed"
