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
echo "  and the report must say WHICH of the four situations the reader is in:"

# --- 8. the report's Genome panel, driven by the PRODUCER's own record ------
# THE FIXTURE IS THIS SCRIPT. scenicplus_09_report.py compared `input.assembly`
# to `chr1_matches` with `==`. The record says "hg38/GRCh38" -- one length, both
# accepted spellings -- and the config says "hg38", so the equality was false,
# the two `not` branches were false, and a REAL report printed no panel at all:
# no confirmation, no warning, a section that reads as if it were never written.
#
# tests/report_render.sh could not see it because it hand-wrote
# `"chr1_matches": "hg38"`, a value this script has never emitted. Consumer and
# fixture were written from the same imagination and agreed with each other and
# with nothing. So the record below comes from the producer, which is the only
# pairing that can disagree usefully -- the same lesson as tests/record_keys.py,
# one level down: that file anchors the key NAMES, this one the VALUES.
panel() {  # panel <chr1 length> <assembly> -> the rendered Genome section
    mk_chromsizes "$W/cs_p.tsv" "$1"
    "$PY" "$PREP" --annotation "$W/ann.tsv" --chromsizes "$W/cs_p.tsv" \
        --atac_regions "$W/peaks.tsv" --out_annotation "$W/pa.tsv" \
        --out_chromsizes "$W/pc.tsv" --record "$W/rec_p.json" \
        --species hsapiens --assembly "$2" >/dev/null 2>&1
    [[ -s "$W/rec_p.json" ]] || { echo "PRODUCER-WROTE-NO-RECORD"; return 0; }
    "$PY" - "$W/rec_p.json" <<'PY'
import os, sys
sys.path.insert(0, os.path.join(os.environ["SCENICPLUS_PATH"], "scripts"))
import scenicplus_09_report as rep
print(rep.sec_assembly(sys.argv[1]))
PY
}
n_panels() { grep -o 'class="\(ok\|miss\|bad\)"' <<<"$1" | wc -l; }

# The regression itself: SOME panel, in every one of the four situations. This
# is the assertion the shipped bug would have failed and the equality-based
# code cannot satisfy for a recognised assembly.
for spec in "248956422 hg38 confirmed" "248956422 hg19 mismatch" \
            "248956422 '' unset" "123456789 hg19 unrecognised"; do
    eval "set -- $spec"
    P="$(panel "$1" "$2")"
    if [[ "$(n_panels "$P")" == 1 ]]; then
        say ok "$3: exactly one panel renders"
    else
        say FAIL "$3: $(n_panels "$P") panels rendered, expected exactly 1"
    fi
done

# ...and it must be the RIGHT panel. Without this, four empty-but-present divs
# would pass the count above.
P="$(panel 248956422 hg38)"
grep -q "Assembly confirmed" <<<"$P" \
    && say ok "the producer's \"hg38/GRCh38\" confirms a configured \"hg38\"" \
    || say FAIL "a correct pair did not produce the confirmation panel"
# Both halves, because the count above cannot separate them: restore the `==`
# and a CORRECT pair falls through to the mismatch panel, which is one panel
# and the wrong one. Measured -- that is what the negative control produced.
grep -q 'class="bad"' <<<"$P" \
    && say FAIL "a correct pair ALSO raised the mismatch panel" \
    || say ok "...and raises no mismatch panel alongside it"

# NEGATIVE CONTROL for the check above: hg19 is absent from R07's EXPECTED_CHR1,
# so R07 does NOT gate it and the disagreement reaches the report. If this
# renders "confirmed", assembly_agrees() is matching anything.
P="$(panel 248956422 hg19)"
grep -q "Assembly MISMATCH" <<<"$P" \
    && say ok "hg38 data configured as hg19 is reported as a MISMATCH" \
    || say FAIL "an hg19/hg38 disagreement was not reported as a mismatch"

P="$(panel 248956422 '')"
grep -q "Assembly NOT confirmed" <<<"$P" \
    && say ok "an unset input.assembly is reported as unchecked" \
    || say FAIL "an unset input.assembly produced the wrong panel"

# chr1 at a length in no reference table: detected is null, configured is not.
P="$(panel 123456789 hg19)"
grep -q "Assembly unrecognised" <<<"$P" \
    && say ok "a chr1 length in no reference table is reported as unrecognised" \
    || say FAIL "an unrecognised chr1 length produced the wrong panel"

echo
[[ "$FAIL" -ne 0 ]] && { echo "FAILED"; exit 1; }
echo "all checks passed"
