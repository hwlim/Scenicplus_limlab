#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Gate for scenicplus_09_report.py (I6): does it RENDER, and does it tell the
# truth about what it found.
#
#   tests/report_render.sh
#
# Needs python3 only -- no snakemake, no cluster, no scenicplus. Builds a
# synthetic workspace with real PNG bytes and real TSVs, renders, and greps the
# page.
#
# WHY A SEPARATE GATE FROM dryrun.sh. A dry run never executes a rule, so it
# proves the report rule is WIRED and nothing about whether the script works.
# The sibling repo learned this the expensive way: its report driver shipped an
# undeclared render param, `snakemake -n` was green, and the failure surfaced at
# the END of a cluster run after every expensive rule had already succeeded.
# This is the check that would have caught it.
#
# THE HALF THAT MATTERS IS THE MISSING-DATA HALF. A report that renders a
# complete run is easy. The failure this must catch is a report that silently
# omits the section whose data never arrived, because that page looks finished
# and is not. So the fixtures below are deliberately PARTIAL, and the assertions
# are about what the page SAYS IS MISSING.
# -----------------------------------------------------------------------------
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$HERE")"
GEN="$ROOT/scripts/scenicplus_09_report.py"
[[ -f "$GEN" ]] || { echo "missing $GEN"; exit 1; }

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
FAIL=0
say() { printf '  %-4s %s\n' "$1" "$2"; [[ "$1" == FAIL ]] && FAIL=1; return 0; }

# --- fixture -----------------------------------------------------------------
WS="$WORK/ws"
mkdir -p "$WS/config" "$WS/5.analysis/tsv" "$WS/5.analysis/plots" \
         "$WS/QC" "$WS/logs/lsf"

cat > "$WS/config/config.yaml" <<'YML'
Pipeline: "ScenicPlus"
input:
  species: hsapiens
  assembly: hg38
  celltype_column: cell_type
  reduction: wnn.umap
  seurat_rds: /nowhere/obj.rds
resources:
  n_cpu: 16
YML

# Real PNG bytes: a 1x1 and a deliberately LARGE one, so the embed budget is
# exercised in both directions rather than asserted about.
python3 - "$WS/5.analysis/plots" <<'PY'
import struct, sys, zlib, os
out = sys.argv[1]
def png(path, w, h):
    raw = b"".join(b"\x00" + bytes((x * 7 + y) % 256 for x in range(w * 3))
                   for y in range(h))
    def chunk(tag, data):
        c = tag + data
        return struct.pack(">I", len(data)) + c + struct.pack(">I", zlib.crc32(c))
    ihdr = struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0)
    open(path, "wb").write(
        b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", ihdr)
        + chunk(b"IDAT", zlib.compress(raw, 0)) + chunk(b"IEND", b""))
for stem in ("01_umap_celltype", "04_heatmap_dotplot_direct",
             "06_TF_target_count", "09_region_overlap_direct",
             "99_unexpected_extra"):
    png(os.path.join(out, stem + ".png"), 8, 8)
# ~6 MB uncompressed, comfortably over the 1 MB budget the test passes.
png(os.path.join(out, "08_eGRN_network_top50.png"), 900, 2200)
# A PDF sibling for exactly one figure, so the "PDF" link is exercised AND its
# absence elsewhere is too.
open(os.path.join(out, "01_umap_celltype.pdf"), "wb").write(b"%PDF-1.4\n%%EOF\n")
PY

printf 'TF\tRegion\tGene\timportance\n' > "$WS/5.analysis/tsv/eRegulons_direct.tsv"
for i in $(seq 1 40); do
    printf 'TF%d\tchr1:%d-%d\tGENE%d\t0.%03d\n' "$i" "$((i*100))" "$((i*100+50))" "$i" "$i" \
        >> "$WS/5.analysis/tsv/eRegulons_direct.tsv"
done
printf 'TF\tn_targets\n' > "$WS/5.analysis/tsv/TF_summary.tsv"
printf 'SPIB\t42\nLEF1\t17\n' >> "$WS/5.analysis/tsv/TF_summary.tsv"
# eRegulons_extended.tsv and eRegulons_combined.tsv are DELIBERATELY ABSENT.

# Model-selection figures, as step 04 now writes them. The combined plot plus
# one per-metric panel, so the section's ordering (combined first) is exercised.
python3 - "$WS/QC" <<'PY'
import struct, sys, zlib, os
out = sys.argv[1]
def png(path, w, h):
    raw = b"".join(b"\x00" + bytes((x * 5 + y) % 256 for x in range(w * 3))
                   for y in range(h))
    def chunk(tag, data):
        c = tag + data
        return struct.pack(">I", len(data)) + c + struct.pack(">I", zlib.crc32(c))
    open(path, "wb").write(
        b"\x89PNG\r\n\x1a\n"
        + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0))
        + chunk(b"IDAT", zlib.compress(raw, 0)) + chunk(b"IEND", b""))
png(os.path.join(out, "topic_model_selection.png"), 12, 8)
png(os.path.join(out, "topic_model_selection_mimno_2011_maximize.png"), 12, 8)
PY
open "$WS/QC/topic_model_selection.pdf" 2>/dev/null || printf '%%PDF-1.4\n%%%%EOF\n' > "$WS/QC/topic_model_selection.pdf"

cat > "$WS/QC/assembly.json" <<'JSN'
{"assembly_detected": "hg38", "assembly_configured": "hg19",
 "n_chromosomes": 24, "n_annotation_rows": 1200,
 "n_annotation_chromosomes": 24, "peak_chromosomes": 25,
 "peak_chromosomes_annotated": 24, "peak_chromosomes_unannotated": ["chrM"],
 "source": {"annotation": {"path": "/ref/ann.tsv", "sha256": "abc123def456789a"},
            "chromsizes": {"path": "/ref/chrom.tsv", "sha256": "0f0f0f0f0f0f0f0f"}}}
JSN

cat > "$WS/logs/lsf/R09_cistarget.11.out" <<'LSF'
Job was executed on host(s) <16*bmi-200m5-04>, in queue <normal>, as user <x>
    Run time :                                   682 sec.
    Max Memory :                                 152566 MB
    Max Processes :                              28
    Max Threads :                                1151
LSF
cat > "$WS/logs/lsf/R04_topic_modeling.3.out" <<'LSF'
Job was executed on host(s) <16*bmi-200m5-03>, in queue <normal>, as user <x>
    Run time :                                   4651 sec.
    Max Memory :                                 21745 MB
    Max Processes :                              31
    Max Threads :                                1332
LSF
# A PREVIOUS run's leftovers, then a fresh marker: exactly what `--forcerun R20`
# leaves behind. Reported from a real forcerun, where the report listed the whole
# previous run and put R21_report in the Compute table with the PREVIOUS run's
# numbers -- a row that looks current and is not.
printf 'from an older run\n' > "$WS/logs/R98_ancient.log"
cat > "$WS/logs/lsf/R97_stale.1.out" <<'LSF'
Job was executed on host(s) <8*a-node-from-last-week>, in queue <normal>, as user <x>
    Run time :                                   999 sec.
    Max Memory :                                 4242 MB
    Max Processes :                              3
    Max Threads :                                4
LSF
cat > "$WS/logs/lsf/R21_report.7.out" <<'LSF'
Job was executed on host(s) <bmi-200m5-07>, in queue <normal>, as user <x>
    Run time :                                   6 sec.
    Max Memory :                                 72 MB
    Max Processes :                              9
    Max Threads :                                11
LSF
# Only the two that belong to the previous run. R04's epilogue stays current --
# backdating it too silently removed the data an unrelated check reads, and that
# check went red rather than the scoping one. A fixture change that breaks
# another assertion is a fixture bug, not a finding.
touch -t 202001010000 "$WS/logs/R98_ancient.log" \
                      "$WS/logs/lsf/R21_report.7.out" \
                      "$WS/logs/lsf/R97_stale.1.out"
sleep 0.1
: > "$WS/logs/.run_started"          # this run begins HERE
sleep 0.1
printf 'ran fine\n' > "$WS/logs/R01_seurat_export.log"
: > "$WS/logs/R20_visualize.log"          # empty on purpose -- a REAL alarm
# The report's own log, empty exactly as it is on a real run: `tee` creates it
# when the job starts and this script's output arrives only after the page is
# written. It must NOT be reported as an empty log.
: > "$WS/logs/R21_report.log"

# --- render ------------------------------------------------------------------
OUT="$WS/report.html"
python3 "$GEN" --workspace "$WS" --config "$WS/config/config.yaml" \
    --out "$OUT" --pipeline "$ROOT" --head-rows 10 --max-embed-mb 1 \
    --self-log "$WS/logs/R21_report.log" \
    >"$WORK/stdout" 2>"$WORK/stderr"
rc=$?
[[ $rc -eq 0 ]] && say ok "the generator exits 0" \
                || { say FAIL "the generator exited $rc"; sed -n '1,10p' "$WORK/stderr"; }
[[ -s "$OUT" ]] && say ok "report.html is non-empty" || say FAIL "report.html is empty"

has() { grep -qF -- "$2" "$OUT" && say ok "$1" || say FAIL "$1"; }
hasnt() { grep -qF -- "$2" "$OUT" && say FAIL "$1" || say ok "$1"; }

# --- 1. it renders the things that ARE there ---------------------------------
has "the workspace config is embedded"            "cell_type"
has "the direct eRegulon table is rendered"       "GENE1"
has "the TF summary is rendered"                  "SPIB"
has "a table says how many rows it is showing of" "of 40 rows"
has "LSF accounting is read from the epilogues"   "bmi-200m5-04"
# SCOPED to this run. An older epilogue must not appear, however plausible.
# A stale epilogue from a rule that is NOT R21, so this is independent of the
# report's own-job case below -- otherwise one fix would appear to satisfy both.
grep -q "R97_stale\|a-node-from-last-week" <<<"$(sed -n '/id=compute/,/id=logs/p' "$OUT")" \
  && say FAIL "a PREVIOUS run's epilogue leaked into Compute" \
  || say ok "an older run's epilogue is NOT in Compute"
grep -q "R98_ancient" "$OUT" \
  && say FAIL "a previous run's log is listed as this run's" \
  || say ok "an older run's log is NOT listed"
has "...and the page says how many it left out" "belong to earlier runs"
has "wall clock is humanised, not raw seconds"    "1:17:31"
has "the compute note explains the -M caveat"     "per process"
# R21's own job cannot be in its own Compute table: LSF appends the accounting
# when the job ENDS. Reported from a real run, where R21 appeared under Logs and
# not under Compute and read as a missing job. BOTH directions -- it must be
# absent from the table AND explained, since suppressing the explanation and
# omitting the row look identical to a reader counting rules.
grep -q "R21_report" <<<"$(sed -n '/id=compute/,/id=logs/p' "$OUT" | grep -o '<td[^>]*>R21_report</td>')" \
  && say FAIL "R21 should not be in its own compute table" \
  || say ok "the report's own job is absent from Compute, as it must be"
has "...and the page explains WHY it is absent"   "own epilogue does not exist"

# --- 2. THE HALF THAT MATTERS: it says what is NOT there ---------------------
# NOT just the filename: the section opens with a summary table that lists all
# four names and their row counts, so the name is present even when the
# per-table block was skipped. Mutation-tested -- dropping the per-table
# `missing()` call left a bare filename check GREEN. Assert the missing-div that
# only the per-table branch emits.
grep -qF 'Not in this run:</strong> eRegulons_extended.tsv' "$OUT" \
  && say ok "an absent table gets its OWN missing block, not just a row" \
  || say FAIL "an absent table gets its OWN missing block, not just a row"
has "the summary row still reports it absent"     "absent"
has "...and marked as missing"                    "Not in this run"
has "an absent figure family is reported"         "figure 02_"
# 09_/10_ are the region-overlap pair. The direct one is in the fixture and the
# extended one is not, so both halves of FIGURE_ORDER's handling are exercised
# in the same render: a family that is present and one that is not.
has "the region-overlap figure is placed, not dumped in Other" "Target-region overlap, direct"
has "...with the caption that says what to read"  "one finding reported twice"
has "...and its missing twin is reported"         "figure 10_"
has "an empty log is called out"                  "R20_visualize.log"

# The report's OWN log is empty on every real run: tee creates it at job start
# and this script prints afterwards. Reporting it would be a false alarm every
# single time, and a check that cries wolf is one people learn to skip.
# BOTH directions, because either alone is satisfiable the wrong way --
# suppressing the alarm entirely passes the first, alarming on everything passes
# the second.
alarm="$(grep -o 'Empty log(s):</strong>[^<]*' "$OUT" || true)"
if grep -q 'R20_visualize.log' <<<"$alarm" && ! grep -q 'R21_report.log' <<<"$alarm"; then
    say ok "the report's own log is excluded from the empty-log alarm"
else
    say FAIL "the empty-log alarm is wrong: [$alarm]"
fi
has "...and the page explains why it reads as 0 bytes" "written after this page"
has "...naming the mechanism, not just asserting it"   "cannot describe that run"

# --- 3. the assembly mismatch, which is the one real alarm -------------------
has "the model-selection section renders"         "Model selection"
has "...showing the combined sweep figure"        "topic_model_selection.png"
has "...explaining that rescaling can mislead"    "flatter a weak optimum"
# ONLY the combined figure is shown. The per-metric panels are written and are
# pointed at, not rendered -- five near-identical line plots is more page than a
# QC detail earns. Both halves asserted: shown-once, and named-not-shown.
qc_imgs="$(grep -o '<img src="[^"]*" alt="topic_model_selection[^"]*"' "$OUT" | wc -l)"
[[ "$qc_imgs" == "1" ]] \
  && say ok "exactly ONE model-selection figure is embedded" \
  || say FAIL "expected 1 embedded model-selection figure, got $qc_imgs"
has "...and the per-metric panels are POINTED AT"  "panel(s) alongside this one"
has "...naming the QC folder they are in"          "QC/"
has "...and the all-panels PDF"                    "topic_model_selection_all_pages.pdf"
has "an assembly mismatch is raised loudly"       "Assembly mismatch"
has "...naming both sides"                        "hg19"
has "unannotated peak chromosomes are named"      "chrM"

# --- 4. the embed budget, in BOTH directions ---------------------------------
has "a small PNG is embedded as a data URI"       "data:image/png;base64,"
has "an oversized PNG is linked instead"          "Linked, not embedded"
has "...and the link is a relative path"          "5.analysis/plots/08_eGRN"
has "a PDF sibling is offered when it exists"     "01_umap_celltype.pdf"
has "an unexpected figure is listed, not dropped" "99_unexpected_extra"

# The oversized figure must NOT also be embedded -- that is the whole budget.
if grep -o 'data:image/png;base64,[A-Za-z0-9+/=]*' "$OUT" \
   | awk '{ if (length($0) > 2000000) exit 1 } END { exit 0 }'; then
    say ok "no single embedded image exceeds the budget"
else
    say FAIL "an oversized image was embedded despite the budget"
fi

# --- 5. escaping, because config and paths are untrusted text ----------------
mkdir -p "$WORK/ws2/config"
printf 'Pipeline: "ScenicPlus"\ninput:\n  species: "<script>x</script>"\n' \
    > "$WORK/ws2/config/config.yaml"
python3 "$GEN" --workspace "$WORK/ws2" --config "$WORK/ws2/config/config.yaml" \
    --out "$WORK/ws2/report.html" >/dev/null 2>&1
if grep -qF "<script>x</script>" "$WORK/ws2/report.html"; then
    say FAIL "a config value is escaped before it reaches the page"
else
    grep -qF "&lt;script&gt;" "$WORK/ws2/report.html" \
        && say ok "a config value is escaped before it reaches the page" \
        || say FAIL "the config value did not reach the page at all"
fi

# --- 6. an EMPTY workspace still renders, and says so ------------------------
mkdir -p "$WORK/ws3/config"
printf 'Pipeline: "ScenicPlus"\ninput: {species: hsapiens}\n' \
    > "$WORK/ws3/config/config.yaml"
python3 "$GEN" --workspace "$WORK/ws3" --config "$WORK/ws3/config/config.yaml" \
    --out "$WORK/ws3/report.html" >"$WORK/o3" 2>"$WORK/e3"
if [[ $? -eq 0 && -s "$WORK/ws3/report.html" ]]; then
    # grep -c counts LINES, and several missing-divs land on one line
    # because the sections are joined without newlines. Count OCCURRENCES.
    n="$(grep -o "Not in this run" "$WORK/ws3/report.html" | wc -l)"
    [[ "$n" -ge 4 ]] \
        && say ok "an empty workspace renders and reports every section missing ($n)" \
        || say FAIL "an empty workspace rendered but only flagged $n missing section(s)"
else
    say FAIL "an empty workspace did not render at all"
fi

echo
if [[ $FAIL -eq 0 ]]; then echo "all checks passed"; else echo "FAILED"; fi
exit $FAIL
