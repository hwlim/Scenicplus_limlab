#!/usr/bin/env python
"""Gate for schemas/config.schema.yaml.

The point of this file is the NEGATIVE cases. A schema that accepts the shipped
template proves only that it is not actively broken; what has to be true is that
it REJECTS the mistakes it exists to catch. Each case below names the mistake it
stands for, and a case that stops failing is a hole in the contract.

Run:  python tests/test_config_schema.py [schema.yaml]

The optional argument points the whole suite at a different schema file, which
is how this file gets checked itself: feed it a deliberately weakened copy and
the matching case must go red. A gate nobody has watched fail is not evidence.
"""
import copy
import sys
from pathlib import Path

import yaml
from snakemake.utils import validate          # the same call the Snakefile makes

ROOT = Path(__file__).resolve().parent.parent
SCHEMA = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else ROOT / "schemas" / "config.schema.yaml"
TEMPLATE = ROOT / "config" / "config.yaml"

FAILURES = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + (f"   [{detail}]" if detail and not ok else ""))
    if not ok:
        FAILURES.append(name)


def accepts(cfg):
    """True if the schema accepts this config."""
    try:
        validate(copy.deepcopy(cfg), str(SCHEMA))
        return True
    except Exception:
        return False


def rejects(name, mutate):
    """Apply `mutate` to a copy of the template and require the schema to refuse."""
    cfg = copy.deepcopy(TEMPLATE_CFG)
    mutate(cfg)
    check(name, not accepts(cfg), "accepted, but should not be")


TEMPLATE_CFG = yaml.safe_load(TEMPLATE.read_text())

print(f"schema:   {SCHEMA}")
print(f"template: {TEMPLATE.relative_to(ROOT)}\n")

# --- the schema's own editing rules -----------------------------------------
raw = yaml.safe_load(SCHEMA.read_text())


def find_defaults(node, path=""):
    hits = []
    if isinstance(node, dict):
        for k, v in node.items():
            if k == "default":
                hits.append(path or "<root>")
            hits += find_defaults(v, f"{path}.{k}" if path else k)
    elif isinstance(node, list):
        for i, v in enumerate(node):
            hits += find_defaults(v, f"{path}[{i}]")
    return hits


defaults = find_defaults(raw)
check("no `default:` anywhere in the schema", not defaults, ", ".join(defaults[:3]))

# --- the positive case ------------------------------------------------------
check("the shipped template validates", accepts(TEMPLATE_CFG))

# --- the negative cases, which are the reason this file exists --------------
print("\n  each of these MUST be refused:")


def set_in(cfg, section, key, value):
    cfg[section][key] = value


rejects("an unknown TOP-LEVEL key (`inputs:` for `input:`)",
        lambda c: c.__setitem__("inputs", {}))

rejects("an unknown key inside input (`celltype_col` for `celltype_column`)",
        lambda c: set_in(c, "input", "celltype_col", "cell_type"))

rejects("an unknown key inside cistopic (`n_topic` for `n_topics`)",
        lambda c: set_in(c, "cistopic", "n_topic", [10]))

rejects("another pipeline's config (wrong Pipeline value)",
        lambda c: c.__setitem__("Pipeline", "scRNA_LimLab_Snake"))

rejects("a missing Pipeline line altogether",
        lambda c: c.pop("Pipeline"))

rejects("a required input dropped (seurat_rds)",
        lambda c: c["input"].pop("seurat_rds"))

rejects("an empty path where one is required",
        lambda c: set_in(c, "input", "seurat_rds", ""))

rejects("a species the pipeline has no vocabulary for",
        lambda c: set_in(c, "input", "species", "mouse"))

rejects("a count given as text (`n_cpu: x`)",
        lambda c: set_in(c, "resources", "n_cpu", "x"))

rejects("zero cores",
        lambda c: set_in(c, "resources", "n_cpu", 0))

rejects("a probability above 1 (dar_adjpval_thr)",
        lambda c: set_in(c, "cistopic", "dar_adjpval_thr", 2))

rejects("an empty topic list, which fits no model at all",
        lambda c: set_in(c, "cistopic", "n_topics", []))

rejects("a duplicated topic count, which fits the same model twice",
        lambda c: set_in(c, "cistopic", "n_topics", [10, 10, 20]))

rejects("a search space with one bound where the CLI takes two",
        lambda c: set_in(c, "scenicplus", "search_space_upstream", "1000"))

rejects("a quantile above 1",
        lambda c: set_in(c, "grn", "quantile_thresholds_region_to_gene", [0.9, 1.5]))

# overlap_top_n is a COUNT OF eREGULONS entering a pairwise matrix, so one is
# not a comparison and a fraction is not a count. `minimum: 2` rather than 1,
# because the figure needs two things to overlap -- the script skips at one and
# the schema should say so rather than leaving it to a runtime message.
rejects("an overlap_top_n of 1, which is not a comparison",
        lambda c: set_in(c, "visualization", "overlap_top_n", 1))
rejects("a fractional overlap_top_n",
        lambda c: set_in(c, "visualization", "overlap_top_n", 12.5))
rejects("a typo'd visualization key",
        lambda c: set_in(c, "visualization", "overlap_topn", 40))

# --- shapes that must STAY legal --------------------------------------------
# Guarding the other direction: a schema that rejects a valid config is a worse
# failure than a loose one, because it blocks a run that would have worked.
print("\n  each of these MUST still be accepted:")


def keeps(name, mutate):
    cfg = copy.deepcopy(TEMPLATE_CFG)
    mutate(cfg)
    check(name, accepts(cfg), "refused, but is valid")


keeps("a search space written as a list, which toks() also accepts",
      lambda c: set_in(c, "scenicplus", "search_space_upstream", [1000, 150000]))

keeps("an empty reduction, meaning `compute a layout`",
      lambda c: set_in(c, "input", "reduction", ""))

keeps("a named reduction",
      lambda c: set_in(c, "input", "reduction", "wnn.umap"))

keeps("empty genome files, meaning `try to fetch them`",
      lambda c: (set_in(c, "input", "genome_annotation", ""),
                 set_in(c, "input", "chromsizes", "")))

keeps("an omitted optional section",
      lambda c: c.pop("visualization"))

keeps("an empty celltype_scope, meaning `every cell`",
      lambda c: set_in(c, "input", "celltype_scope", []))
keeps("a raised overlap_top_n, which is the reason the key exists",
      lambda c: set_in(c, "visualization", "overlap_top_n", 120))

print()
if FAILURES:
    print(f"FAILED: {len(FAILURES)}")
    for f in FAILURES:
        print(f"  - {f}")
    sys.exit(1)
print("all checks passed")
