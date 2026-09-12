#!/usr/bin/env python
"""
Run LDA topic modeling on the CistopicObject and select the best model.
The chosen model is attached to cto.selected_model and the object is
re-pickled for downstream steps.

ALSO EMITS THE MODEL-SELECTION FIGURES. `evaluate_models` already computes four
metrics across every topic count in the sweep and picks the optimum from them;
until now the pipeline passed `plot=False` and threw all of that away, so the
one number that decides the whole run -- how many topics -- arrived with no
evidence behind it. The sweep is the most expensive step in the workflow and it
was the least inspectable.

WHY THE FIGURES ARE CAPTURED RATHER THAN `save=`d. `evaluate_models(save=...)`
writes ONE MULTI-PAGE PDF and no PNG, and this project's output contract is one
plot per file in both formats (CLAUDE.md). It also builds the figure whether or
not `plot` is true -- `plot=False` merely closes it -- so the figures are taken
off pyplot afterwards and written out individually. `plot=True` is passed for
exactly that reason: to stop the library closing what we are about to save.

Names come from each figure's own axes TITLE, not from the argument order.
`plot_metrics` emits Arun_2010, Cao_Juan_2009, Mimno_2011 and loglikelihood in
an order this file would otherwise have to assume, and an assumption there
mislabels a metric rather than failing.
"""
import argparse
import pickle
import re
from pathlib import Path

import yaml



def _slug(text):
    """A filename-safe stem from a figure title."""
    return re.sub(r"[^A-Za-z0-9]+", "_", text.strip()).strip("_").lower()


def save_open_figures(out_dir, stem):
    """Write every figure pyplot currently holds, one plot per file, pdf+png.

    Returns the paths written. The FIRST figure `evaluate_models` builds is the
    combined rescaled-metric plot and carries no title -- it gets `stem`. The
    per-metric figures each set their own title, which becomes the suffix, so a
    reordering upstream renames a file instead of mislabelling one.
    """
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for num in plt.get_fignums():
        fig = plt.figure(num)
        title = ""
        for ax in fig.get_axes():
            if ax.get_title():
                title = ax.get_title()
                break
        name = stem if not title else f"{stem}_{_slug(title)}"
        for ext in ("pdf", "png"):
            path = out_dir / f"{name}.{ext}"
            fig.savefig(path, bbox_inches="tight",
                        **({"dpi": 200} if ext == "png" else {}))
            written.append(path)
        plt.close(fig)
    return written


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--in_pkl", required=True)
    p.add_argument("--out_pkl", required=True)
    p.add_argument("--config", required=True, help="Top-level pipeline config.yaml")
    p.add_argument("--tmp_dir", required=True)
    p.add_argument("--n_cpu", type=int, default=8)
    p.add_argument("--qc_dir", default=None, metavar="DIR",
                   help="Where to write the model-selection figures. Omitted "
                        "= do not draw them (the sweep still runs).")
    p.add_argument("--plot_metrics", action="store_true",
                   help="Also emit one figure per metric. The combined plot "
                        "RESCALES all four onto one axis, which is convenient "
                        "and can flatter a weak optimum; these are the "
                        "unscaled read.")
    args = p.parse_args()

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)
    ct_cfg = cfg["cistopic"]

    with open(args.in_pkl, "rb") as fh:
        cto = pickle.load(fh)

    from pycisTopic.lda_models import evaluate_models, run_cgs_models

    Path(args.tmp_dir).mkdir(parents=True, exist_ok=True)
    models = run_cgs_models(
        cto,
        n_topics=ct_cfg["n_topics"],
        n_cpu=args.n_cpu,
        n_iter=ct_cfg["n_iter"],
        random_state=ct_cfg["random_state"],
        alpha=ct_cfg["alpha"],
        alpha_by_topic=ct_cfg["alpha_by_topic"],
        eta=ct_cfg["eta"],
        eta_by_topic=ct_cfg["eta_by_topic"],
        save_path=args.tmp_dir,
    )

    # Backend before pyplot: compute nodes have no display, and the import
    # order is what decides whether that matters.
    if args.qc_dir:
        import matplotlib
        matplotlib.use("Agg")

    # plot=True does NOT mean "show it" here -- under Agg `plt.show()` is a
    # no-op. It means "do not close the figure", which is the only way to get
    # it out of the library as something we can save per the output contract.
    # `plt.show()` under Agg warns once per figure that the canvas is
    # non-interactive. It is true and irrelevant -- we never wanted a window --
    # and five copies of it in a rule log is noise where a real warning should
    # stand out.
    import warnings
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*non-interactive.*")
        # `save=` ALSO writes pycisTopic's own multi-page PDF, every panel in
        # one file to page through. Kept alongside the per-figure output rather
        # than instead of it: a multi-page PDF cannot be embedded in the report
        # (a browser will not inline one, and there is no single page to show),
        # so on its own it would be invisible there -- verified by rendering a
        # workspace holding only that file, which reported the figures as never
        # drawn. Different name, or it would collide with the combined figure's
        # own single-page PDF.
        all_pages = (str(Path(args.qc_dir) / "topic_model_selection_all_pages.pdf")
                     if args.qc_dir else None)
        if args.qc_dir:
            Path(args.qc_dir).mkdir(parents=True, exist_ok=True)
        best = evaluate_models(
            models, select_model=None, return_model=True,
            plot=bool(args.qc_dir),
            plot_metrics=bool(args.qc_dir and args.plot_metrics),
            save=all_pages,
        )
    cto.add_LDA_model(best)

    if args.qc_dir:
        figs = save_open_figures(args.qc_dir, "topic_model_selection")
        print(f"[topic_modeling] all panels in one file: {all_pages}")
        print(f"[topic_modeling] model-selection figures: "
              f"{', '.join(str(f) for f in figs) or 'NONE -- pycisTopic drew nothing'}")

    out = Path(args.out_pkl)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(out.suffix + ".partial")
    with open(tmp, "wb") as fh:
        pickle.dump(cto, fh)
    tmp.replace(out)
    print(f"[topic_modeling] Selected {best.n_topic} topics; wrote {out}")


if __name__ == "__main__":
    main()
