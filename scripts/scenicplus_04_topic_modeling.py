#!/usr/bin/env python
"""
Run LDA topic modeling on the CistopicObject and select the best model.
The chosen model is attached to cto.selected_model and the object is
re-pickled for downstream steps.
"""
import argparse
import pickle
from pathlib import Path

import yaml


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--in_pkl", required=True)
    p.add_argument("--out_pkl", required=True)
    p.add_argument("--config", required=True, help="Top-level pipeline config.yaml")
    p.add_argument("--tmp_dir", required=True)
    p.add_argument("--n_cpu", type=int, default=8)
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

    best = evaluate_models(
        models, select_model=None, return_model=True, plot=False
    )
    cto.add_LDA_model(best)

    out = Path(args.out_pkl)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(out.suffix + ".partial")
    with open(tmp, "wb") as fh:
        pickle.dump(cto, fh)
    tmp.replace(out)
    print(f"[topic_modeling] Selected {best.n_topic} topics; wrote {out}")


if __name__ == "__main__":
    main()
