#!/usr/bin/env python
"""
Tiny helper used by scenicplus_run_pipeline.sh to slice config.yaml and hash
the slice deterministically. Lets the master script answer "does the relevant
chunk of config.yaml still match what step <N> was last run with?"

Subcommands
-----------
hash <config> <dotted_keys>
    Load <config>, pluck the comma-separated <dotted_keys> (e.g.
    "cistopic.n_topics,cistopic.alpha"), serialize the resulting nested dict to
    canonical JSON (sorted keys, no spaces) and print sha256 of that.

slice <config> <dotted_keys>
    Same selection, but print the canonical JSON. Handy for debugging.

get <config> <dotted_key>
    Print one scalar value (used to read e.g. resources.n_cpu from bash).

getcsv <config> <dotted_key>
    Print a list value as comma-separated. Empty / missing → empty string.
"""
import hashlib
import json
import sys
from pathlib import Path


def _load(path: str) -> dict:
    import yaml
    with open(path) as fh:
        return yaml.safe_load(fh)


def _pluck(cfg: dict, dotted: str):
    cur = cfg
    for part in dotted.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return None
        cur = cur[part]
    return cur


def _slice(cfg: dict, dotted_keys: list[str]) -> dict:
    out: dict = {}
    for k in dotted_keys:
        if not k:
            continue
        cur_in = cfg
        cur_out = out
        parts = k.split(".")
        for part in parts[:-1]:
            cur_in = cur_in.get(part, {}) if isinstance(cur_in, dict) else {}
            cur_out = cur_out.setdefault(part, {})
        leaf = parts[-1]
        if isinstance(cur_in, dict) and leaf in cur_in:
            cur_out[leaf] = cur_in[leaf]
    return out


def _canonical(obj) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), default=str)


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print(__doc__, file=sys.stderr)
        return 2
    cmd = argv[1]

    if cmd == "hash":
        if len(argv) != 4:
            print("usage: scenicplus_helper hash <config> <dotted_keys_csv>", file=sys.stderr)
            return 2
        cfg = _load(argv[2])
        keys = [k.strip() for k in argv[3].split(",") if k.strip()]
        h = hashlib.sha256(_canonical(_slice(cfg, keys)).encode()).hexdigest()
        print(h)
        return 0

    if cmd == "slice":
        if len(argv) != 4:
            print("usage: scenicplus_helper slice <config> <dotted_keys_csv>", file=sys.stderr)
            return 2
        cfg = _load(argv[2])
        keys = [k.strip() for k in argv[3].split(",") if k.strip()]
        print(_canonical(_slice(cfg, keys)))
        return 0

    if cmd == "get":
        if len(argv) != 4:
            print("usage: scenicplus_helper get <config> <dotted_key>", file=sys.stderr)
            return 2
        cfg = _load(argv[2])
        v = _pluck(cfg, argv[3])
        if v is None:
            return 1
        print(v)
        return 0

    if cmd == "getcsv":
        if len(argv) != 4:
            print("usage: scenicplus_helper getcsv <config> <dotted_key>", file=sys.stderr)
            return 2
        cfg = _load(argv[2])
        v = _pluck(cfg, argv[3]) or []
        if not isinstance(v, (list, tuple)):
            print(f"{argv[3]} is not a list (got {type(v).__name__})", file=sys.stderr)
            return 2
        print(",".join(str(x) for x in v))
        return 0

    print(f"unknown subcommand: {cmd}", file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv))
