#!/usr/bin/env python3
"""Emit a frozen EC Bench spec for the curves compiled into gpu/ecc/bench.

The generator and curve parameters come from ecref.py, the same source that
generates curve_*.h. EC Bench validates subgroup membership when planning.
"""

import argparse
import json

from ecref import get_curve


def curve_spec(name):
    if name not in ("toy36j0", "toy40"):
        raise ValueError("matched EC Bench curves: toy36j0, toy40")
    c = get_curve(name)
    return {
        "kind": "prime_explicit",
        "p": c.p,
        "a": c.a,
        "b": c.b,
        "group_order": c.h * c.n,
        "r": c.n,
        "gx": c.G[0],
        "gy": c.G[1],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--curve", choices=("toy36j0", "toy40"), required=True)
    parser.add_argument("--targets", type=int, default=2)
    parser.add_argument("--rounds", type=int, default=5)
    parser.add_argument("--target-seed", type=int, default=20261013)
    args = parser.parse_args()
    if args.targets < 1 or args.rounds < 1:
        parser.error("targets and rounds must be positive")
    doc = {
        "schema": "ecbench.spec/v1",
        "label": f"CUDA matched control: {args.curve}",
        "description": "Frozen CPU controls and workload export for one-target CUDA BSGS receipts",
        "workloads": {
            "curves": [curve_spec(args.curve)],
            "targets_per_curve": args.targets,
            "target_seed": args.target_seed,
        },
        "arms": [
            {"name": "rho", "role": "reference", "method": {"id": "rho.negation"}},
            {"name": "bsgs", "role": "candidate", "method": {"id": "bsgs.negation"}},
            {"name": "bsgs-aa", "role": "control", "method": {"id": "bsgs.negation"}},
        ],
        "measurement": {
            "rounds": args.rounds,
            "warmup": 1,
            "order": "alternate",
            "seed": 20261013,
            "isolation_required": "L2",
            "timeout_seconds": 600,
        },
    }
    print(json.dumps(doc, indent=2))


if __name__ == "__main__":
    main()
