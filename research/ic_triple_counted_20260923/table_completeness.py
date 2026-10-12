#!/usr/bin/env python3
"""Is the triple table complete?  Enumerate every C_0 + P_j + P_l on a real base.

    python3 table_completeness.py TRIPLE_WORKER [--cell n23a1]

Runs the triple arm once (the first fixture of the exposed probe stream,
`random.Random(20260923)`), rebuilds its base from the reported orbit
representatives with `oracle.py`'s own curve arithmetic, forms every sum the
table is meant to hold -- row 0, `j ≤ l` -- and counts distinct Frobenius
orbits.  If that count equals the arm's reported `sum_orbits`, every orbit that
exists is stored, and any shortfall against the multiset count is sums that
coincide, not entries the table lost.  Pure Python: `n23a1` takes about a minute.
"""
import argparse
import json
import random
import subprocess
import sys
from math import comb
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "research/ic_candidate_tournament_20260915"))
from oracle import Curve  # noqa: E402

CELLS = {"n23a1": (23, 1), "n37a0": (37, 0), "n43a1": (43, 1)}
CONFIG = {"batch_trials": 1, "linear_algebra": "sparse", "max_trials": 65536,
          "solver": "triple_table", "summands": 4}


def run(worker, payload):
    out = subprocess.run([worker], input=json.dumps(payload), text=True, capture_output=True, timeout=1800)
    return json.loads(out.stdout)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("worker")
    ap.add_argument("--cell", default="n23a1")
    args = ap.parse_args()
    rng = random.Random(20260923)
    tseed, aseed = rng.getrandbits(64), rng.getrandbits(64)
    n, a = CELLS[args.cell]
    job = {"degree": n, "curve_a": a, "target_seeds": [tseed], "algorithm_seed": aseed, "config": CONFIG,
           "factor_base": {"kind": "subgroup_orbits", "seed": 43, "points": 222}}
    fixture = run(args.worker, dict(job, mode="fixture"))["fixture"]
    report = run(args.worker, dict(job, mode="ic"))
    c = Curve(fixture)
    base = []
    for rep in (c.decode(p) for p in report["factor_base_orbits"]):
        q = rep
        for _ in range(c.n):
            base += [q, c.neg(q)]
            q = c.frob(q)

    def orbit_key(p):
        best, q = None, p
        for _ in range(c.n):
            best = q[0] if best is None or q[0] < best else best
            q = c.frob(q)
        return best

    # The arm skips a row whose pair C_0 + P_j is the identity (P_j = −C_0): its
    # "triples" are bare base points, and a witness through it would carry a
    # cancelling pair.  Count with and without it.
    pairs = [c.add(base[0], x) for x in base]
    keys, kept, entries = set(), set(), 0
    for j in range(len(base)):
        for l in range(j, len(base)):
            s = c.add(pairs[j], base[l])
            entries += 1
            if s is not None:
                keys.add(orbit_key(s))
                if pairs[j] is not None:
                    kept.add(orbit_key(s))
    size, per = len(base), 2 * c.n
    multisets = (comb(size + 2, 3) - comb(size - per + 2, 3)) // per
    print(json.dumps({"cell": args.cell, "target_seed": tseed, "algorithm_seed": aseed,
                      "base_points": size, "table_entries": entries,
                      "cancelling_rows": [j for j, p in enumerate(pairs) if p is None],
                      "distinct_orbits_enumerated": len(keys),
                      "distinct_orbits_without_cancelling_rows": len(kept),
                      "sum_orbits_reported_by_arm": report["pair_table"]["sum_orbits"],
                      "multiset_orbits_counted": multisets}))


if __name__ == "__main__":
    main()
