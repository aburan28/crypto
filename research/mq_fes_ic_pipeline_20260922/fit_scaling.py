#!/usr/bin/env python3
"""Exponent fits for the oracle scaling stage diagnostic.

For each arm: least-squares slope of log2(cost per call) on ℓ over the
distinct ℓ measured (duplicates at one ℓ are averaged in log space), for
word ops and for wall time.  A 2^{c·ℓ} cost has slope c.  The linear split
is also fitted after dividing out its polynomial factor ℓ(ℓ-1)/2 + 3ℓ + 1.

    python3 fit_scaling.py results/oracle_scaling_01
"""
import json
import math
import statistics
import sys
from collections import defaultdict


def slope(xs, ys):
    mx, my = statistics.fmean(xs), statistics.fmean(ys)
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)


def main(directory):
    rows = [json.loads(line) for line in open(f"{directory}/raw.jsonl")]
    assert all(r["agree"] for r in rows), "oracles disagree on a target"
    per = defaultdict(lambda: defaultdict(list))
    for r in rows:
        for arm in r["arms"]:
            per[arm["arm"]][r["ell"]].append(arm)
    fits = {}
    for arm, by_ell in per.items():
        ells = sorted(by_ell)
        ops = [statistics.fmean(math.log2(a["word_ops_per_call"]) for a in by_ell[e]) for e in ells]
        wall = [statistics.fmean(math.log2(a["wall_us_per_call"]) for a in by_ell[e]) for e in ells]
        entry = {"sizes": ells, "word_ops_slope": slope(ells, ops), "wall_slope": slope(ells, wall)}
        if arm == "mqfes-linear":
            poly = [math.log2(e * (e - 1) / 2 + 3 * e + 1) for e in ells]
            entry["word_ops_slope_poly_removed"] = slope(ells, [o - p for o, p in zip(ops, poly)])
        fits[arm] = entry
    big = [r for r in rows if r["ell"] >= 11]
    out = {
        "schema_version": 1,
        "source": f"{directory}/raw.jsonl",
        "regime": "n >= 2*ell (overdetermined; refutations pay the full enumeration)",
        "fits": fits,
        "prediction": {"baseline": "2 * 4^ell word ops (slope 2)",
                       "mqfes-linear": "2^ell * (ell(ell-1)/2 + 3 ell + 1) word ops (slope 1 after the polynomial)"},
        "per_call_wall_ratio_at_ell_ge_11": [
            {"degree": r["degree"], "ell": r["ell"],
             "ratio": r["arms"][0]["wall_us_per_call"] / r["arms"][1]["wall_us_per_call"]} for r in big],
    }
    with open(f"{directory}/summary.json", "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main(sys.argv[1])
