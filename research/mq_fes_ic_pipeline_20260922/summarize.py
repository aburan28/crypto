#!/usr/bin/env python3
"""Summarise a matched round: correctness, oracle counters, paired ratios.

Ratios are baseline / arm, per (instance, split): totals over seeds of the
median-over-repetitions value, plus a paired geometric mean over seeds of
the elapsed-time ratio with a 95% t interval on the log ratios.

    python3 summarize.py results/round_02 [--baseline baseline]
"""
import argparse
import json
import math
import statistics
from collections import defaultdict

T95 = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447, 7: 2.365, 8: 2.306,
       9: 2.262, 10: 2.228, 11: 2.201, 12: 2.179, 13: 2.160, 14: 2.145, 15: 2.131}


def ci(logs):
    mean = statistics.fmean(logs)
    if len(logs) < 2:
        return [math.exp(mean), None, None]
    half = T95.get(len(logs) - 1, 1.96) * statistics.stdev(logs) / math.sqrt(len(logs))
    return [math.exp(mean), math.exp(mean - half), math.exp(mean + half)]


def verified(r):
    rep = r.get("report", {})
    return rep.get("status") == "complete" and rep.get("result", {}).get("verified") is True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("directory")
    ap.add_argument("--baseline", default="baseline")
    args = ap.parse_args()
    runs = [json.loads(line) for line in open(f"{args.directory}/raw.jsonl")]
    arms = list(dict.fromkeys(r["arm"] for r in runs))
    by = defaultdict(list)
    for r in runs:
        by[(r["degree"], r["curve_a"], r["split"], r["seed"], r["arm"])].append(r)
    rows = []
    for d, a, split in sorted({k[:3] for k in by}):
        seeds = sorted({k[3] for k in by if k[:3] == (d, a, split)})
        row = {"degree": d, "curve_a": a, "split": split, "seeds": len(seeds), "arms": {}}
        med_elapsed = defaultdict(dict)
        for arm in arms:
            tot = defaultdict(float)
            for seed in seeds:
                reps = by[(d, a, split, seed, arm)]
                tot["runs"] += len(reps)
                tot["verified_runs"] += sum(verified(r) for r in reps)
                if not all(verified(r) for r in reps):
                    continue
                first = reps[0]["report"]["counts"]
                tot["trials"] += first["trials"]
                tot["relations"] += first["relations"]
                prof = first.get("mq_fes_profile")
                if prof:
                    tot["oracle_calls"] += prof["calls"]
                    tot["gray_points"] += prof["points"]
                    tot["word_ops"] += prof.get("word_ops", 2 * prof["points"])
                    tot["linear_steps"] += prof.get("linear_steps", 0)
                    for key in ("build_ns", "walk_ns", "lift_ns"):
                        tot[key] += statistics.median(r["report"]["counts"]["mq_fes_profile"][key] for r in reps)
                for key in ("relation_collection", "linear_algebra", "pair_table"):
                    tot[key + "_s"] += statistics.median(r["report"]["timing_seconds"][key] for r in reps)
                m = statistics.median(r["report"]["elapsed_seconds"] for r in reps)
                tot["elapsed_s"] += m
                med_elapsed[arm][seed] = m
            calls = max(tot["oracle_calls"], 1)
            tot["word_ops_per_call"] = tot["word_ops"] / calls if tot["oracle_calls"] else None
            tot["oracle_us_per_call"] = ((tot["build_ns"] + tot["walk_ns"] + tot["lift_ns"]) / calls / 1e3
                                         if tot["oracle_calls"] else None)
            row["arms"][arm] = dict(tot)
        base = row["arms"][args.baseline]
        for arm in arms:
            x = row["arms"][arm]
            common = [s for s in seeds if s in med_elapsed[arm] and s in med_elapsed[args.baseline]]
            x["elapsed_ratio_total"] = base["elapsed_s"] / x["elapsed_s"] if x["elapsed_s"] else None
            x["elapsed_paired"] = ci([math.log(med_elapsed[args.baseline][s] / med_elapsed[arm][s])
                                      for s in common]) if common else None
            x["word_ops_ratio"] = (base["word_ops"] / x["word_ops"]) if x.get("word_ops") else None
            x["trials_ratio_arm_over_baseline"] = x["trials"] / base["trials"] if base["trials"] else None
        rows.append(row)
    with open(f"{args.directory}/summary.json", "w") as f:
        json.dump({"schema_version": 2, "source": f"{args.directory}/raw.jsonl",
                   "baseline": args.baseline, "arms": arms, "rows": rows}, f, indent=2)
    for row in rows:
        print(f"\n== n={row['degree']} a={row['curve_a']} {row['split']} ({row['seeds']} seeds)")
        print(f"{'arm':>14} {'ver':>7} {'trials':>7} {'oracle µs/call':>15} {'wordops/call':>13} "
              f"{'elapsed s':>10} {'paired × vs base [95% CI]':>28}")
        for arm, x in row["arms"].items():
            p = x["elapsed_paired"]
            pc = f"{p[0]:.2f} [{p[1]:.2f}, {p[2]:.2f}]" if p and p[1] else "—"
            us = f"{x['oracle_us_per_call']:.1f}" if x.get("oracle_us_per_call") else "—"
            wo = f"{x['word_ops_per_call']:.3g}" if x.get("word_ops_per_call") else "—"
            print(f"{arm:>14} {int(x['verified_runs']):>3}/{int(x['runs']):<3} {int(x['trials']):>7} "
                  f"{us:>15} {wo:>13} {x['elapsed_s']:>10.3f} {pc:>28}")


if __name__ == "__main__":
    main()
