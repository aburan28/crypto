#!/usr/bin/env python3
"""Summarise a matched mq-fes round: correctness, oracle counters, paired ratios.

    python3 summarize.py results/round_01
"""
import json
import math
import statistics
import sys
from collections import defaultdict

T95 = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447, 7: 2.365, 8: 2.306,
       9: 2.262, 10: 2.228, 11: 2.201, 12: 2.179, 13: 2.160, 14: 2.145, 15: 2.131}


def ci_of_log_ratios(logs):
    k = len(logs)
    mean = statistics.fmean(logs)
    if k < 2:
        return math.exp(mean), None, None
    half = T95.get(k - 1, 1.96) * statistics.stdev(logs) / math.sqrt(k)
    return math.exp(mean), math.exp(mean - half), math.exp(mean + half)


def main(directory):
    runs = [json.loads(line) for line in open(f"{directory}/raw.jsonl")]
    by = defaultdict(list)
    for r in runs:
        by[(r["degree"], r["curve_a"], r["split"], r["seed"], r["arm"])].append(r)

    def ok(r):
        rep = r.get("report", {})
        return rep.get("status") == "complete" and rep.get("result", {}).get("verified") is True

    rows = []
    groups = sorted({(d, a, s) for (d, a, s, _, _) in by})
    for d, a, split in groups:
        seeds = sorted({k[3] for k in by if k[:3] == (d, a, split)})
        agg = {arm: defaultdict(float) for arm in ("baseline", "candidate")}
        wall_logs, proc_logs = [], []
        complete = {"baseline": 0, "candidate": 0}
        for seed in seeds:
            med = {}
            for arm in ("baseline", "candidate"):
                reps = by[(d, a, split, seed, arm)]
                complete[arm] += sum(ok(r) for r in reps)
                first = reps[0]["report"]
                c = first.get("counts", {})
                p = c.get("mq_fes_profile", {})
                agg[arm]["trials"] += c.get("trials", 0)
                agg[arm]["relations"] += c.get("relations", 0)
                agg[arm]["points"] += p.get("points", 0)
                agg[arm]["calls"] += p.get("calls", 0)
                for key in ("build_ns", "walk_ns", "lift_ns"):
                    agg[arm][key] += statistics.median(
                        r["report"]["counts"]["mq_fes_profile"][key] for r in reps)
                agg[arm]["collection_s"] += statistics.median(
                    r["report"]["timing_seconds"]["relation_collection"] for r in reps)
                agg[arm]["la_s"] += statistics.median(
                    r["report"]["timing_seconds"]["linear_algebra"] for r in reps)
                med[arm] = (statistics.median(r["report"]["elapsed_seconds"] for r in reps),
                            statistics.median(r["process_seconds"] for r in reps))
                agg[arm]["elapsed_s"] += med[arm][0]
            wall_logs.append(math.log(med["baseline"][0] / med["candidate"][0]))
            proc_logs.append(math.log(med["baseline"][1] / med["candidate"][1]))
        g, lo, hi = ci_of_log_ratios(wall_logs)
        pg, plo, phi = ci_of_log_ratios(proc_logs)
        b, c = agg["baseline"], agg["candidate"]
        rows.append({
            "degree": d, "curve_a": a, "split": split, "seeds": len(seeds),
            "verified_runs": complete, "runs_per_arm": len(seeds) * len(by[(d, a, split, seeds[0], "baseline")]),
            "baseline": dict(b), "candidate": dict(c),
            "points_per_call": {"baseline": b["points"] / max(b["calls"], 1),
                                "candidate": c["points"] / max(c["calls"], 1)},
            "build_us_per_call": {"baseline": b["build_ns"] / max(b["calls"], 1) / 1e3,
                                  "candidate": c["build_ns"] / max(c["calls"], 1) / 1e3},
            "trials_ratio_candidate_over_baseline": c["trials"] / max(b["trials"], 1),
            "oracle_points_ratio_baseline_over_candidate": b["points"] / max(c["points"], 1),
            "elapsed_ratio_total": b["elapsed_s"] / c["elapsed_s"],
            "elapsed_paired_geomean": g, "elapsed_paired_ci95": [lo, hi],
            "process_paired_geomean": pg, "process_paired_ci95": [plo, phi],
        })
    out = {"schema_version": 1, "source": f"{directory}/raw.jsonl", "rows": rows}
    with open(f"{directory}/summary.json", "w") as f:
        json.dump(out, f, indent=2)
    hdr = f"{'inst':>8} {'split':>7} {'ver b/c':>9} {'trials b/c':>11} {'pts/call b→c':>18} {'build µs b→c':>15} {'elapsed b/c s':>14} {'paired wall ×':>22}"
    print(hdr)
    for r in rows:
        b, c = r["baseline"], r["candidate"]
        print(f"{r['degree']:>3}/a={r['curve_a']} {r['split']:>7} "
              f"{r['verified_runs']['baseline']:>4}/{r['verified_runs']['candidate']:<4} "
              f"{int(b['trials']):>5}/{int(c['trials']):<5} "
              f"{r['points_per_call']['baseline']:>9.3g}→{r['points_per_call']['candidate']:<8.3g}"
              f"{r['build_us_per_call']['baseline']:>7.0f}→{r['build_us_per_call']['candidate']:<7.1f}"
              f"{b['elapsed_s']:>7.3f}/{c['elapsed_s']:<7.3f}"
              f"{r['elapsed_paired_geomean']:>7.2f} [{r['elapsed_paired_ci95'][0]:.2f},{r['elapsed_paired_ci95'][1]:.2f}]")


if __name__ == "__main__":
    main(sys.argv[1])
