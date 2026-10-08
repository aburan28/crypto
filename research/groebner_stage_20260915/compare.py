#!/usr/bin/env python3
"""Compare a candidate Gröbner-stage run against the frozen baseline.

    python3 compare.py baseline candidate --output summary.json

Refuses any comparison whose two sides did not decide the same instances
identically: the verdict digest over every target, and the solver's own
counters (reductions, infeasibility certificates, propagations, splits)
must match rung for rung.  A faster stage that answers differently is not
a faster stage.

Ratios are baseline / candidate, so above 1 is an improvement, and are
taken over the median of the repetitions.  The unit is 64-bit word XORs
in the Macaulay elimination; wall time is reported beside it as a
practicality note, never as the metric (AGENTS.md §6).
"""
import argparse, json, pathlib, statistics, sys

IDENTITY = ["verdict_digest", "reductions", "infeasible_branches", "propagations",
            "splits", "oversize", "targets", "decomposed", "f4_calls", "ell",
            "factor_base_points"]

def load(root):
    reps = sorted(pathlib.Path(root).glob("rep*/stage.json"))
    if not reps:
        sys.exit(f"no rep*/stage.json under {root}")
    docs = [json.loads(p.read_text()) for p in reps]
    keys = [tuple((r["curve"], r["m"]) for r in d["rows"]) for d in docs]
    if len(set(keys)) != 1:
        sys.exit(f"{root}: repetitions do not cover the same instances")
    return docs

def median(docs, i, field):
    return statistics.median(d["rows"][i][field] for d in docs)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("baseline"); ap.add_argument("candidate")
    ap.add_argument("--output")
    args = ap.parse_args()
    base, cand = load(args.baseline), load(args.candidate)
    if len(base[0]["rows"]) != len(cand[0]["rows"]):
        sys.exit("instance ladders differ")

    rows, failures = [], []
    for i, (b, c) in enumerate(zip(base[0]["rows"], cand[0]["rows"])):
        for field in IDENTITY:
            if b[field] != c[field]:
                failures.append(f"{b['curve']} m={b['m']}: {field} {b[field]} -> {c[field]}")
        # Deterministic counters must also be stable across repetitions.
        for docs, side in ((base, "baseline"), (cand, "candidate")):
            if len({d["rows"][i]["word_ops"] for d in docs}) != 1:
                failures.append(f"{b['curve']} m={b['m']}: {side} word_ops not reproducible")
        ops_b, ops_c = median(base, i, "word_ops"), median(cand, i, "word_ops")
        rows.append({
            "curve": b["curve"], "m": b["m"], "ell": b["ell"],
            "targets": b["targets"], "decomposed": b["decomposed"],
            "verdict_digest": b["verdict_digest"],
            "baseline_word_ops": ops_b, "candidate_word_ops": ops_c,
            "word_op_ratio": ops_b / ops_c if ops_c else None,
            "baseline_wall_ns": median(base, i, "wall_ns"),
            "candidate_wall_ns": median(cand, i, "wall_ns"),
            "wall_ratio": median(base, i, "wall_ns") / median(cand, i, "wall_ns"),
            "phases_baseline_ns": {p: median(base, i, p + "_ns") for p in ("build", "reduce", "readback")},
            "phases_candidate_ns": {p: median(cand, i, p + "_ns") for p in ("build", "reduce", "readback")},
        })

    if failures:
        print("REFUSED: the two sides did not decide the ladder identically")
        for f in failures:
            print("  " + f)
        sys.exit(1)

    tot_b = sum(r["baseline_word_ops"] for r in rows)
    tot_c = sum(r["candidate_word_ops"] for r in rows)
    wall_b = sum(r["baseline_wall_ns"] for r in rows)
    wall_c = sum(r["candidate_wall_ns"] for r in rows)
    summary = {
        "unit": "64-bit word XORs in the Macaulay elimination",
        "scope": "decomposition-oracle stage only; not an end-to-end ECDLP cost",
        "identity_checked": IDENTITY,
        "repetitions": {"baseline": len(base), "candidate": len(cand)},
        "rows": rows,
        "total_word_op_ratio": tot_b / tot_c,
        "total_wall_ratio": wall_b / wall_c,
        "min_word_op_ratio": min(r["word_op_ratio"] for r in rows),
        "min_wall_ratio": min(r["wall_ratio"] for r in rows),
    }
    width = max(len(f"{r['curve']} m={r['m']}") for r in rows)
    print(f"{'instance':<{width}}  {'word XORs base':>15} {'cand':>13} {'ratio':>6}  {'wall ratio':>10}")
    for r in rows:
        name = f"{r['curve']} m={r['m']}"
        print(f"{name:<{width}}  {r['baseline_word_ops']:>15.0f} {r['candidate_word_ops']:>13.0f} "
              f"{r['word_op_ratio']:>6.2f}  {r['wall_ratio']:>10.2f}")
    print(f"{'total':<{width}}  {tot_b:>15.0f} {tot_c:>13.0f} {tot_b/tot_c:>6.2f}  {wall_b/wall_c:>10.2f}")
    print("\nEvery instance decided identically; ratios are baseline / candidate.")
    if args.output:
        pathlib.Path(args.output).write_text(json.dumps(summary, indent=2) + "\n")
        print(f"Wrote {args.output}")

if __name__ == "__main__":
    main()
