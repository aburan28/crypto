#!/usr/bin/env python3
"""Compare two Gröbner-stage runs that walk *different* splitting trees.

    python3 compare_cross_tree.py reference candidate --output summary.json

`research/groebner_stage_20260915/compare.py` refuses a comparison unless
both sides decided every target identically down to the verdict digest,
which pins the exact decomposition returned and the solver's node counts.
A change of split rule keeps the algebra and the answers but not the tree:
when several decompositions exist a different one is found first, and the
node counts move.  That is a legitimate whole-method comparison under
AGENTS.md §8 — same inputs, correctness preserved, total work compared —
but it is not the same-tree identity the frozen suite enforces, so this
script states exactly what it does check, refuses anything less, and
labels its output as a cross-tree comparison.

Checked, rung for rung:
  - the same instances (curve, m, ell, factor-base size, target count);
  - the same number of targets decomposed and none exhausted — every
    decomposition the harness counts has already been lifted and verified
    against the group identity by `groebner_decompose`;
  - the same oversize count (no side silently skipped a matrix).

Reported: word operations, medians over repetitions, with the node
counts of both sides beside them so the tree change is visible.  Ratios
are reference / candidate.  Wall time is a practicality note.
"""
import argparse, json, pathlib, statistics, sys

SAME_INSTANCE = ["curve", "m", "ell", "factor_base_points", "targets"]
SAME_ANSWERS = ["decomposed", "oversize", "exhausted"]

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
    ap.add_argument("reference"); ap.add_argument("candidate")
    ap.add_argument("--output")
    args = ap.parse_args()
    ref, cand = load(args.reference), load(args.candidate)
    if len(ref[0]["rows"]) != len(cand[0]["rows"]):
        sys.exit("instance ladders differ")

    rows, failures = [], []
    for i, (b, c) in enumerate(zip(ref[0]["rows"], cand[0]["rows"])):
        for field in SAME_INSTANCE + SAME_ANSWERS:
            if b.get(field) != c.get(field):
                failures.append(f"{b['curve']} m={b['m']}: {field} {b.get(field)} -> {c.get(field)}")
        if b.get("exhausted") is None or c.get("exhausted") is None:
            failures.append(f"{b['curve']} m={b['m']}: a side does not record exhausted budgets")
        if c.get("exhausted", 1) != 0:
            failures.append(f"{b['curve']} m={b['m']}: candidate exhausted its budget on {c.get('exhausted')} target(s)")
        for docs, side in ((ref, "reference"), (cand, "candidate")):
            if len({d["rows"][i]["word_ops"] for d in docs}) != 1:
                failures.append(f"{b['curve']} m={b['m']}: {side} word_ops not reproducible")
        ops_b, ops_c = median(ref, i, "word_ops"), median(cand, i, "word_ops")
        rows.append({
            "curve": b["curve"], "m": b["m"], "ell": b["ell"],
            "targets": b["targets"], "decomposed": b["decomposed"],
            "same_tree": b["verdict_digest"] == c["verdict_digest"]
                and all(b[k] == c[k] for k in ("reductions", "splits", "propagations", "infeasible_branches")),
            "reference": {"word_ops": ops_b, "reductions": b["reductions"], "splits": b["splits"],
                          "infeasible_branches": b["infeasible_branches"], "wall_ns": median(ref, i, "wall_ns")},
            "candidate": {"word_ops": ops_c, "reductions": c["reductions"], "splits": c["splits"],
                          "infeasible_branches": c["infeasible_branches"], "wall_ns": median(cand, i, "wall_ns")},
            "word_op_ratio": ops_b / ops_c if ops_c else None,
            "wall_ratio": median(ref, i, "wall_ns") / median(cand, i, "wall_ns"),
        })

    if failures:
        print("REFUSED: the two sides did not answer the ladder identically")
        for f in failures:
            print("  " + f)
        sys.exit(1)

    tot_b = sum(r["reference"]["word_ops"] for r in rows)
    tot_c = sum(r["candidate"]["word_ops"] for r in rows)
    wall_b = sum(r["reference"]["wall_ns"] for r in rows)
    wall_c = sum(r["candidate"]["wall_ns"] for r in rows)
    summary = {
        "kind": "cross-tree comparison: same instances, same verified answers, splitting tree may differ",
        "unit": "64-bit word operations",
        "scope": "decomposition-oracle stage only; not an end-to-end ECDLP cost",
        "identity_checked": SAME_INSTANCE + SAME_ANSWERS,
        "repetitions": {"reference": len(ref), "candidate": len(cand)},
        "rows": rows,
        "total_word_op_ratio": tot_b / tot_c,
        "total_wall_ratio": wall_b / wall_c,
        "min_word_op_ratio": min(r["word_op_ratio"] for r in rows),
    }
    width = max(len(f"{r['curve']} m={r['m']}") for r in rows)
    print(f"{'instance':<{width}}  {'word ops ref':>14} {'cand':>13} {'ratio':>6}  {'nodes ref':>9} {'cand':>6}  same tree")
    for r in rows:
        name = f"{r['curve']} m={r['m']}"
        print(f"{name:<{width}}  {r['reference']['word_ops']:>14.0f} {r['candidate']['word_ops']:>13.0f} "
              f"{r['word_op_ratio']:>6.2f}  {r['reference']['reductions']:>9} {r['candidate']['reductions']:>6}  {r['same_tree']}")
    print(f"{'total':<{width}}  {tot_b:>14.0f} {tot_c:>13.0f} {tot_b/tot_c:>6.2f}")
    print("\nSame instances and the same verified answers on every rung; the tree is not required to match.")
    if args.output:
        pathlib.Path(args.output).write_text(json.dumps(summary, indent=2) + "\n")
        print(f"Wrote {args.output}")

if __name__ == "__main__":
    main()
