#!/usr/bin/env python3
"""The signature engine beside F4 on the same systems (note sections 13.4-13.6).

    python3 research/pkm_tower_round4_20260926/compare_engines.py

1. This round's F4 rows (cells T2-T4) must repeat round 2's rows of the same
   systems on every field but the wall clock, and the trace flag's absence.
2. Each signature row (cells G2-G4) is set beside its F4 row: the solving
   degree, the verdict, the width, the multiply-adds, the rows built and the
   rows reduced to zero, and each process's peak memory from `memory.jsonl`.
3. The predictions of section 13.5 and the adoption rule of section 13.6 are
   evaluated on that table.

Exits non-zero if an F4 row differs from round 2's or a verdict splits.
"""

import glob
import json
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROUND2 = os.path.join(HERE, "..", "pkm_tower_round2_20260925", "runs")
TIMING = {"ms", "wall_s"}


def instance_key(r):
    """The same key as the pilot's analyze.py: the system, whichever engine."""
    return (
        r["p"], r["kind"], r["m"], r["control"], r["t"], r["g"], r["target"],
        r["target_index"], r["x_r"], r["curve"]["a"], r["curve"]["b"],
        json.dumps(r["tower"], sort_keys=True),
    )


def rows_of(pattern):
    out = []
    for path in sorted(glob.glob(pattern)):
        for line in open(path):
            if line.strip():
                r = json.loads(line)
                if "N" in r:
                    r["_file"] = os.path.basename(path)
                    out.append(r)
    return out


def peaks():
    """Peak resident memory per run file stem, in MB."""
    path = os.path.join(HERE, "memory.jsonl")
    if not os.path.exists(path):
        return {}
    recs = [json.loads(l) for l in open(path) if l.strip()]
    return {r["run"]: r for r in recs}


def f4_rows_built(r):
    return r["critical_pairs_reduced"] + r["tower_pairs_reduced"] + r["reducer_rows"]


def sig_rows_built(r):
    return r["sig"]["s_rows"] + r["sig"]["reducer_rows"]


def main():
    bad = 0
    round2 = {instance_key(r): r for r in rows_of(os.path.join(ROUND2, "[KM][134]-kummer-*.jsonl"))
              if r.get("engine") == "f4_fp_tower"}
    tower = rows_of(os.path.join(HERE, "runs", "T*.jsonl"))
    sig = rows_of(os.path.join(HERE, "runs", "G*.jsonl"))
    mem = peaks()

    # 1. F4 against round 2.
    same = 0
    for r in tower:
        old = round2.get(instance_key(r))
        if old is None:
            print(f"F4 row without a round-2 row: {r['_file']} N={r['N']} target {r['target_index']}")
            bad += 1
            continue
        diff = sorted(f for f in (set(r) | set(old)) - TIMING - {"_file"} if r.get(f) != old.get(f))
        if diff:
            print(f"{r['_file']} N={r['N']} target {r['target_index']}: differs from round 2 on {diff}")
            bad += 1
        else:
            same += 1
    print(f"## F4 against round 2\n\n{same} of {len(tower)} rows identical but for the wall clock.\n")

    # 2. The table.
    by_key = {instance_key(r): r for r in tower}
    print("## The signature engine beside F4\n")
    print("| m | N | target | verdict F4 / sig | D | D_sig | width F4 | width sig | "
          "muladds F4 | muladds sig | sig / F4 | rows F4 | rows sig | zero rows sig | "
          "peak MB F4 | peak MB sig |")
    print("|--:|--:|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|")
    finished, timed_out = [], []
    for s in sorted(sig, key=lambda r: (r["m"], r["N"], r["target_index"])):
        f = by_key.get(instance_key(s))
        f_mem = mem.get(f["_file"][:-6], {}).get("peak_rss_mb") if f else None
        s_mem = mem.get(s["_file"][:-6], {}).get("peak_rss_mb")
        if s["timed_out"]:
            timed_out.append(s)
            print(f"| {s['m']} | {s['N']} | {s['target_index']} | timed out | "
                  f"{f['solving_degree_max'] if f else '—'} | ≥ {s['solving_degree_max']} | "
                  f"| {s['max_cols']} | | {s['muladds']:.3g} | | | {sig_rows_built(s)} | "
                  f"{s['sig']['zero_rows']} | {f_mem} | {s_mem} |")
            continue
        if f is None or f["timed_out"]:
            print(f"| {s['m']} | {s['N']} | {s['target_index']} | no F4 row | | "
                  f"{s['solving_degree_max']} | | | | | | | | | | |")
            bad += 1
            continue
        finished.append((f, s))
        if f["inconsistent"] != s["inconsistent"]:
            print(f"VERDICT MISMATCH m={s['m']} N={s['N']} target {s['target_index']}")
            bad += 1
        v = lambda r: "refuted" if r["inconsistent"] else "solvable"
        print(f"| {s['m']} | {s['N']} | {s['target_index']} | {v(f)} / {v(s)} | "
              f"{f['solving_degree_max']} | {s['solving_degree_max']} | {f['max_cols']} | "
              f"{s['max_cols']} | {f['muladds']:.3g} | {s['muladds']:.3g} | "
              f"{s['muladds'] / max(f['muladds'], 1):.2f} | {f4_rows_built(f)} | "
              f"{sig_rows_built(s)} | {s['sig']['zero_rows']} | {f_mem} | {s_mem} |")

    # 3. Predictions and the rule.
    print("\n## Section 13.5 and 13.6\n")
    agree = sum(f["inconsistent"] == s["inconsistent"] for f, s in finished)
    print(f"1. verdicts agree on {agree} of {len(finished)} systems both engines finished "
          f"({len(timed_out)} signature systems timed out); verify.py is run separately.")
    ge = sum(s["solving_degree_max"] >= f["solving_degree_max"] for f, s in finished)
    gt = sum(s["solving_degree_max"] > f["solving_degree_max"] for f, s in finished)
    eq = sum(s["solving_degree_max"] == f["solving_degree_max"] for f, s in finished)
    print(f"2. D_sig >= D on {ge} of {len(finished)}, D_sig > D on {gt}, D_sig = D on {eq}.")
    zero = sum(s["sig"]["zero_rows"] for _, s in finished)
    built = sum(sig_rows_built(s) for _, s in finished)
    worst = max((s["sig"]["zero_rows"] / max(sig_rows_built(s), 1) for _, s in finished), default=0)
    print(f"3. zero rows {zero} of {built} built ({100 * zero / max(built, 1):.3f}%); "
          f"worst system {100 * worst:.3f}%.")
    m4 = [(f, s) for f, s in finished if f["m"] == 4]
    for n in sorted({s["N"] for _, s in m4}):
        pairs = [(f, s) for f, s in m4 if s["N"] == n]
        more = sum(s["muladds"] > f["muladds"] for f, s in pairs)
        print(f"4. m = 4, N = {n}: the signature engine does more multiply-adds on {more} "
              f"of {len(pairs)} systems.")
    for n in (8, 12):
        if not any(s["N"] == n for _, s in m4):
            print(f"4. m = 4, N = {n}: no system both engines finished.")
    t4 = [r for r in mem.values() if r["m"] == 4 and r["N"] == 12]
    if len(t4) == 2:
        f_mb = next(r["peak_rss_mb"] for r in t4 if r["engine"] == "tower")
        s_mb = next(r["peak_rss_mb"] for r in t4 if r["engine"] == "sig")
        print(f"4. m = 4, N = 12: peak memory F4 {f_mb} MB, signature {s_mb} MB.")

    adopt = bool(finished) and eq == len(finished)
    for m in sorted({s["m"] for _, s in finished}):
        top = max(s["N"] for f, s in finished if s["m"] == m)
        pairs = [(f, s) for f, s in finished if s["m"] == m and s["N"] == top]
        cheaper = all(s["muladds"] <= f["muladds"] for f, s in pairs)
        f_mb = mem.get(pairs[0][0]["_file"][:-6], {}).get("peak_rss_mb")
        s_mb = mem.get(pairs[0][1]["_file"][:-6], {}).get("peak_rss_mb")
        lighter = f_mb is not None and s_mb is not None and s_mb <= f_mb
        print(f"   m = {m}, largest finished N = {top}: multiply-adds no more than F4's: {cheaper}; "
              f"peak memory {s_mb} MB against {f_mb} MB, no more: {lighter}.")
        adopt = adopt and cheaper and lighter
    print(f"\nAdoption (section 13.6): {'yes' if adopt else 'no, F4 stays the measuring engine'}.")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
