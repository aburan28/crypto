#!/usr/bin/env python3
"""The gate of note section 14.2, evaluated, with the context of section 14.3.

    python3 research/pkm_tower_round5_20260926/compare_gate.py

Sets every signature row (cells GP, GD, S4) beside round 4's F4 row of the
same system and prints, per system: the verdicts, F4's D, D_sig, D_lm, the
widths, the multiply-adds and each process's peak memory (the F4 peaks are
round 4's). Then it applies the gate: a module order passes when steps by
polynomial degree give D_sig = D on all four m = 3, N = 9 and 12 systems.

Exits non-zero if a verdict splits or a signature row has no F4 row.
"""

import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROUND4 = os.path.join(HERE, "..", "pkm_tower_round4_20260926")
GATE = {(3, 9), (3, 12)}
ORDERS = {"GP": "position over term", "GD": "signature degree first"}


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
                    r["_file"] = os.path.basename(path)[:-len(".jsonl")]
                    out.append(r)
    return out


def peaks(path):
    if not os.path.exists(path):
        return {}
    return {r["run"]: r for r in (json.loads(l) for l in open(path) if l.strip())}


def main():
    bad = 0
    f4 = {instance_key(r): r for r in rows_of(os.path.join(ROUND4, "runs", "T*.jsonl"))}
    f4_mem = peaks(os.path.join(ROUND4, "memory.jsonl"))
    mem = peaks(os.path.join(HERE, "memory.jsonl"))
    sig = rows_of(os.path.join(HERE, "runs", "*.jsonl"))

    print("| cell | m | N | target | verdict F4 / sig | D | D_sig | D_lm | width F4 / sig | "
          "muladds F4 | sig / F4 | zero rows | peak MB F4 / sig |")
    print("|:--|--:|--:|--:|:--|--:|--:|--:|:--|--:|--:|--:|:--|")
    gate = {cell: [] for cell in ORDERS}
    for s in sorted(sig, key=lambda r: (r["_file"][:2], r["m"], r["N"], r["target_index"])):
        cell = s["_file"][:2]
        f = f4.get(instance_key(s))
        if f is None:
            print(f"no F4 row for {s['_file']} target {s['target_index']}")
            bad += 1
            continue
        v = lambda r: "timed out" if r["timed_out"] else ("refuted" if r["inconsistent"] else "solvable")
        if not s["timed_out"] and f["inconsistent"] != s["inconsistent"]:
            print(f"VERDICT MISMATCH {s['_file']} target {s['target_index']}")
            bad += 1
        f_mb = f4_mem.get(f["_file"], {}).get("peak_rss_mb")
        s_mb = mem.get(s["_file"], {}).get("peak_rss_mb")
        ratio = s["muladds"] / max(f["muladds"], 1)
        print(f"| {cell} | {s['m']} | {s['N']} | {s['target_index']} | {v(f)} / {v(s)} | "
              f"{f['solving_degree_max']} | {s['solving_degree_max']} | "
              f"{s['sig'].get('lm_degree_max')} | {f['max_cols']} / {s['max_cols']} | "
              f"{f['muladds']:.3g} | {ratio:.2f} | {s['sig']['zero_rows']} | {f_mb} / {s_mb} |")
        if cell in gate and (s["m"], s["N"]) in GATE:
            ok = not s["timed_out"] and s["solving_degree_max"] == f["solving_degree_max"]
            gate[cell].append(ok)

    lost = [r for r in mem.values() if r["exit"] != 0 or r["rows"] < 2]
    if lost:
        print("\nProcesses that ended without writing both rows:\n")
        for r in sorted(lost, key=lambda r: r["run"]):
            print(f"- {r['run']}: exit {r['exit']}, {r['rows']} rows, peak {r['peak_rss_mb']} MB")
            if r["cell"] in gate and (r["m"], r["N"]) in GATE:
                gate[r["cell"]].extend([False] * (2 - r["rows"]))

    print("\n## The gate (note section 14.2)\n")
    if not any(gate.values()):
        print("No gate rows yet: the gate has not run.")
        return 1 if bad else 0
    passed = []
    for cell, name in ORDERS.items():
        results = gate[cell]
        ok = len(results) == 4 and all(results)
        print(f"- {name} ({cell}): D_sig = D on {sum(results)} of {len(results)} gate systems "
              f"(4 needed): {'passes' if ok else 'fails'}.")
        if ok:
            passed.append(cell)
    if passed:
        print(f"\nThe gate passes for {', '.join(ORDERS[c] for c in passed)}: pre-register "
              "round 5 (section 14.2) before any other run of the variant.")
    else:
        print("\nThe gate fails under both orders: the variant is not worth a round "
              "(sections 13.10 and 14.2), and F4 stays the measuring engine.")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
