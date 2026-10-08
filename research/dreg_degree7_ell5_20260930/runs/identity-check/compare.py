#!/usr/bin/env python3
"""Compare the budgeted dense-finish binary's rows with committed rows.

    python3 compare.py [RUN_DIR]

Each <tag>.jsonl in RUN_DIR (default: this file's directory) is the
budgeted binary run whole-cell (--unsat 4 --controls 0 --ffd-max 5) at the
committed cell's seed with KIC_SPARSE_DENSE_BUDGET_MB set; each
<tag>.dmin<D>.u<K>.jsonl is one --unsat-index K --d-min D replay.  The key
compared is draw index, subspace, target, solution count, outcome and FFD.
"""
import json, sys
from pathlib import Path

R = Path(__file__).resolve().parents[3]  # research/
H = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent
REF = {
    "grid": R / "dreg_ell_grid_20260925/runs",
    "ladder28": R / "dreg_fixed_surplus_20260923/runs",
    "control": R / "dreg_surplus_control_20260925/runs",
}
key = lambda r: (r["draw"], tuple(r["v_basis"]), r["x_r"], r["solutions"],
                 json.dumps(r["outcome"], sort_keys=True), r["ffd"])
load = lambda p: [json.loads(l) for l in p.read_text().splitlines()
                  if l.strip() and '"control"' not in l]
rows = mism = measured = 0
for f in sorted(H.glob("*.jsonl")):
    fam, _, rest = f.stem.partition("-")
    cell, _, replay = rest.partition(".")
    want = load(REF[fam] / f"cell-{cell}.jsonl")
    got = load(f)
    if not got:
        print(f"{f.name}: empty"); continue
    if replay:  # one unsat draw, replayed from a higher degree
        k = int(replay.split(".u")[1])
        ref = [r for r in want if r["outcome"]["kind"] != "satisfiable"][k]
        m = int(key(ref) != key(got[0]))
        print(f"{f.name}: draw {ref['draw']} {ref['outcome']} -> {got[0]['outcome']}"
              f" d_min {got[0].get('d_min')}, mismatches {m}, secs {ref['secs']:.1f} -> {got[0]['secs']:.1f}")
        rows += 1; mism += m; measured += 1
        continue
    m = sum(1 for a, b in zip(want, got) if key(a) != key(b)) + abs(len(want) - len(got))
    u = [(a, b) for a, b in zip(want, got) if a["outcome"]["kind"] != "satisfiable"]
    print(f"{f.name}: {len(want)} rows, {len(u)} measured, mismatches {m}, "
          f"measured secs {sum(a['secs'] for a, _ in u):.1f} -> {sum(b['secs'] for _, b in u):.1f}")
    rows += len(want); mism += m; measured += len(u)
print(f"TOTAL rows {rows} measured {measured} mismatches {mism}")
