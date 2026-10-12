#!/usr/bin/env python3
"""Compare the dense-finish binary's rows with committed rows, draw by draw.

    python3 research/dreg_fixed_surplus_20260923/runs/identity-check-dense/compare.py

Each <tag>.jsonl here is the dense-finish dreg_ladder run whole-cell
(--unsat 4 --controls 0) at the committed cell's seed; the key compared is
draw index, subspace, target, solution count, outcome and FFD.
"""
import json, sys
from pathlib import Path
R = Path(__file__).resolve().parents[3]  # research/
H = Path(__file__).resolve().parent
PAIRS = [(f"grid-{c}", R / "dreg_ell_grid_20260925/runs" / f"cell-{c}.jsonl") for c in ("4-3-7", "5-3-7", "8-2-7", "10-2-7", "12-2-7")]
PAIRS += [(f"ladder28-{c}", R / "dreg_fixed_surplus_20260923/runs" / f"cell-{c}.jsonl") for c in ("5-2-7", "7-2-7", "7-3-7", "9-3-7")]
PAIRS += [(f"control-{c}", R / "dreg_surplus_control_20260925/runs" / f"cell-{c}.jsonl") for c in ("7-4-7", "8-4-7")]
key = lambda r: (r["draw"], tuple(r["v_basis"]), r["x_r"], r["solutions"], json.dumps(r["outcome"], sort_keys=True), r["ffd"])
rows = mism = unsat = 0
t_old = t_new = 0.0
for tag, ref in PAIRS:
    got_f = H / f"{tag}.jsonl"
    if not got_f.exists() or not got_f.read_text().strip():
        print(f"{tag}: not run"); continue
    want = [json.loads(l) for l in ref.read_text().splitlines() if l.strip() and '"control"' not in l]
    got = [json.loads(l) for l in got_f.read_text().splitlines() if l.strip()]
    m = sum(1 for a, b in zip(want, got) if key(a) != key(b)) + abs(len(want) - len(got))
    u = [(a, b) for a, b in zip(want, got) if a["outcome"]["kind"] != "satisfiable"]
    to, tn = sum(a["secs"] for a, _ in u), sum(b["secs"] for _, b in u)
    rows += len(want); mism += m; unsat += len(u); t_old += to; t_new += tn
    print(f"{tag}: {len(want)} rows, {len(u)} measured, mismatches {m}, measured secs {to:.1f} -> {tn:.1f}")
print(f"TOTAL rows {rows} measured {unsat} mismatches {mism} measured secs {t_old:.0f} -> {t_new:.0f}")
