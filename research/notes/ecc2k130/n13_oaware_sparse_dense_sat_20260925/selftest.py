#!/usr/bin/env python3
"""Synthetic ABBA accounting controls; these are not measured solver outcomes."""
from __future__ import annotations

import json
from pathlib import Path

import analyze

HERE = Path(__file__).resolve().parent
TRUTH = {row["id"]: row["point_oracle_positive"] for row in
         json.loads((HERE.parent / "n13_oaware_sat_benchmark_20260925/INPUT.json").read_text())["targets"]}


def fixture(export_walls):
    entries = []
    for q in range(8):
        for t in range(4):
            target = f"Q{q}T{t}"
            for engine in analyze.SOLVERS:
                for rep in analyze.REPS:
                    entries.append({"id": target, "solver": engine,
                                    "representation": rep,
                                    "verdict": "SAT" if TRUTH[target] else "UNSAT",
                                    "wall_seconds": 1.0 if rep == "dense" else .5,
                                    "query_setup_wall_seconds": .01,
                                    "verify_wall_seconds": .01,
                                    "user_cpu_seconds": .1,
                                    "system_cpu_seconds": .01,
                                    "sampled_peak_tree_rss_bytes": 100})
    exports = [{"representation": rep, "wall_seconds": wall}
               for rep, wall in zip(("dense", "sparse", "sparse", "dense"), export_walls)]
    return analyze.summarize(
        {"mode": "smoke", "pass": True},
        {"mode": "panel", "decision": "COMPLETE",
         "preflight_wall_seconds": .1, "process_wall_seconds": 1.0},
        {"mode": "panel", "exports": exports,
         "entries": entries, "base_load_wall_seconds": .1},
        {"dense_replay": 2.0, "sparse_replay": 2.0})


def main():
    robust = fixture((1.0, 1.0, 1.0, 1.0))
    assert robust["decision"] == "SPARSE_TOY_NEXT_RUNG_PREFERENCE"
    assert all(robust["engines"][engine]["sparse_robust_toy_wall_improvement"]
               for engine in analyze.SOLVERS)
    reversed_order = fixture((1.0, 100.0, 1.0, 100.0))
    assert reversed_order["decision"] == "CORRECT_BUT_MIXED_OR_INCONCLUSIVE"
    assert all(reversed_order["engines"][engine]["pair_order_reverses_ranking"]
               for engine in analyze.SOLVERS)
    print("synthetic ABBA accounting controls PASS; no solver ran")


if __name__ == "__main__":
    main()
