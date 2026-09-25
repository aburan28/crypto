#!/usr/bin/env python3
"""Predeclared paired dense/sparse toy-stage cost views."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

SOLVERS = ("cryptominisat5", "kissat", "cadical")
REPS = ("dense", "sparse")


def stage(rows):
    return sum(r["query_setup_wall_seconds"] + r["wall_seconds"] +
               r["verify_wall_seconds"] for r in rows)


def summarize(smoke, receipt, panel, archived_audit):
    assert smoke["mode"] == "smoke" and smoke["pass"]
    assert receipt["mode"] == panel["mode"] == "panel"
    assert receipt["decision"] == "COMPLETE"
    assert len(panel["exports"]) == 4 and len(panel["entries"]) == 192
    exports = panel["exports"]
    assert [r["representation"] for r in exports] == ["dense", "sparse", "sparse", "dense"]
    export_wall = {("dense", 0): exports[0]["wall_seconds"],
                   ("sparse", 0): exports[1]["wall_seconds"],
                   ("sparse", 1): exports[2]["wall_seconds"],
                   ("dense", 1): exports[3]["wall_seconds"]}
    common = receipt["preflight_wall_seconds"] + panel["base_load_wall_seconds"]
    by = {(r["id"], r["solver"], r["representation"]): r for r in panel["entries"]}
    engines = {}
    for engine in SOLVERS:
        engines[engine] = {}
        for rep in REPS:
            rows = [by[(f"Q{q}T{t}", engine, rep)] for q in range(8) for t in range(4)]
            q_rows = []
            for q in range(8):
                block = rows[4*q:4*q+4]
                first = next((t for t, row in enumerate(block) if row["verdict"] == "SAT"), None)
                verdict = ("POSITIVE" if first is not None and
                           all(row["verdict"] == "UNSAT" for row in block[:first]) else
                           "NEGATIVE" if all(row["verdict"] == "UNSAT" for row in block) else
                           "CENSORED")
                q_rows.append({"q": q, "status": verdict,
                               "branch_statuses": [row["verdict"] for row in block],
                               "first_witness_t": first if verdict == "POSITIVE" else None,
                               "first_witness_stage_wall_seconds":
                                   stage(block[:first+1]) if verdict == "POSITIVE" else None,
                               "all_four_stage_wall_seconds": stage(block)})
            complete = all(row["verdict"] in ("SAT", "UNSAT") for row in rows)
            correct = all(q["status"] == ("POSITIVE" if q["q"] < 4 else "NEGATIVE")
                          for q in q_rows)
            stage_wall = stage(rows)
            pair_walls = [common + stage_wall + export_wall[(rep, pair)] for pair in (0, 1)]
            engines[engine][rep] = {
                "complete_correct": complete and correct,
                "branch_sat": sum(row["verdict"] == "SAT" for row in rows),
                "branch_unsat_proof_unchecked": sum(row["verdict"] == "UNSAT" for row in rows),
                "branch_censored_invalid": sum(row["verdict"] not in ("SAT", "UNSAT") for row in rows),
                "solver_child_wall_sum_seconds": sum(row["wall_seconds"] for row in rows),
                "solver_child_cpu_sum_seconds": sum(row["user_cpu_seconds"] + row["system_cpu_seconds"]
                                                     for row in rows),
                "query_setup_wall_sum_seconds": sum(row["query_setup_wall_seconds"] for row in rows),
                "point_verify_wall_sum_seconds": sum(row["verify_wall_seconds"] for row in rows),
                "max_child_wall_seconds": max(row["wall_seconds"] for row in rows),
                "max_sampled_tree_rss_bytes": max(row["sampled_peak_tree_rss_bytes"] for row in rows),
                "query_plus_child_plus_verification_wall_seconds": stage_wall,
                "full_portfolio_wall_seconds_by_export_pair": pair_walls,
                "archived_audit_inclusive_wall_seconds_by_export_pair":
                    [value + archived_audit[f"{rep}_replay"] for value in pair_walls],
                "q": q_rows,
            }
        dense = engines[engine]["dense"]
        sparse = engines[engine]["sparse"]
        ratios = [s / d for s, d in zip(
            sparse["full_portfolio_wall_seconds_by_export_pair"],
            dense["full_portfolio_wall_seconds_by_export_pair"])]
        engines[engine]["paired_sparse_over_dense_full_wall"] = ratios
        engines[engine]["sparse_robust_toy_wall_improvement"] = (
            dense["complete_correct"] and sparse["complete_correct"] and
            all(ratio < .9 for ratio in ratios))
        engines[engine]["pair_order_reverses_ranking"] = (ratios[0] - 1) * (ratios[1] - 1) < 0
    complete_all = all(engines[name][rep]["complete_correct"]
                       for name in SOLVERS for rep in REPS)
    robust_all = all(engines[name]["sparse_robust_toy_wall_improvement"] for name in SOLVERS)
    return {"decision": ("SPARSE_TOY_NEXT_RUNG_PREFERENCE" if complete_all and robust_all
                         else "CORRECT_BUT_MIXED_OR_INCONCLUSIVE" if complete_all
                         else "CENSORED_OR_FAILED"),
            "proof_checked_unsat": False,
            "common_preflight_plus_two_base_load_wall_seconds": common,
            "exports_full_four_panel_wall_seconds": [r["wall_seconds"] for r in exports],
            "actual_campaign_process_wall_seconds": receipt["process_wall_seconds"],
            "archived_semantic_audit_wall_seconds": archived_audit,
            "engines": engines}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("smoke", type=Path)
    parser.add_argument("panel", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    assert not args.output.exists(), "do not overwrite analysis"
    archived = json.loads((Path(__file__).resolve().parent.parent /
                           "rotated_s3_sparse_cnf_20260925/evidence/final/summary.json").read_text())
    audit = {rep + "_replay": archived["cold_children"][rep + "_replay"]["wall_seconds"]
             for rep in ("dense", "sparse")}
    result = summarize(json.loads((args.smoke / "result.json").read_text()),
                       json.loads((args.panel / "receipt.json").read_text()),
                       json.loads((args.panel / "result.json").read_text()), audit)
    args.output.write_text(json.dumps(result, sort_keys=True, separators=(",", ":")) + "\n")
    print(json.dumps({"decision": result["decision"], "output": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
