#!/usr/bin/env python3
"""Fixed-order toy SAT stage costs; no ECDLP/rho extrapolation."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

SOLVERS = ("cryptominisat5", "kissat", "cadical")


def sum_cost(rows):
    return sum(row["wall_seconds"] + row["query_setup_wall_seconds"] +
               row["verifier_wall_seconds"] for row in rows)


def analyze(smoke, panel):
    assert smoke["mode"] == "smoke" and panel["mode"] == "panel"
    assert len(smoke["entries"]) == 6 and len(panel["entries"]) == 96
    by = {(row["id"], row["solver"]): row for row in panel["entries"]}
    common = panel["preflight_wall_seconds"] + panel["base_load_wall_seconds"]
    engines = {}
    for name in SOLVERS:
        arm_smoke = [row for row in smoke["entries"] if row["solver"] == name]
        rows = [by[(f"Q{q}T{t}", name)] for q in range(8) for t in range(4)]
        fixed_q = []
        for q in range(8):
            blocks = rows[4*q:4*q+4]
            statuses = [row["verdict"] for row in blocks]
            first = next((t for t, row in enumerate(blocks) if row["verdict"] == "SAT"), None)
            projected = ("POSITIVE" if first is not None and
                         all(row["verdict"] == "UNSAT" for row in blocks[:first]) else
                         "NEGATIVE" if all(status == "UNSAT" for status in statuses) else
                         "CENSORED")
            fixed_q.append({"q": q, "target_class": "planted" if q < 4 else "negative",
                            "branch_statuses": statuses, "projected_status": projected,
                            "first_witness_t": first if projected == "POSITIVE" else None,
                            "first_witness_wall_seconds":
                                common + sum_cost(blocks[:first+1]) if projected == "POSITIVE" else None,
                            "all_four_wall_seconds": common + sum_cost(blocks),
                            "all_four_child_wall_seconds": sum(row["wall_seconds"] for row in blocks),
                            "max_child_wall_seconds": max(row["wall_seconds"] for row in blocks)})
        complete = (all(row["pass"] for row in arm_smoke) and
                    all(row["verdict"] in ("SAT", "UNSAT") for row in rows) and
                    all(row["projected_status"] == ("POSITIVE" if row["q"] < 4 else "NEGATIVE")
                        for row in fixed_q))
        engines[name] = {"smoke_pass": all(row["pass"] for row in arm_smoke),
                         "complete_and_consistent": complete,
                         "branch_sat": sum(row["verdict"] == "SAT" for row in rows),
                         "branch_unsat_proof_unchecked": sum(row["verdict"] == "UNSAT" for row in rows),
                         "branch_censored_or_invalid": sum(row["verdict"] not in ("SAT", "UNSAT") for row in rows),
                         "child_wall_sum_seconds": sum(row["wall_seconds"] for row in rows),
                         "query_setup_wall_sum_seconds": sum(row["query_setup_wall_seconds"] for row in rows),
                         "model_verify_wall_sum_seconds": sum(row["verifier_wall_seconds"] for row in rows),
                         "charged_portfolio_wall_seconds": common + sum_cost(rows),
                         "max_child_wall_seconds": max(row["wall_seconds"] for row in rows),
                         "q": fixed_q}
    return {"decision": "ADMITTED_TOY_SOLVER_STAGE" if all(
                arm["complete_and_consistent"] for arm in engines.values()) else "CENSORED_OR_BLOCKED",
            "proof_checked_unsat": False,
            "common_preflight_plus_base_load_wall_seconds": common,
            "actual_96_child_process_wall_seconds": panel["full_process_wall_seconds"],
            "engine_costs": engines}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("smoke", type=Path)
    parser.add_argument("panel", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    assert not args.output.exists()
    result = analyze(json.loads((args.smoke / "result.json").read_text()),
                     json.loads((args.panel / "result.json").read_text()))
    args.output.write_text(json.dumps(result, sort_keys=True, separators=(",", ":")) + "\n")
    print(json.dumps({"decision": result["decision"], "output": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
