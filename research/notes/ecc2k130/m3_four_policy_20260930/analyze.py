#!/usr/bin/env python3
"""Recompute exact m3 support and the frozen four-policy decision."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import verify


HERE = Path(__file__).resolve().parent
VARIANTS = ("original", "transported", "descendant_native", "pullback")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def analyze(evidence: Path) -> dict:
    lock = json.loads((HERE / "FROZEN_REPLAY.json").read_text())
    assert sha(evidence / "result.json") == lock["result_sha256"]
    assert sha(evidence / "replay.json") == lock["failed_replay_sha256"]
    result = json.loads((evidence / "result.json").read_text())
    replay = json.loads((evidence / "replay_repaired.json").read_text())
    assert result["status"] == "PASS_PANEL" and replay["status"] == "PASS"
    assert replay["result_sha256"] == lock["result_sha256"]
    assert replay["repair_lock_sha256"] == sha(HERE / "FROZEN_REPLAY.json")
    assert replay["cells_replayed"] == 16 and replay["cases_replayed"] == 8192

    verify.restore_bare_curve()
    F = verify.pilot.FastGF2m(21, verify.pilot.IRR)
    source = verify.pilot.Koblitz(F, 0, 1)
    twist = verify.pilot.Koblitz(F, 1, 1)
    lines, selected, _ = verify.pilot.order_seven_lines(F, twist)
    phi = verify.pilot.BinaryVeluMap.from_generator(
        source, twist, lines[selected]["generator"], 7)
    leaf = phi.codomain
    G = tuple(result["challenge"]["G"])
    reps, by_point, _ = verify.independent_orbits(source, G)
    allowed = {}
    for label, orbits in (("A", set(reps[:5])), ("B", set(reps[5:]))):
        source_points = {P for P, orbit in by_point.items() if orbit in orbits}
        assert len(source_points) == 210
        leaf_points = {phi(P) for P in source_points}
        assert len(leaf_points) == 210
        allowed[label] = {False: source_points, True: leaf_points}

    cells, comparisons = [], []
    for seed in result["variants"]:
        for label in ("A", "B"):
            row_by_name = {}
            for name in VARIANTS:
                is_leaf = name in ("transported", "descendant_native")
                E = leaf if is_leaf else source
                base = [tuple(P) for P in result["bases"][seed][name]]
                triples, _ = verify.brute_triples(E, base)
                support = set(triples)
                assert len(support) == result[
                    "triple_support_cardinality"][seed][label][name]
                eligible = support & allowed[label][is_leaf]
                first_rows = []
                for target in sorted(eligible):
                    k, a, b = triples[target][0]
                    coefficients = [0] * 8
                    for index in (a, b, k):
                        coefficients[index] += 1
                    first_rows.append(tuple(coefficients))
                row_rank, _ = verify.independent_linear_system(
                    [(list(row), 0) for row in first_rows], 8)
                variant = result["variants"][seed][label][name]
                cost = result["cold_cost_to_rank_or_512"][seed][label][name]
                modular = sum(variant["mod_r_ops"].values())
                cell = {
                    "seed": int(seed), "holdout": label, "policy": name,
                    "hits_out_of_512": variant["hits"],
                    "rank": variant["rank"],
                    "first_full_rank_attempt": variant["first_full_rank_attempt"],
                    "distinct_triple_sums_all_421": len(support),
                    "eligible_support_out_of_210": len(eligible),
                    "distinct_first_witness_rows": len(set(first_rows)),
                    "first_witness_base_row_rank_out_of_8": row_rank,
                    "cold_field_mul": cost["mul"],
                    "cold_field_sqr": cost["sqr"],
                    "cold_inversion_calls": cost["inv"],
                    "cold_group_add": cost["group_add"],
                    "modular_row_ops": modular,
                    "cold_cpu_ns_host_diagnostic": cost["cpu_ns"],
                    "verified": variant["verified"],
                }
                assert cell["rank"] == 9 and cell["verified"]
                cells.append(cell)
                row_by_name[name] = cell
            original, transported, native, pullback = (
                row_by_name[name] for name in VARIANTS)
            comparison = {
                "seed": int(seed), "holdout": label,
                "native_minus_transported_cold_field_mul": (
                    native["cold_field_mul"] - transported["cold_field_mul"]),
                "native_div_original_cold_field_mul": (
                    native["cold_field_mul"] / original["cold_field_mul"]),
                "pullback_div_original_cold_field_mul": (
                    pullback["cold_field_mul"] / original["cold_field_mul"]),
                "native_vs_transported_lower_mul": (
                    native["cold_field_mul"] < transported["cold_field_mul"]),
                "native_vs_transported_no_hit_regression": (
                    native["hits_out_of_512"] >= transported["hits_out_of_512"]),
                "native_vs_transported_no_rank_delay": (
                    native["first_full_rank_attempt"] <= transported[
                        "first_full_rank_attempt"]),
                "native_vs_transported_no_modular_burden": (
                    native["modular_row_ops"] <= transported["modular_row_ops"]),
            }
            comparisons.append(comparison)
    assert len(cells) == 16 and len(comparisons) == 4
    decision = "NATIVE_ADVANTAGE" if all(all(value for key, value in row.items()
        if key.startswith("native_vs_transported_")) for row in comparisons) else (
        "NO_NATIVE_ADVANTAGE")
    return {
        "schema": "ecc2k130-degree7-m3-four-policy-analysis-v1",
        "run_url": "https://github.com/aburan28/crypto/actions/runs/36722040881",
        "producer_source_head": result["source_head"],
        "result_sha256": lock["result_sha256"],
        "failed_replay_sha256": lock["failed_replay_sha256"],
        "repaired_replay_sha256": sha(evidence / "replay_repaired.json"),
        "repair_lock_sha256": sha(HERE / "FROZEN_REPLAY.json"),
        "panel_status": result["status"], "replay_status": replay["status"],
        "run_host": {"platform": result["platform"],
                     "python": result["python"],
                     "wall_seconds": result["wall_seconds"],
                     "cpu_seconds": result["cpu_seconds"],
                     "peak_rss_bytes": result["peak_rss_bytes"]},
        "cells": cells, "native_comparisons": comparisons,
        "decision": decision,
        "ECC2K_130_PDP_yield": None,
        "full_ECDLP_cost": None,
        "rho_crossover": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists(), "never overwrite a saved analysis"
    data = analyze(args.evidence)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"decision": data["decision"],
                      "cells": len(data["cells"])}, sort_keys=True))


if __name__ == "__main__":
    main()
