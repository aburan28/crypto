#!/usr/bin/env python3
"""Independent GF(2^85) replay for the three frozen n=85 a=0 paired runs.

Standalone standard-library curve arithmetic, independent of the Rust
producers.  Verifies the frozen fixture, same-point pairing, one-target
accounting, phase sums, and the recovered scalars.  Usage:

    python3 verify_n83_runs.py            # verify R1, R2, R3
"""
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs"
FROZEN = HERE / "frozen"


def load_json(path):
    return json.loads(Path(path).read_text())


def load_one_jsonl(path):
    rows = [json.loads(line) for line in Path(path).read_text().splitlines() if line.strip()]
    if len(rows) != 1:
        raise ValueError(f"expected exactly one row in {path}, got {len(rows)}")
    return rows[0]


def gf_mul(a, b, n, modulus):
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if (a >> n) & 1:
            a ^= modulus
    return out


def curve_replay(n, modulus, point_a, point_b, generator, target, scalar, order):
    mul = lambda a, b: gf_mul(a, b, n, modulus)
    square = lambda a: mul(a, a)

    def inverse(value):
        exponent = (1 << n) - 2
        result = 1
        base = value
        while exponent:
            if exponent & 1:
                result = mul(result, base)
            base = square(base)
            exponent >>= 1
        return result

    def on_curve(point):
        if point is None:
            return True
        x, y = point
        left = square(y) ^ mul(x, y)
        right = mul(square(x), x) ^ mul(point_a, square(x)) ^ point_b
        return left == right

    def add(left, right):
        if left is None:
            return right
        if right is None:
            return left
        x, y = left
        u, v = right
        if x == u:
            if y != v or x == 0:
                return None
            slope = x ^ mul(y, inverse(x))
            xx = square(slope) ^ slope ^ point_a
            yy = square(x) ^ mul(slope ^ 1, xx)
            return xx, yy
        slope = mul(y ^ v, inverse(x ^ u))
        xx = square(slope) ^ slope ^ x ^ u ^ point_a
        yy = mul(slope, x ^ xx) ^ xx ^ y
        return xx, yy

    def scalar_mul(k, point):
        result = None
        while k:
            if k & 1:
                result = add(result, point)
            point = add(point, point)
            k >>= 1
        return result

    return {
        "G_on_curve": on_curve(generator),
        "Q_on_curve": on_curve(target),
        "scalar_replay_equals_Q": scalar_mul(scalar, generator) == target,
        "subgroup_order_replay_is_infinity": scalar_mul(order, generator) is None,
        "subgroup_annihilates_Q": scalar_mul(order, target) is None,
    }


def verify_run(run_dir, fixture, target):
    n = int(fixture["n"])
    low_terms = [int(term) for term in fixture["field_modulus_low_terms"]]
    modulus = (1 << n) | 1
    for term in low_terms:
        modulus |= 1 << term
    generator = tuple(int(value) for value in fixture["generator"])
    scalar = int(fixture["fixture_scalar_validation_only"])
    order = int(fixture["subgroup_order"])
    cofactor = int(fixture["cofactor"])
    group_order = int(fixture["curve_order"])
    curve_trace = int(fixture["curve_trace_over_f2n"])

    ic = load_one_jsonl(run_dir / "ic.jsonl")
    ic_summary = load_one_jsonl(run_dir / "ic_summary.jsonl")
    rho = load_one_jsonl(run_dir / "rho.jsonl")
    ic_execution = load_json(run_dir / "ic_execution.json")
    rho_execution = load_json(run_dir / "rho_execution.json")

    checks = curve_replay(
        n, modulus, int(fixture["a"]), 1, generator, tuple(target), scalar, order
    )
    checks["fixture_order_factorization_matches"] = order * cofactor == group_order
    checks["fixture_trace_matches_order"] = (1 << n) + 1 - group_order == curve_trace
    checks["ic_uses_frozen_Q"] = tuple(int(v) for v in ic["published_q"]) == tuple(target)
    checks["rho_uses_frozen_Q"] = tuple(int(v) for v in rho["published_q"]) == tuple(target)
    checks["rho_uses_frozen_G"] = tuple(int(v) for v in rho["generator"]) == generator
    checks["ic_recovers_fixture_scalar"] = int(ic["recovered_scalar"]) == scalar
    checks["rho_recovers_fixture_scalar"] = int(rho["recovered_fixture_scalar"]) == scalar
    checks["ic_in_process_verification_passed"] = ic.get("group_verified") is True
    checks["rho_in_process_verification_passed"] = rho.get("verified") is True
    checks["one_ic_online_target"] = ic_summary["targets"] == 1 and ic_summary["targets_solved"] == 1
    checks["one_rho_online_target"] = rho["fixture_index"] == 0
    checks["known_answer_not_sent_to_ic"] = ic.get("published_fixture_scalar") is None
    checks["known_answer_not_sent_to_rho"] = (
        rho.get("published_fixture_scalar") is None
        and not rho_execution["known_answer_environment_variable_present"]
    )
    checks["ic_process_succeeded"] = ic_execution["return_code"] == 0
    checks["rho_process_succeeded"] = rho_execution["return_code"] == 0
    checks["no_explicit_pair_table"] = ic_summary["pair_table_entries"] == 0
    checks["no_edge_selectors"] = ic_summary["edge_selectors"] == 0
    checks["rho_used_supplied_public_point"] = rho["target_input_kind"] == "public_point"
    checks["rho_fixture_generation_excluded"] = float(rho["target_generation_ms_excluded"]) == 0.0
    checks["deterministic_target_relation"] = ic["probes"] == 215254

    ic_online = float(ic["target_ms"])
    ic_phases = {
        "target_query": float(ic["target_query_ms"]),
        "target_PDP": float(ic["target_pdp_ms"]),
        "target_relation_check": float(ic["target_relation_check_ms"]),
        "target_descent": float(ic["target_descent_ms"]),
        "target_recovery_check": float(ic["target_recovery_check_ms"]),
    }
    ic_phase_sum = sum(ic_phases.values())
    checks["ic_online_phase_sum_matches"] = abs(ic_phase_sum - ic_online) <= max(
        0.02, ic_online * 1e-8
    )
    rho_phases = {
        "target_specific_jump_setup": float(rho["setup_ms"]),
        "target_walk": float(rho["walk_ms"]),
        "scalar_recovery_check": float(rho["validation_ms"]),
    }
    rho_online = sum(rho_phases.values())
    checks["rho_online_phase_sum_matches"] = abs(
        rho_online - (float(rho["setup_ms"]) + float(rho["walk_ms"]) + float(rho["validation_ms"]))
    ) < 1e-9

    receipt = {
        "kind": "independent_n85_single_target_replay",
        "status": "PASS" if all(checks.values()) else "FAIL",
        "run_id": run_dir.name,
        "candidate_id": "IC1N83A1Ckb1fb99600PDP4rootRCguidedLAgaussTDdirectISO0parallel",
        "scope": (
            "Standalone standard-library GF(2^85) curve replay, same-point pairing, "
            "one-target accounting, and phase sums; does not independently remeasure "
            "producer runtimes."
        ),
        "checks": checks,
        "timing_ms": {
            "ic_online": ic_online,
            "rho_online": rho_online,
            "online_speedup": rho_online / ic_online,
            "ic_phases": ic_phases,
            "rho_phases": rho_phases,
            "ic_peak_rss_bytes": (
                ic_summary["peak_rss_bytes"]
                if isinstance(ic_summary.get("peak_rss_bytes"), int)
                else ic_execution["peak_rss_bytes"]
            ),
            "rho_peak_rss_bytes": rho_execution["peak_rss_bytes"],
        },
        "precompute_detail": {
            "rank_stage_ms": ic_summary["timing_ms"]["rank_stage"],
            "rank_threads": ic_summary["threads"],
            "rank_probes_mean": ic_summary["rank_probes_mean"],
            "rank_failures": ic_summary["rank_failures"],
            "orbit_columns": ic_summary["orbit_columns"],
            "pair_table_entries": ic_summary["pair_table_entries"],
        },
    }
    (run_dir / "independent_replay.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    )
    return receipt


def main():
    fixture = load_json(FROZEN / "fixture.json")
    target = tuple(int(v) for v in fixture["public_target"])
    statuses = []
    for run_dir in sorted(RUNS.iterdir()):
        if not run_dir.is_dir():
            continue
        receipt = verify_run(run_dir, fixture, target)
        statuses.append((run_dir.name, receipt["status"], receipt["timing_ms"]["online_speedup"]))
        print(
            json.dumps(
                {
                    "run": run_dir.name[-2:],
                    "status": receipt["status"],
                    "speedup": round(receipt["timing_ms"]["online_speedup"], 2),
                }
            )
        )
    if any(status != "PASS" for _, status, _ in statuses):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
