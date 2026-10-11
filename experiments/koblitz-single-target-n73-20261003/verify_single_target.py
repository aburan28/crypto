#!/usr/bin/env python3
"""Independent GF(2^73) replay and one-target accounting check."""
import hashlib
import json
import platform
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
FROZEN = ROOT / "frozen"
RUNS = ROOT / "runs"
def load_json(path):
    return json.loads(path.read_text())


def load_one_jsonl(path):
    rows = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    if len(rows) != 1:
        raise ValueError("expected exactly one row in %s, got %d" % (path, len(rows)))
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
    infinity = None

    def inverse(value):
        if not value:
            raise ZeroDivisionError
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
                return infinity
            slope = x ^ mul(y, inverse(x))
            xx = square(slope) ^ slope ^ point_a
            yy = square(x) ^ mul(slope ^ 1, xx)
            return xx, yy
        slope = mul(y ^ v, inverse(x ^ u))
        xx = square(slope) ^ slope ^ x ^ u ^ point_a
        yy = mul(slope, x ^ xx) ^ xx ^ y
        return xx, yy

    def scalar_mul(k, point):
        result = infinity
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
        "subgroup_order_replay_is_infinity": scalar_mul(order, generator) is infinity,
    }


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    protocol = load_json(ROOT / "protocol.json")
    base_run_id = protocol["run_id"].removesuffix("R1")
    if len(sys.argv) == 1:
        run_id = protocol["run_id"]
    elif len(sys.argv) == 2 and sys.argv[1] in {"R1", "R2", "R3"}:
        run_id = base_run_id + sys.argv[1]
    else:
        raise SystemExit("usage: verify_single_target.py [R1|R2|R3]")
    run_dir = RUNS / run_id
    fixture = load_json(FROZEN / "fixture.json")
    target_rows = [line for line in (FROZEN / "target_points.jsonl").read_text().splitlines() if line.strip()]
    if len(target_rows) != 1:
        raise ValueError("frozen target input must contain exactly one point")
    ic = load_one_jsonl(run_dir / "ic.jsonl")
    ic_summary = load_one_jsonl(run_dir / "ic_summary.jsonl")
    rho = load_one_jsonl(run_dir / "rho.jsonl")
    ic_execution = load_json(run_dir / "ic_execution.json")
    rho_execution = load_json(run_dir / "rho_execution.json")
    rho_peak_rss_bytes = int(rho_execution["peak_rss_bytes"])
    ic_peak_rss_bytes = (
        int(ic_summary["peak_rss_bytes"])
        if isinstance(ic_summary.get("peak_rss_bytes"), int)
        else int(ic_execution["peak_rss_bytes"])
    )

    n = int(fixture["n"])
    low_terms = [int(term) for term in fixture["field_modulus_low_terms"]]
    modulus = (1 << n) | 1
    for term in low_terms:
        modulus |= 1 << term
    generator = tuple(int(value) for value in fixture["generator"])
    target = tuple(int(value) for value in fixture["public_target"])
    scalar = int(fixture["fixture_scalar_validation_only"])
    order = int(fixture["subgroup_order"])
    cofactor = int(fixture["cofactor"])
    group_order = int(fixture["curve_order"])
    curve_trace = int(fixture["curve_trace_over_f2n"])

    checks = curve_replay(n, modulus, int(fixture["a"]), 1, generator, target, scalar, order)
    checks["fixture_order_factorization_matches"] = order * cofactor == group_order
    checks["fixture_trace_matches_order"] = (1 << n) + 1 - group_order == curve_trace
    checks["ic_uses_frozen_Q"] = tuple(int(v) for v in ic["published_q"]) == target
    checks["rho_uses_frozen_Q"] = tuple(int(v) for v in rho["published_q"]) == target
    checks["rho_uses_frozen_G"] = tuple(int(v) for v in rho["generator"]) == generator
    checks["ic_recovers_fixture_scalar"] = int(ic["recovered_scalar"]) == scalar
    checks["rho_recovers_fixture_scalar"] = int(rho["recovered_fixture_scalar"]) == scalar
    checks["ic_in_process_verification_passed"] = ic.get("group_verified") is True
    checks["rho_in_process_verification_passed"] = rho.get("verified") is True
    checks["one_ic_online_target"] = ic_summary["targets"] == 1 and ic_summary["targets_solved"] == 1
    run_record = load_json(run_dir / "run.json")
    checks["one_rho_online_target"] = (
        run_id in {base_run_id + tag for tag in ("R1", "R2", "R3")}
        and run_record["run_id"] == run_id
        and run_record["target"] == protocol["paired_target"]
        and run_record["target_count"] == 1
        and protocol["target_count"] == 1
        and rho["fixture_index"] == 0
    )
    checks["known_answer_not_sent_to_ic"] = ic.get("published_fixture_scalar") is None
    checks["known_answer_not_sent_to_rho"] = (
        rho.get("published_fixture_scalar") is None
        and not rho_execution["known_answer_environment_variable_present"]
    )
    checks["ic_process_succeeded"] = ic_execution["return_code"] == 0
    checks["rho_process_succeeded"] = rho_execution["return_code"] == 0
    checks["no_explicit_pair_table"] = ic_summary["pair_table_entries"] == 0
    checks["no_edge_selectors"] = ic_summary["edge_selectors"] == 0
    declared_cap = protocol["resource_cap_bytes_per_arm"]
    checks["ic_peak_rss_within_declared_cap"] = declared_cap is None or ic_peak_rss_bytes <= int(declared_cap)
    checks["rho_peak_rss_within_declared_cap"] = declared_cap is None or rho_peak_rss_bytes <= int(declared_cap)

    ic_phases = {
        "target_query": float(ic["target_query_ms"]),
        "target_PDP": float(ic["target_pdp_ms"]),
        "target_relation_check": float(ic["target_relation_check_ms"]),
        "target_descent": float(ic["target_descent_ms"]),
        "target_recovery_check": float(ic["target_recovery_check_ms"]),
    }
    ic_online = float(ic["target_ms"])
    ic_phase_sum = sum(ic_phases.values())
    rho_phases = {
        "target_specific_jump_setup": float(rho["setup_ms"]),
        "target_walk": float(rho["walk_ms"]),
        "scalar_recovery_check": float(rho["validation_ms"]),
    }
    rho_online = sum(rho_phases.values())
    checks["ic_online_phase_sum_matches"] = abs(ic_phase_sum - ic_online) <= max(0.02, ic_online * 1e-8)
    checks["rho_online_phase_sum_matches"] = abs(rho_online - float(rho["setup_ms"] + rho["walk_ms"] + rho["validation_ms"])) < 1e-9
    checks["rho_used_supplied_public_point"] = rho["target_input_kind"] == "public_point"
    checks["rho_fixture_generation_excluded"] = float(rho["target_generation_ms_excluded"]) == 0.0

    online_speedup = rho_online / ic_online
    passed = all(checks.values())
    receipt = {
        "kind": "independent_n73_single_target_replay",
        "status": "PASS" if passed else "FAIL",
        "scope": "Standalone standard-library GF(2^73) curve replay, same-point pairing, one-target accounting, phase sums, and resource checks; does not independently remeasure producer runtimes.",
        "candidate_id": protocol["candidate_id"],
        "workload_id": protocol["workload_id"],
        "run_id": run_id,
        "curve": {
            "n": n,
            "a": int(fixture["a"]),
            "b": 1,
            "field_modulus_low_terms": low_terms,
            "curve_order": str(group_order),
            "curve_trace_over_f2n": str(curve_trace),
            "subgroup_order": str(order),
            "cofactor": str(cofactor),
            "frobenius_eigenvalue_mod_subgroup_order": fixture["frobenius_eigenvalue_mod_subgroup_order"],
        },
        "public_input": {
            "G": [str(value) for value in generator],
            "Q": [str(value) for value in target],
            "recovered_scalar": str(scalar),
            "target_count": 1,
        },
        "checks": checks,
        "timing_ms": {
            "ic_online": ic_online,
            "ic_phase_sum": ic_phase_sum,
            "ic_phases": ic_phases,
            "ic_peak_rss_bytes": ic_peak_rss_bytes,
            "rho_online": rho_online,
            "rho_peak_rss_bytes": rho_peak_rss_bytes,
            "rho_phase_sum": sum(rho_phases.values()),
            "rho_phases": rho_phases,
            "online_speedup": online_speedup,
        },
        "resource_policy": {
            "finite_memory_cap_bytes_per_arm": declared_cap,
            "peak_rss_source": "wait4(2) per-process ru_maxrss; IC summary may be null when host process inspection is restricted",
            "claim": "same host, sequential arms; peak RSS measured; no finite cap",
        },
        "inputs_sha256": {
            "fixture": sha256(FROZEN / "fixture.json"),
            "target_points": sha256(FROZEN / "target_points.jsonl"),
            "ic": sha256(run_dir / "ic.jsonl"),
            "ic_summary": sha256(run_dir / "ic_summary.jsonl"),
            "ic_execution": sha256(run_dir / "ic_execution.json"),
            "rho": sha256(run_dir / "rho.jsonl"),
            "rho_execution": sha256(run_dir / "rho_execution.json"),
        },
        "environment": {"python": sys.version.split()[0], "platform": platform.platform()},
    }
    out = run_dir / "independent_replay.json"
    out.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
