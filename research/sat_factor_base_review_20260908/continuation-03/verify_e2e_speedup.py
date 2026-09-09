"""Verify matched end-to-end runs before computing the speedup aggregate."""
from pathlib import Path
from statistics import median
import json
import math

ROOT = Path(__file__).resolve().parent
PAIRS = [(1, 4242), (2, 9001), (3, 17001), (4, 30001)]
CAP = 10_000


def load(name):
    lines = (ROOT / name).read_text().splitlines()
    assert len(lines) == 1
    return json.loads(lines[0])


def external_real_seconds(name):
    fields = (ROOT / name).read_text().splitlines()[0].split()
    assert fields[0] == "real" and len(fields) == 2
    return float(fields[1])


rows = []
for seed, secret in PAIRS:
    optimized = load(f"optimized-seed-{seed}-cap-{CAP}.jsonl")
    control = load(f"frobenius-seed-{seed}-cap-{CAP}.jsonl")
    for result, arm, collapse, stop, unknowns, relations in [
        (optimized, "optimized", True, True, 4, 5),
        (control, "frobenius", False, False, 8, 10),
    ]:
        assert result["kind"] == "koblitz_selected_factor_base_e2e"
        assert result["arm"] == arm
        assert result["n"] == 19 and result["a"] == 1
        assert result["secret"] == secret
        assert result["sampler_seed"] == 0xE2E0000000000000 ^ seed
        assert result["per_target_conflict_budget"] == CAP
        assert result["factor_base_points"] == 152
        assert result["frobenius_orbits"] == 8
        assert result["signed_frobenius_orbits"] == 4
        assert result["collapse_negation"] is collapse
        assert result["stop_on_verified_rank"] is stop
        assert result["relation_unknowns"] == unknowns
        assert result["relations"] == relations
        assert result["verified_unknown_scalar_recovery"] is True
        assert result["direct_relation"] is False
        assert result["sat_invalid_models"] == 0
        assert result["linear_solve_attempts"] == 1
        assert result["end_to_end_ns"] >= result["driver_ns"]
        assert result["driver_ns"] >= (
            result["relation_collection_ns"] + result["linear_algebra_ns"]
        )
    assert optimized["sampler_seed"] == control["sampler_seed"]
    assert optimized["secret"] == control["secret"]
    optimized_seconds = optimized["end_to_end_ns"] / 1e9
    control_seconds = control["end_to_end_ns"] / 1e9
    optimized_external = external_real_seconds(f"optimized-seed-{seed}-cap-{CAP}-time.txt")
    control_external = external_real_seconds(f"frobenius-seed-{seed}-cap-{CAP}-time.txt")
    assert abs(optimized_external - optimized_seconds) / optimized_seconds < 0.02
    assert abs(control_external - control_seconds) / control_seconds < 0.02
    speedup = control_seconds / optimized_seconds
    assert speedup > 1.0
    rows.append(
        {
            "run_seed": seed,
            "secret": secret,
            "execution_order": "control_then_optimized" if seed == 4 else "optimized_then_control",
            "optimized_seconds": optimized_seconds,
            "control_seconds": control_seconds,
            "optimized_external_real_seconds": optimized_external,
            "control_external_real_seconds": control_external,
            "paired_speedup": speedup,
            "optimized_trials": optimized["trials"],
            "control_trials": control["trials"],
            "optimized_sat_conflicts": optimized["sat_conflicts"],
            "control_sat_conflicts": control["sat_conflicts"],
            "optimized_sat_unknowns": optimized["sat_unknowns"],
            "control_sat_unknowns": control["sat_unknowns"],
        }
    )

optimized_total = sum(row["optimized_seconds"] for row in rows)
control_total = sum(row["control_seconds"] for row in rows)
paired = [row["paired_speedup"] for row in rows]
summary = {
    "kind": "koblitz_selected_factor_base_e2e_speedup",
    "curve": "K_1/F_(2^19)",
    "factor_base_points": 152,
    "control": "eight Frobenius relation columns, fixed two-relation surplus",
    "optimized": "four signed Frobenius columns, verified rank-aware stopping",
    "shared_per_target_conflict_budget": CAP,
    "completed_pairs": len(rows),
    "all_unknown_scalars_recovered_and_verified": True,
    "optimized_total_seconds": optimized_total,
    "control_total_seconds": control_total,
    "aggregate_speedup": control_total / optimized_total,
    "minimum_paired_speedup": min(paired),
    "median_paired_speedup": median(paired),
    "geometric_mean_paired_speedup": math.prod(paired) ** (1 / len(paired)),
    "optimized_total_trials": sum(row["optimized_trials"] for row in rows),
    "control_total_trials": sum(row["control_trials"] for row in rows),
    "optimized_total_sat_conflicts": sum(row["optimized_sat_conflicts"] for row in rows),
    "control_total_sat_conflicts": sum(row["control_sat_conflicts"] for row in rows),
    "pairs": rows,
    "scope": "paired complete toy n=19 index-calculus runs including factor-base construction, SAT relation collection, modular linear algebra, and verified unknown-scalar recovery",
}
assert summary["aggregate_speedup"] > 2.0
assert summary["minimum_paired_speedup"] > 1.0
(ROOT / "e2e-speedup-summary.json").write_text(json.dumps(summary, indent=2) + "\n")

cap_rows = [
    load("optimized-seed-1-cap-5000.jsonl"),
    load("optimized-seed-1-cap-10000.jsonl"),
    load("optimized-seed-1-cap-20000.jsonl"),
]
assert all(row["verified_unknown_scalar_recovery"] for row in cap_rows)
assert min(cap_rows, key=lambda row: row["end_to_end_ns"])["per_target_conflict_budget"] == CAP

messages = [
    "Four matched pairs recovered and independently reverified unknown scalars 4242, 9001, 17001, and 30001; seed 4 ran control first.",
    f"Optimized total {optimized_total:.6f}s versus control {control_total:.6f}s: {summary['aggregate_speedup']:.6f}x aggregate end-to-end speedup.",
    f"Every pair improved: minimum {min(paired):.6f}x, median {median(paired):.6f}x, geometric mean {summary['geometric_mean_paired_speedup']:.6f}x.",
    f"The optimized arm used {summary['optimized_total_trials']} trials and {summary['optimized_total_sat_conflicts']} conflicts; the control used {summary['control_total_trials']} trials and {summary['control_total_sat_conflicts']} conflicts.",
    "The seed-1 cutoff sweep completed at 5000, 10000, and 20000 conflicts; 10000 was the fastest completing operating point.",
]
(ROOT / "verification.txt").write_text("\n".join(messages) + "\n")
print("\n".join(messages))
