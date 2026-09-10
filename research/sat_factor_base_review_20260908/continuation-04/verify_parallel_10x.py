"""Fail-closed verification of the fixed batch-10 end-to-end speedup."""
from pathlib import Path
from statistics import median
import json
import math

ROOT = Path(__file__).resolve().parent
CONTROL_ROOT = ROOT.parent / "continuation-03"
PAIRS = [(1, 4242), (2, 9001), (3, 17001), (4, 30001)]
CAP = 10_000
BATCH = 10
THREADS = 10


def load(path):
    lines = path.read_text().splitlines()
    assert len(lines) == 1
    return json.loads(lines[0])


def external_real_seconds(path):
    fields = path.read_text().splitlines()[0].split()
    assert fields[0] == "real" and len(fields) == 2
    return float(fields[1])


pairs = []
for seed, secret in PAIRS:
    optimized = load(ROOT / f"optimized-seed-{seed}-batch-{BATCH}.jsonl")
    control = load(CONTROL_ROOT / f"frobenius-seed-{seed}-cap-{CAP}.jsonl")
    for result, arm, collapse, stop, relation_unknowns in [
        (optimized, "optimized", True, True, 4),
        (control, "frobenius", False, False, 8),
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
        assert result["relation_unknowns"] == relation_unknowns
        assert result["verified_unknown_scalar_recovery"] is True
        assert result["direct_relation"] is False
        assert result["sat_invalid_models"] == 0
        assert result["linear_solve_attempts"] == 1
        assert result["end_to_end_ns"] >= result["driver_ns"]
        assert result["driver_ns"] >= result["relation_collection_ns"]
    assert optimized["relation_batch_size"] == BATCH
    assert optimized["rayon_threads"] == THREADS
    assert optimized["relation_batches"] in (2, 3)
    assert optimized["relations"] >= 5
    # Continuation-03 predates the batch-size report field; its driver
    # was serial by construction. A present field must still read one.
    assert control.get("relation_batch_size", 1) == 1
    assert control["relations"] == 10

    optimized_seconds = optimized["end_to_end_ns"] / 1e9
    control_seconds = control["end_to_end_ns"] / 1e9
    optimized_external = external_real_seconds(
        ROOT / f"optimized-seed-{seed}-batch-{BATCH}-time.txt"
    )
    control_external = external_real_seconds(
        CONTROL_ROOT / f"frobenius-seed-{seed}-cap-{CAP}-time.txt"
    )
    assert abs(optimized_external - optimized_seconds) / optimized_seconds < 0.02
    assert abs(control_external - control_seconds) / control_seconds < 0.02
    speedup = control_seconds / optimized_seconds
    assert speedup > 1.0
    pairs.append(
        {
            "run_seed": seed,
            "secret": secret,
            "optimized_seconds": optimized_seconds,
            "control_seconds": control_seconds,
            "paired_speedup": speedup,
            "optimized_external_real_seconds": optimized_external,
            "control_external_real_seconds": control_external,
            "optimized_trials": optimized["trials"],
            "control_trials": control["trials"],
            "optimized_relations": optimized["relations"],
            "control_relations": control["relations"],
            "optimized_sat_conflicts": optimized["sat_conflicts"],
            "control_sat_conflicts": control["sat_conflicts"],
            "optimized_sat_unknowns": optimized["sat_unknowns"],
            "control_sat_unknowns": control["sat_unknowns"],
        }
    )

optimized_total = sum(pair["optimized_seconds"] for pair in pairs)
control_total = sum(pair["control_seconds"] for pair in pairs)
speedups = [pair["paired_speedup"] for pair in pairs]
holdout_pairs = [pair for pair in pairs if pair["run_seed"] != 3]
holdout_optimized_total = sum(pair["optimized_seconds"] for pair in holdout_pairs)
holdout_control_total = sum(pair["control_seconds"] for pair in holdout_pairs)
summary = {
    "kind": "koblitz_parallel_e2e_10x_speedup",
    "curve": "K_1/F_(2^19)",
    "factor_base_points": 152,
    "control": "serial, eight Frobenius columns, fixed two-relation surplus",
    "optimized": "ten-thread batch-10 SAT collection, four signed columns, verified rank stopping",
    "shared_per_target_conflict_budget": CAP,
    "relation_batch_size": BATCH,
    "rayon_threads": THREADS,
    "completed_pairs": len(pairs),
    "batch_selection_seed": 3,
    "holdout_seeds": [1, 2, 4],
    "all_unknown_scalars_recovered_and_verified": True,
    "optimized_total_seconds": optimized_total,
    "control_total_seconds": control_total,
    "aggregate_speedup": control_total / optimized_total,
    "holdout_optimized_total_seconds": holdout_optimized_total,
    "holdout_control_total_seconds": holdout_control_total,
    "holdout_aggregate_speedup": holdout_control_total / holdout_optimized_total,
    "minimum_paired_speedup": min(speedups),
    "median_paired_speedup": median(speedups),
    "geometric_mean_paired_speedup": math.prod(speedups) ** (1 / len(speedups)),
    "optimized_total_trials": sum(pair["optimized_trials"] for pair in pairs),
    "control_total_trials": sum(pair["control_trials"] for pair in pairs),
    "optimized_total_relations": sum(pair["optimized_relations"] for pair in pairs),
    "control_total_relations": sum(pair["control_relations"] for pair in pairs),
    "optimized_total_sat_conflicts": sum(pair["optimized_sat_conflicts"] for pair in pairs),
    "control_total_sat_conflicts": sum(pair["control_sat_conflicts"] for pair in pairs),
    "pairs": pairs,
    "scope": "paired complete toy n=19 index-calculus wall time including factor-base construction, parallel SAT relation collection, modular linear algebra, and verified unknown-scalar recovery",
}
assert summary["aggregate_speedup"] >= 10.0
assert summary["holdout_aggregate_speedup"] >= 10.0
assert summary["geometric_mean_paired_speedup"] >= 10.0
assert summary["minimum_paired_speedup"] >= 5.0
(ROOT / "parallel-10x-summary.json").write_text(json.dumps(summary, indent=2) + "\n")

# The fixed policy is selected from a bounded batch-size screen on seed 3.
batch_rows = {
    size: load(ROOT / f"optimized-seed-3-batch-{size}.jsonl")
    for size in (8, 10, 12, 16)
}
assert all(row["verified_unknown_scalar_recovery"] for row in batch_rows.values())
assert min(batch_rows, key=lambda size: batch_rows[size]["end_to_end_ns"]) == BATCH

# Preserve the failed lower-cutoff arm as operational evidence, not a speed result.
failed = load(ROOT / "optimized-seed-4-cap-5000-batch-10.jsonl")
assert failed["verified_unknown_scalar_recovery"] is False
assert failed["relations"] == 4 and failed["trials"] == 64
assert failed["sat_invalid_models"] == 0

messages = [
    "Four fixed-policy parallel runs and four serial controls recovered and independently reverified their unknown scalars.",
    f"Parallel optimized total {optimized_total:.6f}s versus control {control_total:.6f}s: {summary['aggregate_speedup']:.6f}x aggregate end-to-end speedup.",
    f"Excluding batch-selection seed 3, holdout seeds 1, 2, and 4 give {summary['holdout_aggregate_speedup']:.6f}x aggregate speedup.",
    f"Geometric mean {summary['geometric_mean_paired_speedup']:.6f}x, median {summary['median_paired_speedup']:.6f}x, minimum {summary['minimum_paired_speedup']:.6f}x.",
    f"Batch 10 was fastest among completing seed-3 batches 8, 10, 12, and 16; every admitted run used ten Rayon threads and a 10000-conflict cutoff.",
    "The seed-4 5000-conflict arm stopped at the 64-trial bound with four relations and is retained only as an operational failure.",
]
(ROOT / "verification.txt").write_text("\n".join(messages) + "\n")
print("\n".join(messages))
