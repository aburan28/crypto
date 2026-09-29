"""Scaling screen for compiled reservoir-guided current-plus-addend schedules.

The compiler uses public coefficient geometry only.  Its selected hot-addend
indices form a target-independent schedule that can be replayed without online
planning for any public target in the same prime-order group.
"""

from hashlib import sha256
from math import comb, isqrt
from pathlib import Path
import json

from check_hiaa_collision_scaling import is_prime_64


SPECS = (
    {"label": "fixed", "order": 65_537, "budget": 427, "trials": 16, "reservoir_size": 32, "seed_offset": 0},
    {"label": "fixed", "order": 1_048_583, "budget": 1_536, "trials": 16, "reservoir_size": 32, "seed_offset": 1_000},
    {"label": "scaled", "order": 1_048_583, "budget": 1_536, "trials": 4, "reservoir_size": 128, "seed_offset": 1_000},
    {"label": "fixed", "order": 4_194_301, "budget": 3_072, "trials": 8, "reservoir_size": 32, "seed_offset": 2_000},
    {"label": "scaled", "order": 4_194_301, "budget": 3_072, "trials": 2, "reservoir_size": 256, "seed_offset": 2_000},
)
PARENT_WINDOW = 64
CHOICES = 4
BASE_SEED = 20261500


def enc(value):
    return int(value).to_bytes(16, "big")


def prf(domain, *values):
    return int.from_bytes(
        sha256(domain.encode() + b"".join(enc(v) for v in values)).digest(),
        "big",
    )


def add(left, right, order):
    return ((left[0] + right[0]) % order, (left[1] + right[1]) % order)


def slope(left, right, order):
    dx = (left[0] - right[0]) % order
    if not dx:
        return None
    return ((left[1] - right[1]) * pow(dx, -1, order)) % order


def prefix_metrics(events, budget, order):
    covered = set()
    curve = []
    for index, point in enumerate(events):
        for prior in events[:index]:
            value = slope(point, prior, order)
            if value is not None:
                covered.add(value)
        if index >= 2:
            curve.append(len(covered))
    assert len(curve) == budget + 1
    numerator = sum(order - curve[time] for time in range(budget))
    return {
        "covered_slopes": curve[-1],
        "censored_mean_numerator": numerator,
        "censored_mean_additions": numerator / order,
        "coverage_curve_sha256": sha256(
            json.dumps(curve, separators=(",", ":")).encode()
        ).hexdigest(),
    }


def addend_indices(event_count, choices, seed, step):
    low = max(0, event_count - PARENT_WINDOW)
    width = event_count - low
    values = set()
    lane = 0
    while len(values) < min(choices, width):
        values.add(low + prf("compiled-addend", seed, step, lane) % width)
        lane += 1
    return tuple(sorted(values))


def bottom_k_indices(events, size, seed):
    ranked = sorted(
        (prf("compiled-reservoir", seed, index), index)
        for index in range(len(events))
    )
    return tuple(sorted(index for _, index in ranked[:size]))


def compile_or_control(order, budget, reservoir_size, seed, guided):
    events = [(0, 0), (1, 0), (0, 1)]
    materialized = set(events)
    sketch = set()
    addends = []
    slope_queries = 0
    insert_attempts = 0
    candidate_evaluations = 0
    for step in range(budget):
        current = events[-1]
        witnesses = [
            events[index]
            for index in bottom_k_indices(events, reservoir_size, seed)
        ]
        ranked = []
        for addend_index in addend_indices(len(events), CHOICES, seed, step):
            candidate = add(current, events[addend_index], order)
            slopes = set()
            if guided:
                candidate_evaluations += 1
                for witness in witnesses:
                    value = slope(candidate, witness, order)
                    slope_queries += 1
                    if value is not None:
                        slopes.add(value)
                gain = len(slopes - sketch)
            else:
                gain = 0
            exploration = prf(
                "compiled-explore", seed, step, candidate[0], candidate[1]
            )
            ranked.append((
                gain,
                candidate not in materialized,
                exploration,
                -addend_index,
                addend_index,
                candidate,
                slopes,
            ))
        best = max(ranked)
        addend_index, candidate, selected_slopes = best[-3], best[-2], best[-1]
        assert candidate == add(events[-1], events[addend_index], order)
        events.append(candidate)
        materialized.add(candidate)
        addends.append(addend_index)
        if not guided:
            selected_slopes = set()
            for witness in witnesses:
                value = slope(candidate, witness, order)
                slope_queries += 1
                if value is not None:
                    selected_slopes.add(value)
        sketch.update(selected_slopes)
        insert_attempts += len(selected_slopes)

    metrics = prefix_metrics(events, budget, order)
    offsets = [(3 + step) - index for step, index in enumerate(addends)]
    assert max(offsets) <= PARENT_WINDOW
    return {
        "guided": guided,
        "seed": seed,
        "metrics": metrics,
        "compiler": {
            "candidate_evaluations": candidate_evaluations,
            "sampled_slope_queries": slope_queries,
            "sampled_slope_insert_attempts": insert_attempts,
            "final_exact_sampled_sketch_entries": len(sketch),
            "planning_archive_reads": 0,
        },
        "compiled_schedule": {
            "addend_indices": addends,
            "backward_offsets": offsets,
            "maximum_backward_offset": max(offsets),
            "bits_per_offset": (PARENT_WINDOW - 1).bit_length(),
            "packed_schedule_bytes_ceiling": (
                len(offsets) * (PARENT_WINDOW - 1).bit_length() + 7
            ) // 8,
            "sha256": sha256(
                json.dumps(offsets, separators=(",", ":")).encode()
            ).hexdigest(),
        },
        "events_sha256": sha256(
            json.dumps(events, separators=(",", ":")).encode()
        ).hexdigest(),
    }


def sign_p(wins, losses):
    n = wins + losses
    if not n:
        return 1.0
    tail = sum(comb(n, k) for k in range(min(wins, losses) + 1)) / 2**n
    return min(1.0, 2 * tail)


def main():
    all_rows = []
    summaries = []
    for spec_index, spec in enumerate(SPECS):
        order, budget = spec["order"], spec["budget"]
        assert is_prime_64(order)
        rows = []
        for trial in range(spec["trials"]):
            seed = BASE_SEED + spec["seed_offset"] + trial
            guided = compile_or_control(order, budget, spec["reservoir_size"], seed, True)
            control = compile_or_control(order, budget, spec["reservoir_size"], seed, False)
            row = {
                "order": order,
                "sqrt_floor": isqrt(order),
                "budget": budget,
                "label": spec["label"],
                "reservoir_size": spec["reservoir_size"],
                "trial": trial,
                "seed": seed,
                "guided": guided,
                "control": control,
                "endpoint_delta_guided_minus_control": (
                    guided["metrics"]["covered_slopes"]
                    - control["metrics"]["covered_slopes"]
                ),
                "censored_addition_delta_guided_minus_control": (
                    guided["metrics"]["censored_mean_additions"]
                    - control["metrics"]["censored_mean_additions"]
                ),
            }
            rows.append(row)
            all_rows.append(row)
        ew = sum(r["endpoint_delta_guided_minus_control"] > 0 for r in rows)
        et = sum(r["endpoint_delta_guided_minus_control"] == 0 for r in rows)
        el = len(rows) - ew - et
        pw = sum(r["censored_addition_delta_guided_minus_control"] < 0 for r in rows)
        pt = sum(r["censored_addition_delta_guided_minus_control"] == 0 for r in rows)
        pl = len(rows) - pw - pt
        summaries.append({
            **spec,
            "guided_endpoint_higher_tie_lower": [ew, et, el],
            "guided_prefix_better_tie_worse": [pw, pt, pl],
            "mean_endpoint_delta": sum(r["endpoint_delta_guided_minus_control"] for r in rows) / len(rows),
            "mean_censored_addition_delta": sum(r["censored_addition_delta_guided_minus_control"] for r in rows) / len(rows),
            "endpoint_sign_p_descriptive": sign_p(ew, el),
            "prefix_sign_p_descriptive": sign_p(pw, pl),
            "guided_mean_censored_additions": sum(r["guided"]["metrics"]["censored_mean_additions"] for r in rows) / len(rows),
            "control_mean_censored_additions": sum(r["control"]["metrics"]["censored_mean_additions"] for r in rows) / len(rows),
        })
    report = {
        "evidence_type": "finite_compiled_reservoir_guided_walk_scaling_screen",
        "specifications": SPECS,
        "parent_window": PARENT_WINDOW,
        "choices": CHOICES,
        "reservoir_size_is_per_specification": True,
        "target_independent_compilation": True,
        "online_planning_after_compilation": False,
        "summaries": summaries,
        "rows": all_rows,
        "limitations": [
            "Compiled schedules are finite per-order constants and do not prove an asymptotic family.",
            "The matched control is not an optimized fixed-r rho implementation.",
            "Compilation cost is target-independent but must be charged or amortized over a declared target count.",
            "The projected transition remains time/history dependent and lacks rho coalescence.",
            "No novelty, hardware, private-target, or ECDLP speedup claim follows.",
        ],
        "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    output = Path(__file__).with_name("compiled_walk_scaling_checks.json")
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"summaries": summaries, "output": str(output)}, indent=2))


if __name__ == "__main__":
    main()
