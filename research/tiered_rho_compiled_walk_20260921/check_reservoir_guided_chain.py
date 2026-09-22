"""RAM-reservoir coverage-guided legal addition-chain screen.

The planner keeps a bounded witness reservoir and hot parent window in RAM.
It never scans or randomly samples the disk archive for planning.  A sampled
slope sketch guides selection; exact full coverage is offline evaluation only.
"""

from fractions import Fraction
from hashlib import sha256
from pathlib import Path
from statistics import median
import json

from check_charged_polynomial_construction import P, Transcript, add, direction_profile
from check_multiline_stripes import prefix_metrics, strongest_three_line_control


BUDGET = 427
PARENT_WINDOW = 64
CHOICE_COUNTS = (4, 8, 16)
RESERVOIR_SIZES = (16, 32, 64)
RESERVOIR_MODES = ("hot", "bottom_k")
SEEDS = tuple(20261000 + i for i in range(8))
COEFFICIENT_RECORD_BYTES = 64


def prf(domain, *values):
    material = domain.encode() + b"".join(
        int(v).to_bytes(8, "big") for v in values
    )
    return int.from_bytes(sha256(material).digest(), "big")


def candidate_pairs(event_count, seed, step, choices):
    low = max(0, event_count - PARENT_WINDOW)
    width = event_count - low
    pairs = set()
    last = event_count - 1
    for choice in range(choices):
        left = low + prf("reservoir-left", seed, step, choice) % width
        right = low + prf("reservoir-right", seed, step, choice) % width
        pairs.add(tuple(sorted((left, right))))
        if choice % 8 == 0:
            pairs.add(tuple(sorted((last, right))))
    return tuple(sorted(pairs))


def slope(candidate, witness):
    dx = (candidate[0] - witness[0]) % P
    if not dx:
        return None
    return ((candidate[1] - witness[1]) * pow(dx, -1, P)) % P


def reservoir_indices(events, mode, size, seed):
    if mode == "hot":
        return tuple(range(max(0, len(events) - size), len(events)))
    assert mode == "bottom_k"
    ranked = sorted(
        (prf("reservoir-priority", seed, index), index)
        for index in range(len(events))
    )
    return tuple(sorted(index for _, index in ranked[:size]))


def run(mode, choices, reservoir_size, seed):
    transcript = Transcript()
    materialized = set(transcript.events)
    sketch = set()
    selected_parents = []
    candidate_pairs_scored = 0
    slope_queries = 0
    slope_insert_attempts = 0
    reservoir_changes = 0
    previous_reservoir = ()
    for step in range(BUDGET):
        witness_indices = reservoir_indices(
            transcript.events, mode, reservoir_size, seed
        )
        reservoir_changes += int(witness_indices != previous_reservoir)
        previous_reservoir = witness_indices
        witnesses = [transcript.events[index] for index in witness_indices]
        pairs = candidate_pairs(len(transcript.events), seed, step, choices)
        candidate_pairs_scored += len(pairs)
        ranked = []
        seen = set()
        for left, right in pairs:
            candidate = add(transcript.events[left], transcript.events[right])
            if candidate in seen:
                continue
            seen.add(candidate)
            slopes = set()
            for witness in witnesses:
                value = slope(candidate, witness)
                slope_queries += 1
                if value is not None:
                    slopes.add(value)
            novel = slopes - sketch
            exploration = prf(
                "reservoir-explore", seed, step, candidate[0], candidate[1]
            )
            ranked.append((
                len(novel),
                candidate not in materialized,
                exploration,
                min(left, right),
                -max(left, right),
                -left,
                -right,
                left,
                right,
                candidate,
                slopes,
            ))
        best = max(ranked)
        left, right, candidate, slopes = best[-4], best[-3], best[-2], best[-1]
        result = transcript.charge_add(
            transcript.events[left], transcript.events[right],
            "reservoir_coverage_guided", f"step_{step}_parents_{left}_{right}",
        )
        assert result == candidate
        materialized.add(result)
        selected_parents.append((left, right))
        sketch.update(slopes)
        slope_insert_attempts += len(slopes)

    profile = direction_profile(transcript.events)
    prefix = prefix_metrics(transcript.events)
    assert profile["unsigned_useful_finite_slopes"] == prefix["covered_slopes"]
    assert all(
        (3 + step - min(pair)) <= PARENT_WINDOW
        for step, pair in enumerate(selected_parents)
    )
    return {
        "reservoir_mode": mode,
        "choices": choices,
        "reservoir_size": reservoir_size,
        "seed": seed,
        "charged_additions": BUDGET,
        "profile": profile,
        "prefix_metrics": prefix,
        "planner_cost": {
            "candidate_pairs_scored": candidate_pairs_scored,
            "sampled_slope_queries": slope_queries,
            "sampled_slope_insert_attempts": slope_insert_attempts,
            "final_exact_sampled_sketch_entries": len(sketch),
            "reservoir_changes": reservoir_changes,
            "reservoir_ram_bytes": reservoir_size * COEFFICIENT_RECORD_BYTES,
            "hot_parent_ram_bytes": PARENT_WINDOW * COEFFICIENT_RECORD_BYTES,
            "hypothetical_16_bit_filter_bytes": 2 * slope_insert_attempts,
            "planning_archive_reads": 0,
            "batch_inversion_groups_if_one_per_step": BUDGET,
        },
        "selected_parent_indices": selected_parents,
        "events_sha256": sha256(
            json.dumps(transcript.events, separators=(",", ":")).encode()
        ).hexdigest(),
    }


def aggregate(rows):
    endpoints = [row["profile"]["unsigned_useful_finite_slopes"] for row in rows]
    means = [Fraction(row["prefix_metrics"]["exact_censored_mean_min_T_budget"]) for row in rows]
    return {
        "trials": len(rows),
        "endpoint_min_median_max": [min(endpoints), median(endpoints), max(endpoints)],
        "censored_mean_min_median_max": [str(min(means)), str(median(means)), str(max(means))],
        "best_endpoint_seed": max(rows, key=lambda r: r["profile"]["unsigned_useful_finite_slopes"])["seed"],
        "best_prefix_seed": min(rows, key=lambda r: Fraction(r["prefix_metrics"]["exact_censored_mean_min_T_budget"]))["seed"],
    }


def main():
    rows = []
    summaries = []
    for mode in RESERVOIR_MODES:
        for choices in CHOICE_COUNTS:
            for size in RESERVOIR_SIZES:
                current = [run(mode, choices, size, seed) for seed in SEEDS]
                rows.extend(current)
                summaries.append({
                    "reservoir_mode": mode,
                    "choices": choices,
                    "reservoir_size": size,
                    **aggregate(current),
                })
    three_final, three_prefix = strongest_three_line_control(BUDGET)
    best_endpoint = max(rows, key=lambda r: r["profile"]["unsigned_useful_finite_slopes"])
    best_prefix = min(rows, key=lambda r: Fraction(r["prefix_metrics"]["exact_censored_mean_min_T_budget"]))
    report = {
        "evidence_type": "finite_ram_reservoir_coverage_guided_chain_screen",
        "field_prime": P,
        "budget": BUDGET,
        "planner_grid": {
            "parent_window": PARENT_WINDOW,
            "choice_counts": CHOICE_COUNTS,
            "reservoir_sizes": RESERVOIR_SIZES,
            "reservoir_modes": RESERVOIR_MODES,
            "public_seeds": SEEDS,
            "planning_archive_reads": 0,
            "target_independent": True,
            "planner_uses_full_slope_coverage": False,
        },
        "summaries": summaries,
        "best_endpoint": best_endpoint,
        "best_prefix": best_prefix,
        "three_line_final_control": three_final,
        "three_line_prefix_control": three_prefix,
        "comparisons": {
            "best_endpoint_minus_three_line": (
                best_endpoint["profile"]["unsigned_useful_finite_slopes"]
                - three_final["all_materialized_profile"]["unsigned_useful_finite_slopes"]
            ),
            "best_prefix_minus_three_line_numerator": (
                best_prefix["prefix_metrics"]["censored_mean_numerator"]
                - three_prefix["prefix_metrics"]["censored_mean_numerator"]
            ),
        },
        "trials": rows,
        "limitations": [
            "The finite planner uses an exact sampled-slope set; approximate filter false positives can change selections.",
            "The slope reservoir is bounded but supplies no asymptotic candidate-ordering guarantee.",
            "Coefficient arithmetic and sketch costs are counted but not converted to wall time.",
            "The experiment covers one small public field.",
            "No novelty, ECDLP speedup, private target, hardware, or production claim follows.",
        ],
        "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    output = Path(__file__).with_name("reservoir_guided_chain_checks.json")
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({
        "summaries": summaries,
        "best_endpoint": best_endpoint["profile"]["unsigned_useful_finite_slopes"],
        "best_prefix": best_prefix["prefix_metrics"]["exact_censored_mean_min_T_budget"],
        "three_line_endpoint": three_final["all_materialized_profile"]["unsigned_useful_finite_slopes"],
        "three_line_prefix": three_prefix["prefix_metrics"]["exact_censored_mean_min_T_budget"],
        "output": str(output),
    }, indent=2))


if __name__ == "__main__":
    main()
