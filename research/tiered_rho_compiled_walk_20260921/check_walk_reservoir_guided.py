"""Single-register reservoir-guided history walk with paired controls.

Every online step is current + archived_addend.  The active register therefore
forms one history-dependent walk, while a bounded RAM reservoir and sampled
slope sketch choose among public hot-addend proposals.  Public-synthetic only.
"""

from fractions import Fraction
from hashlib import sha256
from math import comb
from pathlib import Path
import json

from check_charged_polynomial_construction import Transcript, add, direction_profile
from check_multiline_stripes import prefix_metrics, strongest_three_line_control
from check_reservoir_guided_chain import BUDGET, PARENT_WINDOW, prf, slope


CONFIGS = (
    {"choices": 4, "reservoir_size": 32},
    {"choices": 4, "reservoir_size": 64},
    {"choices": 8, "reservoir_size": 32},
    {"choices": 8, "reservoir_size": 64},
    {"choices": 16, "reservoir_size": 64},
)
SEEDS = tuple(20261300 + i for i in range(32))


def addend_indices(event_count, choices, seed, step):
    low = max(0, event_count - PARENT_WINDOW)
    width = event_count - low
    values = set()
    lane = 0
    while len(values) < min(choices, width):
        values.add(low + prf("walk-addend", seed, step, lane) % width)
        lane += 1
    return tuple(sorted(values))


def bottom_k_indices(events, size, seed):
    ranked = sorted(
        (prf("walk-reservoir", seed, index), index)
        for index in range(len(events))
    )
    return tuple(sorted(index for _, index in ranked[:size]))


def run(guided, choices, reservoir_size, seed):
    transcript = Transcript()
    materialized = set(transcript.events)
    sketch = set()
    selected_addends = []
    slope_queries = 0
    insert_attempts = 0
    candidate_evaluations = 0
    for step in range(BUDGET):
        current_index = len(transcript.events) - 1
        witnesses = [
            transcript.events[index]
            for index in bottom_k_indices(transcript.events, reservoir_size, seed)
        ]
        ranked = []
        for addend_index in addend_indices(
            len(transcript.events), choices, seed, step
        ):
            candidate = add(
                transcript.events[current_index], transcript.events[addend_index]
            )
            slopes = set()
            if guided:
                candidate_evaluations += 1
                for witness in witnesses:
                    value = slope(candidate, witness)
                    slope_queries += 1
                    if value is not None:
                        slopes.add(value)
                gain = len(slopes - sketch)
            else:
                gain = 0
            exploration = prf(
                "walk-explore", seed, step, candidate[0], candidate[1]
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
        result = transcript.charge_add(
            transcript.events[current_index], transcript.events[addend_index],
            "walk_reservoir_guided" if guided else "walk_exploration_control",
            f"step_{step}_addend_{addend_index}",
        )
        assert result == candidate
        materialized.add(result)
        selected_addends.append(addend_index)
        if not guided:
            selected_slopes = set()
            for witness in witnesses:
                value = slope(result, witness)
                slope_queries += 1
                if value is not None:
                    selected_slopes.add(value)
        sketch.update(selected_slopes)
        insert_attempts += len(selected_slopes)

    profile = direction_profile(transcript.events)
    prefix = prefix_metrics(transcript.events)
    assert all(
        (3 + step - index) <= PARENT_WINDOW
        for step, index in enumerate(selected_addends)
    )
    return {
        "guided": guided,
        "choices": choices,
        "reservoir_size": reservoir_size,
        "seed": seed,
        "profile": profile,
        "prefix_metrics": prefix,
        "planner_cost": {
            "candidate_evaluations": candidate_evaluations,
            "sampled_slope_queries": slope_queries,
            "slope_insert_attempts": insert_attempts,
            "final_exact_sampled_sketch_entries": len(sketch),
            "planning_archive_reads": 0,
            "active_walk_registers": 1,
            "hot_addend_records": PARENT_WINDOW,
            "reservoir_records": reservoir_size,
        },
        "selected_addend_indices": selected_addends,
        "events_sha256": sha256(
            json.dumps(transcript.events, separators=(",", ":")).encode()
        ).hexdigest(),
    }


def sign_p(wins, losses):
    n = wins + losses
    if not n:
        return 1.0
    tail = sum(comb(n, k) for k in range(min(wins, losses) + 1)) / 2**n
    return min(1.0, 2 * tail)


def summarize(rows):
    ew = sum(r["guided_endpoint"] > r["control_endpoint"] for r in rows)
    et = sum(r["guided_endpoint"] == r["control_endpoint"] for r in rows)
    el = len(rows) - ew - et
    pw = sum(r["guided_prefix_numerator"] < r["control_prefix_numerator"] for r in rows)
    pt = sum(r["guided_prefix_numerator"] == r["control_prefix_numerator"] for r in rows)
    pl = len(rows) - pw - pt
    ed = sum(r["guided_endpoint"] - r["control_endpoint"] for r in rows)
    pd = sum(r["guided_prefix_numerator"] - r["control_prefix_numerator"] for r in rows)
    return {
        "paired_trials": len(rows),
        "guided_endpoint_higher_tie_lower": [ew, et, el],
        "guided_prefix_better_tie_worse": [pw, pt, pl],
        "mean_endpoint_difference": ed / len(rows),
        "mean_censored_addition_difference": pd / (len(rows) * 65537),
        "endpoint_sign_p_descriptive_unadjusted": sign_p(ew, el),
        "prefix_sign_p_descriptive_unadjusted": sign_p(pw, pl),
    }


def main():
    all_rows = []
    summaries = []
    for config in CONFIGS:
        rows = []
        for seed in SEEDS:
            guided = run(True, config["choices"], config["reservoir_size"], seed)
            control = run(False, config["choices"], config["reservoir_size"], seed)
            row = {
                **config,
                "seed": seed,
                "guided_endpoint": guided["profile"]["unsigned_useful_finite_slopes"],
                "control_endpoint": control["profile"]["unsigned_useful_finite_slopes"],
                "guided_prefix_numerator": guided["prefix_metrics"]["censored_mean_numerator"],
                "control_prefix_numerator": control["prefix_metrics"]["censored_mean_numerator"],
                "guided_prefix": guided["prefix_metrics"]["exact_censored_mean_min_T_budget"],
                "control_prefix": control["prefix_metrics"]["exact_censored_mean_min_T_budget"],
                "guided_planner_cost": guided["planner_cost"],
                "guided_events_sha256": guided["events_sha256"],
                "control_events_sha256": control["events_sha256"],
            }
            rows.append(row)
            all_rows.append(row)
        summaries.append({**config, **summarize(rows)})

    three_final, three_prefix = strongest_three_line_control(BUDGET)
    report = {
        "evidence_type": "finite_single_register_walk_reservoir_guidance_control",
        "field_prime": 65537,
        "budget": BUDGET,
        "configs": CONFIGS,
        "public_seeds": SEEDS,
        "transition": "current plus one selected archived hot addend",
        "projected_state_is_not_fixed_function_of_current_group_element": True,
        "summaries": summaries,
        "rows": all_rows,
        "three_line_final_control": three_final,
        "three_line_prefix_control": three_prefix,
        "limitations": [
            "The growing archive and coverage planner destroy fixed projected-state coalescence.",
            "The paired exploration control is not a fully optimized fixed-r rho implementation.",
            "Exact sampled membership is used; Bloom false positives require separate audit.",
            "Coefficient planning costs are counted but not converted to wall time.",
            "One public field supplies no asymptotic, novelty, or ECDLP speedup result.",
        ],
        "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    output = Path(__file__).with_name("walk_reservoir_guided_checks.json")
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"summaries": summaries, "output": str(output)}, indent=2))


if __name__ == "__main__":
    main()
