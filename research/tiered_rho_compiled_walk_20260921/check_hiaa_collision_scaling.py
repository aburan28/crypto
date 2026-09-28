"""Collision-only scaling screen for the history-indexed archive accumulator.

This is a deterministic public-synthetic experiment in additive prime-order
cyclic groups.  It accepts no external target and makes no speedup claim.
"""

from hashlib import sha256
from math import comb, isqrt
from pathlib import Path
import json


SPECS = (
    {"order": 65_537, "budget": 768, "trials": 512, "r_values": (16, 64)},
    {"order": 1_048_583, "budget": 3_072, "trials": 256, "r_values": (16, 64, 256)},
    {"order": 16_777_259, "budget": 12_288, "trials": 128, "r_values": (64, 256)},
)
EXPERIMENT_SEED = 20260921


def enc(value):
    return int(value).to_bytes(16, "big", signed=False)


def prf(domain, *values):
    material = domain.encode() + b"".join(enc(v) for v in values)
    return int.from_bytes(sha256(material).digest(), "big")


def is_prime_64(n):
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for p in small:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 325, 9375, 28178, 450775, 9780504, 1795265022):
        if a % n == 0:
            continue
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def target_for(order, spec_index, trial):
    # Avoid 0, 1, and -1 so the three public initial records are distinct.
    return 2 + prf("public-target", EXPERIMENT_SEED, spec_index, trial) % (order - 3)


def seed_for(spec_index, r, trial):
    return prf("public-trace-seed", EXPERIMENT_SEED, spec_index, r, trial) % (1 << 64)


def add(left, right, order):
    return ((left[0] + right[0]) % order, (left[1] + right[1]) % order)


def scalar(record, target, order):
    return (record[0] + target * record[1]) % order


def observe(record, index, additions, target, order, first_by_scalar):
    value = scalar(record, target, order)
    prior = first_by_scalar.get(value)
    if prior is None:
        first_by_scalar[value] = (index, record)
        return None, False
    prior_index, prior_record = prior
    delta_a = (record[0] - prior_record[0]) % order
    delta_b = (record[1] - prior_record[1]) % order
    if delta_b == 0:
        assert delta_a == 0
        return None, True
    recovered = -delta_a * pow(delta_b, -1, order) % order
    assert recovered == target
    assert scalar(prior_record, recovered, order) == value
    assert scalar(record, recovered, order) == value
    return {
        "addition": additions,
        "prior_index": prior_index,
        "current_index": index,
    }, False


def run_trace(kind, order, target, budget, r, seed):
    assert kind in ("history", "fixed_r")
    assert budget >= r - 2
    events = [(0, 0), (1, 0), (0, 1)]
    step_indices = [1, 2]
    first_by_scalar = {}
    degenerate = 0
    for index, record in enumerate(events):
        collision, is_degenerate = observe(
            record, index, 0, target, order, first_by_scalar
        )
        assert collision is None
        degenerate += int(is_degenerate)
    assert degenerate == 0

    for build_step in range(r - 2):
        left = step_indices[prf("build-left", seed, build_step) % len(step_indices)]
        right = step_indices[prf("build-right", seed, build_step) % len(step_indices)]
        record = add(events[left], events[right], order)
        events.append(record)
        index = len(events) - 1
        step_indices.append(index)
        additions = build_step + 1
        collision, is_degenerate = observe(
            record, index, additions, target, order, first_by_scalar
        )
        degenerate += int(is_degenerate)
        if collision is not None:
            return collision["addition"], additions, degenerate

    fixed_indices = tuple(step_indices)
    available = list(step_indices)
    current = events[step_indices[-1]]
    additions = r - 2
    for _online_step in range(budget - additions):
        raw = prf("select-state", seed, scalar(current, target, order))
        choices = available if kind == "history" else fixed_indices
        selected = choices[raw % len(choices)]
        record = add(current, events[selected], order)
        events.append(record)
        index = len(events) - 1
        additions += 1
        collision, is_degenerate = observe(
            record, index, additions, target, order, first_by_scalar
        )
        degenerate += int(is_degenerate)
        if collision is not None:
            return collision["addition"], additions, degenerate
        current = record
        if kind == "history":
            available.append(index)
    return None, additions, degenerate


def exact_two_sided_sign_p(wins, losses):
    n = wins + losses
    if n == 0:
        return 1.0
    tail = sum(comb(n, k) for k in range(min(wins, losses) + 1)) / (2**n)
    return min(1.0, 2.0 * tail)


def summarize(rows, budget):
    wins = sum(row["history_first"] < row["fixed_first"] for row in rows)
    ties = sum(row["history_first"] == row["fixed_first"] for row in rows)
    losses = len(rows) - wins - ties
    history_censored = sum(row["history_collision"] is None for row in rows)
    fixed_censored = sum(row["fixed_collision"] is None for row in rows)
    history_restricted_total = sum(row["history_first"] for row in rows)
    fixed_restricted_total = sum(row["fixed_first"] for row in rows)
    return {
        "paired_trials": len(rows),
        "history_earlier_tie_later": [wins, ties, losses],
        "history_censored": history_censored,
        "fixed_censored": fixed_censored,
        "restricted_mean_history_additions": history_restricted_total / len(rows),
        "restricted_mean_fixed_additions": fixed_restricted_total / len(rows),
        "restricted_mean_ratio_history_over_fixed": (
            history_restricted_total / fixed_restricted_total
        ),
        "paired_restricted_total_difference_history_minus_fixed": (
            history_restricted_total - fixed_restricted_total
        ),
        "two_sided_exact_sign_test_p_descriptive_unadjusted": exact_two_sided_sign_p(
            wins, losses
        ),
        "censor_value": budget + 1,
    }


def main():
    for spec in SPECS:
        assert is_prime_64(spec["order"])
        assert spec["budget"] > 2 * isqrt(spec["order"])

    all_rows = []
    summaries = []
    for spec_index, spec in enumerate(SPECS):
        order = spec["order"]
        budget = spec["budget"]
        for r in spec["r_values"]:
            rows = []
            for trial in range(spec["trials"]):
                target = target_for(order, spec_index, trial)
                seed = seed_for(spec_index, r, trial)
                history, history_done, history_degenerate = run_trace(
                    "history", order, target, budget, r, seed
                )
                fixed, fixed_done, fixed_degenerate = run_trace(
                    "fixed_r", order, target, budget, r, seed
                )
                row = {
                    "order": order,
                    "sqrt_floor": isqrt(order),
                    "budget": budget,
                    "r": r,
                    "trial": trial,
                    "public_target": target,
                    "seed": seed,
                    "history_collision": history,
                    "fixed_collision": fixed,
                    "history_first": budget + 1 if history is None else history,
                    "fixed_first": budget + 1 if fixed is None else fixed,
                    "history_additions_executed": history_done,
                    "fixed_additions_executed": fixed_done,
                    "history_degenerate_repeats_before_stop": history_degenerate,
                    "fixed_degenerate_repeats_before_stop": fixed_degenerate,
                }
                rows.append(row)
                all_rows.append(row)
            summaries.append({
                "order": order,
                "sqrt_floor": isqrt(order),
                "budget": budget,
                "r": r,
                **summarize(rows, budget),
            })

    report = {
        "evidence_type": "deterministic_public_synthetic_collision_scaling_screen",
        "scope": {
            "group_model": "additive_prime_order_cyclic_group",
            "selector": "state_hash",
            "external_target_interface": False,
            "experiment_seed": EXPERIMENT_SEED,
            "specifications": SPECS,
            "all_group_additions_table_builds_archive_probes_and_selectors_charged": True,
            "candidate_and_control_share_table_build_rule_and_exact_archive": True,
        },
        "interpretation_limits": [
            "The deterministic suite is finite and was not sampled from a claimed deployment distribution.",
            "The unadjusted sign-test value is descriptive and does not prove a speedup.",
            "Restricted means replace censoring by budget plus one and are not uncensored expected times.",
            "The experiment models group additions and exact records, not NVMe latency or wall time.",
            "No novelty, asymptotic improvement, ECDLP break, or production claim follows.",
        ],
        "summaries": summaries,
        "rows": all_rows,
        "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    output = Path(__file__).with_name("hiaa_collision_scaling.json")
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({
        "paired_trials": len(all_rows),
        "summary_rows": len(summaries),
        "output": str(output),
    }, indent=2))


if __name__ == "__main__":
    main()
