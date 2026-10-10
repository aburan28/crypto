"""Charged finite multi-stripe generic-DLP screen.

All computations are public and take place in the coefficient plane F_65537^2.
The checker charges binary scalar setup, line starts, and every sequential line
addition.  It excludes the vertical direction from useful slopes.  It is not an
elliptic-curve solver and does not accept an external target.
"""

from decimal import Decimal, getcontext
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
from statistics import stdev
import json
import platform

from check_charged_polynomial_construction import (
    A,
    B,
    O,
    P,
    Transcript,
    add,
    direction,
    direction_profile,
    random_control,
)


RANDOM_TRIALS = 64

# Explicit finite four-stripe candidate from a limited deterministic offset
# screen. Its displayed points improve the uncharged three-line count; this
# checker tests whether the lead survives legal setup accounting.
# Each tuple is (x-coordinate, y-offset, y-stride), all modulo P.
CANDIDATE_LINES = (
    (0, 0, 1),
    (1, 96, 96),
    (2, (-97) % P, (-97) % P),
    (3, 19969, (-146) % P),
)
CANDIDATE_POINTS_PER_LINE = 96


def phase_counts(transcript):
    result = {}
    for operation in transcript.operations:
        phase = operation["phase"]
        result[phase] = result.get(phase, 0) + 1
    return result


def prefix_metrics(events):
    """Exact discovery curve for a uniform finite slope target."""
    assert events[:3] == [O, A, B]
    slopes = set()
    curve = []
    for event_index, point in enumerate(events):
        for prior in events[:event_index]:
            dx = (point[0] - prior[0]) % P
            dy = (point[1] - prior[1]) % P
            if dx:
                slopes.add(direction((dx, dy)))
        if event_index >= 2:
            curve.append(len(slopes))
    # curve[t] is coverage after t additions; t=0 is the three public inputs.
    additions = len(events) - 3
    assert len(curve) == additions + 1
    discovered_weight = 0
    previous = 0
    for time, covered in enumerate(curve):
        discovered_weight += time * (covered - previous)
        previous = covered
    censored_numerator = discovered_weight + additions * (P - curve[-1])
    survival_numerator = sum(P - curve[time] for time in range(additions))
    assert survival_numerator == censored_numerator
    return {
        "coverage_curve": curve,
        "coverage_curve_sha256": sha256(
            json.dumps(curve, separators=(",", ":")).encode()
        ).hexdigest(),
        "covered_slopes": curve[-1],
        "covered_slope_discovery_weight": discovered_weight,
        "conditional_mean_additions_for_covered_slopes": str(
            Fraction(discovered_weight, curve[-1])
        ),
        "exact_censored_mean_min_T_budget": str(Fraction(censored_numerator, P)),
        "censored_mean_numerator": censored_numerator,
    }


def build_y_scalars(transcript, values):
    """Build requested multiples of B with one shared binary-power table."""
    values = sorted({value % P for value in values})
    maximum_bits = max(1, max(values).bit_length())
    powers = {0: B}
    for bit in range(1, maximum_bits):
        powers[bit] = transcript.charge_add(
            powers[bit - 1], powers[bit - 1], "y_power_setup", f"2^{bit}B"
        )
    scalars = {0: O, 1: B}
    for value in values:
        if value in scalars:
            continue
        bits = [bit for bit in range(maximum_bits) if (value >> bit) & 1]
        assert bits
        current = powers[bits[0]]
        for bit in bits[1:]:
            current = transcript.charge_add(
                current, powers[bit], "y_scalar_setup", f"{value}B:add_bit_{bit}"
            )
        scalars[value] = current
        assert current == (0, value)
    return scalars


def build_stripes(lines, points_per_line):
    assert points_per_line >= 1
    assert tuple(line[0] for line in lines) == tuple(range(len(lines)))
    transcript = Transcript()

    x_values = {0: O, 1: A}
    for x in range(2, len(lines)):
        x_values[x] = transcript.charge_add(
            x_values[x - 1], A, "x_setup", f"{x}A"
        )
        assert x_values[x] == (x, 0)

    requested_y = []
    for _, offset, stride in lines:
        requested_y.extend((offset, stride))
    y_values = build_y_scalars(transcript, requested_y)

    displayed = []
    currents = []
    strides = []
    line_lengths = []
    for x, offset, stride in lines:
        offset %= P
        stride %= P
        x_value = x_values[x]
        y_value = y_values[offset]
        if x_value == O:
            current = y_value
        elif y_value == O:
            current = x_value
        else:
            current = transcript.charge_add(
                x_value, y_value, "line_start", f"line_{x}_start"
            )
        assert current == (x, offset)
        displayed.append(current)
        currents.append(current)
        strides.append(y_values[stride])
        line_lengths.append(points_per_line)

    # Interleave the lines.  Each individual line is still sequential, while
    # cross-line slopes become available as early as the charged schedule allows.
    for step in range(1, points_per_line):
        for line_index, (x, offset, stride) in enumerate(lines):
            currents[line_index] = transcript.charge_add(
                currents[line_index], strides[line_index], "line_step",
                f"line_{x}_point_{step}",
            )
            assert currents[line_index] == (x, (offset + step * stride) % P)
            displayed.append(currents[line_index])

    for operation_index, operation in enumerate(transcript.operations):
        previous = transcript.events[: operation_index + 3]
        assert tuple(operation["left"]) in previous
        assert tuple(operation["right"]) in previous
        assert tuple(operation["result"]) == add(
            tuple(operation["left"]), tuple(operation["right"])
        )
    return transcript, displayed, currents, strides, line_lengths


def extend_round_robin(transcript, displayed, currents, strides, line_lengths, budget):
    line = 0
    while len(transcript.operations) < budget:
        currents[line] = transcript.charge_add(
            currents[line], strides[line], "budget_extension",
            f"line_{line}_extra_{line_lengths[line]}",
        )
        displayed.append(currents[line])
        line_lengths[line] += 1
        line = (line + 1) % len(currents)
    assert len(transcript.operations) == budget


def result_record(name, lines, points_per_line, target_budget=None):
    transcript, displayed, currents, strides, lengths = build_stripes(lines, points_per_line)
    base_additions = len(transcript.operations)
    if target_budget is not None:
        assert base_additions <= target_budget
        extend_round_robin(
            transcript, displayed, currents, strides, lengths, target_budget
        )
    all_profile = direction_profile(transcript.events)
    displayed_profile = direction_profile(displayed)
    additions = len(transcript.operations)
    assert all_profile["comparison_events"] == additions + 3
    assert not all_profile["unsigned_vertical_direction_present"] or (
        all_profile["unsigned_projective_directions"]
        == all_profile["unsigned_useful_finite_slopes"] + 1
    )
    return {
        "name": name,
        "lines": [list(line) for line in lines],
        "initial_points_per_line": points_per_line,
        "final_line_lengths": lengths,
        "base_additions_before_extension": base_additions,
        "charged_additions": additions,
        "phase_addition_counts": phase_counts(transcript),
        "displayed_profile": displayed_profile,
        "all_materialized_profile": all_profile,
        "prefix_metrics": prefix_metrics(transcript.events),
        "exact_generic_slope_success": str(
            Fraction(all_profile["unsigned_useful_finite_slopes"], P)
        ),
        "slopes_per_charged_addition": str(
            Fraction(all_profile["unsigned_useful_finite_slopes"], additions)
        ),
        "event_values_sha256": sha256(
            json.dumps(transcript.events, separators=(",", ":")).encode()
        ).hexdigest(),
        "operation_transcript_sha256": sha256(
            json.dumps(transcript.operations, separators=(",", ":"), sort_keys=True).encode()
        ).hexdigest(),
    }


def strongest_three_line_control(budget):
    best_final = None
    best_censored = None
    rows = []
    for n in range(2, budget + 1):
        lines = ((0, 0, 1), (1, n, n), (2, (-(n + 1)) % P, (-(n + 1)) % P))
        transcript, displayed, currents, strides, lengths = build_stripes(lines, n)
        if len(transcript.operations) > budget:
            continue
        extend_round_robin(transcript, displayed, currents, strides, lengths, budget)
        profile = direction_profile(transcript.events)
        prefix = prefix_metrics(transcript.events)
        row = {
            "n": n,
            "base_additions": sum(phase_counts(transcript).get(key, 0) for key in (
                "x_setup", "y_power_setup", "y_scalar_setup", "line_start", "line_step"
            )),
            "final_line_lengths": lengths,
            "useful_finite_slopes": profile["unsigned_useful_finite_slopes"],
            "censored_mean_numerator": prefix["censored_mean_numerator"],
        }
        rows.append(row)
        if best_final is None or row["useful_finite_slopes"] > best_final["useful_finite_slopes"]:
            best_final = row
        if best_censored is None or row["censored_mean_numerator"] < best_censored["censored_mean_numerator"]:
            best_censored = row
    assert best_final is not None and best_censored is not None

    def rebuild(best, name):
        n = best["n"]
        lines = ((0, 0, 1), (1, n, n), (2, (-(n + 1)) % P, (-(n + 1)) % P))
        record = result_record(name, lines, n, budget)
        assert record["all_materialized_profile"]["unsigned_useful_finite_slopes"] == best["useful_finite_slopes"]
        assert record["prefix_metrics"]["censored_mean_numerator"] == best["censored_mean_numerator"]
        record["searched_n_values"] = len(rows)
        return record

    return (
        rebuild(best_final, "charged_three_line_final_coverage_control"),
        rebuild(best_censored, "charged_three_line_censored_mean_control"),
    )


def summary(values, candidate):
    mean = Fraction(sum(values), len(values))
    return {
        "mean_exact": str(mean),
        "mean": float(mean),
        "sample_stdev": stdev(values),
        "minimum": min(values),
        "maximum": max(values),
        "controls_greater_than_candidate": sum(value > candidate for value in values),
        "controls_equal_to_candidate": sum(value == candidate for value in values),
        "controls_less_than_candidate": sum(value < candidate for value in values),
    }


def compare_curves(left, right):
    assert len(left) == len(right)
    return {
        "left_ahead": sum(a > b for a, b in zip(left, right)),
        "ties": sum(a == b for a, b in zip(left, right)),
        "left_behind": sum(a < b for a, b in zip(left, right)),
        "maximum_left_advantage": max(a - b for a, b in zip(left, right)),
        "maximum_left_deficit": max(b - a for a, b in zip(left, right)),
    }


def ideal_random_references(additions):
    getcontext().prec = 60
    q = additions + 3
    survival = Decimal(1)
    decimal_p = Decimal(P)
    for index in range(q):
        survival *= Decimal(P - index) / decimal_p
    birthday_success = Decimal(1) - survival
    # Bernstein--Lange's usual favorable Poisson-scale random-walk reference.
    poisson_success = Decimal(1) - (-(Decimal(additions) ** 2) / (2 * decimal_p)).exp()
    return {
        "independent_uniform_group_samples": q,
        "exact_birthday_success_decimal": str(birthday_success),
        "poisson_random_walk_reference_decimal": str(poisson_success),
        "references_are_not_fixed_r_rho_measurements": True,
    }


def main():
    candidate = result_record(
        "explicit_4_stripe_candidate", CANDIDATE_LINES, CANDIDATE_POINTS_PER_LINE
    )
    budget = candidate["charged_additions"]
    assert budget == 427
    assert candidate["displayed_profile"]["unsigned_useful_finite_slopes"] == 48629
    assert candidate["all_materialized_profile"]["unsigned_useful_finite_slopes"] == 50147

    three_line, three_line_censored = strongest_three_line_control(budget)
    candidate_slopes = candidate["all_materialized_profile"]["unsigned_useful_finite_slopes"]
    three_line_slopes = three_line["all_materialized_profile"]["unsigned_useful_finite_slopes"]

    random_rows = []
    random_counts = []
    for trial in range(RANDOM_TRIALS):
        seed = 20260921 + 10_000 * budget + trial
        events = random_control("random_accumulator", budget, seed)
        profile = direction_profile(events)
        prefix = prefix_metrics(events)
        count = profile["unsigned_useful_finite_slopes"]
        random_counts.append(count)
        random_rows.append({
            "trial": trial,
            "seed": seed,
            "useful_finite_slopes": count,
            "exact_generic_slope_success": str(Fraction(count, P)),
            "prefix_metrics": {
                "covered_slope_discovery_weight": prefix["covered_slope_discovery_weight"],
                "conditional_mean_additions_for_covered_slopes": prefix["conditional_mean_additions_for_covered_slopes"],
                "exact_censored_mean_min_T_budget": prefix["exact_censored_mean_min_T_budget"],
                "censored_mean_numerator": prefix["censored_mean_numerator"],
                "coverage_curve_sha256": prefix["coverage_curve_sha256"],
            },
            "event_values_sha256": sha256(
                json.dumps(events, separators=(",", ":")).encode()
            ).hexdigest(),
        })

    report = {
        "evidence_type": "finite_charged_multiline_stripe_screen",
        "field_prime": P,
        "python_version": platform.python_version(),
        "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "metric": {
            "useful": "Distinct finite slopes dy/dx with dx nonzero among every materialized coefficient state.",
            "generic_success": "For a target-independent generic addition chain, d useful finite slopes give exact collision success d/p before any final guess.",
            "vertical": "The vertical projective direction is excluded.",
        },
        "candidate": candidate,
        "three_line_control": three_line,
        "three_line_censored_mean_control": three_line_censored,
        "random_accumulator_summary": summary(random_counts, candidate_slopes),
        "random_accumulator_censored_mean_numerator_summary": summary(
            [row["prefix_metrics"]["censored_mean_numerator"] for row in random_rows],
            candidate["prefix_metrics"]["censored_mean_numerator"],
        ),
        "random_accumulator_trials": random_rows,
        "ideal_random_references": ideal_random_references(budget),
        "comparisons": {
            "candidate_minus_three_line_slopes": candidate_slopes - three_line_slopes,
            "candidate_beats_three_line_at_equal_addition_budget": candidate_slopes > three_line_slopes,
            "candidate_beats_every_random_accumulator_control": candidate_slopes > max(random_counts),
            "candidate_beats_any_random_accumulator_control": candidate_slopes > min(random_counts),
            "candidate_censored_mean_minus_best_three_line_numerator": (
                candidate["prefix_metrics"]["censored_mean_numerator"]
                - three_line_censored["prefix_metrics"]["censored_mean_numerator"]
            ),
            "candidate_has_lower_censored_mean_than_best_three_line": (
                candidate["prefix_metrics"]["censored_mean_numerator"]
                < three_line_censored["prefix_metrics"]["censored_mean_numerator"]
            ),
            "candidate_has_lower_censored_mean_than_any_random_accumulator": (
                candidate["prefix_metrics"]["censored_mean_numerator"]
                < min(row["prefix_metrics"]["censored_mean_numerator"] for row in random_rows)
            ),
            "candidate_curve_vs_best_three_line_censored_mean": compare_curves(
                candidate["prefix_metrics"]["coverage_curve"],
                three_line_censored["prefix_metrics"]["coverage_curve"],
            ),
            "candidate_curve_vs_best_three_line_final_coverage": compare_curves(
                candidate["prefix_metrics"]["coverage_curve"],
                three_line["prefix_metrics"]["coverage_curve"],
            ),
        },
        "storage_mapping": {
            "exact_records": budget + 3,
            "hypothetical_record_slot_bytes_from_capacity_envelope": 128,
            "final_record_payload_bytes_excluding_directory_and_filter": (budget + 3) * 128,
            "external_sort_role": "Sort canonical group encodings to expose exact equalities; coefficient records supply certificates.",
            "external_sort_io_measured": False,
            "large_capacity_transfer_proved": False,
        },
        "verification": {
            "candidate_addition_equations_checked": budget,
            "candidate_expected_cost_and_slope_count_asserted": True,
            "three_line_n_values_exhaustively_screened_under_budget": three_line["searched_n_values"],
            "matched_addition_budget": candidate["charged_additions"] == three_line["charged_additions"] == budget,
            "all_random_controls_match_addition_and_event_budget": True,
        },
        "limitations": [
            "The explicit offsets and strides were selected by a finite seeded screen; no optimality or asymptotic construction is proved.",
            "The random-accumulator ranges are sampled controls, not a theorem about random chains.",
            "The ideal-random references are favorable random-mapping comparators, not measured fixed-r rho runs.",
            "The primary source already reports a four-giant variant and speculates about increasing giants, so the multi-line concept is not a novelty claim.",
            "No external-sort runtime, I/O volume, hardware feasibility, ECDLP speedup, or 1 TB / 100 TB throughput result is established.",
        ],
    }
    output = Path(__file__).with_name("multiline_stripe_checks.json")
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({
        "budget": budget,
        "candidate_slopes": candidate_slopes,
        "three_line_slopes": three_line_slopes,
        "random_accumulator_min": min(random_counts),
        "random_accumulator_mean": float(Fraction(sum(random_counts), len(random_counts))),
        "random_accumulator_max": max(random_counts),
        "output": str(output),
    }, indent=2))


if __name__ == "__main__":
    main()
