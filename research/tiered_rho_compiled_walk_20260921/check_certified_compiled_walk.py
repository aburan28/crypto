"""Exact all-target certificate for one compiled reservoir-guided archive walk."""

from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import json

from check_charged_polynomial_construction import P, add
from check_walk_reservoir_guided import run
from check_walk_reservoir_bloom import run_bloom


CHOICES = 4
RESERVOIR_SIZE = 32
SEED = 20261305
BUDGET = 427


def offsets(selected_addends):
    result = []
    for step, index in enumerate(selected_addends):
        distance = 3 + step - index
        assert 1 <= distance <= 64
        result.append(distance)
    return result


def pack_six_bit(values):
    accumulator = 0
    bits = 0
    output = bytearray()
    for value in values:
        accumulator |= (value - 1) << bits
        bits += 6
        while bits >= 8:
            output.append(accumulator & 0xFF)
            accumulator >>= 8
            bits -= 8
    if bits:
        output.append(accumulator & 0xFF)
    return bytes(output)


def unpack_six_bit(data, count):
    accumulator = 0
    bits = 0
    cursor = 0
    values = []
    while len(values) < count:
        while bits < 6:
            accumulator |= data[cursor] << bits
            cursor += 1
            bits += 8
        values.append((accumulator & 0x3F) + 1)
        accumulator >>= 6
        bits -= 6
    assert cursor == len(data)
    return values


def replay_coefficients(schedule):
    events = [(0, 0), (1, 0), (0, 1)]
    for distance in schedule:
        new_index = len(events)
        addend_index = new_index - distance
        assert 0 <= addend_index < new_index
        events.append(add(events[-1], events[addend_index]))
    assert len(events) == BUDGET + 3
    return events


def first_useful_collision(events, target):
    first_by_scalar = {}
    for index, coefficients in enumerate(events):
        value = (coefficients[0] * target + coefficients[1]) % P
        prior_index = first_by_scalar.get(value)
        if prior_index is None:
            first_by_scalar[value] = index
            continue
        prior = events[prior_index]
        delta_a = (coefficients[0] - prior[0]) % P
        delta_b = (coefficients[1] - prior[1]) % P
        if not delta_a:
            assert delta_b == 0
            continue
        recovered = -delta_b * pow(delta_a, -1, P) % P
        assert recovered == target
        return max(0, index - 2)
    return BUDGET


def all_target_audit(events):
    histogram = Counter()
    for target in range(P):
        histogram[first_useful_collision(events, target)] += 1
    assert sum(histogram.values()) == P
    numerator = sum(time * count for time, count in histogram.items())
    return {
        "censored_mean_numerator": numerator,
        "censored_mean": str(Fraction(numerator, P)),
        "censored_mean_float": numerator / P,
        "targets_solved_before_budget": P - histogram[BUDGET],
        "targets_censored_at_budget": histogram[BUDGET],
        "stopping_histogram": dict(sorted(histogram.items())),
        "histogram_sha256": sha256(
            json.dumps(dict(sorted(histogram.items())), separators=(",", ":")).encode()
        ).hexdigest(),
    }


def main():
    guided = run(True, CHOICES, RESERVOIR_SIZE, SEED)
    control = run(False, CHOICES, RESERVOIR_SIZE, SEED)
    bloom = run_bloom(CHOICES, RESERVOIR_SIZE, SEED)
    assert bloom["events_sha256"] == guided["events_sha256"]
    assert bloom["bloom"]["false_negatives"] == 0
    assert bloom["bloom"]["changed_choices_vs_exact_online_score"] == 0

    guided_offsets = offsets(guided["selected_addend_indices"])
    control_offsets = offsets(control["selected_addend_indices"])
    packed = pack_six_bit(guided_offsets)
    assert len(packed) == 321
    assert unpack_six_bit(packed, BUDGET) == guided_offsets
    guided_events = replay_coefficients(guided_offsets)
    control_events = replay_coefficients(control_offsets)
    assert sha256(json.dumps(guided_events, separators=(",", ":")).encode()).hexdigest() == guided["events_sha256"]
    assert sha256(json.dumps(control_events, separators=(",", ":")).encode()).hexdigest() == control["events_sha256"]

    guided_audit = all_target_audit(guided_events)
    control_audit = all_target_audit(control_events)
    assert guided_audit["censored_mean_numerator"] == guided["prefix_metrics"]["censored_mean_numerator"]
    assert control_audit["censored_mean_numerator"] == control["prefix_metrics"]["censored_mean_numerator"]
    assert guided_audit["censored_mean_numerator"] < control_audit["censored_mean_numerator"]

    schedule_path = Path(__file__).with_name("certified_compiled_walk_schedule.bin")
    schedule_path.write_bytes(packed)
    report = {
        "evidence_type": "exact_all_target_compiled_archive_walk_certificate",
        "field_prime": P,
        "budget": BUDGET,
        "compiler_parameters": {
            "choices": CHOICES,
            "reservoir_size": RESERVOIR_SIZE,
            "public_seed": SEED,
        },
        "schedule": {
            "packed_bytes": len(packed),
            "bits_per_offset": 6,
            "maximum_backward_offset": max(guided_offsets),
            "sha256": sha256(packed).hexdigest(),
            "path": schedule_path.name,
        },
        "guided": guided_audit,
        "control": control_audit,
        "exact_improvement": {
            "censored_mean_numerator_control_minus_guided": (
                control_audit["censored_mean_numerator"]
                - guided_audit["censored_mean_numerator"]
            ),
            "censored_mean_additions_control_minus_guided": str(
                Fraction(
                    control_audit["censored_mean_numerator"]
                    - guided_audit["censored_mean_numerator"],
                    P,
                )
            ),
            "relative_censored_mean_reduction": (
                1
                - guided_audit["censored_mean_numerator"]
                / control_audit["censored_mean_numerator"]
            ),
        },
        "bloom_equivalence": {
            "events_match_exact": True,
            "false_positives": bloom["bloom"]["false_positives"],
            "false_negatives": 0,
            "changed_choices": 0,
        },
        "limitations": [
            "The theorem is finite for one public order, budget, seed, and matched control.",
            "It proves online group-addition stopping improvement after compilation, not amortized wall time.",
            "The schedule is time/history dependent and is not a fixed coalescing rho function.",
            "No asymptotic, novelty, private-target, or cryptographic-scale claim follows.",
        ],
        "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    output = Path(__file__).with_name("certified_compiled_walk.json")
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({
        "schedule_sha256": report["schedule"]["sha256"],
        "guided_mean": guided_audit["censored_mean"],
        "control_mean": control_audit["censored_mean"],
        "improvement": report["exact_improvement"],
        "output": str(output),
    }, indent=2))


if __name__ == "__main__":
    main()
