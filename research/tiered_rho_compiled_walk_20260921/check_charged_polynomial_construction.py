"""Exact finite charged-addition screen for polynomial coefficient families.

This is a coefficient-plane experiment, not an elliptic-curve implementation
or a discrete-log solver.  The public inputs O=(0,0), A=(1,0), and B=(0,1)
represent generic-group coefficient pairs.  Every vector addition is charged
as one group addition and every addition output is an eligible comparison
event.  The signed view applies negation once to every comparison event for
all constructions, including repeated values, so the operation counts match.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from pathlib import Path
from random import Random
from statistics import stdev
import json
import platform


P = 65537
DEGREES = (2, 3, 4)
SIZES = (16, 64, 256)
EQUAL_ADDITION_BUDGETS = (64, 256, 1024)
CONTROL_FAMILIES = ("random_pair_dag", "random_accumulator")
CONTROL_TRIALS = 16
O = (0, 0)
A = (1, 0)
B = (0, 1)
INVERSES = [0] + [pow(value, -1, P) for value in range(1, P)]


def add(left, right):
    return ((left[0] + right[0]) % P, (left[1] + right[1]) % P)


def negate(value):
    return ((-value[0]) % P, (-value[1]) % P)


def direction(delta):
    x, y = delta[0] % P, delta[1] % P
    assert x or y
    return -1 if x == 0 else y * INVERSES[x] % P


def direction_profile(events):
    """Count directions of the event values and their explicit signed view.

    Events may repeat.  Direction sets operate on distinct group values, while
    event_count and sign_applications preserve all paid/materialized events.
    O is always present, so signed directions are unsigned difference
    directions plus pair-sum directions.
    """
    points = sorted(set(events))
    assert O in points
    unsigned = set()
    for left, right in combinations(points, 2):
        unsigned.add(direction((left[0] - right[0], left[1] - right[1])))
    nonzero = [value for value in points if value != O]
    signed = set(unsigned)
    for index, left in enumerate(nonzero):
        for right in nonzero[index:]:
            total = add(left, right)
            if total != O:
                signed.add(direction(total))
    closed = set(points)
    closed.update(negate(value) for value in points)
    return {
        "comparison_events": len(events),
        "distinct_comparison_states": len(points),
        "unsigned_projective_directions": len(unsigned),
        "unsigned_useful_finite_slopes": len(unsigned - {-1}),
        "unsigned_vertical_direction_present": -1 in unsigned,
        "sign_applications": len(events),
        "distinct_signed_states": len(closed),
        "signed_projective_directions": len(signed),
        "signed_useful_finite_slopes": len(signed - {-1}),
        "signed_vertical_direction_present": -1 in signed,
    }


def crosscheck_signed_closure(events, profile):
    """Independently enumerate every pair after explicit negation closure."""
    points = set(events)
    points.update(negate(value) for value in events)
    values = set()
    for left, right in combinations(sorted(points), 2):
        values.add(direction((left[0] - right[0], left[1] - right[1])))
    assert len(values) == profile["signed_projective_directions"]
    assert len(values - {-1}) == profile["signed_useful_finite_slopes"]
    assert (-1 in values) == profile["signed_vertical_direction_present"]


class Transcript:
    def __init__(self):
        self.events = [O, A, B]
        self.operations = []

    def charge_add(self, left, right, phase, label):
        result = add(left, right)
        self.events.append(result)
        self.operations.append({
            "phase": phase,
            "label": label,
            "left": list(left),
            "right": list(right),
            "result": list(result),
        })
        return result


def initialize_polynomial(transcript, degree):
    """Build forward-difference registers from A and B with charged chains."""
    x1 = transcript.charge_add(A, B, "initialization", "x1=A+B")
    b2 = transcript.charge_add(B, B, "initialization", "b2=B+B")
    if degree == 2:
        return [O, x1, b2]

    b3 = transcript.charge_add(b2, B, "initialization", "b3=b2+B")
    b6 = transcript.charge_add(b3, b3, "initialization", "b6=b3+b3")
    if degree == 3:
        # Delta^2 t^3 at zero and Delta^3 t^3 are both 6.
        return [O, x1, b6, b6]

    assert degree == 4
    b12 = transcript.charge_add(b6, b6, "initialization", "b12=b6+b6")
    b24 = transcript.charge_add(b12, b12, "initialization", "b24=b12+b12")
    b36 = transcript.charge_add(b24, b12, "initialization", "b36=b24+b12")
    b7 = transcript.charge_add(b6, B, "initialization", "b7=b6+B")
    b14 = transcript.charge_add(b7, b7, "initialization", "b14=b7+b7")
    # Delta^k t^4 at zero for k=1..4 is 1,14,36,24.
    return [O, x1, b14, b36, b24]


def polynomial_transcript(degree, size):
    assert degree in DEGREES and size in SIZES
    transcript = Transcript()
    registers = initialize_polynomial(transcript, degree)
    initialization_additions = len(transcript.operations)
    vertices = [O, registers[1]]
    registers[0] = registers[1]

    # X_1 is already available as Delta X_0.  Advance the other forward
    # differences from t=0 to t=1 without paying the redundant O+X_1 add.
    for level in range(1, degree):
        registers[level] = transcript.charge_add(
            registers[level], registers[level + 1],
            "difference_advance", f"t=1,level={level}",
        )

    for t in range(2, size):
        registers[0] = transcript.charge_add(
            registers[0], registers[1], "vertex", f"x{t}",
        )
        vertices.append(registers[0])
        if t < size - 1:
            for level in range(1, degree):
                registers[level] = transcript.charge_add(
                    registers[level], registers[level + 1],
                    "difference_advance", f"t={t},level={level}",
                )

    expected_vertices = [(t % P, pow(t, degree, P)) for t in range(size)]
    assert vertices == expected_vertices
    expected_initialization = {2: 2, 3: 4, 4: 9}[degree]
    expected_additions = expected_initialization + degree * (size - 2)
    assert initialization_additions == expected_initialization
    assert len(transcript.operations) == expected_additions
    assert len(transcript.events) == expected_additions + 3
    for operation in transcript.operations:
        assert tuple(operation["result"]) == add(
            tuple(operation["left"]), tuple(operation["right"])
        )
    return transcript, vertices, initialization_additions


def polynomial_budget_transcript(degree, budget):
    """Run the conventional update schedule to exactly ``budget`` additions.

    A budget may stop during a difference-register advance.  Such paid partial
    work remains in the comparison event stream, but it is not mislabeled as a
    completed graph vertex.
    """
    assert degree in DEGREES
    transcript = Transcript()
    registers = initialize_polynomial(transcript, degree)
    initialization_additions = len(transcript.operations)
    assert budget >= initialization_additions
    registers[0] = registers[1]
    vertices = [O, registers[0]]
    remaining = budget - initialization_additions
    next_t = 2
    partial_levels = 0
    while remaining:
        partial_levels = 0
        for level in range(1, degree):
            if not remaining:
                break
            registers[level] = transcript.charge_add(
                registers[level], registers[level + 1],
                "difference_advance", f"before_x{next_t},level={level}",
            )
            remaining -= 1
            partial_levels += 1
        if not remaining:
            break
        registers[0] = transcript.charge_add(
            registers[0], registers[1], "vertex", f"x{next_t}",
        )
        vertices.append(registers[0])
        next_t += 1
        remaining -= 1
        partial_levels = 0
    expected_vertices = [(t % P, pow(t, degree, P)) for t in range(len(vertices))]
    assert vertices == expected_vertices
    assert len(transcript.operations) == budget
    assert len(transcript.events) == budget + 3
    completed_cycles, remainder = divmod(budget - initialization_additions, degree)
    assert len(vertices) == 2 + completed_cycles
    assert partial_levels == remainder
    for operation in transcript.operations:
        assert tuple(operation["result"]) == add(
            tuple(operation["left"]), tuple(operation["right"])
        )
    return transcript, vertices, initialization_additions, remainder


def random_control(family, additions, seed):
    assert family in CONTROL_FAMILIES
    rng = Random(seed)
    events = [O, A, B]
    selectable = [1, 2]
    current = None
    for operation_index in range(additions):
        if family == "random_pair_dag":
            left = events[selectable[rng.randrange(len(selectable))]]
            right = events[selectable[rng.randrange(len(selectable))]]
        else:
            if operation_index == 0:
                left, right = A, B
            else:
                left = current
                right = events[selectable[rng.randrange(len(selectable))]]
        current = add(left, right)
        events.append(current)
        if current != O:
            selectable.append(len(events) - 1)
    assert len(events) == additions + 3
    return events


def fibonacci_control(additions):
    events = [O, A, B]
    left, right = A, B
    for _ in range(additions):
        result = add(left, right)
        events.append(result)
        left, right = right, result
    assert len(events) == additions + 3
    return events


def seed_for(family, additions, trial):
    offset = {"random_pair_dag": 0, "random_accumulator": 10_000_000}[family]
    return 20260921 + 10_000 * additions + offset + trial


def compact_hash(value):
    encoded = json.dumps(value, separators=(",", ":"), sort_keys=True).encode()
    return sha256(encoded).hexdigest()


def sample_summary(values):
    mean = Fraction(sum(values), len(values))
    return {
        "mean_exact": str(mean),
        "mean": float(mean),
        "sample_stdev": stdev(values),
        "minimum": min(values),
        "maximum": max(values),
    }


def metric_summary(values, candidate):
    result = sample_summary(values)
    result.update({
        "controls_greater_than_candidate": sum(value > candidate for value in values),
        "controls_equal_to_candidate": sum(value == candidate for value in values),
        "controls_less_than_candidate": sum(value < candidate for value in values),
    })
    return result


def main():
    assert all(value * INVERSES[value] % P == 1 for value in range(1, P))
    candidate_rows = []
    control_rows = []
    control_summaries = []
    fibonacci_rows = []
    equal_budget_candidate_rows = []
    equal_budget_control_rows = []
    equal_budget_control_summaries = []
    equal_budget_fibonacci_controls = []
    operation_checks = 0
    explicit_signed_closure_crosschecks = 0

    for degree in DEGREES:
        for size in SIZES:
            transcript, vertices, init_additions = polynomial_transcript(degree, size)
            operation_checks += len(transcript.operations)
            all_profile = direction_profile(transcript.events)
            vertex_profile = direction_profile(vertices)
            crosscheck_signed_closure(transcript.events, all_profile)
            crosscheck_signed_closure(vertices, vertex_profile)
            explicit_signed_closure_crosschecks += 2
            additions = len(transcript.operations)
            phase_counts = {}
            for operation in transcript.operations:
                phase = operation["phase"]
                phase_counts[phase] = phase_counts.get(phase, 0) + 1
            candidate = {
                "degree": degree,
                "displayed_vertices": size,
                "charged_additions": additions,
                "initialization_additions": init_additions,
                "phase_addition_counts": phase_counts,
                "comparison_events_expected": additions + 3,
                "sign_applications_expected": additions + 3,
                "displayed_vertex_profile": vertex_profile,
                "all_materialized_state_profile": all_profile,
                "transcript_sha256": compact_hash(transcript.operations),
                "event_values_sha256": compact_hash(transcript.events),
            }
            assert all_profile["comparison_events"] == additions + 3
            assert all_profile["sign_applications"] == additions + 3
            candidate_rows.append(candidate)

            fib_events = fibonacci_control(additions)
            fib_profile = direction_profile(fib_events)
            crosscheck_signed_closure(fib_events, fib_profile)
            explicit_signed_closure_crosschecks += 1
            fib_row = {
                "degree": degree,
                "displayed_vertices_for_candidate": size,
                "charged_additions": additions,
                "profile": fib_profile,
                "event_values_sha256": compact_hash(fib_events),
                "weakly_dominates_candidate_both_direction_metrics": (
                    fib_profile["unsigned_useful_finite_slopes"] >= all_profile["unsigned_useful_finite_slopes"]
                    and fib_profile["signed_useful_finite_slopes"] >= all_profile["signed_useful_finite_slopes"]
                ),
                "strictly_dominates_candidate_both_direction_metrics": (
                    fib_profile["unsigned_useful_finite_slopes"] > all_profile["unsigned_useful_finite_slopes"]
                    and fib_profile["signed_useful_finite_slopes"] > all_profile["signed_useful_finite_slopes"]
                ),
            }
            fibonacci_rows.append(fib_row)

            for family in CONTROL_FAMILIES:
                unsigned_counts = []
                signed_counts = []
                weak_dominance = 0
                strict_dominance = 0
                for trial in range(CONTROL_TRIALS):
                    seed = seed_for(family, additions, trial)
                    events = random_control(family, additions, seed)
                    profile = direction_profile(events)
                    if trial == 0:
                        crosscheck_signed_closure(events, profile)
                        explicit_signed_closure_crosschecks += 1
                    assert profile["comparison_events"] == additions + 3
                    assert profile["sign_applications"] == additions + 3
                    unsigned_counts.append(profile["unsigned_useful_finite_slopes"])
                    signed_counts.append(profile["signed_useful_finite_slopes"])
                    weak = (
                        profile["unsigned_useful_finite_slopes"] >= all_profile["unsigned_useful_finite_slopes"]
                        and profile["signed_useful_finite_slopes"] >= all_profile["signed_useful_finite_slopes"]
                    )
                    strict = (
                        profile["unsigned_useful_finite_slopes"] > all_profile["unsigned_useful_finite_slopes"]
                        and profile["signed_useful_finite_slopes"] > all_profile["signed_useful_finite_slopes"]
                    )
                    weak_dominance += int(weak)
                    strict_dominance += int(strict)
                    control_rows.append({
                        "degree": degree,
                        "displayed_vertices_for_candidate": size,
                        "charged_additions": additions,
                        "family": family,
                        "trial": trial,
                        "seed": seed,
                        "profile": profile,
                        "event_values_sha256": compact_hash(events),
                        "weakly_dominates_candidate_both_direction_metrics": weak,
                        "strictly_dominates_candidate_both_direction_metrics": strict,
                    })
                control_summaries.append({
                    "degree": degree,
                    "displayed_vertices_for_candidate": size,
                    "charged_additions": additions,
                    "family": family,
                    "trials": CONTROL_TRIALS,
                    "candidate_all_state_unsigned_useful_finite_slopes": all_profile["unsigned_useful_finite_slopes"],
                    "candidate_all_state_signed_useful_finite_slopes": all_profile["signed_useful_finite_slopes"],
                    "unsigned_useful_finite_slope_summary": metric_summary(unsigned_counts, all_profile["unsigned_useful_finite_slopes"]),
                    "signed_useful_finite_slope_summary": metric_summary(signed_counts, all_profile["signed_useful_finite_slopes"]),
                    "controls_weakly_dominating_both_metrics": weak_dominance,
                    "controls_strictly_dominating_both_metrics": strict_dominance,
                })

    for additions in EQUAL_ADDITION_BUDGETS:
        candidates_at_budget = {}
        for degree in DEGREES:
            transcript, vertices, init_additions, partial_updates = polynomial_budget_transcript(
                degree, additions
            )
            operation_checks += len(transcript.operations)
            profile = direction_profile(transcript.events)
            vertex_profile = direction_profile(vertices)
            crosscheck_signed_closure(transcript.events, profile)
            explicit_signed_closure_crosschecks += 1
            candidates_at_budget[degree] = profile
            equal_budget_candidate_rows.append({
                "degree": degree,
                "charged_additions": additions,
                "initialization_additions": init_additions,
                "completed_displayed_vertices": len(vertices),
                "paid_partial_difference_updates_after_last_vertex": partial_updates,
                "displayed_vertex_profile": vertex_profile,
                "all_materialized_state_profile": profile,
                "transcript_sha256": compact_hash(transcript.operations),
                "event_values_sha256": compact_hash(transcript.events),
            })

        fib_events = fibonacci_control(additions)
        fib_profile = direction_profile(fib_events)
        equal_budget_fibonacci_controls.append({
            "charged_additions": additions,
            "profile": fib_profile,
            "event_values_sha256": compact_hash(fib_events),
        })

        for family in CONTROL_FAMILIES:
            unsigned_counts = []
            signed_counts = []
            strict_by_degree = {str(degree): 0 for degree in DEGREES}
            for trial in range(CONTROL_TRIALS):
                seed = seed_for(family, additions, trial)
                events = random_control(family, additions, seed)
                profile = direction_profile(events)
                if trial == 0:
                    crosscheck_signed_closure(events, profile)
                    explicit_signed_closure_crosschecks += 1
                unsigned_counts.append(profile["unsigned_useful_finite_slopes"])
                signed_counts.append(profile["signed_useful_finite_slopes"])
                strict_flags = {}
                for degree, candidate_profile in candidates_at_budget.items():
                    strict = (
                        profile["unsigned_useful_finite_slopes"]
                        > candidate_profile["unsigned_useful_finite_slopes"]
                        and profile["signed_useful_finite_slopes"]
                        > candidate_profile["signed_useful_finite_slopes"]
                    )
                    strict_flags[str(degree)] = strict
                    strict_by_degree[str(degree)] += int(strict)
                equal_budget_control_rows.append({
                    "charged_additions": additions,
                    "family": family,
                    "trial": trial,
                    "seed": seed,
                    "profile": profile,
                    "strictly_dominates_candidate_both_finite_slope_metrics_by_degree": strict_flags,
                    "event_values_sha256": compact_hash(events),
                })
            equal_budget_control_summaries.append({
                "charged_additions": additions,
                "family": family,
                "trials": CONTROL_TRIALS,
                "unsigned_useful_finite_slope_summary": sample_summary(unsigned_counts),
                "signed_useful_finite_slope_summary": sample_summary(signed_counts),
                "controls_strictly_dominating_both_metrics_by_candidate_degree": strict_by_degree,
            })

    report = {
        "evidence_type": "exact_finite_charged_addition_screen_with_seeded_controls",
        "field_prime": P,
        "python_version": platform.python_version(),
        "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "degrees": DEGREES,
        "displayed_vertex_sizes": SIZES,
        "equal_charged_addition_budgets": EQUAL_ADDITION_BUDGETS,
        "random_trials_per_family_and_row": CONTROL_TRIALS,
        "charged_model": {
            "public_inputs": {"O": list(O), "A": list(A), "B": list(B)},
            "addition": "One coefficient-vector addition modulo p models one charged generic-group addition.",
            "initialization": "All non-input finite-difference constants are built by the transcript's explicit addition chain.",
            "comparison": "O, A, B, and every addition output are comparison events; repeated values remain events but direction sets deduplicate values.",
            "negation": "Negation/inversion is not free: the signed view applies one separately charged abstract sign map to every comparison event, including repeated values and O; thus sign applications equal additions+3 for every matched construction. A sign map is not equated to a group addition.",
            "useful_slope": "A useful finite slope has nonzero A-coordinate difference. The vertical projective direction (encoded internally as -1) is reported separately and excluded from all candidate/control rankings.",
            "unpriced": "Coefficient/control bookkeeping, comparison lookup, memory traffic, and field-level sign-map cost are not converted into group additions.",
        },
        "candidate_rows": candidate_rows,
        "fibonacci_controls": fibonacci_rows,
        "random_control_summaries": control_summaries,
        "random_control_trials": control_rows,
        "equal_budget_candidate_rows": equal_budget_candidate_rows,
        "equal_budget_fibonacci_controls": equal_budget_fibonacci_controls,
        "equal_budget_random_control_summaries": equal_budget_control_summaries,
        "equal_budget_random_control_trials": equal_budget_control_rows,
        "verification": {
            "candidate_operation_equations_checked": operation_checks,
            "candidate_rows_matching_exact_vertices_t_td": len(candidate_rows),
            "explicit_signed_closure_crosschecks": explicit_signed_closure_crosschecks,
            "all_rows_have_equal_addition_comparison_and_sign_budgets_against_controls": True,
            "equal_budget_rows_use_all_budget_and_retain_paid_partial_updates": True,
            "inverse_table_checked": True,
        },
        "limitations": [
            "Coefficient pairs are a finite generic-group model, not elliptic-curve points or a DLP solver.",
            "Direction count is a geometric diagnostic, not a success probability or end-to-end search cost.",
            "Random control ranges are seeded finite samples, not population or asymptotic theorems.",
            "The signed accounting charges one abstract sign application per event; it does not assert that negation equals one group addition on any curve.",
            "Lookup, storage, comparison, control-flow, and I/O costs are not priced.",
            "No ECDLP speedup, novelty, asymptotic, or 1 TB RAM / 100 TB disk feasibility claim.",
        ],
    }
    output = Path(__file__).with_name("charged_polynomial_construction_checks.json")
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({
        "candidate_rows": len(candidate_rows),
        "candidate_operation_equations_checked": operation_checks,
        "random_control_trials": len(control_rows),
        "fibonacci_controls": len(fibonacci_rows),
        "output": str(output),
    }, indent=2))


if __name__ == "__main__":
    main()
