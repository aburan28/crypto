"""Target-independent coverage-guided legal addition-chain experiment.

The planner scores coefficient sums of archived parents with a Bloom sketch of
all already covered finite slopes.  Only the selected sum consumes one charged
group addition.  All intermediate records are retained.  This is a public
finite coefficient-plane calculation, not an elliptic-curve solver.
"""

from collections import OrderedDict
from fractions import Fraction
from hashlib import sha256
from math import exp, log, sqrt
from pathlib import Path
from statistics import stdev
import json
import platform

from check_charged_polynomial_construction import A, B, O, P, add, direction, direction_profile, random_control
from check_multiline_stripes import prefix_metrics, strongest_three_line_control


BUDGET = 427
CANDIDATES_PER_STEP = 32
HOT_WINDOW = 64
TRIALS = 8
RECORD_BYTES = 128
PAGE_RECORDS = 64
PARENT_CACHE_PAGES = 8
BLOOM_BITS = 1 << 20
BLOOM_HASHES = 7
LARGE_SKETCH_BYTES = 128_000_000_000
RAM_LIMIT = 1_000_000_000_000
DISK_LIMIT = 100_000_000_000_000
MODE_B_FIXED_RAM_WITHOUT_SLOPE_SKETCH = 369_310_378_496
MODE_B_RECORD_RAM_BYTES = 2


def mix64(value):
    value &= (1 << 64) - 1
    value ^= value >> 30
    value = (value * 0xBF58476D1CE4E5B9) & ((1 << 64) - 1)
    value ^= value >> 27
    value = (value * 0x94D049BB133111EB) & ((1 << 64) - 1)
    return value ^ (value >> 31)


class SlopeBloom:
    def __init__(self):
        assert BLOOM_BITS & (BLOOM_BITS - 1) == 0
        self.data = bytearray(BLOOM_BITS // 8)
        self.queries = 0
        self.hash_probes = 0
        self.insertions = 0
        self.bit_sets = 0

    def positions(self, slope):
        first = mix64(slope + 0x9E3779B97F4A7C15)
        second = mix64(slope + 0xD1B54A32D192ED03) | 1
        mask = BLOOM_BITS - 1
        for index in range(BLOOM_HASHES):
            yield (first + index * second) & mask

    def contains(self, slope):
        self.queries += 1
        for position in self.positions(slope):
            self.hash_probes += 1
            if not (self.data[position >> 3] & (1 << (position & 7))):
                return False
        return True

    def insert(self, slope):
        self.insertions += 1
        for position in self.positions(slope):
            self.hash_probes += 1
            mask = 1 << (position & 7)
            byte = position >> 3
            if not self.data[byte] & mask:
                self.data[byte] |= mask
                self.bit_sets += 1


class ParentPageCache:
    def __init__(self):
        self.pages = OrderedDict()
        self.page_hits = 0
        self.page_misses = 0
        self.parent_queries = 0

    def append(self, record_index):
        page = record_index // PAGE_RECORDS
        self.pages[page] = None
        self.pages.move_to_end(page)
        while len(self.pages) > PARENT_CACHE_PAGES:
            self.pages.popitem(last=False)

    def batch_query(self, record_indices):
        self.parent_queries += len(record_indices)
        for page in sorted({index // PAGE_RECORDS for index in record_indices}):
            if page in self.pages:
                self.page_hits += 1
                self.pages.move_to_end(page)
            else:
                self.page_misses += 1
                self.pages[page] = None
                while len(self.pages) > PARENT_CACHE_PAGES:
                    self.pages.popitem(last=False)


def candidate_pairs(record_count, step, seed):
    all_pairs = record_count * (record_count + 1) // 2
    if all_pairs <= CANDIDATES_PER_STEP:
        return [(left, right) for left in range(record_count) for right in range(left, record_count)]
    pairs = set()
    attempt = 0
    recent_start = max(0, record_count - HOT_WINDOW)
    while len(pairs) < CANDIDATES_PER_STEP:
        raw_left = mix64(seed ^ (step << 32) ^ (attempt * 2 + 1))
        raw_right = mix64(seed ^ (step << 32) ^ (attempt * 2 + 2))
        lane = attempt % 4
        if lane == 0:
            left = record_count - 1
            right = raw_right % record_count
        elif lane == 1:
            left = recent_start + raw_left % (record_count - recent_start)
            right = raw_right % record_count
        elif lane == 2:
            left = raw_left % record_count
            right = raw_right % record_count
        else:
            left = recent_start + raw_left % (record_count - recent_start)
            right = recent_start + raw_right % (record_count - recent_start)
        pairs.add(tuple(sorted((left, right))))
        attempt += 1
    return sorted(pairs)


def slopes_to_archive(candidate, events):
    result = set()
    for point in set(events):
        dx = (candidate[0] - point[0]) % P
        if dx:
            result.add(direction((dx, (candidate[1] - point[1]) % P)))
    return result


def initialize_slopes(events):
    slopes = set()
    for index, point in enumerate(events):
        slopes.update(slopes_to_archive(point, events[:index]))
    return slopes


def run_guided(mode, seed):
    assert mode in ("bloom", "exact")
    events = [O, A, B]
    exact_slopes = initialize_slopes(events)
    bloom = SlopeBloom()
    for slope in exact_slopes:
        bloom.insert(slope)
    cache = ParentPageCache()
    for index in range(3):
        cache.append(index)

    estimated_gain_sum = 0
    exact_gain_sum = 0
    false_positive_queries = 0
    false_negative_queries = 0
    candidate_evaluations = 0
    archive_scan_records = 0
    parent_pair_records = 0
    chosen_pairs = []

    for step in range(BUDGET):
        pairs = candidate_pairs(len(events), step, seed)
        parent_indices = [index for pair in pairs for index in pair]
        parent_pair_records += len(parent_indices)
        cache.batch_query(parent_indices)
        archive_scan_records += len(events)
        scored = []
        for left, right in pairs:
            candidate = add(events[left], events[right])
            slopes = slopes_to_archive(candidate, events)
            estimated_gain = 0
            for slope in slopes:
                present = bloom.contains(slope)
                actually_present = slope in exact_slopes
                false_positive_queries += int(present and not actually_present)
                false_negative_queries += int((not present) and actually_present)
                estimated_gain += int(not present)
            exact_gain = len(slopes - exact_slopes)
            score = exact_gain if mode == "exact" else estimated_gain
            scored.append((score, estimated_gain, exact_gain, -left, -right, candidate, slopes, left, right))
            candidate_evaluations += 1
        chosen = max(scored)
        _, estimated_gain, exact_gain, _, _, candidate, slopes, left, right = chosen
        assert candidate == add(events[left], events[right])
        events.append(candidate)
        cache.append(len(events) - 1)
        chosen_pairs.append((left, right))
        estimated_gain_sum += estimated_gain
        exact_gain_sum += exact_gain
        new_slopes = slopes - exact_slopes
        exact_slopes.update(slopes)
        for slope in new_slopes:
            bloom.insert(slope)

    assert false_negative_queries == 0
    profile = direction_profile(events)
    assert len(exact_slopes) == profile["unsigned_useful_finite_slopes"]
    prefix = prefix_metrics(events)
    return {
        "mode": mode,
        "seed": seed,
        "charged_additions": BUDGET,
        "profile": profile,
        "prefix_metrics": {
            "covered_slope_discovery_weight": prefix["covered_slope_discovery_weight"],
            "conditional_mean_additions_for_covered_slopes": prefix["conditional_mean_additions_for_covered_slopes"],
            "exact_censored_mean_min_T_budget": prefix["exact_censored_mean_min_T_budget"],
            "censored_mean_numerator": prefix["censored_mean_numerator"],
            "coverage_curve_sha256": prefix["coverage_curve_sha256"],
        },
        "planner": {
            "candidate_evaluations": candidate_evaluations,
            "estimated_gain_sum_for_chosen_candidates": estimated_gain_sum,
            "exact_gain_sum_for_chosen_candidates": exact_gain_sum,
            "bloom_queries": bloom.queries,
            "bloom_hash_probes": bloom.hash_probes,
            "bloom_insertions": bloom.insertions,
            "bloom_bits_set": bloom.bit_sets,
            "false_positive_queries": false_positive_queries,
            "false_negative_queries": false_negative_queries,
            "archive_scan_records": archive_scan_records,
            "archive_scan_bytes": archive_scan_records * RECORD_BYTES,
            "parent_pair_record_queries": parent_pair_records,
            "parent_page_hits": cache.page_hits,
            "parent_page_misses": cache.page_misses,
        },
        "event_values_sha256": sha256(json.dumps(events, separators=(",", ":")).encode()).hexdigest(),
        "chosen_pairs_sha256": sha256(json.dumps(chosen_pairs, separators=(",", ":")).encode()).hexdigest(),
    }


def summary(values, candidate=None):
    mean = Fraction(sum(values), len(values))
    result = {
        "mean_exact": str(mean),
        "mean": float(mean),
        "sample_stdev": stdev(values),
        "minimum": min(values),
        "maximum": max(values),
    }
    if candidate is not None:
        result.update({
            "greater_than_comparator": sum(value > candidate for value in values),
            "equal_to_comparator": sum(value == candidate for value in values),
            "less_than_comparator": sum(value < candidate for value in values),
        })
    return result


def main():
    three_final, three_censored = strongest_three_line_control(BUDGET)
    guided_rows = []
    for trial in range(TRIALS):
        seed = 20260921 + trial
        guided_rows.append(run_guided("bloom", seed))
        guided_rows.append(run_guided("exact", seed))

    random_rows = []
    for trial in range(TRIALS):
        seed = 20260921 + 10_000 * BUDGET + trial
        events = random_control("random_accumulator", BUDGET, seed)
        profile = direction_profile(events)
        prefix = prefix_metrics(events)
        random_rows.append({
            "trial": trial,
            "seed": seed,
            "useful_finite_slopes": profile["unsigned_useful_finite_slopes"],
            "censored_mean_numerator": prefix["censored_mean_numerator"],
            "exact_censored_mean_min_T_budget": prefix["exact_censored_mean_min_T_budget"],
        })

    bloom_rows = [row for row in guided_rows if row["mode"] == "bloom"]
    exact_rows = [row for row in guided_rows if row["mode"] == "exact"]
    bloom_slopes = [row["profile"]["unsigned_useful_finite_slopes"] for row in bloom_rows]
    exact_slopes = [row["profile"]["unsigned_useful_finite_slopes"] for row in exact_rows]
    bloom_censored = [row["prefix_metrics"]["censored_mean_numerator"] for row in bloom_rows]
    exact_censored = [row["prefix_metrics"]["censored_mean_numerator"] for row in exact_rows]
    random_slopes = [row["useful_finite_slopes"] for row in random_rows]
    random_censored = [row["censored_mean_numerator"] for row in random_rows]

    large_sketch_bits = LARGE_SKETCH_BYTES * 8
    target_fpr = 0.01
    max_insertions_at_target_fpr = int(
        -(large_sketch_bits / BLOOM_HASHES)
        * log(1 - target_fpr ** (1 / BLOOM_HASHES))
    )
    approximate_records_before_quadratic_unique_slopes_reach_capacity = int(
        (1 + sqrt(1 + 8 * max_insertions_at_target_fpr)) / 2
    )
    fixed_ram = MODE_B_FIXED_RAM_WITHOUT_SLOPE_SKETCH + LARGE_SKETCH_BYTES
    max_records_by_ram = (RAM_LIMIT - fixed_ram) // MODE_B_RECORD_RAM_BYTES
    max_records_by_peak_disk = (DISK_LIMIT - 65_000_000_000) // 324
    max_records = min(max_records_by_ram, max_records_by_peak_disk)

    report = {
        "evidence_type": "finite_target_independent_coverage_guided_chain",
        "field_prime": P,
        "python_version": platform.python_version(),
        "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "parameters": {
            "charged_addition_budget": BUDGET,
            "candidate_pairs_per_step": CANDIDATES_PER_STEP,
            "hot_parent_window": HOT_WINDOW,
            "trials": TRIALS,
            "bloom_bits": BLOOM_BITS,
            "bloom_hashes": BLOOM_HASHES,
            "record_bytes": RECORD_BYTES,
            "page_records": PAGE_RECORDS,
            "parent_cache_pages": PARENT_CACHE_PAGES,
        },
        "mechanism": {
            "legal_step": "Propose archived parent pairs from a deterministic hot/global schedule; score coefficient sums; spend one group addition only on the chosen pair.",
            "target_independent": True,
            "coverage_state": "Bloom sketch contains every exactly covered finite slope; no false negatives in the verified finite runs.",
            "planning_scan": "One logical sequential scan of all archived coefficient records per step evaluates all candidates in a batch.",
            "difference_steps_enabled": False,
            "reason_difference_disabled": "The primary experiment stays in the non-inverting generic model; no free negation is assumed.",
        },
        "guided_trials": guided_rows,
        "random_accumulator_trials": random_rows,
        "summaries": {
            "bloom_final_slopes": summary(bloom_slopes, three_final["all_materialized_profile"]["unsigned_useful_finite_slopes"]),
            "exact_planner_final_slopes": summary(exact_slopes, three_final["all_materialized_profile"]["unsigned_useful_finite_slopes"]),
            "random_final_slopes": summary(random_slopes, three_final["all_materialized_profile"]["unsigned_useful_finite_slopes"]),
            "bloom_censored_mean_numerator": summary(bloom_censored, three_censored["prefix_metrics"]["censored_mean_numerator"]),
            "exact_planner_censored_mean_numerator": summary(exact_censored, three_censored["prefix_metrics"]["censored_mean_numerator"]),
            "random_censored_mean_numerator": summary(random_censored, three_censored["prefix_metrics"]["censored_mean_numerator"]),
        },
        "three_line_controls": {
            "final_coverage": three_final,
            "censored_mean": three_censored,
        },
        "large_capacity_envelope": {
            "slope_sketch_bytes": LARGE_SKETCH_BYTES,
            "fixed_ram_including_sketch_bytes": fixed_ram,
            "record_ram_bytes": MODE_B_RECORD_RAM_BYTES,
            "maximum_records_by_ram": max_records_by_ram,
            "maximum_records_by_peak_disk": max_records_by_peak_disk,
            "maximum_records_under_both_capacity_inequalities": max_records,
            "bloom_target_false_positive_probability_assumption": target_fpr,
            "bloom_hashes": BLOOM_HASHES,
            "maximum_slope_insertions_at_target_fpr_under_independent_uniform_hash_assumption": max_insertions_at_target_fpr,
            "approximate_records_if_every_pair_slope_is_unique_at_that_insertion_limit": approximate_records_before_quadratic_unique_slopes_reach_capacity,
            "quadratic_scan_bytes_for_N_records_formula": "64*N*(N+5) bytes for 128-byte records and N-3 planned additions.",
            "capacity_does_not_establish_planning_feasibility": True,
        },
        "verification": {
            "all_guided_steps_use_archived_parents": True,
            "all_selected_outputs_recomputed_as_parent_sums": True,
            "bloom_false_negatives": sum(row["planner"]["false_negative_queries"] for row in bloom_rows),
            "all_runs_match_addition_and_record_budget": True,
        },
        "limitations": [
            "Candidate-pair generation and Bloom parameters are finite design choices with no optimality proof.",
            "Bloom false-positive formulas assume independent uniform hashes; finite counters are logical, not hardware timings.",
            "Maintaining a no-false-negative all-slope sketch requires quadratic slope generation and archive scanning.",
            "The 1 TB / 100 TB inequalities bound resident records only and do not make quadratic planning traffic feasible.",
            "No novelty, ECDLP speedup, fixed-r rho improvement, hardware feasibility, or production claim.",
        ],
    }
    output = Path(__file__).with_name("coverage_guided_chain_checks.json")
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({
        "bloom_slopes": summary(bloom_slopes),
        "exact_planner_slopes": summary(exact_slopes),
        "random_slopes": summary(random_slopes),
        "three_line_slopes": three_final["all_materialized_profile"]["unsigned_useful_finite_slopes"],
        "output": str(output),
    }, indent=2))


if __name__ == "__main__":
    main()
