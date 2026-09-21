"""Finite public-synthetic multi-target collision/rank experiment.

Records represent X=aP+sum_j b_j Q_j in the additive group F_65537.  A
time/seed history-indexed legal addition chain is compared with K independent
single-target chains and an uncharged iid-coefficient sample reference.
"""

from fractions import Fraction
from hashlib import sha256
from math import ceil, sqrt
from pathlib import Path
from random import Random
from statistics import stdev
import json
import platform


P = 65537
K_VALUES = (1, 2, 4, 8, 16)
TRIALS = 32
RAM_LIMIT = 1_000_000_000_000
DISK_LIMIT = 100_000_000_000_000
FIXED_RAM = 369_310_378_496
PEAK_STATIC_DISK = 65_000_000_000


def mix64(value):
    value &= (1 << 64) - 1
    value ^= value >> 30
    value = (value * 0xBF58476D1CE4E5B9) & ((1 << 64) - 1)
    value ^= value >> 27
    value = (value * 0x94D049BB133111EB) & ((1 << 64) - 1)
    return value ^ (value >> 31)


def selector_index(seed, step, length):
    return mix64(seed ^ (step * 0x9E3779B97F4A7C15)) % length


def group_scalar(coefficients, logs):
    return (coefficients[0] + sum(b * x for b, x in zip(coefficients[1:], logs))) % P


class RankSystem:
    def __init__(self, k):
        self.k = k
        self.rows = {}
        self.field_additions = 0
        self.field_multiplications = 0
        self.field_inversions = 0
        self.dependent = 0

    def add_equation(self, coefficients, rhs):
        row = [value % P for value in coefficients]
        rhs %= P
        for pivot in sorted(self.rows):
            if not row[pivot]:
                continue
            existing, existing_rhs = self.rows[pivot]
            factor = row[pivot]
            for column in range(pivot, self.k):
                row[column] = (row[column] - factor * existing[column]) % P
                self.field_multiplications += 1
                self.field_additions += 1
            rhs = (rhs - factor * existing_rhs) % P
            self.field_multiplications += 1
            self.field_additions += 1
        pivot = next((column for column, value in enumerate(row) if value), None)
        if pivot is None:
            assert rhs == 0
            self.dependent += 1
            return False
        inverse = pow(row[pivot], -1, P)
        self.field_inversions += 1
        for column in range(pivot, self.k):
            row[column] = row[column] * inverse % P
            self.field_multiplications += 1
        rhs = rhs * inverse % P
        self.field_multiplications += 1
        self.rows[pivot] = (row, rhs)
        return True

    @property
    def rank(self):
        return len(self.rows)

    def solve(self):
        assert self.rank == self.k
        solution = [0] * self.k
        for pivot in sorted(self.rows, reverse=True):
            row, rhs = self.rows[pivot]
            value = rhs
            for column in range(pivot + 1, self.k):
                value = (value - row[column] * solution[column]) % P
                self.field_multiplications += 1
                self.field_additions += 1
            solution[pivot] = value
        return solution

    def counters(self):
        return {
            "field_additions": self.field_additions,
            "field_multiplications": self.field_multiplications,
            "field_inversions": self.field_inversions,
            "dependent_rows": self.dependent,
        }


def public_logs(k, seed):
    rng = Random(seed)
    values = []
    seen = {0, 1}
    while len(values) < k:
        value = rng.randrange(2, P)
        if value not in seen:
            seen.add(value)
            values.append(value)
    return values


def run_chain(logs, seed, cap):
    k = len(logs)
    zero = (0,) * (k + 1)
    generator = (1,) + (0,) * k
    targets = [
        (0,) + tuple(1 if column == target else 0 for column in range(k))
        for target in range(k)
    ]
    records = [zero, generator] + targets
    scalars = [group_scalar(record, logs) for record in records]
    directory = {}
    for index, scalar in enumerate(scalars):
        assert scalar not in directory
        directory[scalar] = index
    current = generator
    rank = RankSystem(k)
    collision_probes = len(records)
    exact_verification_reads = 0
    useful_rows = 0
    zero_coefficient_repeats = 0
    parent_reads = 0
    first_rank_times = []
    lineage = []
    current_index = 1

    for step in range(cap):
        parent_index = selector_index(seed, step, len(records))
        parent_reads += 1
        parent = records[parent_index]
        current = tuple((left + right) % P for left, right in zip(current, parent))
        scalar = group_scalar(current, logs)
        records.append(current)
        scalars.append(scalar)
        result_index = len(records) - 1
        lineage.append((current_index, parent_index, result_index))
        current_index = result_index
        collision_probes += 1
        if scalar in directory:
            prior_index = directory[scalar]
            exact_verification_reads += 1
            prior = records[prior_index]
            delta = tuple((left - right) % P for left, right in zip(current, prior))
            assert group_scalar(delta, logs) == 0
            if all(value == 0 for value in delta):
                zero_coefficient_repeats += 1
            else:
                useful_rows += 1
                previous_rank = rank.rank
                rank.add_equation(delta[1:], (-delta[0]) % P)
                if rank.rank > previous_rank:
                    first_rank_times.append(step + 1)
                if rank.rank == k:
                    solution = rank.solve()
                    assert solution == logs
                    additions = step + 1
                    break
        else:
            directory[scalar] = len(records) - 1
    else:
        additions = None

    replay = list(records[: k + 2])
    for parent_index, addend_index, result_index in lineage:
        assert result_index == len(replay)
        replay.append(tuple(
            (left + right) % P
            for left, right in zip(replay[parent_index], replay[addend_index])
        ))
    assert replay == records

    return {
        "k": k,
        "cap": cap,
        "completed": additions is not None,
        "group_additions": additions,
        "records_at_completion": None if additions is None else additions + k + 2,
        "collision_directory_probes": collision_probes,
        "exact_verification_reads": exact_verification_reads,
        "useful_collision_rows": useful_rows,
        "independent_rows": rank.rank,
        "pivot_columns": sorted(rank.rows),
        "free_columns": [column for column in range(k) if column not in rank.rows],
        "dependent_collision_rows": rank.dependent,
        "zero_coefficient_repeats": zero_coefficient_repeats,
        "rank_arrival_additions": first_rank_times,
        "parent_record_reads": parent_reads,
        "full_audit_replay_group_additions": len(lineage),
        "full_audit_replay_coefficient_field_additions": len(lineage) * (k + 1),
        "rank_field_operations": rank.counters(),
        "coefficient_vector_sha256": sha256(
            json.dumps(records, separators=(",", ":")).encode()
        ).hexdigest(),
        "lineage_sha256": sha256(
            json.dumps(lineage, separators=(",", ":")).encode()
        ).hexdigest(),
    }


def run_iid_reference(logs, seed, cap):
    k = len(logs)
    rng = Random(seed)
    directory = {}
    coefficients = []
    rank = RankSystem(k)
    useful = 0
    for sample in range(cap):
        record = tuple(rng.randrange(P) for _ in range(k + 1))
        scalar = group_scalar(record, logs)
        coefficients.append(record)
        if scalar in directory:
            prior = coefficients[directory[scalar]]
            delta = tuple((left - right) % P for left, right in zip(record, prior))
            if any(delta):
                useful += 1
                rank.add_equation(delta[1:], (-delta[0]) % P)
                if rank.rank == k:
                    assert rank.solve() == logs
                    return {
                        "samples": sample + 1,
                        "useful_collision_rows": useful,
                        "dependent_collision_rows": rank.dependent,
                    }
        else:
            directory[scalar] = sample
    return {"samples": None, "useful_collision_rows": useful, "dependent_collision_rows": rank.dependent}


def aligned_slot_bytes(k):
    raw = 33 + 32 * (k + 1) + 8 + 8 + 3 + 4
    return ((raw + 63) // 64) * 64


def capacity(k):
    slot = aligned_slot_bytes(k)
    directory_and_filter = 34
    final_per_record = slot + directory_and_filter
    peak_per_record = 2 * final_per_record
    rank_ram = 32 * k * (k + 1)
    active_coefficient_ram = 65_536 * 32 * (k + 1)
    fixed_ram = FIXED_RAM + rank_ram + active_coefficient_ram
    max_ram = (RAM_LIMIT - fixed_ram) // 2
    max_disk = (DISK_LIMIT - PEAK_STATIC_DISK) // peak_per_record
    maximum = min(max_ram, max_disk)
    return {
        "coefficient_bytes_per_record": 32 * (k + 1),
        "aligned_exact_record_slot_bytes": slot,
        "external_directory_and_filter_bytes_per_record": directory_and_filter,
        "final_bytes_per_record": final_per_record,
        "peak_bytes_per_record": peak_per_record,
        "rank_workspace_ram_bytes": rank_ram,
        "active_walker_coefficient_ram_bytes": active_coefficient_ram,
        "fixed_ram_bytes": fixed_ram,
        "maximum_records_by_ram": max_ram,
        "maximum_records_by_peak_disk": max_disk,
        "maximum_records": maximum,
        "binding_resource": "ram" if max_ram <= max_disk else "peak_disk",
        "final_disk_at_maximum_bytes": 1_000_000_000 + final_per_record * maximum,
        "peak_disk_at_maximum_bytes": PEAK_STATIC_DISK + peak_per_record * maximum,
        "coefficients_stored_exactly_no_replay_required": True,
        "parent_addend_indices_retained_for_audit_replay": True,
    }


def summarize(values):
    mean = Fraction(sum(values), len(values))
    return {
        "mean_exact": str(mean),
        "mean": float(mean),
        "sample_stdev": 0.0 if len(values) == 1 else stdev(values),
        "minimum": min(values),
        "maximum": max(values),
    }


def main():
    trial_rows = []
    summaries = []
    for k in K_VALUES:
        multi_values = []
        independent_values = []
        iid_values = []
        fallback_values = []
        rank_deficient_at_cutoff = 0
        scale = sqrt(k * P)
        cap = ceil(24 * scale)
        for trial in range(TRIALS):
            seed = 20260921 + 100_000 * k + trial
            logs = public_logs(k, seed)
            multi = run_chain(logs, seed ^ 0xA5A5A5A5, cap)
            assert multi["completed"]
            independent_runs = []
            for target, value in enumerate(logs):
                single = run_chain([value], seed ^ (0xC3C3C3C3 + target), cap)
                assert single["completed"]
                independent_runs.append(single)
            independent_total = sum(run["group_additions"] for run in independent_runs)
            cutoff = ceil(4 * scale)
            partial = run_chain(logs, seed ^ 0xA5A5A5A5, cutoff)
            if partial["completed"]:
                shared_plus_fallback = partial["group_additions"]
                fallback_target_indices = []
            else:
                rank_deficient_at_cutoff += 1
                fallback_target_indices = partial["free_columns"]
                shared_plus_fallback = cutoff + sum(
                    independent_runs[index]["group_additions"]
                    for index in fallback_target_indices
                )
            iid = run_iid_reference(logs, seed ^ 0x5A5A5A5A, cap)
            assert iid["samples"] is not None
            multi_values.append(multi["group_additions"])
            independent_values.append(independent_total)
            iid_values.append(iid["samples"])
            fallback_values.append(shared_plus_fallback)
            trial_rows.append({
                "k": k,
                "trial": trial,
                "seed": seed,
                "public_logs_sha256": sha256(json.dumps(logs).encode()).hexdigest(),
                "multi_target": multi,
                "independent_single_target_total_additions": independent_total,
                "independent_single_target_runs": independent_runs,
                "shared_cutoff_multiple_sqrt_kN": 4,
                "shared_cutoff_run": partial,
                "fallback_target_indices": fallback_target_indices,
                "shared_plus_charged_residual_fallback_additions": shared_plus_fallback,
                "iid_coefficient_sample_reference": iid,
            })
        multi_summary = summarize(multi_values)
        independent_summary = summarize(independent_values)
        iid_summary = summarize(iid_values)
        fallback_summary = summarize(fallback_values)
        summaries.append({
            "k": k,
            "trials": TRIALS,
            "sqrt_kN": scale,
            "sqrt_2kN": sqrt(2 * k * P),
            "cap_multiple_of_sqrt_kN": 24,
            "multi_target_additions": multi_summary,
            "independent_total_additions": independent_summary,
            "iid_coefficient_samples_not_additions": iid_summary,
            "shared_plus_charged_residual_fallback_additions": fallback_summary,
            "rank_deficient_shared_runs_at_4_sqrt_kN_cutoff": rank_deficient_at_cutoff,
            "mean_multi_over_sqrt_kN": multi_summary["mean"] / scale,
            "mean_independent_over_k_sqrt_N": independent_summary["mean"] / (k * sqrt(P)),
            "mean_speedup_independent_total_over_multi": independent_summary["mean"] / multi_summary["mean"],
            "mean_fallback_over_independent_total": fallback_summary["mean"] / independent_summary["mean"],
            "capacity_256_bit_hypothesis": capacity(k),
        })

    report = {
        "evidence_type": "finite_public_synthetic_multi_target_collision_rank_control",
        "field_prime": P,
        "k_values": K_VALUES,
        "trials_per_k": TRIALS,
        "python_version": platform.python_version(),
        "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "distribution": {
            "logs": "K distinct uniform public synthetic nonzero values in F_p, selected by fixed seeds.",
            "legal_chain": "Inputs O,P,Q_1,...,Q_K; current record adds one seed/time-selected prior archive record per charged group addition.",
            "iid_reference": "Independent uniform coefficient vectors in F_p^(K+1); samples are not legal one-addition outputs and are not counted as group additions.",
        },
        "summaries": summaries,
        "trials": trial_rows,
        "verification": {
            "every_collision_equation_checked_against_public_logs": True,
            "every_full_rank_solution_equals_public_logs": True,
            "all_multi_and_independent_runs_completed_within_cap": True,
            "all_iid_references_completed_within_cap": True,
        },
        "limitations": [
            "The history-indexed chain is time-dependent and memory-full, not a coalescing Pollard-rho function.",
            "Finite scaling does not prove the legal chain has iid collision or rank behavior asymptotically.",
            "The iid sample reference hides coefficient-materialization cost and is not a group-addition algorithm.",
            "Capacity arithmetic does not measure directory latency, archive bandwidth, coefficient-update cost, rank throughput, or replay performance.",
            "No private targets, key recovery, novelty, concrete ECDLP speedup, or hardware feasibility claim.",
        ],
    }
    output = Path(__file__).with_name("multi_target_linear_collision_checks.json")
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({
        "summaries": summaries,
        "output": str(output),
    }, indent=2))


if __name__ == "__main__":
    main()
