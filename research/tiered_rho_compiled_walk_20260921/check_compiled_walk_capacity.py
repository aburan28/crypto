"""Exact integer capacity envelope for compiled reservoir-guided schedules."""

from hashlib import sha256
from math import exp, isqrt, log
from pathlib import Path
import json


RAM_LIMIT = 10**12
DISK_LIMIT = 10**14
RAM_SKETCH_BYTES = 800_000_000_000
RAM_COMPILER_OTHER_BYTES = 200_000_000_000
EXTERNAL_SKETCH_BYTES = 100_000_000_000_000
BLOOM_HASHES = 7
TARGET_FPR = 0.01
ONLINE_RECORD_CAPACITY = 308_441_358_024
ONLINE_PEAK_STATIC_BYTES = 65_000_000_000
ONLINE_PEAK_PER_RECORD_BYTES = 324
PARENT_WINDOW = 64
FIXED_RESERVOIR = 32


def ceil_sqrt(value):
    root = isqrt(value)
    return root if root * root == value else root + 1


def schedule_steps(order):
    return (3 * ceil_sqrt(order) + 1) // 2


def scaled_reservoir(order):
    return max(FIXED_RESERVOIR, (ceil_sqrt(order) + 7) // 8)


def bloom_max_insertions(bytes_available):
    bits = bytes_available * 8
    return int(
        -(bits / BLOOM_HASHES)
        * log(1 - TARGET_FPR ** (1 / BLOOM_HASHES))
    )


def compiler_insertions(order, mode):
    reservoir = FIXED_RESERVOIR if mode == "fixed" else scaled_reservoir(order)
    return reservoir * schedule_steps(order)


def fits(order, mode, max_insertions):
    records = schedule_steps(order) + 3
    return (
        compiler_insertions(order, mode) <= max_insertions
        and records <= ONLINE_RECORD_CAPACITY
        and ONLINE_PEAK_STATIC_BYTES + ONLINE_PEAK_PER_RECORD_BYTES * records <= DISK_LIMIT
    )


def maximum_order(mode, max_insertions):
    low, high = 1, 1
    while fits(high, mode, max_insertions):
        low, high = high, high * 2
    while low + 1 < high:
        mid = (low + high) // 2
        if fits(mid, mode, max_insertions):
            low = mid
        else:
            high = mid
    return low


def row(mode, placement, sketch_bytes):
    max_insertions = bloom_max_insertions(sketch_bytes)
    order = maximum_order(mode, max_insertions)
    steps = schedule_steps(order)
    reservoir = FIXED_RESERVOIR if mode == "fixed" else scaled_reservoir(order)
    insertions = compiler_insertions(order, mode)
    schedule_bytes = (steps * (PARENT_WINDOW - 1).bit_length() + 7) // 8
    online_records = steps + 3
    peak_disk = ONLINE_PEAK_STATIC_BYTES + ONLINE_PEAK_PER_RECORD_BYTES * online_records
    return {
        "compiler_mode": mode,
        "sketch_placement": placement,
        "sketch_bytes": sketch_bytes,
        "bloom_hashes": BLOOM_HASHES,
        "target_fpr_assumption": TARGET_FPR,
        "maximum_insertions_at_target_fpr": max_insertions,
        "maximum_order": order,
        "maximum_order_log2_floor": order.bit_length() - 1,
        "schedule_steps": steps,
        "reservoir_records": reservoir,
        "compiler_insertions": insertions,
        "compiler_insertion_slack": max_insertions - insertions,
        "packed_six_bit_schedule_bytes": schedule_bytes,
        "online_exact_records": online_records,
        "online_peak_disk_bytes": peak_disk,
        "maximum_plus_one_fails_some_bound": not fits(order + 1, mode, max_insertions),
    }


def main():
    assert RAM_SKETCH_BYTES + RAM_COMPILER_OTHER_BYTES == RAM_LIMIT
    placements = []
    for mode in ("fixed", "scaled"):
        placements.append(row(mode, "800GB_RAM", RAM_SKETCH_BYTES))
        placements.append(row(mode, "100TB_external_capacity_only", EXTERNAL_SKETCH_BYTES))

    order_256 = 1 << 256
    hypothetical = {}
    for mode in ("fixed", "scaled"):
        needed = compiler_insertions(order_256, mode)
        hypothetical[mode] = {
            "reservoir_records": FIXED_RESERVOIR if mode == "fixed" else scaled_reservoir(order_256),
            "schedule_steps": schedule_steps(order_256),
            "compiler_insertions": needed,
            "over_800GB_ram_sketch_capacity_ratio": needed / bloom_max_insertions(RAM_SKETCH_BYTES),
            "over_100TB_external_sketch_capacity_ratio": needed / bloom_max_insertions(EXTERNAL_SKETCH_BYTES),
            "online_records_fit": schedule_steps(order_256) + 3 <= ONLINE_RECORD_CAPACITY,
        }

    report = {
        "evidence_type": "compiled_reservoir_guided_walk_capacity_envelope",
        "decimal_limits": {"ram_bytes": RAM_LIMIT, "peak_disk_bytes": DISK_LIMIT},
        "compiler_ram_placement": {
            "slope_sketch_bytes": RAM_SKETCH_BYTES,
            "other_runtime_and_buffers_bytes": RAM_COMPILER_OTHER_BYTES,
        },
        "online_record_model": {
            "record_capacity": ONLINE_RECORD_CAPACITY,
            "peak_static_bytes": ONLINE_PEAK_STATIC_BYTES,
            "peak_per_record_bytes": ONLINE_PEAK_PER_RECORD_BYTES,
            "schedule_bits_per_step": (PARENT_WINDOW - 1).bit_length(),
        },
        "placements": placements,
        "hypothetical_256_bit_order": hypothetical,
        "limitations": [
            "Bloom capacities assume independent uniform hashes and target occupancy.",
            "The external-sketch rows are capacity only; random bit probes over NVMe have no service guarantee.",
            "Compilation coefficient work, filter probes, and schedule search are not wall-clock measurements.",
            "The fixed reservoir fits larger orders but its measured relative gain shrinks with order.",
            "The scaled reservoir restores finite absolute gain by using linear-in-order compiler insertions.",
            "Capacity does not establish novelty, asymptotic speedup, or cryptographic-scale feasibility.",
        ],
        "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    output = Path(__file__).with_name("compiled_walk_capacity.json")
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"placements": placements, "output": str(output)}, indent=2))


if __name__ == "__main__":
    main()
