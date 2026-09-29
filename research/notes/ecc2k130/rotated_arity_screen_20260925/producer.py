#!/usr/bin/env python3
"""Exact Gray-code F0 census for the frozen n131 rotated-arity arms."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import resource
import signal
import sys
import time
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent / "rotated_subspace_support_20260925/gate.py"
spec = importlib.util.spec_from_file_location("arity_parent_gate", PARENT)
assert spec is not None and spec.loader is not None
gate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gate)


def save(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n")


def peak_rss() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def run(arm: dict, data: dict, out: Path) -> None:
    m, d = arm["m"], arm["d"]
    started, cpu_started = time.perf_counter(), time.process_time()
    f = gate.Field(131, [0, 1, 2, 13])
    assert f.poly == data["field_poly"]
    f.rabin_prime_degree()
    assert gate.source_group_order(131) == data["group_order"] == 4 * data["q"]
    conjugates = gate.normal_conjugates(f, data["beta"])
    assert gate.rank(conjugates) == 131
    slots = gate.subspace_basis(conjugates, m, d)
    basis = slots[0]
    assert all(f.trace(b) == 1 for b in basis)
    assert gate.rank(basis + [1]) == d + 1
    trace_mask = gate.trace_mask(f)
    assert all(gate.trace_fast(trace_mask, b) == 1 for b in basis)
    setup_wall, setup_cpu = time.perf_counter() - started, time.process_time() - cpu_started
    setup_ops = dict(f.operations)

    columns: dict[int, int] = {}
    rows_hash = hashlib.sha256()
    chunk_hash = hashlib.sha256()
    x = 0
    zero_x = one_x = liftable = 0
    chunk_liftable = chunk_index = 0
    limit, chunk = 1 << d, data["chunk_rows"]
    with (out / "chunks.jsonl").open("w") as stream:
        for ordinal in range(limit):
            if ordinal:
                changed = ordinal & -ordinal
                x ^= basis[changed.bit_length() - 1]
            mask = ordinal ^ (ordinal >> 1)
            if x == 0:
                assert ordinal == 0 and mask == 0
                zero_x += 1
                lift, projected = 1, -1
            else:
                assert x != 1
                inverse_x = f.inverse(x)
                trace = (mask.bit_count() ^ (trace_mask & inverse_x).bit_count()) & 1
                if trace:
                    lift, projected = 0, -2
                else:
                    liftable += 1
                    chunk_liftable += 1
                    u = f.square(x ^ inverse_x)
                    assert u != 0
                    projected = f.square(u ^ f.inverse(u))
                    assert projected != 0
                    columns[projected] = columns.get(projected, 0) + 1
                    assert columns[projected] <= 4
                    lift = 2
            row = f"{mask},{x},{lift},{projected}\n".encode("ascii")
            rows_hash.update(row)
            chunk_hash.update(row)
            if (ordinal + 1) % chunk == 0 or ordinal + 1 == limit:
                assert peak_rss() <= data["caps"]["rss_bytes"]
                save_row = {"start_ordinal": chunk_index * chunk,
                            "stop_ordinal": ordinal + 1,
                            "row_sha256": chunk_hash.hexdigest(),
                            "liftable_nonzero_x": chunk_liftable,
                            "column_count_so_far": len(columns)}
                stream.write(json.dumps(save_row, sort_keys=True, separators=(",", ":")) + "\n")
                stream.flush()
                chunk_index += 1
                chunk_hash = hashlib.sha256()
                chunk_liftable = 0
    assert zero_x == 1 and one_x == 0
    assert sum(columns.values()) == liftable
    multiplicities = Counter(columns.values())
    assert not set(multiplicities) - {1, 2, 3, 4}
    c = len(columns)
    physical = 1 + 2 * liftable
    ordered_tuples = physical ** m
    threshold = data["threshold"]
    followup = (ordered_tuples * threshold["minimum_ordered_tuple_count_over_q_den"] >=
                data["q"] * threshold["minimum_ordered_tuple_count_over_q_num"] and
                c <= threshold["max_nonzero_signed_columns"])
    result = {"arm": arm, "domain": data["domain"], "beta": data["beta"],
              "field_poly": f.poly, "q": data["q"],
              "total_masks": limit, "zero_x_count": zero_x, "one_x_count": one_x,
              "liftable_nonzero_x": liftable, "nonzero_signed_columns": c,
              "physical_f0_points": physical,
              "distinct_projected_f0_points": 1 + 2 * c,
              "projected_signed_columns_including_O": 1 + c,
              "column_x_multiplicity_histogram": {str(i): multiplicities[i] for i in range(1, 5)},
              "ordered_physical_tuples": ordered_tuples,
              "necessary_support_upper_num": min(ordered_tuples, data["q"]),
              "necessary_support_upper_den": data["q"],
              "followup_screen_pass": followup,
              "row_sha256": rows_hash.hexdigest(), "chunk_count": chunk_index,
              "setup_field_operations": setup_ops,
              "total_field_operations": dict(f.operations),
              "setup_wall_seconds": setup_wall, "setup_cpu_seconds": setup_cpu,
              "total_wall_seconds": time.perf_counter() - started,
              "total_cpu_seconds": time.process_time() - cpu_started,
              "peak_rss_bytes": peak_rss()}
    assert result["total_wall_seconds"] <= data["caps"]["producer_wall_seconds"]
    assert result["peak_rss_bytes"] <= data["caps"]["rss_bytes"]
    save(out / "result.json", result)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--m", type=int, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    data = json.loads((HERE / "INPUT.json").read_text())
    arm = next((row for row in data["arms"] if row["m"] == args.m), None)
    assert arm is not None
    args.out.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter()
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(TimeoutError("producer wall cap")))
    signal.alarm(data["caps"]["producer_wall_seconds"])
    try:
        run(arm, data, args.out)
    except BaseException as error:
        save(args.out / "failure.json", {"arm": arm, "error": repr(error),
                                          "wall_seconds": time.perf_counter() - started,
                                          "peak_rss_bytes": peak_rss()})
        raise
    finally:
        signal.alarm(0)


if __name__ == "__main__":
    main()
