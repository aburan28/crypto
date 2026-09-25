#!/usr/bin/env python3
"""Exact Gray-code census of rational F0 x and canonical signed [4] columns.

The protocol and frozen hashes in this directory govern measured use.
"""
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
PARENT = HERE.parent / "rotated_subspace_support_20260925" / "gate.py"
spec = importlib.util.spec_from_file_location("rotated_parent_gate", PARENT)
assert spec is not None and spec.loader is not None
parent = importlib.util.module_from_spec(spec)
spec.loader.exec_module(parent)

CAP_WALL = {"pilot": 90, "full": 1800}
CAP_RSS = 768 * 1024 * 1024
CHUNK = 8192


def compact(obj: object) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"))


def save(path: Path, obj: object) -> None:
    path.write_text(compact(obj) + "\n")


def rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def coordinates(basis: list[int]):
    """Return exact F2 coordinates in the independent F0 basis, or None."""
    rows: dict[int, tuple[int, int]] = {}
    for j, value in enumerate(basis):
        v, coeff = value, 1 << j
        while v:
            pivot = v.bit_length() - 1
            if pivot in rows:
                row, row_coeff = rows[pivot]
                v ^= row
                coeff ^= row_coeff
            else:
                rows[pivot] = (v, coeff)
                break
        else:
            raise AssertionError("dependent F0 basis")

    def decode(value: int) -> int | None:
        coeff = 0
        while value:
            pivot = value.bit_length() - 1
            if pivot not in rows:
                return None
            row, row_coeff = rows[pivot]
            value ^= row
            coeff ^= row_coeff
        return coeff

    return decode


def gray_ordinal(mask: int) -> int:
    """Invert reflected Gray code to its unique traversal ordinal."""
    ordinal = 0
    while mask:
        ordinal ^= mask
        mask >>= 1
    return ordinal


def setup(data: dict):
    model = next(m for m in data["models"] if m["name"] == "n131-full")
    assert (model["n"], model["m"], model["d"], model["beta"]) == (131, 6, 21, 3)
    f = parent.Field(model["n"], [0, 1, 2, 13])
    assert f.poly == model["poly"]
    f.rabin_prime_degree()
    assert parent.source_group_order(131) == 4 * parent.Q131
    conjugates = parent.normal_conjugates(f, 3)
    assert len(conjugates) == 131 and parent.rank(conjugates) == 131
    assert parent.rank([conjugates[6 * j] for j in range(21)]) == 21
    assert parent.rank([conjugates[6 * j + i] for i in range(6) for j in range(21)]) == 126
    assert f.trace(3) == 1 and parent.rank(conjugates) == 131
    assert 1 == __import__("functools").reduce(int.__xor__, conjugates, 0)
    basis = [conjugates[6 * j] for j in range(21)]
    decode = coordinates(basis)
    # The all-one normal coordinates of 1 cannot fit in the 21 selected slots.
    assert decode(1) is None
    tmask = parent.trace_mask(f)
    assert all(parent.trace_fast(tmask, b) == 1 for b in basis)
    curve = parent.Curve(f)
    anchor = parent.n131_projected_lambda_point(curve, tmask)
    assert anchor["group_order"] == 4 * parent.Q131
    return f, curve, basis, decode, tmask, anchor


def run(data: dict, stage: str, out: Path) -> dict:
    limit = data["pilot_masks"] if stage == "pilot" else data["full_masks"]
    assert limit == (1 << 15 if stage == "pilot" else 1 << 21)
    started = time.perf_counter()
    cpu_started = time.process_time()
    f, curve, basis, decode, tmask, anchor = setup(data)
    setup_done = time.perf_counter()
    setup_cpu_done = time.process_time()
    setup_field_ops = f.operations.copy()
    setup_curve_ops = curve.operations.copy()
    columns: dict[int, int] = {}
    row_hash = hashlib.sha256()
    chunk_hash = hashlib.sha256()
    chunk_lifts = 0
    chunks = 0
    x = 0
    zero_x = one_x = liftable = inverse_pairs = 0
    chunk_path = out / "chunks.jsonl"
    with chunk_path.open("w") as stream:
        for ordinal in range(limit):
            if ordinal:
                changed = ordinal & -ordinal
                x ^= basis[changed.bit_length() - 1]
            mask = ordinal ^ (ordinal >> 1)
            if x == 0:
                assert ordinal == 0 and mask == 0
                zero_x += 1
                lift, column = 1, -1
            elif x == 1:
                # Explicit torsion branch: forbidden by the normal-basis proof.
                one_x += 1
                raise AssertionError("x=1 entered F0")
            else:
                inverse_x = f.inverse(x)
                # Tr(x + x^-2) = Tr(x + x^-1); Tr(x) is mask parity.
                solvable = ((mask.bit_count() ^ (tmask & inverse_x).bit_count()) & 1) == 0
                if not solvable:
                    lift, column = 0, -2
                else:
                    liftable += 1
                    chunk_lifts += 1
                    u = f.square(x ^ inverse_x)
                    assert u != 0
                    column = f.square(u ^ f.inverse(u))
                    assert column != 0
                    columns[column] = columns.get(column, 0) + 1
                    assert columns[column] <= 4, (mask, x, column, columns[column])
                    lift = 2
                    inverse_mask = decode(inverse_x)
                    if (x < inverse_x and inverse_mask is not None
                            and gray_ordinal(inverse_mask) < limit):
                        inverse_pairs += 1
            row = f"{mask},{x},{lift},{column}\n".encode("ascii")
            row_hash.update(row)
            chunk_hash.update(row)
            if (ordinal + 1) % CHUNK == 0 or ordinal + 1 == limit:
                if rss_bytes() > CAP_RSS:
                    raise MemoryError(f"peak RSS exceeded {CAP_RSS} bytes")
                save_row = {"ordinal_start": chunks * CHUNK,
                            "ordinal_stop": ordinal + 1,
                            "sha256": chunk_hash.hexdigest(),
                            "liftable_nonzero_x": chunk_lifts,
                            "distinct_columns_so_far": len(columns),
                            "wall_seconds_so_far": time.perf_counter() - started,
                            "peak_rss_bytes_so_far": rss_bytes()}
                stream.write(compact(save_row) + "\n")
                stream.flush()
                chunks += 1
                chunk_hash = hashlib.sha256()
                chunk_lifts = 0
    assert zero_x == 1 and one_x == 0
    assert sum(columns.values()) == liftable
    histogram = Counter(columns.values())
    assert not set(histogram) - {1, 2, 3, 4}
    assert len(columns) >= (liftable + 3) // 4
    assert len(columns) <= liftable
    assert inverse_pairs <= liftable - len(columns)
    result = {
        "stage": stage,
        "model": {"n": 131, "m": 6, "d": 21, "beta": 3, "poly": f.poly},
        "total_masks": limit,
        "zero_x_count": zero_x,
        "one_x_count": one_x,
        "liftable_nonzero_x": liftable,
        "nonzero_signed_columns": len(columns),
        "zero_column_count": 1,
        "physical_factor_points": 1 + 2 * liftable,
        "distinct_projected_points": 1 + 2 * len(columns),
        "total_signed_columns_including_O": 1 + len(columns),
        "multiplicity_histogram": {str(i): histogram[i] for i in range(1, 5)},
        "max_column_preimages": max(histogram, default=0),
        "inverse_pair_collisions": inverse_pairs,
        "row_sha256": row_hash.hexdigest(),
        "chunk_count": chunks,
        "setup_anchor": anchor,
        "field_operations_setup": dict(setup_field_ops),
        "curve_operations_setup": dict(setup_curve_ops),
        "field_operations_enumeration": {k: f.operations[k] - setup_field_ops[k]
                                         for k in sorted(set(f.operations) | set(setup_field_ops))},
        "curve_operations_enumeration": {k: curve.operations[k] - setup_curve_ops[k]
                                         for k in sorted(set(curve.operations) | set(setup_curve_ops))},
        "setup_wall_seconds": setup_done - started,
        "setup_cpu_seconds": setup_cpu_done - cpu_started,
        "enumeration_wall_seconds": time.perf_counter() - setup_done,
        "enumeration_cpu_seconds": time.process_time() - setup_cpu_done,
        "total_wall_seconds": time.perf_counter() - started,
        "total_cpu_seconds": time.process_time() - cpu_started,
        "peak_rss_bytes": rss_bytes(),
    }
    assert result["total_wall_seconds"] <= CAP_WALL[stage]
    assert result["peak_rss_bytes"] <= CAP_RSS
    save(out / "summary.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--stage", choices=["pilot", "full"], required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    data = json.loads(args.input.read_text())
    assert data["domain"] == "ECC2K130-N131-BETA3-F0-CENSUS-20260925-v1"
    args.out.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter()
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(TimeoutError("wall cap")))
    signal.alarm(CAP_WALL[args.stage])
    try:
        result = run(data, args.stage, args.out)
        print(compact({"stage": args.stage,
                       "L": result["liftable_nonzero_x"],
                       "C": result["nonzero_signed_columns"],
                       "wall_seconds": result["total_wall_seconds"]}), flush=True)
    except BaseException as error:
        save(args.out / "failure.json", {"stage": args.stage, "error": repr(error),
                                          "elapsed_wall_seconds": time.perf_counter() - started,
                                          "peak_rss_bytes": rss_bytes()})
        raise
    finally:
        signal.alarm(0)


if __name__ == "__main__":
    main()
