#!/usr/bin/env python3
"""Exact per-slot physical and [4]-projected support on two fixed 263-leaves."""
from __future__ import annotations

import argparse
from collections import Counter, deque
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import resource
import signal
import subprocess
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
CONFIG = HERE / "CONFIG.json"
GATE = ROOT / "research/notes/ecc2k130/rotated_subspace_support_20260925/gate.py"
spec = importlib.util.spec_from_file_location("leaf_support_gate", GATE)
assert spec is not None and spec.loader is not None
gate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gate)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()


def peak_rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def checked_inputs(require_lock: bool = True) -> dict:
    config = json.loads(CONFIG.read_text())
    assert config["schema"] == "ecc2k130-leaf-m10-support-v1"
    for relative, expected in config["inputs_sha256"].items():
        assert sha(ROOT / relative) == expected, relative
    parent = json.loads((ROOT / "research/notes/ecc2k130/m10_export_capacity_20260925/INPUT.json").read_text())
    assert config["field_degree"] == parent["field_degree"] == 131
    assert config["field_modulus_hex"] == parent["field_modulus_hex"]
    assert config["normal_beta"] == parent["normal_beta"] == 3
    assert [arm["slot_dimensions"] for arm in config["arms"]] == [
        arm["dimensions"] for arm in parent["arms"]]
    certificate = json.loads((ROOT / "research/notes/ecc2k130/dual_all_lines_20260930/RESULT.json").read_text())
    assert certificate["status"] == "PASS_ALL_LINES"
    rows = {tuple(row["line"]): row["forward_b"] for row in certificate["runs"][0]["rows"]}
    assert config["curves"][0]["b_hex"] == "0x1"
    for curve in config["curves"][1:]:
        assert curve["b_hex"] == rows[tuple(curve["line"])]
    if require_lock:
        lock = json.loads((HERE / "FROZEN.json").read_text())
        assert lock["schema"] == "ecc2k130-leaf-m10-source-lock-v1"
        for relative, expected in lock["sha256"].items():
            assert sha(ROOT / relative) == expected, relative
    return config


def slot_bases(field, beta: int) -> tuple[dict[str, list[int]], str]:
    conjugates = gate.normal_conjugates(field, beta)
    assert gate.rank(conjugates) == 131
    bases = {f"low_{i}": [conjugates[10 * j + i] for j in range(13)]
             for i in range(10)}
    bases["high_0"] = [conjugates[10 * j] for j in range(14)]
    assert all(gate.rank(basis) == len(basis) for basis in bases.values())
    assert all(gate.rank(basis + [1]) == len(basis) + 1
               for basis in bases.values())
    assert gate.rank([value for i in range(10) for value in bases[f"low_{i}"]]) == 130
    assert gate.rank(bases["high_0"] + [value for i in range(1, 10)
                                         for value in bases[f"low_{i}"]]) == 131
    basis_bytes = "".join(f"{value}\n" for value in conjugates).encode("ascii")
    return bases, hashlib.sha256(basis_bytes).hexdigest()


def lift_and_project(field, trace_mask: int, a2: int, b: int, x: int) -> tuple[int, int]:
    if x == 0:
        return 1, -1
    inv_x = field.inverse(x)
    b_over_x2 = field.mul(b, field.square(inv_x))
    if gate.trace_fast(trace_mask, x ^ a2 ^ b_over_x2):
        return 0, -2
    twice_x = field.square(x) ^ b_over_x2
    if twice_x == 0:
        return 2, -1
    projected = field.square(twice_x) ^ field.mul(
        b, field.square(field.inverse(twice_x)))
    return 2, projected


def scan(field, trace_mask: int, curve: dict, slot: str, basis: list[int],
         output: Path, config: dict) -> tuple[dict, set[int]]:
    assert not output.exists()
    output.mkdir(parents=True)
    a2, b = int(curve["a2_hex"], 16), int(curve["b_hex"], 16)
    started, cpu_started = time.perf_counter(), time.process_time()
    before = Counter(field.operations)
    limit, chunk_rows = 1 << len(basis), config["gray_checkpoint_rows"]
    rows_hash, chunk_hash = hashlib.sha256(), hashlib.sha256()
    columns: Counter[int] = Counter()
    first_samples: list[dict] = []
    last_samples: deque[dict] = deque(maxlen=2)
    x = liftable = projected_infinity = 0
    chunk_liftable = 0
    chunks = []
    raw_path = output / "rows.csv.gz"
    with raw_path.open("wb") as raw, gzip.GzipFile(
            filename="", fileobj=raw, mode="wb", mtime=0,
            compresslevel=9) as compressed:
        for ordinal in range(limit):
            if ordinal:
                bit = (ordinal & -ordinal).bit_length() - 1
                x ^= basis[bit]
            mask = ordinal ^ (ordinal >> 1)
            assert x != 1
            lift, projected = lift_and_project(field, trace_mask, a2, b, x)
            if x:
                if lift == 2:
                    liftable += 1
                    chunk_liftable += 1
                    sample = {"mask": mask, "x": x, "projected_x": projected}
                    if len(first_samples) < 2:
                        first_samples.append(sample)
                    last_samples.append(sample)
                    if projected == -1:
                        projected_infinity += 1
                    else:
                        assert projected > 0, "[4] image cannot be 2-torsion in order-4q group"
                        columns[projected] += 1
                        assert columns[projected] <= 4
            else:
                assert ordinal == 0 and mask == 0 and lift == 1
            row = f"{mask},{x},{lift},{projected}\n".encode("ascii")
            compressed.write(row)
            rows_hash.update(row)
            chunk_hash.update(row)
            if (ordinal + 1) % chunk_rows == 0 or ordinal + 1 == limit:
                compressed.flush()
                assert peak_rss_bytes() <= config["process_rss_cap_bytes"]
                chunks.append({"start_ordinal": ordinal + 1 - (
                                   chunk_rows if (ordinal + 1) % chunk_rows == 0
                                   else (ordinal + 1) % chunk_rows),
                               "stop_ordinal": ordinal + 1,
                               "sha256": chunk_hash.hexdigest(),
                               "liftable_nonzero_x": chunk_liftable})
                save(output / "chunks.json", chunks)
                chunk_hash = hashlib.sha256()
                chunk_liftable = 0
    assert len(first_samples) == len(last_samples) == 2
    column_bytes = "".join(f"{value}\n" for value in sorted(columns)).encode("ascii")
    histogram = Counter(columns.values())
    summary = {"schema": "ecc2k130-leaf-m10-slot-v1",
               "curve": curve["id"], "slot": slot, "dimension": len(basis),
               "basis_sha256": hashlib.sha256("".join(
                   f"{value}\n" for value in basis).encode("ascii")).hexdigest(),
               "total_masks": limit, "zero_x_count": 1,
               "liftable_nonzero_x": liftable,
               "nonliftable_nonzero_x": limit - 1 - liftable,
               "physical_points": 1 + 2 * liftable,
               "projected_infinity_nonzero_x": projected_infinity,
               "projected_sign_classes": len(columns),
               "column_x_multiplicity_histogram": {
                   str(i): histogram[i] for i in range(1, 5)},
               "column_set_sha256": hashlib.sha256(column_bytes).hexdigest(),
               "rows_sha256": rows_hash.hexdigest(),
               "rows_gzip_sha256": sha(raw_path),
               "chunks_sha256": sha(output / "chunks.json"),
               "chunk_count": len(chunks),
               "sample_lifts": first_samples + list(last_samples),
               "field_operations": dict(field.operations - before),
               "wall_seconds": time.perf_counter() - started,
               "cpu_seconds": time.process_time() - cpu_started,
               "peak_rss_bytes": peak_rss_bytes()}
    save(output / "summary.json", summary)
    return summary, set(columns)


def run(output: Path) -> dict:
    config = checked_inputs()
    field = gate.Field(131, [0, 1, 2, 13])
    assert field.poly == int(config["field_modulus_hex"], 16)
    field.rabin_prime_degree()
    bases, normal_hash = slot_bases(field, config["normal_beta"])
    trace_mask = gate.trace_mask(field)
    scans, column_sets, arms = {}, {}, {}
    for curve in config["curves"]:
        curve_id = curve["id"]
        scans[curve_id], column_sets[curve_id] = {}, {}
        for slot in [f"low_{i}" for i in range(10)] + ["high_0"]:
            summary, columns = scan(field, trace_mask, curve, slot, bases[slot],
                                    output / curve_id / slot, config)
            scans[curve_id][slot] = summary
            column_sets[curve_id][slot] = columns
            if curve_id == "source":
                expected = (16125, 8062) if slot == "high_0" else (7977, 3988)
                assert (summary["physical_points"], summary["projected_sign_classes"]) == expected
        arms[curve_id] = {}
        for arm in config["arms"]:
            slots = (["high_0"] if arm["id"] == "unequal" else ["low_0"])
            slots += [f"low_{i}" for i in range(1, 10)]
            counts = [scans[curve_id][slot]["physical_points"] for slot in slots]
            union = set().union(*(column_sets[curve_id][slot] for slot in slots))
            product = 1
            for count in counts:
                product *= count
            arms[curve_id][arm["id"]] = {
                "slots": slots, "physical_point_counts": counts,
                "physical_tuple_product": str(product),
                "tuple_product_over_q_num": str(product),
                "tuple_product_over_q_den": config["subgroup_order"],
                "uncompressed_projected_sign_union": len(union),
                "uncompressed_projected_union_sha256": hashlib.sha256(
                    "".join(f"{value}\n" for value in sorted(union)).encode("ascii")
                ).hexdigest()}
    result = {"schema": "ecc2k130-leaf-m10-support-result-v1",
              "status": "PASS_CENSUS", "source_head": subprocess.check_output(
                  ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
              "config_sha256": sha(CONFIG),
              "frozen_sha256": sha(HERE / "FROZEN.json"),
              "normal_basis_sha256": normal_hash,
              "trace_mask_hex": hex(trace_mask),
              "scans": scans, "arms": arms,
              "PDP_yield": None, "full_ECDLP_cost": None,
              "method_crossover": None,
              "peak_rss_bytes": peak_rss_bytes()}
    save(output / "result.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists(), "never overwrite a support census"
    args.out.mkdir(parents=True)
    config = json.loads(CONFIG.read_text())
    started = time.perf_counter()
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(
        TimeoutError("frozen producer wall cap")))
    signal.alarm(config["producer_wall_cap_seconds"])
    try:
        result = run(args.out)
        assert peak_rss_bytes() <= config["process_rss_cap_bytes"]
        print(json.dumps({"status": result["status"],
                          "scans": sum(map(len, result["scans"].values()))},
                         sort_keys=True))
    except BaseException as error:
        save(args.out / "failure.json", {"status": "CENSORED" if isinstance(
            error, (TimeoutError, MemoryError)) else "FAIL",
            "error_type": type(error).__name__, "error": str(error),
            "traceback": traceback.format_exc(),
            "elapsed_seconds": time.perf_counter() - started,
            "peak_rss_bytes": peak_rss_bytes()})
        raise
    finally:
        signal.alarm(0)


if __name__ == "__main__":
    main()
