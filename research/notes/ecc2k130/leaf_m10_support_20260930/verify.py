#!/usr/bin/env python3
"""Independent full-row arithmetic and sampled point-law replay of leaf support."""
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
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
CONFIG = HERE / "CONFIG.json"
REF_PATH = ROOT / "research/notes/ecc2k130/rotated_unequal_arity_20260925/verify.py"
spec = importlib.util.spec_from_file_location("leaf_support_independent_field", REF_PATH)
assert spec is not None and spec.loader is not None
ref = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ref)
sys.path.insert(0, str(ROOT / "research/ecc2k130_relations"))
from fastfield import FastGF2m  # noqa: E402
from relations import Koblitz  # noqa: E402


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def peak_rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def checked_config() -> dict:
    config = json.loads(CONFIG.read_text())
    assert config["schema"] == "ecc2k130-leaf-m10-support-v1"
    for relative, expected in config["inputs_sha256"].items():
        assert sha(ROOT / relative) == expected, relative
    lock = json.loads((HERE / "FROZEN.json").read_text())
    assert lock["schema"] == "ecc2k130-leaf-m10-source-lock-v1"
    for relative, expected in lock["sha256"].items():
        assert sha(ROOT / relative) == expected, relative
    parent = json.loads((ROOT / "research/notes/ecc2k130/m10_export_capacity_20260925/INPUT.json").read_text())
    assert [arm["slot_dimensions"] for arm in config["arms"]] == [
        arm["dimensions"] for arm in parent["arms"]]
    certificate = json.loads((ROOT / "research/notes/ecc2k130/dual_all_lines_20260930/RESULT.json").read_text())
    assert certificate["status"] == "PASS_ALL_LINES"
    lines = {tuple(row["line"]): row["forward_b"]
             for row in certificate["runs"][0]["rows"]}
    for curve in config["curves"][1:]:
        assert curve["b_hex"] == lines[tuple(curve["line"])]
    return config


def linear_rank(values: list[int]) -> int:
    pivots = {}
    for value in values:
        while value:
            bit = value.bit_length() - 1
            if bit in pivots:
                value ^= pivots[bit]
            else:
                pivots[bit] = value
                break
    return len(pivots)


def independent_bases(field: FastGF2m, beta: int) -> tuple[dict[str, list[int]], str]:
    conjugates = []
    value = beta
    for _ in range(131):
        conjugates.append(value)
        value = field.sqr(value)
    assert value == beta and linear_rank(conjugates) == 131
    bases = {f"low_{i}": [conjugates[10 * j + i] for j in range(13)]
             for i in range(10)}
    bases["high_0"] = [conjugates[10 * j] for j in range(14)]
    assert all(linear_rank(basis) == len(basis) for basis in bases.values())
    assert linear_rank([v for i in range(10) for v in bases[f"low_{i}"]]) == 130
    assert linear_rank(bases["high_0"] + [v for i in range(1, 10)
                                            for v in bases[f"low_{i}"]]) == 131
    digest = hashlib.sha256("".join(f"{v}\n" for v in conjugates).encode("ascii")).hexdigest()
    return bases, digest


def x_from_mask(basis: list[int], mask: int) -> int:
    x = 0
    for bit, value in enumerate(basis):
        if mask & (1 << bit):
            x ^= value
    return x


def row_result(field: FastGF2m, a2: int, b: int, x: int) -> tuple[int, int]:
    if x == 0:
        return 1, -1
    inv_x = field.inv(x)
    inv_x2 = field.sqr(inv_x)
    b_over_x2 = field.mul(b, inv_x2)
    if field.trace(x ^ a2 ^ b_over_x2):
        return 0, -2
    twice_x = field.sqr(x) ^ b_over_x2
    if twice_x == 0:
        return 2, -1
    return 2, field.sqr(twice_x) ^ field.mul(
        b, field.sqr(field.inv(twice_x)))


def replay_scan(field: FastGF2m, alternate: ref.Field, curve: dict,
                slot: str, basis: list[int], folder: Path,
                config: dict) -> tuple[dict, set[int], int]:
    summary = json.loads((folder / "summary.json").read_text())
    chunks = json.loads((folder / "chunks.json").read_text())
    raw = folder / "rows.csv.gz"
    assert sha(raw) == summary["rows_gzip_sha256"]
    assert sha(folder / "chunks.json") == summary["chunks_sha256"]
    a2, b = int(curve["a2_hex"], 16), int(curve["b_hex"], 16)
    E = Koblitz(field, a=a2, b=b)
    columns = Counter()
    first: list[dict] = []
    last: deque[dict] = deque(maxlen=2)
    limit, chunk_rows = 1 << len(basis), config["gray_checkpoint_rows"]
    rows_hash, chunk_hash = hashlib.sha256(), hashlib.sha256()
    liftable = projected_infinity = sample_checks = 0
    chunk_liftable = 0
    expected_chunks = []
    with gzip.open(raw, "rb") as stream:
        for ordinal in range(limit):
            mask = ordinal ^ (ordinal >> 1)
            x = x_from_mask(basis, mask)
            assert x != 1
            lift, projected = row_result(field, a2, b, x)
            line = f"{mask},{x},{lift},{projected}\n".encode("ascii")
            assert stream.readline() == line, (curve["id"], slot, ordinal)
            rows_hash.update(line)
            chunk_hash.update(line)
            if x:
                if lift == 2:
                    liftable += 1
                    chunk_liftable += 1
                    sample = {"mask": mask, "x": x, "projected_x": projected}
                    if len(first) < 2:
                        first.append(sample)
                    last.append(sample)
                    if projected == -1:
                        projected_infinity += 1
                    else:
                        assert projected > 0
                        columns[projected] += 1
                        assert columns[projected] <= 4
            else:
                assert ordinal == mask == 0 and lift == 1 and projected == -1
            if (ordinal + 1) % chunk_rows == 0 or ordinal + 1 == limit:
                expected_chunks.append({
                    "start_ordinal": ordinal + 1 - (chunk_rows if (
                        ordinal + 1) % chunk_rows == 0 else (ordinal + 1) % chunk_rows),
                    "stop_ordinal": ordinal + 1,
                    "sha256": chunk_hash.hexdigest(),
                    "liftable_nonzero_x": chunk_liftable})
                chunk_hash = hashlib.sha256()
                chunk_liftable = 0
                assert peak_rss_bytes() <= config["process_rss_cap_bytes"]
        assert stream.read(1) == b""
    assert chunks == expected_chunks
    samples = first + list(last)
    assert samples == summary["sample_lifts"] and len(samples) == 4
    zero_points = E.points_over(0)
    assert len(zero_points) == 1 and E.on_curve(zero_points[0])
    assert E.mul(zero_points[0], 4) is None
    for sample in samples:
        x, projected = sample["x"], sample["projected_x"]
        ref_inverse = alternate.inverse(x)
        assert ref_inverse == field.inv(x)
        rhs = x ^ a2 ^ alternate.mul(b, alternate.square(ref_inverse))
        assert alternate.trace(rhs) == 0
        points = E.points_over(x)
        assert len(points) == 2 and points[0] != points[1]
        for point in points:
            assert E.on_curve(point)
            image = E.mul(point, 4)
            assert (-1 if image is None else image[0]) == projected
            sample_checks += 1
    hist = Counter(columns.values())
    digest = hashlib.sha256("".join(f"{value}\n" for value in sorted(columns)).encode("ascii")).hexdigest()
    expected = {"schema": "ecc2k130-leaf-m10-slot-v1",
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
                    str(i): hist[i] for i in range(1, 5)},
                "column_set_sha256": digest,
                "rows_sha256": rows_hash.hexdigest(),
                "rows_gzip_sha256": sha(raw),
                "chunks_sha256": sha(folder / "chunks.json"),
                "chunk_count": len(chunks), "sample_lifts": samples}
    for key, value in expected.items():
        assert summary[key] == value, (curve["id"], slot, key)
    return summary, set(columns), sample_checks


def replay(folder: Path) -> dict:
    config = checked_config()
    report = json.loads((folder / "result.json").read_text())
    assert report["schema"] == "ecc2k130-leaf-m10-support-result-v1"
    assert report["status"] == "PASS_CENSUS"
    assert report["config_sha256"] == sha(CONFIG)
    assert report["frozen_sha256"] == sha(HERE / "FROZEN.json")
    field = FastGF2m(131, int(config["field_modulus_hex"], 16))
    alternate = ref.Field(131, int(config["field_modulus_hex"], 16))
    bases, normal_hash = independent_bases(field, config["normal_beta"])
    assert report["normal_basis_sha256"] == normal_hash
    assert report["trace_mask_hex"] == hex(field.trace_bits)
    scans, sets = {}, {}
    checks = 0
    for curve in config["curves"]:
        curve_id = curve["id"]
        scans[curve_id], sets[curve_id] = {}, {}
        for slot in [f"low_{i}" for i in range(10)] + ["high_0"]:
            summary, columns, nchecks = replay_scan(
                field, alternate, curve, slot, bases[slot],
                folder / curve_id / slot, config)
            assert summary == report["scans"][curve_id][slot]
            scans[curve_id][slot] = summary
            sets[curve_id][slot] = columns
            checks += nchecks
            if curve_id == "source":
                expected = (16125, 8062) if slot == "high_0" else (7977, 3988)
                assert (summary["physical_points"], summary["projected_sign_classes"]) == expected
    assert set(report["scans"]) == set(scans)
    for curve in config["curves"]:
        curve_id = curve["id"]
        for arm in config["arms"]:
            slots = (["high_0"] if arm["id"] == "unequal" else ["low_0"])
            slots += [f"low_{i}" for i in range(1, 10)]
            counts = [scans[curve_id][slot]["physical_points"] for slot in slots]
            product = 1
            for count in counts:
                product *= count
            union = set().union(*(sets[curve_id][slot] for slot in slots))
            expected = {"slots": slots, "physical_point_counts": counts,
                        "physical_tuple_product": str(product),
                        "tuple_product_over_q_num": str(product),
                        "tuple_product_over_q_den": config["subgroup_order"],
                        "uncompressed_projected_sign_union": len(union),
                        "uncompressed_projected_union_sha256": hashlib.sha256(
                            "".join(f"{value}\n" for value in sorted(union)).encode("ascii")
                        ).hexdigest()}
            assert report["arms"][curve_id][arm["id"]] == expected
    assert report["PDP_yield"] is report["full_ECDLP_cost"] is report["method_crossover"] is None
    return {"schema": "ecc2k130-leaf-m10-support-replay-v1",
            "status": "PASS", "result_sha256": sha(folder / "result.json"),
            "config_sha256": sha(CONFIG), "frozen_sha256": sha(HERE / "FROZEN.json"),
            "scans_replayed": sum(map(len, scans.values())),
            "raw_masks_replayed": sum(summary["total_masks"] for cases in scans.values()
                                      for summary in cases.values()),
            "sample_point_images_checked": checks,
            "source_positive_controls": "PASS"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists(), "never overwrite a replay receipt"
    config = json.loads(CONFIG.read_text())
    started = time.perf_counter()
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(
        TimeoutError("frozen replay wall cap")))
    signal.alarm(config["verifier_wall_cap_seconds"])
    try:
        receipt = replay(args.evidence.resolve())
        assert peak_rss_bytes() <= config["process_rss_cap_bytes"]
        save(args.out, receipt)
        print(json.dumps(receipt, sort_keys=True))
    except BaseException as error:
        save(args.out, {"status": "CENSORED" if isinstance(error, (
            TimeoutError, MemoryError)) else "FAIL", "error_type": type(error).__name__,
            "error": str(error), "traceback": traceback.format_exc(),
            "elapsed_seconds": time.perf_counter() - started,
            "peak_rss_bytes": peak_rss_bytes()})
        raise
    finally:
        signal.alarm(0)


if __name__ == "__main__":
    main()
