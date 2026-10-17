"""Audit the complete quotient reduction after independent GP replay."""

import collections
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verify_adaptive_private_columns import cell, mul, on_curve  # noqa: E402


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


points_path = HERE / "type_i_line_full_new_quotient_points.json"
summary_path = HERE / "type_i_line_new_all_reduction.json"
prepared_path = HERE / "type_i_line_new_all_fixed_origin_input.json"
runtime_path = HERE / "type_i_line_new64_runtime_info.json"
density_path = HERE / "type_i_character_density.json"
points = json.loads(points_path.read_text())
summary = json.loads(summary_path.read_text())
prepared = json.loads(prepared_path.read_text())
runtime = json.loads(runtime_path.read_text())
density = json.loads(density_path.read_text())
assert len(points) == summary["input_classes"] == summary["reduced_classes"] == 9784
assert summary["status"] == "all_9784_new_quotient_classes_reduced"
assert summary["input_sha256"] == sha(points_path)
assert summary["source_sha256"] == sha(HERE / "line_new_quotient_all_driver.py")
assert summary["worker_source_sha256"] == sha(HERE / "line_new_quotient_reduction.py")
assert summary["runtime_info_sha256"] == sha(runtime_path)
assert runtime["status"] == "verified"
assert prepared["cases"] == prepared["nonspecial_representatives"] == 9784
assert prepared["source_sha256"] == sha(HERE / "prepare_line_new_all_fixed_origin_gp.py")
assert prepared["summary_sha256"] == sha(summary_path)
assert prepared["gp_input_sha256"] == sha(HERE / "type_i_line_new_all_fixed_origin_input.gp")
assert prepared["points_sha256"] == sha(points_path)
assert prepared["plane_sha256"] == sha(HERE / "type_i_plane_quartic.json")
native_log_path = HERE / "type_i_line_new_all_fixed_origin_native_verification.log"
native_error_path = HERE / "type_i_line_new_all_fixed_origin_native_verification.stderr"
native_log = native_log_path.read_text()
assert native_log.splitlines()[-1] == "NEW_ALL_LINE_FIXED_ORIGIN_NATIVE_PASS cases=9784 prime_bits=192"
assert "***" not in native_log and not native_error_path.read_text()
assert "NEW_ALL_REDUCTION_PASS" in (HERE / "type_i_line_new_all_driver.log").read_text()

p = int(summary["field_prime"])
r = int("6277101735386680763835789423176059013767194773182842284081")
coefficients = [int(value) for value in density["sextic_coefficients_decimal"]]
sampled_indices = set(range(0, 9784, 153)) | {0, 9783}
cells = collections.Counter()
next_index = 0
for chunk in summary["chunks"]:
    start, count = chunk["start"], chunk["count"]
    assert start == next_index
    path = HERE / "type_i_line_new_all_chunks" / f"{start:05d}-{start+count-1:05d}.json"
    assert chunk["receipt_sha256"] == sha(path)
    log_path = path.with_suffix(".log")
    assert "NEW_QUOTIENT_REDUCTION_PASS" in log_path.read_text()
    assert not path.with_suffix(".stderr").read_text()
    data = json.loads(path.read_text())
    assert data["start_index"] == start
    assert data["input_classes"] == data["reduced_classes"] == count
    assert data["unreduced_classes"] == 0
    assert data["source_sha256"] == sha(HERE / "line_new_quotient_reduction.py")
    assert data["input_sha256"] == sha(points_path)
    assert data["plane_receipt_sha256"] == sha(HERE / "type_i_plane_quartic.json")
    assert data["density_receipt_sha256"] == sha(density_path)
    assert data["quotient_constructor_sha256"] == sha(HERE / "line_quotient_cubic_probe.py")
    assert data["reduction_constructor_sha256"] == sha(HERE / "line_quotient_fixed_origin.py")
    assert len(data["rows"]) == len(data["representatives"]) == count
    for offset, (row, representative) in enumerate(zip(
            data["rows"], data["representatives"])):
        index = start + offset
        assert row["index"] == representative["index"] == index
        assert row["quotient"] == representative["quotient"]
        q = tuple(map(int, row["quotient"]))
        assert list(q) == points[index]
        assert row["status"] == "verified_prime_field_cubic_pair"
        assert representative["status"] == "verified_fixed_origin_degree_three"
        if index in sampled_indices:
            assert on_curve(q, p) and mul(r, q, p) is None
        half = mul(pow(2, -1, r), q, p) if index in sampled_indices else None
        if half is not None:
            assert row["cell"] == cell(half, coefficients, p)[::-1]
        cells[row["cell"]] += 1
        y = representative["representative_y_polynomial"]
        assert len(y) == 3 and int(y[2][0]) != 0
    next_index += count
assert next_index == 9784 and sum(cells.values()) == 9784
assert dict(sorted(cells.items())) == prepared["cells"]

out = {
    "status": "verified_all_new_quotient_reduction_receipt",
    "classes": 9784,
    "nonspecial_representatives": 9784,
    "sampled_subgroup_and_cell_checks": len(sampled_indices),
    "cells": dict(sorted(cells.items())),
    "sage_summary_sha256": sha(summary_path),
    "gp_input_sha256": sha(HERE / "type_i_line_new_all_fixed_origin_input.gp"),
    "gp_verification_sha256": sha(native_log_path),
    "runtime_info_sha256": sha(runtime_path),
    "source_sha256": sha(Path(__file__)),
}
(HERE / "type_i_line_new_all_receipt_verification.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n")
print("NEW_ALL_RECEIPT_VERIFICATION_PASS", json.dumps(out, sort_keys=True))
