"""Audit custody and nonspeciality for the new quotient reduction slice."""

import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verify_adaptive_private_columns import cell, mul, on_curve  # noqa: E402


def sha(name):
    return hashlib.sha256((HERE / name).read_bytes()).hexdigest()


new = json.loads((HERE / "type_i_line_full_new_quotient_points.json").read_text())
reduction = json.loads((HERE / "type_i_line_new64_fixed_origin.json").read_text())
prepared = json.loads((HERE / "type_i_line_new64_fixed_origin_input.json").read_text())
runtime = json.loads((HERE / "type_i_line_new64_runtime_info.json").read_text())
density = json.loads((HERE / "type_i_character_density.json").read_text())
p, r = int(reduction["field_prime"]), int(reduction["subgroup_order"])
coefficients = [int(value) for value in density["sextic_coefficients_decimal"]]
assert runtime["status"] == "verified"
assert reduction["status"] == "bounded_new_quotient_fixed_origin_reduction"
assert reduction["start_index"] == 0
assert reduction["input_classes"] == reduction["reduced_classes"] == 64
assert reduction["unreduced_classes"] == 0
assert len(reduction["rows"]) == len(reduction["representatives"]) == 64
assert reduction["source_sha256"] == sha("line_new_quotient_reduction.py")
assert reduction["input_sha256"] == sha("type_i_line_full_new_quotient_points.json")
assert reduction["plane_receipt_sha256"] == sha("type_i_plane_quartic.json")
assert reduction["density_receipt_sha256"] == sha("type_i_character_density.json")
assert reduction["quotient_constructor_sha256"] == sha("line_quotient_cubic_probe.py")
assert reduction["reduction_constructor_sha256"] == sha("line_quotient_fixed_origin.py")
assert prepared["reduction_sha256"] == sha("type_i_line_new64_fixed_origin.json")
assert prepared["gp_input_sha256"] == sha("type_i_line_new64_fixed_origin_input.gp")
assert prepared["plane_sha256"] == sha("type_i_plane_quartic.json")
assert prepared["source_sha256"] == sha("prepare_line_new_fixed_origin_gp.py")
native_log = (HERE / "type_i_line_new64_fixed_origin_native_verification.log").read_text()
assert "NEW_LINE_FIXED_ORIGIN_NATIVE_PASS cases=64 prime_bits=192" in native_log
assert "***" not in native_log
assert "NEW_QUOTIENT_REDUCTION_PASS" in (HERE / "type_i_line_new64_fixed_origin.log").read_text()

cells = Counter()
for i, (input_row, representative) in enumerate(zip(
        reduction["rows"], reduction["representatives"])):
    assert input_row["index"] == representative["index"] == i
    assert input_row["quotient"] == representative["quotient"]
    assert list(map(int, input_row["quotient"])) == new[i]
    assert input_row["status"] == "verified_prime_field_cubic_pair"
    assert representative["status"] == "verified_fixed_origin_degree_three"
    Q = tuple(map(int, input_row["quotient"]))
    assert on_curve(Q, p) and mul(r, Q, p) is None
    H = mul(pow(2, -1, r), Q, p)
    assert H is not None
    assert input_row["cell"] == cell(H, coefficients, p)[::-1]
    cells[input_row["cell"]] += 1
    y = representative["representative_y_polynomial"]
    assert len(y) == 3 and int(y[2][0]) != 0
assert sum(cells.values()) == 64

out = {
    "status": "verified_new64_quotient_reduction_receipt",
    "classes": 64,
    "nonspecial_representatives": 64,
    "cells": dict(sorted(cells.items())),
    "sage_receipt_sha256": sha("type_i_line_new64_fixed_origin.json"),
    "gp_input_sha256": sha("type_i_line_new64_fixed_origin_input.gp"),
    "gp_verification_sha256": sha("type_i_line_new64_fixed_origin_native_verification.log"),
    "runtime_info_sha256": sha("type_i_line_new64_runtime_info.json"),
    "source_sha256": sha("verify_line_new64_receipt.py"),
}
(HERE / "type_i_line_new64_receipt_verification.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n")
print("NEW64_RECEIPT_VERIFICATION_PASS", json.dumps(out, sort_keys=True))
