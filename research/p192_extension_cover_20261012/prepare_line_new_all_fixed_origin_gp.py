"""Serialize all new quotient cubic cuts for independent PARI replay."""

import hashlib
import json
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
summary_path = HERE / "type_i_line_new_all_reduction.json"
summary = json.loads(summary_path.read_text())
plane_path = HERE / "type_i_plane_quartic.json"
plane = json.loads(plane_path.read_text())
point_path = HERE / "type_i_line_full_new_quotient_points.json"
points = json.loads(point_path.read_text())
assert summary["status"] == "all_9784_new_quotient_classes_reduced"
assert len(points) == summary["input_classes"] == summary["reduced_classes"] == 9784


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coeff(encoded):
    assert all(len(item) == 1 for item in encoded)
    return [int(item[0]) for item in encoded]


terms = [[*map(int, label.split(",")[:2]), int(value)]
         for label, value in plane["quartic_coefficients"].items()]
anchors = [[int(v) for v in item["coordinates_Fp_decimal"]]
           for item in plane["known_points_verified"] if item["source"] == "base"][:3]
exponents = [[3, 0], [0, 3], [0, 0], [2, 1], [2, 0],
             [1, 2], [0, 2], [1, 0], [0, 1], [1, 1]]
cases = []
cells = Counter()
for chunk in summary["chunks"]:
    start, count = chunk["start"], chunk["count"]
    name = f"{start:05d}-{start + count - 1:05d}.json"
    path = HERE / "type_i_line_new_all_chunks" / name
    assert sha(path) == chunk["receipt_sha256"]
    source = json.loads(path.read_text())
    assert source["start_index"] == start
    assert source["input_classes"] == source["reduced_classes"] == count
    assert source["unreduced_classes"] == 0
    assert len(source["rows"]) == len(source["representatives"]) == count
    for offset, (row, representative) in enumerate(zip(
            source["rows"], source["representatives"])):
        index = start + offset
        assert len(cases) == index
        assert row["index"] == representative["index"] == index
        assert list(map(int, row["quotient"])) == points[index]
        assert row["status"] == "verified_prime_field_cubic_pair"
        assert representative["status"] == "verified_fixed_origin_degree_three"
        y = representative["representative_y_polynomial"]
        assert len(y) == 3 and int(y[2][0]) != 0
        cells[row["cell"]] += 1
        A = row["positive_input_residual"]
        B = row["negative_input_residual"]
        cases.append([
            coeff(A["residual_x_polynomial"]), coeff(A["residual_y_polynomial"]),
            coeff(B["residual_x_polynomial"]), coeff(B["residual_y_polynomial"]),
            list(map(int, representative["denominator_cubic_coefficients"])),
            coeff(representative["intermediate_x_polynomial"]),
            coeff(representative["intermediate_y_polynomial"]),
            list(map(int, representative["numerator_cubic_coefficients"])),
            coeff(representative["representative_x_polynomial"]),
            coeff(representative["representative_y_polynomial"]),
        ])
assert len(cases) == 9784 and sum(cells.values()) == 9784
gp_path = HERE / "type_i_line_new_all_fixed_origin_input.gp"
with gp_path.open("w") as out:
    out.write(f"p={plane['p']};\n")
    out.write(f"terms={terms};\n")
    out.write(f"exponents={exponents};\n")
    out.write(f"anchors={anchors};\n")
    out.write(f"cases={cases};\n")
receipt = {
    "status": "prepared_all_new_line_fixed_origin_native_replay",
    "cases": len(cases), "nonspecial_representatives": len(cases),
    "cells": dict(sorted(cells.items())),
    "source_sha256": sha(Path(__file__)),
    "summary_sha256": sha(summary_path),
    "points_sha256": sha(point_path),
    "plane_sha256": sha(plane_path),
    "gp_input_sha256": sha(gp_path),
}
(HERE / "type_i_line_new_all_fixed_origin_input.json").write_text(
    json.dumps(receipt, indent=2, sort_keys=True) + "\n")
print("NEW_ALL_FIXED_ORIGIN_INPUT_PASS", len(cases), dict(sorted(cells.items())))
