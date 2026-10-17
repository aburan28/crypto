"""Freeze decimal PARI input for independent 141-class cubic-cut replay."""

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
source_path = HERE / "type_i_line_quotient_cubic_probe.json"
reduced_path = HERE / "type_i_line_fixed_origin_reduction.json"
plane_path = HERE / "type_i_plane_quartic.json"
source = json.loads(source_path.read_text())
reduced = json.loads(reduced_path.read_text())
plane = json.loads(plane_path.read_text())
assert reduced["reduced_classes"] == 141 and reduced["unreduced_classes"] == 0
assert len(source["rows"]) == len(reduced["representatives"]) == 141


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
for input_row, output_row in zip(source["rows"], reduced["representatives"]):
    assert input_row["index"] == output_row["index"]
    A = input_row["positive_input_residual"]
    B = input_row["negative_input_residual"]
    cases.append([
        coeff(A["residual_x_polynomial"]), coeff(A["residual_y_polynomial"]),
        coeff(B["residual_x_polynomial"]), coeff(B["residual_y_polynomial"]),
        list(map(int, output_row["denominator_cubic_coefficients"])),
        coeff(output_row["intermediate_x_polynomial"]),
        coeff(output_row["intermediate_y_polynomial"]),
        list(map(int, output_row["numerator_cubic_coefficients"])),
        coeff(output_row["representative_x_polynomial"]),
        coeff(output_row["representative_y_polynomial"]),
    ])
assert len(cases) == 141
gp_path = HERE / "type_i_line_fixed_origin_native_input.gp"
with gp_path.open("w") as out:
    out.write(f"p={plane['p']};\n")
    out.write(f"terms={terms};\n")
    out.write(f"exponents={exponents};\n")
    out.write(f"anchors={anchors};\n")
    out.write(f"cases={cases};\n")
receipt = {
    "status": "prepared_line_fixed_origin_native_replay",
    "cases": 141,
    "source_sha256": sha(__file__),
    "cubic_receipt_sha256": sha(source_path),
    "reduction_receipt_sha256": sha(reduced_path),
    "plane_receipt_sha256": sha(plane_path),
    "gp_input_sha256": sha(gp_path),
}
(HERE / "type_i_line_fixed_origin_native_input.json").write_text(
    json.dumps(receipt, indent=2, sort_keys=True) + "\n")
print("LINE_FIXED_ORIGIN_NATIVE_INPUT_PASS", len(cases))
