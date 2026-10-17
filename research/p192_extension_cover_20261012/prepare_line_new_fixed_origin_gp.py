"""Serialize the bounded new quotient reduction for independent PARI replay."""

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
source_path = HERE / "type_i_line_new64_fixed_origin.json"
plane_path = HERE / "type_i_plane_quartic.json"
source = json.loads(source_path.read_text())
plane = json.loads(plane_path.read_text())
assert source["start_index"] == 0
assert source["input_classes"] == source["reduced_classes"] == 64
assert source["unreduced_classes"] == 0
assert len(source["rows"]) == len(source["representatives"]) == 64


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
for row, representative in zip(source["rows"], source["representatives"]):
    assert row["index"] == representative["index"]
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
gp_path = HERE / "type_i_line_new64_fixed_origin_input.gp"
with gp_path.open("w") as out:
    out.write(f"p={plane['p']};\n")
    out.write(f"terms={terms};\n")
    out.write(f"exponents={exponents};\n")
    out.write(f"anchors={anchors};\n")
    out.write(f"cases={cases};\n")
result = {
    "status": "prepared_new_line_fixed_origin_native_replay",
    "cases": len(cases),
    "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "reduction_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
    "plane_sha256": hashlib.sha256(plane_path.read_bytes()).hexdigest(),
    "gp_input_sha256": hashlib.sha256(gp_path.read_bytes()).hexdigest(),
}
(HERE / "type_i_line_new64_fixed_origin_input.json").write_text(
    json.dumps(result, indent=2, sort_keys=True) + "\n")
print("NEW_LINE_FIXED_ORIGIN_INPUT_PASS", len(cases))
