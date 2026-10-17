"""Independent sorted-list replay of the full two-fan capacity screen."""

import hashlib
import json
from itertools import groupby
from pathlib import Path

HERE = Path(__file__).resolve().parent
receipt = json.loads((HERE / "type_i_two_fan_capacity.json").read_text())
paths = [HERE / name for name in receipt["input_sha256"]]
for path in paths:
    assert hashlib.sha256(path.read_bytes()).hexdigest() == (
        receipt["input_sha256"][path.name])
assert receipt["source_sha256"] == hashlib.sha256(
    (HERE / "line_two_fan_capacity.py").read_bytes()).hexdigest()
old = json.loads((HERE / "type_i_line_guided_collector.json").read_text())
new = json.loads((HERE / "type_i_line_full_new_quotient_points.json").read_text())
fans = [json.loads((HERE / f"type_i_residual_pivot{i:03d}_{'three' if i == 0 else 'cross'}.json").read_text())
        for i in (0, 1)]
density = json.loads((HERE / "type_i_character_density.json").read_text())


def distinct(values):
    return [key for key, _ in groupby(sorted(tuple(map(int, point))
                                        for point in values))]


def overlap(a, b):
    i = j = total = 0
    while i < len(a) and j < len(b):
        if a[i] == b[j]:
            total += 1
            i += 1
            j += 1
        elif a[i] < b[j]:
            i += 1
        else:
            j += 1
    return total


selected = distinct(old["base_points_fp"] + new)
fan_points = [distinct(point for row in fan["rows"]
                       for point in row["projected_residuals"])
              for fan in fans]
assert [len(selected), *map(len, fan_points)] == [9925, 10073, 9941]
assert [overlap(selected, fan) for fan in fan_points] == [2, 2]
assert overlap(*fan_points) == 0
union = distinct(selected + fan_points[0] + fan_points[1])
assert len(union) == receipt["union_distinct_usable_upper_base"] == 29935
assert (len(union)*(len(union)+1)*(len(union)+2)//6 ==
        receipt["unordered_three_sum_capacity"] == 4471261391320)
lower = int(density["accepted_nonzero_scalar_lower_bound"])
assert lower == int(receipt["accepted_scalar_set_lower_bound"])
assert receipt["unordered_three_sum_capacity"]*(1 << 147) < lower
output = {
    "status": "verified_independent_two_fan_capacity",
    "selected_classes": len(selected),
    "fan_classes": list(map(len, fan_points)),
    "overlaps": [2, 2], "union_classes": len(union),
    "three_sum_capacity": receipt["unordered_three_sum_capacity"],
    "producer_sha256": hashlib.sha256(
        (HERE / "type_i_two_fan_capacity.json").read_bytes()).hexdigest(),
    "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
}
(HERE / "type_i_two_fan_capacity_verification.json").write_text(
    json.dumps(output, indent=2, sort_keys=True) + "\n")
print("TWO_FAN_CAPACITY_VERIFY_PASS", len(union), output["three_sum_capacity"])
