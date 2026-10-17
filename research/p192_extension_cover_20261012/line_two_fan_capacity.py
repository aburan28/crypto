"""Exact two-fan point union and three-sum capacity screen."""

import hashlib
import json
from math import comb
from pathlib import Path

HERE = Path(__file__).resolve().parent
names = (
    "type_i_line_guided_collector.json",
    "type_i_line_full_new_quotient_points.json",
    "type_i_residual_pivot000_three.json",
    "type_i_residual_pivot001_cross.json",
    "type_i_character_density.json",
)
paths = [HERE / name for name in names]
old, new, fan0, fan1, density = [json.loads(path.read_text()) for path in paths]
p = int(old["field_prime"])
r = int(old["subgroup_order"])
assert p == int(fan0["field_prime"]) == int(fan1["field_prime"])
assert r == int(fan0["subgroup_order"]) == int(fan1["subgroup_order"])


def points(values):
    return {tuple(map(int, value)) for value in values}


selected = points(old["base_points_fp"] + new)
fans = [points(value for row in fan["rows"]
               for value in row["projected_residuals"])
        for fan in (fan0, fan1)]
assert len(selected) == 9925
assert [len(fan) for fan in fans] == [10073, 9941]
overlaps = [len(selected & fan) for fan in fans]
assert overlaps == [2, 2] and not (fans[0] & fans[1])
union = selected | fans[0] | fans[1]
assert len(union) == 29935
assert all((y*y-x*x*x+3*x-int("64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1",16)) % p == 0
           for x, y in union)
capacity = comb(len(union)+2, 3)
assert capacity == 4471261391320
lower = 1569275433846670190958947353987675545133315652035359680560
assert int(density["accepted_nonzero_scalar_lower_bound"]) == lower
assert capacity*(1 << 147) < lower
output = {
    "status": "verified_two_fan_full_base_capacity",
    "candidate_id": None,
    "field_prime": str(p), "subgroup_order": str(r),
    "selected_classes": len(selected), "fan_classes": [len(fan) for fan in fans],
    "selected_fan_overlaps": overlaps, "cross_fan_overlap": 0,
    "union_distinct_usable_upper_base": len(union),
    "unordered_three_sum_capacity": capacity,
    "accepted_scalar_set_lower_bound": str(lower),
    "ideal_uniform_accepted_hit_upper_exponent": 147,
    "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "input_sha256": {name: hashlib.sha256(path.read_bytes()).hexdigest()
                     for name, path in zip(names, paths)},
    "boundary": "Capacity assumes all exact quotient columns become usable fixed-base classes. The accepted-set lower bound is supplied by the previously checked character-density analysis; exact relation yield and rank require a collector.",
}
path = HERE / "type_i_two_fan_capacity.json"
path.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
print("TWO_FAN_CAPACITY_PASS", len(union), capacity)
