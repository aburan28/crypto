"""Independent integer replay of the all-pivot fixed-base support criterion."""

import hashlib
import json
from math import factorial, prod
from pathlib import Path

HERE = Path(__file__).resolve().parent
cap_path = HERE / "type_i_all_pivot_capacity.json"
screen_path = HERE / "type_i_pivot_geometry_screen_141.json"
verification_path = HERE / "type_i_pivot_geometry_screen_141_verification.json"
density_path = HERE / "type_i_character_density.json"
cap = json.loads(cap_path.read_text())
screen = json.loads(screen_path.read_text())
geometry_check = json.loads(verification_path.read_text())
density = json.loads(density_path.read_text())
assert cap["status"] == "all_141_pivot_optimistic_base_capacity_bound"
assert geometry_check["status"] == "independent_pivot_geometry_polynomial_controls_pass"
assert geometry_check["screen_sha256"] == hashlib.sha256(screen_path.read_bytes()).hexdigest()
for key, path in (("screen_sha256", screen_path),
                  ("verification_sha256", verification_path),
                  ("density_sha256", density_path)):
    assert cap[key] == hashlib.sha256(path.read_bytes()).hexdigest()
assert cap["source_sha256"] == hashlib.sha256(
    (HERE / "line_all_pivot_capacity.py").read_bytes()).hexdigest()

selected = screen["selected_classes"]
residual = screen["distinct_residual_points"]
base_upper = selected+residual
accepted_lower = int(density["accepted_nonzero_scalar_lower_bound"])
assert (selected, residual, base_upper, accepted_lower) == (
    cap["selected_classes"], cap["distinct_residual_geometry"],
    cap["base_size_upper_bound"], int(cap["accepted_scalar_set_lower_bound"]))


def unordered_capacity(base, arity):
    return prod(base+i for i in range(arity)) // factorial(arity)


three = unordered_capacity(base_upper, 3)
exponent = cap["ideal_uniform_hit_probability_less_than_two_to_minus"]
m = cap["minimum_arity_whose_capacity_reaches_accepted_set"]
assert three == cap["unordered_three_sum_capacity_upper_bound"]
assert three * 2**exponent < accepted_lower <= three * 2**(exponent+1)
assert m > 1
assert unordered_capacity(base_upper, m-1) == cap["capacity_at_previous_arity"]
assert unordered_capacity(base_upper, m) == cap["capacity_at_threshold_arity"]
assert cap["capacity_at_previous_arity"] < accepted_lower <= cap["capacity_at_threshold_arity"]

out = {
    "status": "independent_all_pivot_capacity_replay_pass",
    "base_size_upper_bound": base_upper,
    "three_sum_capacity_upper_bound": three,
    "hit_exponent": exponent,
    "full_capacity_arity_threshold": m,
    "capacity_receipt_sha256": hashlib.sha256(cap_path.read_bytes()).hexdigest(),
    "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
}
path = HERE / "type_i_all_pivot_capacity_verification.json"
path.write_text(json.dumps(out, indent=2, sort_keys=True)+"\n")
print("ALL_PIVOT_CAPACITY_REPLAY_PASS", base_upper, exponent, m, flush=True)
