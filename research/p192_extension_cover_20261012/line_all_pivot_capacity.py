"""Optimistic fixed-base support bounds from a complete pivot geometry screen."""

import hashlib
import json
from math import comb
from pathlib import Path

HERE = Path(__file__).resolve().parent
screen_path = HERE / "type_i_pivot_geometry_screen_141.json"
verification_path = HERE / "type_i_pivot_geometry_screen_141_verification.json"
density_path = HERE / "type_i_character_density.json"
screen = json.loads(screen_path.read_text())
verification = json.loads(verification_path.read_text())
density = json.loads(density_path.read_text())
assert screen["status"] == "complete_declared_pivot_residual_geometry_screen"
assert screen["pivot_fans"] == verification["pivot_fans"] == 141
assert verification["status"] == "independent_pivot_geometry_polynomial_controls_pass"
assert verification["screen_sha256"] == hashlib.sha256(screen_path.read_bytes()).hexdigest()
assert screen["counts"] == verification["counts"]

# The quotient-point map is a function of a quartic point. Consequently,
# its image has at most as many elements as the geometric input union.
# Adding the two geometric-set cardinalities is an intentionally generous
# upper bound that remains valid without computing their intersection.
base_upper = screen["selected_classes"] + screen["distinct_residual_points"]
accepted_lower = int(density["accepted_nonzero_scalar_lower_bound"])
assert accepted_lower > 0
three_sum_capacity = comb(base_upper+2, 3)
exponent = 0
while three_sum_capacity * (2**(exponent+1)) < accepted_lower:
    exponent += 1
minimum_capacity_arity = next(m for m in range(1, 64)
                              if comb(base_upper+m-1, m) >= accepted_lower)
assert comb(base_upper+minimum_capacity_arity-2,
            minimum_capacity_arity-1) < accepted_lower

out = {
    "status": "all_141_pivot_optimistic_base_capacity_bound",
    "selected_classes": screen["selected_classes"],
    "distinct_residual_geometry": screen["distinct_residual_points"],
    "base_size_upper_bound": base_upper,
    "accepted_scalar_set_lower_bound": str(accepted_lower),
    "unordered_three_sum_capacity_upper_bound": three_sum_capacity,
    "ideal_uniform_hit_probability_less_than_two_to_minus": exponent,
    "minimum_arity_whose_capacity_reaches_accepted_set": minimum_capacity_arity,
    "capacity_at_previous_arity": comb(base_upper+minimum_capacity_arity-2,
                                       minimum_capacity_arity-1),
    "capacity_at_threshold_arity": comb(base_upper+minimum_capacity_arity-1,
                                        minimum_capacity_arity),
    "screen_sha256": hashlib.sha256(screen_path.read_bytes()).hexdigest(),
    "verification_sha256": hashlib.sha256(verification_path.read_bytes()).hexdigest(),
    "density_sha256": hashlib.sha256(density_path.read_bytes()).hexdigest(),
    "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "boundary": "Upper support capacity assumes every selected and residual geometric point yields a distinct usable quotient class and that unordered m-sums are distinct. The arity threshold is only a cardinality screen, not a relation-yield or cost estimate. Ideal-uniform probability is conditional on sampling injectively from the certified accepted scalar set.",
}
path = HERE / "type_i_all_pivot_capacity.json"
path.write_text(json.dumps(out, indent=2, sort_keys=True)+"\n")
print("ALL_PIVOT_CAPACITY_PASS", base_upper, exponent,
      minimum_capacity_arity, flush=True)
