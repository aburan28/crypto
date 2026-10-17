"""Independently check the fixed-base capacity receipt against the density certificate."""

import hashlib
import json
from fractions import Fraction
from math import comb
from pathlib import Path


HERE = Path(__file__).resolve().parent
density_path = HERE / "type_i_character_density.json"
verification_path = HERE / "type_i_character_density_verification.json"
threshold_path = HERE / "type_i_factor_base_threshold.json"
density = json.loads(density_path.read_text())
density_verification = json.loads(verification_path.read_text())
threshold = json.loads(threshold_path.read_text())

assert density_verification["status"] == "verified_type_i_character_density_inputs"
lower = int(density["accepted_nonzero_scalar_lower_bound"])
assert lower == int(density_verification["accepted_nonzero_lower_bound"])
assert lower == int(threshold["accepted_nonzero_scalar_lower_bound"])
assert threshold["status"] == "exact_type_i_fixed_three_sum_capacity_threshold"
assert threshold["density_receipt_sha256"] == hashlib.sha256(density_path.read_bytes()).hexdigest()
assert threshold["source_sha256"] == hashlib.sha256(
    (HERE / "type_i_factor_base_threshold.py").read_bytes()
).hexdigest()

base = int(threshold["minimum_actual_distinct_base_classes_for_capacity"])
previous_capacity = comb(base + 1, 3)
capacity = comb(base + 2, 3)
assert base > 0 and 2 * previous_capacity < lower <= 2 * capacity
assert previous_capacity == int(threshold["previous_base_three_sum_capacity"])
assert capacity == int(threshold["threshold_base_three_sum_capacity"])
assert comb(base + 1, 2) == int(threshold["full_unordered_pair_index_entries_at_threshold"])
assert int(threshold["current_base_classes"]) == 16
assert comb(18, 3) == int(threshold["current_three_sum_support_cap"])
assert Fraction(threshold["current_per_accepted_query_hit_upper_bound"]) == Fraction(816, lower)
assert Fraction(threshold["expected_accepted_queries_lower_bound_for_one_hit"]) == Fraction(lower, 816)
assert Fraction(threshold["expected_accepted_queries_lower_bound_for_sixteen_hits"]) == Fraction(16 * lower, 816)

result = {
    "status": "verified_type_i_fixed_base_capacity_threshold",
    "base_classes": str(base),
    "previous_twice_capacity": str(2 * previous_capacity),
    "threshold_twice_capacity": str(2 * capacity),
    "accepted_lower_bound": str(lower),
    "pair_index_entries": str(comb(base + 1, 2)),
    "threshold_receipt_sha256": hashlib.sha256(threshold_path.read_bytes()).hexdigest(),
    "density_verification_sha256": hashlib.sha256(verification_path.read_bytes()).hexdigest(),
    "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "boundary": "Exact capacity test for fixed-base unordered three-sum support; geometry and density are supplied by separately verified input receipts.",
}
(HERE / "type_i_factor_base_threshold_verification.json").write_text(
    json.dumps(result, indent=2, sort_keys=True) + "\n"
)
print(json.dumps(result, indent=2, sort_keys=True))
