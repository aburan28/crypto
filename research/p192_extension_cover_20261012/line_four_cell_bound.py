"""Arithmetic certificate for the four-cell P-192 character-sum bound."""

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
density_path = HERE / "type_i_character_density.json"
support_path = HERE / "type_i_line_new64_character_collector.json"
density = json.loads(density_path.read_text())
support = json.loads(support_path.read_text())
r = int(density["r"])
p = int(density["p"])
trace = int(density["trace"])
U = int(density["ceil_sqrt_p"])
assert int(support["subgroup_order"]) == r
assert int(support["field_prime"]) == p
assert U * U > p and (U - 1) * (U - 1) < p
assert density["cover_genera_by_riemann_hurwitz"] == [7, 31]
B1 = 14 * U + abs(trace)
B12 = 62 * U + abs(trace)
# The same first-character sum occurs twice: [2] permutes E(Fp), as r is odd.
# For any signs eps, eta, the absolute error is at most 2*B1+B12.
error = 2 * B1 + B12
lower = (r - error + 3) // 4 - 1
assert lower == int(density["accepted_nonzero_scalar_lower_bound"])
assert (r - error) % 4 == 0
cells = support["support_by_Q_H_cell"]
assert set(cells) == {"++", "+-", "-+", "--"}
assert sum(cells.values()) == support["distinct_three_sum_targets"] == 1456904
max_cell = max(cells.values())
assert max_cell == 364824
assert max_cell * (1 << 171) < lower
assert 64 * max_cell * (1 << 165) < lower
out = {
    "status": "p192_four_cell_character_bound_arithmetic_certificate",
    "field_prime": str(p), "subgroup_order": str(r),
    "ceil_sqrt_prime": str(U), "trace": str(trace),
    "first_character_sum_absolute_bound": str(B1),
    "product_character_sum_absolute_bound": str(B12),
    "all_cell_error_bound": str(error),
    "accepted_nonzero_scalar_lower_bound_each_cell": str(lower),
    "support_by_cell": cells,
    "max_cell_support": max_cell,
    "ideal_uniform_conditioned_single_query_hit_lt_two_pow_minus_171": True,
    "ideal_uniform_conditioned_64_query_union_bound_lt_two_pow_minus_165": True,
    "proof_hypotheses": "The same Hasse-Weil bounds for the genus-7 and genus-31 character covers and the identity exception at most one apply to every sign pair; doubling permutes the odd-prime-order elliptic group.",
    "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "input_sha256": {
        "density": hashlib.sha256(density_path.read_bytes()).hexdigest(),
        "support": hashlib.sha256(support_path.read_bytes()).hexdigest(),
    },
}
(HERE / "type_i_line_four_cell_bound.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n")
print("FOUR_CELL_BOUND_PASS", json.dumps({
    "lower_each_cell": str(lower), "max_cell_support": max_cell,
    "single_query_exponent": 171,
}, sort_keys=True))
