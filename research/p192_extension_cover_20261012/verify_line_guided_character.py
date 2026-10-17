"""Independent Sage replay of the two-fiber support and query admissions."""

import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, Integer, prod

import type_i_divisor_transfer as dt

HERE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
receipt_path = HERE / "type_i_line_guided_collector.json"
receipt = json.loads(receipt_path.read_text())
p, r, L = dt.p, dt.r, dt.L
E = EllipticCurve(L, [L(dt.gm.a), L(dt.gm.b)])
G = E(L(dt.gm.gx), L(dt.gm.gy))
O = E(0)
base = [E(L(Integer(x)), L(Integer(y))) for x, y in receipt["base_points_fp"]]
assert len(base) == 141 and r * G == O


def admitted_target(T):
    # T=-2P on E[r]. This inversion is unique because r is odd.
    P = -(Integer(pow(2, -1, int(r)))) * T
    for Q in (P, 2 * P):
        if Q == O or Q[0] == dt.c0:
            return False
        X = dt.inf + dt.d0 / (Q[0] - dt.c0)
        f = prod(X - root ** p for root in dt.quartic_roots)
        assert f ** (p ** 3) == f
        if f == 0 or f ** ((p ** 3 - 1) // 2) != 1:
            return False
    return True


targets = set()
for i in range(16):
    for j in range(i, 16):
        for k in range(j, 16):
            targets.add(base[i] + base[j] + base[k])
assert len(targets) == receipt["first16_distinct_three_sum_targets"] == 816
accepted = sum(admitted_target(T) for T in targets)
assert accepted == receipt["first16_character_admitted_three_sum_targets"]
uniform = [Integer(s) for s in receipt["uniform_scalars"]]
guided = [Integer(s) for s in receipt["guided_scalars"]]
uniform_accept = sum(admitted_target(-2 * scalar * G) for scalar in uniform)
assert uniform_accept == receipt["uniform_character_admitted_queries"]
assert all(admitted_target(-2 * scalar * G) for scalar in guided)
out = {
    "status": "LINE_GUIDED_CHARACTER_INDEPENDENT_PASS",
    "support_targets": len(targets), "support_admitted": accepted,
    "uniform_admitted": uniform_accept, "guided_admitted": len(guided),
    "producer_receipt_sha256": sha(receipt_path),
    "source_sha256": sha(__file__),
    "runtime_sha256": sha(HERE / "type_i_line_character_verify_runtime_info.json"),
}
(HERE / "type_i_line_guided_character_verification.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n")
print("LINE_GUIDED_CHARACTER_INDEPENDENT_PASS", json.dumps(out, sort_keys=True))
