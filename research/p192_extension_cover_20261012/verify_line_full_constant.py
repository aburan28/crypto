"""Independent prime-field replay of all 4,892 split-line quotient sums."""

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
receipt_path = HERE / "type_i_line_full_constant.json"
receipt = json.loads(receipt_path.read_text())
prior_path = HERE / "type_i_line_guided_collector.json"
projection_path = HERE / "type_i_line_projection.json"
leaf_path = HERE / "type_i_line_leaf_relation_probe.json"
new_path = HERE / "type_i_line_full_new_quotient_points.json"
prior = json.loads(prior_path.read_text())
projection = json.loads(projection_path.read_text())
p = int(prior["field_prime"])
r = int(prior["subgroup_order"])
a = p - 3
b = int("64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1", 16)
assert receipt["field_prime"] == str(p) and receipt["subgroup_order"] == str(r)
assert receipt["source_sha256"] == sha(HERE / "line_full_constant_replay.py")
assert receipt["runtime_sha256"] == sha(HERE / "type_i_line_full_constant_runtime_info.json")
assert receipt["prior_sha256"] == sha(prior_path)
assert receipt["projection_sha256"] == sha(projection_path)
assert receipt["new_coordinate_sha256"] == sha(new_path)
for name, digest in receipt["input_sha256"].items():
    assert digest == sha(HERE / name)


def add(P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    x, y = P
    u, v = Q
    if x == u and (y + v) % p == 0:
        return None
    slope = ((3 * x * x + a) * pow(2 * y, -1, p) if P == Q
             else (v - y) * pow((u - x) % p, -1, p)) % p
    z = (slope * slope - x - u) % p
    return z, (slope * (x - z) - y) % p


def neg(P):
    return None if P is None else (P[0], (-P[1]) % p)


def mul(n, P):
    out = None
    while n:
        if n & 1:
            out = add(out, P)
        P = add(P, P)
        n >>= 1
    return out


def point(xy):
    x, y = map(int, xy)
    assert 0 <= x < p and 0 <= y < p
    assert (y * y - x * x * x - a * x - b) % p == 0
    return x, y


old = [point(xy) for xy in prior["base_points_fp"]]
new = [point(xy) for xy in json.loads(new_path.read_text())]
assert len(old) == 141 and len(new) == 9784
assert len(set(old)) == 141 and len(set(new)) == 9784
assert len({min(Q, neg(Q)) for Q in new}) == 9784
assert {min(Q, neg(Q)) for Q in new}.isdisjoint(
    {min(Q, neg(Q)) for Q in old})
key_digest = hashlib.sha256("\n".join(f"{x}:{y}" for x, y in new).encode()).hexdigest()
assert key_digest == receipt["new_key_sha256"] == projection["new_key_sha256"]
pair_rows = [tuple(map(int, line.split())) for line in
             (HERE / "type_i_line_new_point_pairs.tsv").read_text().splitlines()]
assert len(pair_rows) == 9784
C = point(receipt["constant_fp"])
D = point(receipt["center_fp"])
leaf = json.loads(leaf_path.read_text())
assert leaf["constant_key"] == f"{C[0]}:{C[1]}"
assert C is not None and D is not None
assert mul(4, D) == C and mul(r, C) is None and mul(r, D) is None
assert D == mul(pow(4, -1, r), C)
assert receipt["center_fourfold_is_constant"]
for k in range(4892):
    i, j = pair_rows[2 * k]
    assert pair_rows[2 * k + 1] == (i, j) and 1 <= i < j <= 141
    row = (old[i - 1], old[j - 1], new[2 * k], new[2 * k + 1])
    total = None
    centered = None
    for Q in row:
        total = add(total, Q)
        centered = add(centered, add(Q, neg(D)))
    assert total == C and centered is None
assert receipt["line_cuts"] == 4892 and receipt["all_centered_line_sums_zero"]
out = {
    "status": "LINE_FULL_CONSTANT_INDEPENDENT_PASS",
    "line_cuts": 4892,
    "new_points": len(new),
    "constant_fp": list(C),
    "center_fp": list(D),
    "center_in_frozen_base": D in set(old + new),
    "constant_in_frozen_base": C in set(old + new),
    "source_sha256": sha(__file__),
    "receipt_sha256": sha(receipt_path),
    "leaf_sample_sha256": sha(leaf_path),
    "new_coordinate_sha256": sha(new_path),
}
(HERE / "type_i_line_full_constant_verification.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n")
print("LINE_FULL_CONSTANT_INDEPENDENT_PASS", 4892, len(new), flush=True)
