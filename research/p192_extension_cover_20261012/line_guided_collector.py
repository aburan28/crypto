"""Compare frozen known-scalar query laws on the line-derived Type I base."""

import hashlib
import json
import random
import time
from pathlib import Path

from sage.all import EllipticCurve, GF, Integer, matrix, prod, vector

import type_i_divisor_transfer as dt
from plane_point_inverse import PlanePointInverse, elliptic_quotient_image

HERE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
prior = json.loads((HERE / "type_i_geom_character_collector.json").read_text())
projection = json.loads((HERE / "type_i_line_projection.json").read_text())
base_receipt = json.loads((HERE / "type_i_quotient_base_probe.json").read_text())
assert prior["accepted_queries"] == 64 and prior["seed"] == 0x19220261023
assert projection["base_rational_points"] == 141

p, r, L = dt.p, dt.r, dt.L
F = dt.gm.F
E = EllipticCurve(L, [L(dt.gm.a), L(dt.gm.b)])
G = E(L(dt.gm.gx), L(dt.gm.gy))
O = E(0)
h = Integer(base_receipt["candidate_cofactor"])
assert r * G == O and r.is_prime()


def cover_from_x(x):
    X = L(x)
    U = (X - dt.zeta) * (X - dt.beta)
    y = U ** ((dt.q + 1) // 2)
    N = prod(X - root ** (p ** j)
             for root in (dt.zeta, dt.beta) for j in range(3))
    if y * y ** p * y ** (p * p) == -N:
        y = -y
    assert y * y ** p * y ** (p * p) == N
    return (X, *dt.Vandermonde.solve_right(
        vector(L, [y, y ** p, y ** (p * p)])))


def character_test(P):
    if P == O or P[0] == dt.c0:
        return False
    X = dt.inf + dt.d0 / (P[0] - dt.c0)
    rhs = prod(X - z ** p for z in dt.quartic_roots)
    assert rhs ** (p ** 3) == rhs
    return rhs != 0 and rhs ** ((p ** 3 - 1) // 2) == 1


def key(P):
    if P == O:
        return "identity"
    assert P[0] ** p == P[0] and P[1] ** p == P[1]
    return f"{int(P[0].polynomial()[0])}:{int(P[1].polynomial()[0])}"


start = time.monotonic_ns()
rng = random.Random(base_receipt["pilot_seed"])
anchor_x = [Integer(rng.randrange(int(p)))
            for _ in range(base_receipt["anchor_pilot_index"] + 1)][-1]
anchor = cover_from_x(anchor_x)
anchor_images = [elliptic_quotient_image(dt, E, anchor, j) for j in range(3)]
inverse = PlanePointInverse(dt)
rows = [tuple(map(int, line.split())) for line in
        (HERE / "type_i_line_base_points.tsv").read_text().splitlines()]
assert len(rows) == 141
base = []
for _, u, v in rows:
    point = inverse.invert_affine(F(u), F(v))
    images = [elliptic_quotient_image(dt, E, point, j) for j in range(3)]
    projected = [h * (images[j] - anchor_images[j]) for j in range(3)]
    assert projected[0] == projected[1] == projected[2]
    assert projected[0] != O and r * projected[0] == O
    base.append(projected[0])
base_keys = [key(P) for P in base]
assert len(set(base_keys)) == 141
assert hashlib.sha256("\n".join(base_keys).encode()).hexdigest() == projection["base_key_sha256"]

pairs = {}
for i in range(len(base)):
    for j in range(i, len(base)):
        pairs.setdefault(key(base[i] + base[j]), []).append((i, j))
assert sum(map(len, pairs.values())) == 10011

# The first 16 line-derived classes give an exact, tractable support census.
# In the r-torsion group T=-2P uniquely recovers P by multiplying by
# -(2^-1 mod r), so the geometric admission predicate is well-defined on T.
inv2 = Integer(pow(2, -1, int(r)))
support = {}
for i in range(16):
    for j in range(i, 16):
        for k in range(j, 16):
            T = base[i] + base[j] + base[k]
            support.setdefault(key(T), T)
support_admitted = 0
for T in support.values():
    P = -inv2 * T
    if character_test(P) and character_test(2 * P):
        support_admitted += 1

seeded = random.Random(prior["seed"])
uniform = [Integer(seeded.randrange(1, int(r))) for _ in range(64)]
conditioned = [Integer(a) for a in prior["accepted_scalars_decimal"]]
assert len(set(uniform)) == len(set(conditioned)) == 64
assert all(character_test(a * G) and character_test(2 * a * G)
           for a in conditioned)
uniform_admitted = sum(character_test(a * G) and character_test(2 * a * G)
                       for a in uniform)


def collect(scalars):
    hits = []
    for query_index, a in enumerate(scalars):
        T = -2 * a * G
        assert r * T == O
        for k in range(len(base)):
            for i, j in pairs.get(key(T - base[k]), []):
                if j > k:
                    continue
                assert base[i] + base[j] + base[k] == T
                hits.append({"query_index": query_index, "scalar": str(a),
                             "indices": [i, j, k]})
    rows = []
    for hit in hits:
        row = [0] * len(base)
        for index in hit["indices"]:
            row[index] += 1
        rows.append(row)
    rank = int(matrix(GF(r), rows).rank()) if rows else 0
    return hits, rank


uniform_hits, uniform_rank = collect(uniform)
guided_hits, guided_rank = collect(conditioned)
elapsed = time.monotonic_ns() - start
out = {
    "status": "verified_line_base_known_scalar_query_law_comparison",
    "candidate_id": None,
    "field_prime": str(p), "subgroup_order": str(r),
    "base_classes": len(base), "pair_index_entries": 10011,
    "first16_distinct_three_sum_targets": len(support),
    "first16_character_admitted_three_sum_targets": support_admitted,
    "seed": prior["seed"], "uniform_queries": 64,
    "uniform_character_admitted_queries": uniform_admitted,
    "guided_queries": 64, "guided_trial_count": prior["attempted_scalars"],
    "uniform_scalars": list(map(str, uniform)),
    "guided_scalars": list(map(str, conditioned)),
    "uniform_hits": uniform_hits, "guided_hits": guided_hits,
    "uniform_incidence_rank": uniform_rank, "guided_incidence_rank": guided_rank,
    "base_points_fp": [[int(P[0].polynomial()[0]), int(P[1].polynomial()[0])]
                       for P in base],
    "elapsed_ns_unisolated": elapsed,
    "source_sha256": sha(__file__),
    "line_projection_sha256": sha(HERE / "type_i_line_projection.json"),
    "guided_query_sha256": sha(HERE / "type_i_geom_character_collector.json"),
    "base_points_sha256": sha(HERE / "type_i_line_base_points.tsv"),
    "runtime_sha256": sha(HERE / "type_i_line_collector_runtime_info.json"),
    "boundary": "Exact first-quotient three-sum lookup on 141 line-derived classes under two 64-query known-scalar laws; first 16 classes have exact character-conditioned support. Local timing is a correctness diagnostic.",
}
(HERE / "type_i_line_guided_collector.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n")
print("LINE_GUIDED_COLLECTOR_PASS", json.dumps({
    k: out[k] for k in ("base_classes", "first16_distinct_three_sum_targets",
                      "first16_character_admitted_three_sum_targets",
                      "uniform_character_admitted_queries", "uniform_incidence_rank",
                      "guided_incidence_rank")}, sort_keys=True), flush=True)
