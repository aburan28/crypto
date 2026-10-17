"""Independent prime-field group replay of the line-base collector receipt."""

import hashlib
import json
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
receipt = json.loads((HERE / "type_i_line_guided_collector.json").read_text())
prior = json.loads((HERE / "type_i_geom_character_collector.json").read_text())
projection = json.loads((HERE / "type_i_line_projection.json").read_text())
assert receipt["status"] == "verified_line_base_known_scalar_query_law_comparison"
assert receipt["source_sha256"] == sha(HERE / "line_guided_collector.py")
assert receipt["runtime_sha256"] == sha(HERE / "type_i_line_collector_runtime_info.json")
assert receipt["line_projection_sha256"] == sha(HERE / "type_i_line_projection.json")
assert receipt["guided_query_sha256"] == sha(HERE / "type_i_geom_character_collector.json")
assert receipt["base_points_sha256"] == sha(HERE / "type_i_line_base_points.tsv")
character = json.loads((HERE / "type_i_line_guided_character_verification.json").read_text())
assert character["status"] == "LINE_GUIDED_CHARACTER_INDEPENDENT_PASS"
assert character["producer_receipt_sha256"] == sha(HERE / "type_i_line_guided_collector.json")
assert character["source_sha256"] == sha(HERE / "verify_line_guided_character.py")
assert character["runtime_sha256"] == sha(HERE / "type_i_line_character_verify_runtime_info.json")
assert (character["support_targets"], character["support_admitted"],
        character["uniform_admitted"], character["guided_admitted"]) == (
            receipt["first16_distinct_three_sum_targets"],
            receipt["first16_character_admitted_three_sum_targets"],
            receipt["uniform_character_admitted_queries"], 64)

p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF
r = 0xFFFFFFFFFFFFFFFFFFFFFFFF99DEF836146BC9B1B4D22831
a = p - 3
b = 0x64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1
G = (0x188DA80EB03090F67CBF20EB43A18800F4FF0AFD82FF1012,
     0x07192B95FFC8DA78631011ED6B24CDD573F977A11E794811)
assert int(receipt["field_prime"]) == p and int(receipt["subgroup_order"]) == r


def on_curve(P):
    return P is None or (P[1] * P[1] - P[0] ** 3 - a * P[0] - b) % p == 0


def add(P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    if P[0] == Q[0] and (P[1] + Q[1]) % p == 0:
        return None
    if P == Q:
        m = ((3 * P[0] * P[0] + a) * pow(2 * P[1], -1, p)) % p
    else:
        m = ((Q[1] - P[1]) * pow(Q[0] - P[0], -1, p)) % p
    x = (m * m - P[0] - Q[0]) % p
    return (x, (m * (P[0] - x) - P[1]) % p)


def mul(n, P):
    Q = None
    while n:
        if n & 1:
            Q = add(Q, P)
        P = add(P, P)
        n >>= 1
    return Q


def key(P):
    return "identity" if P is None else f"{P[0]}:{P[1]}"


assert on_curve(G) and mul(r, G) is None
base = [tuple(map(int, pair)) for pair in receipt["base_points_fp"]]
assert len(base) == receipt["base_classes"] == 141
assert all(on_curve(P) and P is not None and mul(r, P) is None for P in base)
assert len(set(base)) == 141
assert hashlib.sha256("\n".join(map(key, base)).encode()).hexdigest() == projection["base_key_sha256"]

pairs = {}
for i in range(141):
    for j in range(i, 141):
        pairs.setdefault(key(add(base[i], base[j])), []).append((i, j))
assert sum(map(len, pairs.values())) == receipt["pair_index_entries"] == 10011

support = {key(add(add(base[i], base[j]), base[k]))
           for i in range(16) for j in range(i, 16) for k in range(j, 16)}
assert len(support) == receipt["first16_distinct_three_sum_targets"]
assert 0 <= receipt["first16_character_admitted_three_sum_targets"] <= len(support)

rng = random.Random(receipt["seed"])
uniform = [rng.randrange(1, r) for _ in range(64)]
guided = list(map(int, prior["accepted_scalars_decimal"]))
assert uniform == list(map(int, receipt["uniform_scalars"]))
assert guided == list(map(int, receipt["guided_scalars"]))
assert receipt["uniform_queries"] == receipt["guided_queries"] == 64
assert receipt["guided_trial_count"] == prior["attempted_scalars"]


def replay(scalars):
    result = []
    for query_index, scalar in enumerate(scalars):
        T = mul((-2 * scalar) % r, G)
        for k, B in enumerate(base):
            neg_B = None if B is None else (B[0], (-B[1]) % p)
            for i, j in pairs.get(key(add(T, neg_B)), []):
                if j > k:
                    continue
                assert add(add(base[i], base[j]), base[k]) == T
                result.append({"query_index": query_index, "scalar": str(scalar),
                               "indices": [i, j, k]})
    return result


assert replay(uniform) == receipt["uniform_hits"]
assert replay(guided) == receipt["guided_hits"]


def rank(hits):
    pivots = {}
    for hit in hits:
        row = [0] * 141
        for index in hit["indices"]:
            row[index] += 1
        for column, pivot in sorted(pivots.items()):
            if row[column]:
                factor = row[column]
                row = [(x - factor * y) % r for x, y in zip(row, pivot)]
        column = next((i for i, value in enumerate(row) if value), None)
        if column is not None:
            inv = pow(row[column], -1, r)
            pivots[column] = [value * inv % r for value in row]
    return len(pivots)


assert rank(receipt["uniform_hits"]) == receipt["uniform_incidence_rank"]
assert rank(receipt["guided_hits"]) == receipt["guided_incidence_rank"]
print("LINE_GUIDED_COLLECTOR_INDEPENDENT_GROUP_PASS", json.dumps({
    "base_classes": len(base), "pair_entries": 10011,
    "uniform_hits": len(receipt["uniform_hits"]),
    "guided_hits": len(receipt["guided_hits"]),
}, sort_keys=True))
