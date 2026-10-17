"""Exact additive-pair collision rank on 141 frozen P-192 quotient classes."""

import hashlib
import json
import time
from collections import Counter
from pathlib import Path

from sage.all import EllipticCurve, GF, Integer, matrix

HERE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
source_path = HERE / "type_i_line_guided_collector.json"
source = json.loads(source_path.read_text())
p = Integer(source["field_prime"])
r = Integer(source["subgroup_order"])
assert p == Integer("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF", 16)
assert source["base_classes"] == 141
F = GF(p)
E = EllipticCurve(F, [F(p-3), F(Integer("64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1", 16))])
O = E(0)
base = [E(F(x), F(y)) for x, y in source["base_points_fp"]]
assert len(base) == 141 and len(set(base)) == 141
assert all(P != O and r * P == O for P in base)
start = time.monotonic_ns()

first = {}
multiplicities = Counter()
rows = []
witnesses = []
for i in range(141):
    for j in range(i, 141):
        T = base[i] + base[j]
        multiplicities[T] += 1
        if T in first:
            a, b = first[T]
            row = [0] * 141
            for k in (i, j):
                row[k] += 1
            for k in (a, b):
                row[k] -= 1
            assert sum(row) == 0
            assert base[i] + base[j] == base[a] + base[b]
            rows.append(row)
            witnesses.append({"new": [i, j], "first": [a, b]})
        else:
            first[T] = (i, j)
assert len(first) + len(rows) == 10011
rank = int(matrix(GF(r), rows).rank()) if rows else 0
hist = Counter(multiplicities.values())
out = {
    "status": "verified_p192_line_base_pair_collision_rank",
    "candidate_id": None,
    "field_prime": str(p), "subgroup_order": str(r),
    "base_classes": 141, "unordered_pairs_with_repetition": 10011,
    "distinct_pair_targets": len(first), "collision_rows": len(rows),
    "relation_rank": rank, "rank_deficiency": 141-rank,
    "multiplicity_histogram": {str(k): v for k, v in sorted(hist.items())},
    "witnesses": witnesses,
    "elapsed_ns_unisolated": time.monotonic_ns()-start,
    "source_sha256": sha(__file__),
    "input_sha256": sha(source_path),
    "runtime_sha256": sha(HERE/"type_i_line_pair_rank_runtime_info.json"),
    "boundary": "Homogeneous additive relations among 141 exact r-torsion quotient classes; group-equation rank does not by itself give individual logarithms.",
}
(HERE/"type_i_line_pair_collision_rank.json").write_text(
    json.dumps(out,indent=2,sort_keys=True)+"\n")
print("LINE_PAIR_COLLISION_RANK_PASS", json.dumps({
    k: out[k] for k in ("collision_rows","relation_rank","rank_deficiency",
                        "multiplicity_histogram")},sort_keys=True))
