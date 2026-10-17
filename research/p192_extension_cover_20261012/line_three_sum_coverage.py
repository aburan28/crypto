"""Exact 141-class P-192 three-sum support and uniform-scalar coverage."""

import hashlib
import json
import math
import time
from pathlib import Path

from sage.all import EllipticCurve, GF, Integer, matrix

HERE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
source_path = HERE/"type_i_line_guided_collector.json"
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
assert all(P != O and r*P == O for P in base)

def key(P):
    return None if P == O else (int(P[0]), int(P[1]))

def support_digest(keys):
    h = hashlib.sha256()
    for P in sorted(keys, key=lambda v: (-1,-1) if v is None else v):
        h.update(b"\xff"*48 if P is None else
                 P[0].to_bytes(24,"big") + P[1].to_bytes(24,"big"))
    return h.hexdigest()

start = time.monotonic_ns()
first = {}
duplicate_count = 0
witnesses = []
relation_rows = []
for i in range(141):
    for j in range(i,141):
        S = base[i]+base[j]
        for k in range(j,141):
            target = key(S+base[k])
            if target in first:
                duplicate_count += 1
                a,b,c = first[target]
                if len(witnesses) < 2000:
                    witnesses.append({"first":[a,b,c],"new":[i,j,k]})
                    row = [0]*141
                    for t in (i,j,k): row[t] += 1
                    for t in (a,b,c): row[t] -= 1
                    relation_rows.append(row)
            else:
                first[target] = (i,j,k)
    if (i+1)%16 == 0:
        print(f"LINE_THREE_SUM_PROGRESS {i+1}/141 support={len(first)}",flush=True)
assert len(first)+duplicate_count == math.comb(143,3)
first16 = set()
for i in range(16):
    for j in range(i,16):
        S = base[i]+base[j]
        for k in range(j,16):
            first16.add(key(S+base[k]))
assert len(first16) == source["first16_distinct_three_sum_targets"] == 816
nonzero_count = len(first)-(None in first)
rank = int(matrix(GF(r),relation_rows).rank()) if relation_rows else 0
out = {
    "status":"verified_p192_line_base_exact_three_sum_support",
    "candidate_id":None,
    "field_prime":str(p),"subgroup_order":str(r),
    "base_classes":141,
    "unordered_triples_with_repetition":math.comb(143,3),
    "distinct_three_sum_targets":len(first),
    "nonzero_three_sum_targets":nonzero_count,
    "identity_in_support":None in first,
    "duplicate_decomposition_count":duplicate_count,
    "collision_witness_count":len(witnesses),
    "collision_witness_rank":rank,
    "collision_witness_rank_scope":"first 2000 duplicate rows only" if duplicate_count>2000 else "all duplicate rows",
    "collision_witnesses":witnesses,
    "first16_distinct_three_sum_targets":len(first16),
    "support_sha256_sorted_xy_be24":support_digest(first),
    "ideal_uniform_scalar_hit_probability_numerator":nonzero_count,
    "ideal_uniform_scalar_hit_probability_denominator":int(r-1),
    "ideal_uniform_64_query_union_bound_numerator":64*nonzero_count,
    "ideal_uniform_64_query_union_bound_denominator":int(r-1),
    "elapsed_ns_unisolated":time.monotonic_ns()-start,
    "source_sha256":sha(__file__),
    "input_sha256":sha(source_path),
    "runtime_sha256":sha(HERE/"type_i_line_three_sum_runtime_info.json"),
    "boundary":"Exact support for the frozen 141-class first quotient base. The rational hit probability applies to an ideal uniform nonzero scalar, not a geometry-conditioned law. Local time is diagnostic.",
}
(HERE/"type_i_line_three_sum_coverage.json").write_text(
    json.dumps(out,indent=2,sort_keys=True)+"\n")
print("LINE_THREE_SUM_COVERAGE_PASS",json.dumps({
    k:out[k] for k in ("distinct_three_sum_targets","nonzero_three_sum_targets",
                       "duplicate_decomposition_count","collision_witness_rank",
                       "identity_in_support")},sort_keys=True),flush=True)
