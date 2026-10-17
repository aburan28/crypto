"""Independent integer-field replay of P-192 pair and three-sum support."""

import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
source_path = HERE/"type_i_line_guided_collector.json"
pair_path = HERE/"type_i_line_pair_collision_rank.json"
triple_path = HERE/"type_i_line_three_sum_coverage.json"
source = json.loads(source_path.read_text())
pair = json.loads(pair_path.read_text())
triple = json.loads(triple_path.read_text())
p = int(source["field_prime"])
r = int(source["subgroup_order"])
a = p-3
b = int("64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1",16)
assert p == int("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF",16)
base = [tuple(map(int,xy)) for xy in source["base_points_fp"]]
assert len(base) == len(set(base)) == 141

def plus(P,Q):
    if P is None: return Q
    if Q is None: return P
    x,y = P
    u,v = Q
    if x == u and (y+v)%p == 0: return None
    if P == Q:
        slope = ((3*x*x+a)*pow(2*y,-1,p))%p
    else:
        slope = ((v-y)*pow((u-x)%p,-1,p))%p
    z = (slope*slope-x-u)%p
    return (z,(slope*(x-z)-y)%p)

def mul(n,P):
    R = None
    while n:
        if n&1: R=plus(R,P)
        P=plus(P,P)
        n >>= 1
    return R

assert all((y*y-x*x*x-a*x-b)%p == 0 and mul(r,(x,y)) is None
           for x,y in base)
G = (int("188DA80EB03090F67CBF20EB43A18800F4FF0AFD82FF1012",16),
     int("07192B95FFC8DA78631011ED6B24CDD573F977A11E794811",16))
assert (G[1]*G[1]-G[0]*G[0]*G[0]-a*G[0]-b)%p == 0
assert mul(r,G) is None and G is not None and plus(G,G) is not None
pair_sums = set()
for i in range(141):
    for j in range(i,141):
        pair_sums.add(plus(base[i],base[j]))
assert len(pair_sums) == pair["distinct_pair_targets"] == 10011
first16 = set()
for i in range(16):
    for j in range(i,16):
        S = plus(base[i],base[j])
        for k in range(j,16):
            first16.add(plus(S,base[k]))
assert len(first16) == source["first16_distinct_three_sum_targets"] == 816
support = set()
for i in range(141):
    for j in range(i,141):
        S = plus(base[i],base[j])
        for k in range(j,141):
            support.add(plus(S,base[k]))
    if (i+1)%32 == 0:
        print(f"LINE_SUM_INDEPENDENT_PROGRESS {i+1}/141 support={len(support)}",flush=True)
assert len(support) == triple["distinct_three_sum_targets"] == math.comb(143,3)
assert None not in support and triple["identity_in_support"] is False
digest = hashlib.sha256()
for x,y in sorted(support):
    digest.update(x.to_bytes(24,"big")+y.to_bytes(24,"big"))
assert digest.hexdigest() == triple["support_sha256_sorted_xy_be24"]
assert pair["source_sha256"] == sha(HERE/"line_pair_collision_rank.py")
assert pair["input_sha256"] == sha(source_path)
assert pair["runtime_sha256"] == sha(HERE/"type_i_line_pair_rank_runtime_info.json")
assert triple["source_sha256"] == sha(HERE/"line_three_sum_coverage.py")
assert triple["input_sha256"] == sha(source_path)
assert triple["runtime_sha256"] == sha(HERE/"type_i_line_three_sum_runtime_info.json")
assert pair["collision_rows"] == triple["duplicate_decomposition_count"] == 0
assert r%2 == 1
assert triple["ideal_uniform_scalar_hit_probability_numerator"] == len(support)
assert triple["ideal_uniform_scalar_hit_probability_denominator"] == r-1
assert triple["ideal_uniform_64_query_union_bound_numerator"] == 64*len(support)
assert triple["ideal_uniform_64_query_union_bound_denominator"] == r-1
capacity_threshold = 8379989631760464847
assert 64*math.comb(capacity_threshold+2,3) >= r-1
assert 64*math.comb(capacity_threshold+1,3) < r-1
print("LINE_SUM_INDEPENDENT_PASS",json.dumps({
    "pair_targets":len(pair_sums),"triple_targets":len(support),
    "first16_targets":len(first16),"support_sha256":digest.hexdigest()
},sort_keys=True),flush=True)
