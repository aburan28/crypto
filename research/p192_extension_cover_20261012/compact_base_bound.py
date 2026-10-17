"""Exact coverage bound for the frozen P-192 unguided three-sum policy."""
import hashlib
import json
from math import comb
from pathlib import Path

HERE = Path(__file__).resolve().parent
run = json.loads((HERE/'scholten_compact_base.json').read_text())
assert run['status'] == 'bounded_exact_compact_base_relation_probe'
B = run['base_classes']
r = int(run['r'])
queries = run['ordinary_queries']
assert B == 16 and queries == 64 and r > 1
triples = comb(B+2,3)
assert triples == 816

# Any unordered three-sum with repetition is one target class.  For uniform
# nonzero a*N, at most `triples` of the r-1 possible targets can hit a fixed
# base.  A union bound covers any number of independently sampled queries.
numerator, denominator = queries*triples, r-1
lo,hi = 0,1<<70
while lo < hi:
    mid = (lo+hi)//2
    if comb(mid+2,3)*2 >= r-1:
        hi = mid
    else:
        lo = mid+1
half_coverage_base_min = lo
assert comb(lo+2,3)*2 >= r-1
assert comb(lo+1,3)*2 < r-1

out = {
    'status':'exact_three_sum_coverage_bound',
    'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'probe_source_sha256':run['source_sha256'],
    'base_classes':B,'queries':queries,'subgroup_order':str(r),
    'unordered_three_sums_with_repetition':triples,
    'per_query_probability_upper_bound':f'{triples}/{denominator}',
    'any_of_64_probability_union_upper_bound':f'{numerator}/{denominator}',
    'minimum_base_for_half_single_query_coverage_by_count':half_coverage_base_min,
    'corresponding_pair_index_entries':half_coverage_base_min*(half_coverage_base_min+1)//2,
    'scope':'fixed base, exact unguided three-sum, uniform nonzero subgroup target; no bound on structured decomposition algorithms',
}
(HERE/'compact_base_bound.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps(out,indent=2,sort_keys=True))
