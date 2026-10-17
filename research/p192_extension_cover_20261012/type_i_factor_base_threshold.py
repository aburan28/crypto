"""Exact base-size and pair-index thresholds for the P-192 Type I law."""
import hashlib
import json
from fractions import Fraction
from pathlib import Path

HERE=Path(__file__).resolve().parent
prior=json.loads((HERE/'type_i_character_density.json').read_text())
assert prior['status']=='checked_type_i_character_density_inputs'
L=int(prior['accepted_nonzero_scalar_lower_bound'])
assert L>0

def triples(B):
    return B*(B+1)*(B+2)//6

low,high=0,1
while 2*triples(high)<L:
    high*=2
while low+1<high:
    mid=(low+high)//2
    if 2*triples(mid)>=L:
        high=mid
    else:
        low=mid
B=high
assert B>0 and 2*triples(B-1)<L<=2*triples(B)
pairs=B*(B+1)//2
assert pairs>0
expected_one=Fraction(L,816)
expected_sixteen=16*expected_one
out={
    'status':'exact_type_i_fixed_three_sum_capacity_threshold',
    'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'density_receipt_sha256':hashlib.sha256((HERE/'type_i_character_density.json').read_bytes()).hexdigest(),
    'accepted_nonzero_scalar_lower_bound':str(L),
    'policy':'uniform nonzero scalars conditioned on the certified Fp3 character predicate; fixed B-class unordered three-sum lookup',
    'coverage_target':'at least one half of accepted targets',
    'minimum_actual_distinct_base_classes_for_capacity':str(B),
    'previous_base_three_sum_capacity':str(triples(B-1)),
    'threshold_base_three_sum_capacity':str(triples(B)),
    'full_unordered_pair_index_entries_at_threshold':str(pairs),
    'current_base_classes':16,
    'current_three_sum_support_cap':816,
    'current_per_accepted_query_hit_upper_bound':str(Fraction(816,L)),
    'expected_accepted_queries_lower_bound_for_one_hit':str(expected_one),
    'expected_accepted_queries_lower_bound_for_sixteen_hits':str(expected_sixteen),
    'boundary':'necessary support and direct pair-index capacities for this fixed-base policy; not a complete index-calculus cost or a bound on structured relation methods',
}
(HERE/'type_i_factor_base_threshold.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps(out,indent=2,sort_keys=True))
