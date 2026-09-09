"""Independent integer/combinatorial read-back of the mathematical artifacts."""
from functools import lru_cache
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent

def rows(name):
    return [json.loads(line) for line in (ROOT / name).read_text().splitlines()]

@lru_cache(None)
def partition(n, a, r):
    t = 1 if a == 1 else -1
    previous, current = 2, t
    for _ in range(2, n + 1):
        previous, current = current, t * current - 2 * previous
    cardinality = (1 << n) + 1 - current
    assert cardinality % r == 0
    # Public roots of the characteristic polynomial; no unknown point scalars.
    candidates = [x for x in range(1, r) if (x*x-t*x+2) % r == 0 and pow(x,n,r) == 1]
    assert len(candidates) == 1
    eigenvalue = candidates[0]
    visited = {0}
    representatives = []
    size = None
    for k in range(1, r):
        if k in visited:
            continue
        orbit = {(sign * k * pow(eigenvalue, j, r)) % r for sign in (-1,1) for j in range(n)}
        if size is None:
            size = len(orbit)
        assert len(orbit) == size and not orbit & visited
        visited |= orbit
        representatives.append(k)
    assert visited == set(range(r))
    return representatives, size, cardinality

messages = []
exact_count = 0
for path in sorted(ROOT.glob('exact-n*.jsonl')):
    for d in rows(path.name):
        assert d['exhaustive_via_symmetry']
        representatives, weight, cardinality = partition(d['n'], d['a'], d['r'])
        assert d['target_scalars'] == representatives
        assert d['target_orbit_size'] == weight
        counts = d['representative_multiplicities']
        assert len(counts) == len(representatives)
        assert sum(counts)*weight + d['zero_sum_unordered_triples'] == d['cofactor_admissible_unordered_triples']
        assert sum(x > 0 for x in counts) == d['hits']
        assert d['covered_nonzero_targets'] == d['hits'] * weight
        assert d['represented_nonzero_targets'] == d['r'] - 1
        assert d['coefficient_matrix_rank'] == d['projected_signed_orbits']
        assert cardinality == d['r'] * int(d['cofactor'])
        exact_count += 1
messages.append(f'{exact_count} lifted-base exact records: complete independent target partitions, weighted Burnside totals, support and full recorded coefficient ranks verified.')

quotient_names = ['affine-n17.jsonl','affine-n19.jsonl','affine-conditioned-n19.jsonl',
                  'affine-conditioned-remaining-n19.jsonl','quotient-old-best-n19.jsonl','orbit-exchange-n19.jsonl']
quotient_exact = []
for name in quotient_names:
    for d in rows(name):
        if not d.get('exhaustive_via_symmetry'):
            continue
        representatives, weight, cardinality = partition(d['n'], d['a'], d['r'])
        assert cardinality == 2*d['r']
        assert representatives == d['target_scalars'] and weight == d['target_orbit_size']
        counts = d['representative_multiplicities']
        assert len(counts) == len(representatives)
        assert sum(counts)*weight+d['zero_sum_triples_including_identity'] == d['symmetric_cube_total']
        assert sum(x > 0 for x in counts) == d['hits_at_most_three']
        assert d['hits_exactly_three'] <= d['hits_at_most_three']
        assert d['quotient_points'] == 2*d['n']*d['projected_signed_orbits']
        quotient_exact.append(d)
messages.append(f'{len(quotient_exact)} quotient-base exact records: independent target partitions and complete symmetric-cube counts verified; two/three-summand support kept distinct.')

matched = [d for d in quotient_exact if d['identity']['type'] == 'cardinality_conditioned_affine_u']
assert len(matched) == 8 and len({d['identity']['attempt'] for d in matched}) == 8
assert all(d['quotient_points'] == 152 and d['projected_signed_orbits'] == 4 for d in matched)
assert max(d['hits_at_most_three'] for d in matched) == 6101
messages.append('All eight distinct cardinality-matched affine seed domains have exact results; best coverage is 6101/6909.')

control = rows('quotient-old-best-n19.jsonl')[0]
old = rows('exact-n19-drop-0.jsonl')[0]
assert control['hits_at_most_three'] == old['hits'] == 6165
assert control['quotient_points'] == 152 and old['points'] == 305
messages.append('Independent quotient-side oracle reproduces the earlier lifted-base support: 6165/6909.')
exchange = rows('orbit-exchange-n19.jsonl')
summary = [d for d in exchange if d['kind'] == 'orbit_exchange_summary'][0]
assert summary['pool_orbits'] == 37 and summary['candidates_including_incumbent'] == 1+4*(37-4) == 133
finalists = [d for d in exchange if d.get('exhaustive_via_symmetry')]
assert len(finalists) == 5
winner = max(finalists, key=lambda d:d['hits_at_most_three'])
assert winner['hits_at_most_three'] == 6184
assert (winner['hits_at_most_three'] - control['hits_at_most_three'])*38 == 722
messages.append('Exchange search: 133 fixed-size screened candidates, five exact finalists; selected result 6184/6909, a gain of 722 nonzero target points over the earlier base.')

for name in ['quotient-geometry.jsonl','quotient-general-binary.jsonl']:
    for d in rows(name):
        admitted = d.get('admissible_nonzero_u', d.get('admitted_nonzero_u'))
        assert d['affine_points'] - 1 == 4*admitted + d['u_zero_fiber_size']
messages.append('Eight exhaustive rational-fiber audits account for every affine point, including the exceptional u=0 fiber and rational 2-torsion.')
global_rows = rows('orbit-global-n19.jsonl')
optimum = next(d for d in global_rows if d['kind'] == 'global_pool_optimum')
selected = next(d for d in global_rows if d.get('exhaustive_via_symmetry'))
assert optimum['four_orbit_subsets_examined'] == 66045
assert optimum['triplet_support_types'] == 9139
assert sum(optimum['coverage_histogram'].values()) == 66045
assert max(map(int,optimum['coverage_histogram'])) == optimum['best_covered_target_orbits'] == 6257
assert optimum['coverage_histogram']['6257'] == 1
assert optimum['maximizing_subsets'] == [[6,9,22,28]]
assert len(set(tuple(p) for p in optimum['pool_representatives'])) == 37
bits = optimum['selected_support_words']
assert sum(word.bit_count() for word in bits) == 6257
assert selected['hits_at_most_three'] == 6257 and selected['hits_exactly_three'] == 6224
assert selected['quotient_points'] == 152 and selected['projected_signed_orbits'] == 4
reps,weight,cardinality = partition(selected['n'],selected['a'],selected['r'])
assert reps == selected['target_scalars'] and weight == 38
counts=selected['representative_multiplicities']
assert sum(counts)*weight+selected['zero_sum_triples_including_identity'] == selected['symmetric_cube_total']
assert all(bool((bits[i//64]>>(i%64))&1) == (count>0) for i,count in enumerate(counts))
assert all(word==0 for word in bits[(len(counts)+63)//64:])
assert bits[-1] >> (len(counts)%64) == 0
messages.append('All 66045 four-orbit subsets of the fixed 37-orbit pool were compared: unique maximum 6257/6909; selected support bit-for-bit matches the second exact enumeration method.')
(ROOT / 'verification.txt').write_text('\n'.join(messages) + '\n')
print('\n'.join(messages))
