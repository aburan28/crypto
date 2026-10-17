"""Exact P-192 Scholten three-sum probe with fresh geometric queries.

The 16-class base is reconstructed from the frozen compact-base receipt.
Queries use subsequent eligible F_p abscissae under the same point-pair
construction, with a disjoint abscissa and no base-class reuse.
"""
import hashlib
import json
import random
import time
from pathlib import Path

from sage.all import GF, HyperellipticCurve, Integer, PolynomialRing, matrix
from scholten_twist import K3, descended_K, p, zeta

HERE = Path(__file__).resolve().parent
prior = json.loads((HERE / 'scholten_compact_base.json').read_text())
lower = json.loads((HERE / 'scholten_lower_field.json').read_text())
assert prior['status'] == 'bounded_exact_compact_base_relation_probe'
assert lower['status'] == 'verified_scholten_norm_class_over_Fp3'
r = Integer(prior['r'])
cofactor = Integer(prior['candidate_cofactor'])
R = PolynomialRing(K3, 'T')
T = R.gen()
J = HyperellipticCurve(descended_K).jacobian()

seed = prior['pilot_seed']
rng = random.Random(seed)
xs = [Integer(rng.randrange(int(p))) for _ in range(256)]
assert hashlib.sha256(','.join(map(str, xs)).encode()).hexdigest() == prior['pilot_x_sha256_decimal_csv']

def fiber(x):
    rhs = descended_K(K3(x))
    return None if not rhs.is_square() else rhs.sqrt()

anchor_index = prior['anchor_pilot_index']
assert all(fiber(x) is None for x in xs[:anchor_index])
x0 = K3(xs[anchor_index])
w0 = fiber(xs[anchor_index])
assert w0 is not None

def projected_class(x, w):
    assert K3(x) != x0
    U = (T - x0) * (T - K3(x))
    V = R(w0 + (w - w0) * (T - x0) / (K3(x) - x0))
    assert (V**2 - descended_K) % U == 0
    Q = cofactor * J([U, V])
    assert r * Q == J(0)
    return Q

start = time.monotonic_ns()
base = []
for item in prior['base_selection']:
    index = item['pilot_index']
    x = xs[index]
    assert str(x) == item['x']
    w = fiber(x)
    assert w is not None
    assert [str(c) for c in w.polynomial().list()] == item['w_coefficients']
    Q = projected_class(x, w)
    assert Q != J(0) and Q not in base
    base.append(Q)
assert len(base) == 16
base_ns = time.monotonic_ns() - start

start = time.monotonic_ns()
pairs = {}
for i in range(16):
    for j in range(i, 16):
        pairs.setdefault(str(base[i] + base[j]), []).append((i, j))
assert sum(map(len, pairs.values())) == 136
pair_ns = time.monotonic_ns() - start

start = time.monotonic_ns()
query_indices = []
hits = []
zero_classes = 0
base_collisions = 0
duplicate_classes = 0
queries = []
last_base_index = prior['base_selection'][-1]['pilot_index']
for index in range(last_base_index + 1, len(xs)):
    if len(queries) == 64:
        break
    w = fiber(xs[index])
    if w is None:
        continue
    Q = projected_class(xs[index], w)
    if Q == J(0):
        zero_classes += 1
        continue
    if Q in base:
        base_collisions += 1
        continue
    if Q in queries:
        duplicate_classes += 1
        continue
    query_index = len(queries)
    queries.append(Q)
    query_indices.append(index)
    for k in range(16):
        for i, j in pairs.get(str(Q - base[k]), []):
            assert base[i] + base[j] + base[k] == Q
            hits.append({'query_index': query_index, 'pilot_index': index,
                         'base_indices': [i, j, k]})
    if len(queries) % 8 == 0:
        print(f'checked {len(queries)}/64 geometric queries', flush=True)
assert len(queries) == 64
query_ns = time.monotonic_ns() - start
rows = []
for hit in hits:
    row = [0] * 16
    for index in hit['base_indices']:
        row[index] += 1
    rows.append(row)
rank = int(matrix(GF(r), rows).rank()) if rows else 0

out = {
    'status': 'bounded_exact_structured_query_probe',
    'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'base_receipt_sha256': hashlib.sha256((HERE / 'scholten_compact_base.json').read_bytes()).hexdigest(),
    'lower_field_receipt_sha256': hashlib.sha256((HERE / 'scholten_lower_field.json').read_bytes()).hexdigest(),
    'proposal_id': 'Q-P192-Scholten-Fp3-compact-16-geometric-query',
    'candidate_id': None,
    'p': str(p), 'r': str(r),
    'base_classes': len(base),
    'pair_index_entries': 136,
    'query_policy': 'first 64 distinct nonzero projected classes from subsequent eligible frozen Fp abscissae, excluding base classes; fixed anchor and square-root convention',
    'pilot_seed': seed,
    'pilot_x_sha256_decimal_csv': prior['pilot_x_sha256_decimal_csv'],
    'query_pilot_indices': query_indices,
    'query_x_sha256_decimal_csv': hashlib.sha256(','.join(str(xs[i]) for i in query_indices).encode()).hexdigest(),
    'ordinary_geometric_queries': len(queries),
    'excluded_zero_classes': zero_classes,
    'excluded_base_collisions': base_collisions,
    'excluded_duplicate_classes': duplicate_classes,
    'verified_three_sum_relations': len(hits),
    'relation_witnesses': hits,
    'useful_rank': rank,
    'base_rebuild_ns': base_ns,
    'pair_index_ns': pair_ns,
    'query_projection_and_match_ns': query_ns,
    'boundary': 'stage-level geometric-query yield at a fixed 64-query cap; no transfer to uniform subgroup distribution or full P-192 cost',
}
(HERE / 'scholten_structured_query.json').write_text(json.dumps(out, indent=2, sort_keys=True) + '\n')
print(json.dumps(out, indent=2, sort_keys=True))
