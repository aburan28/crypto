"""Known-scalar shifted relation probe on the frozen P-192 Scholten base.

Each query has the form a*N + B_i + B_j. An exact three-sum match gives
the usable relation a*N = B_u+B_v+B_w-B_i-B_j. The scalar and shifts are
frozen independently of the match search.
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

def decode(coefficients):
    return R([sum((K3(Integer(v)) * zeta**i for i, v in enumerate(row)), K3(0))
              for row in coefficients])

N = J([decode(lower['u_coefficients_Fp3']), decode(lower['v_coefficients_Fp3'])])
assert N != J(0) and r * N == J(0)
rng = random.Random(prior['pilot_seed'])
xs = [Integer(rng.randrange(int(p))) for _ in range(256)]
assert hashlib.sha256(','.join(map(str, xs)).encode()).hexdigest() == prior['pilot_x_sha256_decimal_csv']
anchor_index = prior['anchor_pilot_index']
x0 = K3(xs[anchor_index])
w0 = descended_K(x0).sqrt()
base = []
start = time.monotonic_ns()
for item in prior['base_selection']:
    x = xs[item['pilot_index']]
    assert str(x) == item['x']
    w = descended_K(K3(x)).sqrt()
    assert [str(c) for c in w.polynomial().list()] == item['w_coefficients']
    U = (T - x0) * (T - K3(x))
    V = R(w0 + (w - w0) * (T - x0) / (K3(x) - x0))
    assert (V**2 - descended_K) % U == 0
    Q = cofactor * J([U, V])
    assert Q != J(0) and r * Q == J(0) and Q not in base
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

seed = 0x19220261017
query_rng = random.Random(seed)
scalars = [Integer(query_rng.randrange(1, int(r))) for _ in range(64)]
scalar_digest = hashlib.sha256(','.join(map(str, scalars)).encode()).hexdigest()
hits = []
start = time.monotonic_ns()
for q, a in enumerate(scalars):
    shift_i, shift_j = q % 16, (7*q + 3) % 16
    target = a*N + base[shift_i] + base[shift_j]
    for k in range(16):
        for i, j in pairs.get(str(target - base[k]), []):
            assert target == base[i] + base[j] + base[k]
            assert a*N == base[i] + base[j] + base[k] - base[shift_i] - base[shift_j]
            hits.append({'query_index': q, 'scalar': str(a),
                         'shift_indices': [shift_i, shift_j],
                         'three_sum_indices': [i, j, k]})
    if (q + 1) % 8 == 0:
        print(f'checked {q + 1}/64 shifted queries', flush=True)
query_ns = time.monotonic_ns() - start

rows = []
for hit in hits:
    row = [0] * 16
    for i in hit['three_sum_indices']:
        row[i] += 1
    for i in hit['shift_indices']:
        row[i] -= 1
    rows.append(row)
rank = int(matrix(GF(r), rows).rank()) if rows else 0

out = {
    'status': 'bounded_exact_known_scalar_shifted_relation_probe',
    'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'base_receipt_sha256': hashlib.sha256((HERE / 'scholten_compact_base.json').read_bytes()).hexdigest(),
    'lower_field_receipt_sha256': hashlib.sha256((HERE / 'scholten_lower_field.json').read_bytes()).hexdigest(),
    'proposal_id': 'Q-P192-Scholten-Fp3-compact-16-shifted',
    'candidate_id': None,
    'p': str(p), 'r': str(r),
    'base_classes': 16,
    'pair_index_entries': 136,
    'query_policy': '64 frozen nonzero uniform scalar multiples of N, shifted by B_(q mod 16)+B_((7q+3) mod 16); known formal scalar and base coefficients',
    'query_seed': seed,
    'query_scalar_sha256_decimal_csv': scalar_digest,
    'ordinary_queries': 64,
    'verified_relation_equations': len(hits),
    'relation_witnesses': hits,
    'useful_rank': rank,
    'base_rebuild_ns': base_ns,
    'pair_index_ns': pair_ns,
    'query_and_exact_match_ns': query_ns,
    'bound': 'each shift translates an at-most-816-image support; with uniform nonzero a, each query has match probability at most 816/(r-1)',
    'boundary': 'fixed 64-query stage probe; rank zero gates final LA and target descent',
}
(HERE / 'scholten_shifted_relation.json').write_text(json.dumps(out, indent=2, sort_keys=True) + '\n')
print(json.dumps(out, indent=2, sort_keys=True))
