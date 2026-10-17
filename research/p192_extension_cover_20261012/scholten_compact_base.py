"""Bounded exact compact-base relation probe on the P-192 Scholten Jacobian.

The queried subgroup and every selected base class are checked over F_(p^3).
The pair index tests all three-summand combinations of the frozen base for
each seeded ordinary query.  This is a capped stage probe,
not a complete IC candidate.
"""
import hashlib
import json
import random
import time
from pathlib import Path

from sage.all import GF, HyperellipticCurve, Integer, PolynomialRing, matrix, vector
from scholten_twist import K3, descended_K, p, zeta

HERE = Path(__file__).resolve().parent
proof = json.loads((HERE / 'scholten_lower_field.json').read_text())
assert proof['status'] == 'verified_scholten_norm_class_over_Fp3'
r = Integer('FFFFFFFFFFFFFFFFFFFFFFFF99DEF836146BC9B1B4D22831', 16)
R = PolynomialRing(K3, 'T')
T = R.gen()
J = HyperellipticCurve(descended_K).jacobian()

def decode(coefficients):
    return R([sum((K3(Integer(v)) * zeta**i for i, v in enumerate(row)), K3(0))
              for row in coefficients])

N = J([decode(proof['u_coefficients_Fp3']),
       decode(proof['v_coefficients_Fp3'])])
assert N != J(0) and r * N == J(0)

# #E(F_(p^6)) follows from its Frobenius trace recurrence.  We use its
# quotient by r as a projection *candidate*, and check each output directly
# rather than assume the cover Jacobian's order from this formula.
trace = p + 1 - r
s0, s1 = Integer(2), trace
for _ in range(2, 7):
    s0, s1 = s1, trace * s1 - p * s0
extension_order = p**6 + 1 - s1
assert extension_order % r == 0
cofactor = extension_order // r
assert cofactor % r != 0

SEED = 0x19220261013
rng = random.Random(SEED)
xs = [Integer(rng.randrange(int(p))) for _ in range(256)]
pilot_digest = hashlib.sha256(','.join(map(str, xs)).encode()).hexdigest()
assert pilot_digest == '92166e6477a7738fdae1381531c864499523d404f0ac6e3ba3cadb0cdc8cfaf8'
base = []
selection = []
anchor = None
start = time.monotonic_ns()
for index, x in enumerate(xs):
    if len(base) >= 16:
        break
    rhs = descended_K(K3(x))
    if not rhs.is_square():
        continue
    w = rhs.sqrt()
    if anchor is None:
        anchor = (K3(x), w, index)
        continue
    x0, w0, _ = anchor
    assert K3(x) != x0
    U = (T - x0) * (T - K3(x))
    V = R(w0 + (w - w0) * (T - x0) / (K3(x) - x0))
    assert (V**2 - descended_K) % U == 0
    P = J([U, V])
    Q = cofactor * P
    assert r * Q == J(0)
    if Q == J(0) or any(Q == prior for prior in base):
        continue
    base.append(Q)
    selection.append({'pilot_index': index, 'x': str(x),
                      'w_coefficients': [str(c) for c in w.polynomial().list()]})
    print(f'base {len(base)}/16, pilot index {index}', flush=True)
base_ns = time.monotonic_ns() - start
assert len(base) == 16

print('building complete pair index', flush=True)
start = time.monotonic_ns()
pairs = {}
for i in range(len(base)):
    for j in range(i, len(base)):
        pairs.setdefault(str(base[i] + base[j]), []).append((i, j))
pair_ns = time.monotonic_ns() - start

query_rng = random.Random(0x19220261015)
query_scalars = [Integer(query_rng.randrange(1, int(r))) for _ in range(64)]
query_digest = hashlib.sha256(','.join(map(str, query_scalars)).encode()).hexdigest()
print('checking 64 seeded ordinary queries', flush=True)
start = time.monotonic_ns()
hits = []
for query_index, scalar in enumerate(query_scalars):
    target = scalar * N
    assert target != J(0)
    for k in range(len(base)):
        for i, j in pairs.get(str(target - base[k]), []):
            assert base[i] + base[j] + base[k] == target
            hits.append({'query_index': query_index, 'indices': [i, j, k],
                         'scalar': str(scalar)})
query_ns = time.monotonic_ns() - start
rows = []
for hit in hits:
    row = [0] * len(base)
    for index in hit['indices']:
        row[index] += 1
    rows.append(row)
rank = matrix(GF(r), rows).rank() if rows else 0
assert rank <= len(base)

out = {
    'status': 'bounded_exact_compact_base_relation_probe',
    'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'lower_field_source_sha256': proof['source_sha256'],
    'proposal_id': 'Q-P192-Scholten-Fp3-compact-16',
    'candidate_id': None,
    'p': str(p), 'r': str(r),
    'extension_curve_order': str(extension_order),
    'candidate_cofactor': str(cofactor),
    'base_policy': 'fixed first eligible H(Fp3) anchor plus next 16 distinct degree-two divisors from 256 frozen Fp abscissae; one square root per abscissa; project by candidate cofactor and verify r-torsion',
    'anchor_pilot_index': anchor[2],
    'pilot_seed': SEED, 'pilot_x_sha256_decimal_csv': pilot_digest,
    'base_classes': len(base), 'base_selection': selection,
    'pair_index_entries': sum(len(v) for v in pairs.values()),
    'query_policy': '64 fixed Python MT19937 randrange(1,r) outputs multiplied by the exact lower-field norm class',
    'query_seed': 0x19220261015,
    'query_scalar_sha256_decimal_csv': query_digest,
    'ordinary_queries': len(query_scalars),
    'verified_three_sum_relations': len(hits),
    'relation_witnesses': hits,
    'useful_rank': int(rank),
    'base_setup_ns': base_ns,
    'pair_index_ns': pair_ns,
    'query_and_exact_match_ns': query_ns,
    'boundary': 'bounded 16-point three-sum probe; no final LA or unseen-target descent if no full-rank relation system',
}
(HERE / 'scholten_compact_base.json').write_text(json.dumps(out, indent=2, sort_keys=True) + '\n')
print(json.dumps(out, indent=2, sort_keys=True))
