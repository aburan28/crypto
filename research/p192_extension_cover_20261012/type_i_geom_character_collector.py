"""Frozen P-192 Type I quadratic-character query law on exact quotient classes.

The two inverse quotient fibers above aG and 2aG must each have a nonzero
square quartic RHS already in F_(p^3). Every trial scalar remains charged.
"""
import hashlib
import json
import random
import time
from pathlib import Path

from sage.all import Integer, GF, matrix, prod
import type_i_divisor_transfer as dt
import verify_type_i_quotient_base as vb

HERE = Path(__file__).resolve().parent
p, r, L = dt.p, dt.r, dt.L
G, O, base = vb.G, vb.O, vb.base
assert len(base) == 16 and all(v[0] == v[1] == v[2] for v in base)
N = (-2*G, -2*G, -2*G)
assert all(r*Q == O and Q != O for Q in N)

def character_test(P):
    if P == O or P[0] == dt.c0:
        return False
    X = dt.inf + dt.d0/(P[0]-dt.c0)
    rhs = prod(X-z**p for z in dt.quartic_roots)
    assert rhs**(p**3) == rhs
    return rhs != 0 and rhs**((p**3-1)//2) == 1

pairs = {}
for i in range(16):
    for j in range(i, 16):
        key = str(base[i][0] + base[j][0])
        pairs.setdefault(key, []).append((i, j))
assert sum(map(len, pairs.values())) == 136

seed = 0x19220261023
rng = random.Random(seed)
cap = 512
target_count = 64
attempted = []
accepted = []
first_reject = second_reject = 0
hits = []
start = time.monotonic_ns()
for attempt_index in range(cap):
    a = Integer(rng.randrange(1, int(r)))
    attempted.append(a)
    P = a*G
    if not character_test(P):
        first_reject += 1
        continue
    if not character_test(2*P):
        second_reject += 1
        continue
    accepted.append(a)
    target = tuple(a*Q for Q in N)
    for k in range(16):
        for i, j in pairs.get(str(target[0]-base[k][0]), []):
            assert tuple(base[i][m]+base[j][m]+base[k][m] for m in range(3)) == target
            hits.append({'query_index': len(accepted)-1, 'attempt_index': attempt_index,
                         'scalar': str(a), 'indices': [i,j,k]})
    if len(accepted) == target_count:
        break
elapsed = time.monotonic_ns()-start
assert len(accepted) == target_count
assert len(attempted) == first_reject + second_reject + len(accepted)
rows = []
for hit in hits:
    row = [0]*16
    for i in hit['indices']:
        row[i] += 1
    rows.append(row)
rank = int(matrix(GF(r), rows).rank()) if rows else 0
digest = lambda seq: hashlib.sha256(','.join(map(str,seq)).encode()).hexdigest()
out = {
    'status': 'bounded_type_i_f3_character_query_probe',
    'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'base_verifier_source_sha256': hashlib.sha256((HERE/'verify_type_i_quotient_base.py').read_bytes()).hexdigest(),
    'base_receipt_sha256': hashlib.sha256((HERE/'type_i_quotient_base_probe.json').read_bytes()).hexdigest(),
    'proposal_id': 'Q-P192-TypeI-Fp-quotient-16-character-Fp3', 'candidate_id': None,
    'p': str(p), 'r': str(r), 'seed': seed, 'attempt_cap': cap,
    'query_law': 'first 64 seeded nonzero a<r for which quartic RHS at inverse images of aG and 2aG are nonzero squares in F_(p^3)',
    'attempted_scalar_sha256_decimal_csv': digest(attempted),
    'accepted_scalar_sha256_decimal_csv': digest(accepted),
    'accepted_scalars_decimal': list(map(str,accepted)),
    'attempted_scalars': len(attempted), 'first_character_rejects': first_reject,
    'second_character_rejects': second_reject, 'accepted_queries': len(accepted),
    'base_classes': len(base), 'pair_index_entries': 136,
    'exact_three_sum_matches': len(hits), 'relation_witnesses': hits,
    'quotient_incidence_rank': rank, 'attempt_filter_match_ns': elapsed,
    'boundary': 'exact quotient-coordinate classes; geometric Jacobian reduction and inverse isogeny remain unimplemented',
}
(HERE/'type_i_geom_character_collector.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps(out,indent=2,sort_keys=True))
