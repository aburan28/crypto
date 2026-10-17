"""Independent frozen-law replay and exhaustive triple check for Type I."""
import hashlib
import itertools
import json
import random
from pathlib import Path

from sage.all import Integer, prod
import type_i_divisor_transfer as dt
import verify_type_i_quotient_base as vb

HERE = Path(__file__).resolve().parent
receipt = json.loads((HERE/'type_i_geom_character_collector.json').read_text())
assert receipt['source_sha256'] == hashlib.sha256((HERE/'type_i_geom_character_collector.py').read_bytes()).hexdigest()
assert receipt['base_receipt_sha256'] == hashlib.sha256((HERE/'type_i_quotient_base_probe.json').read_bytes()).hexdigest()
p, r, L = dt.p, dt.r, dt.L
G, O, base = vb.G, vb.O, vb.base
assert str(p) == receipt['p'] and str(r) == receipt['r']

def is_selected(P):
    if P == O or P[0] == dt.c0:
        return False
    x = dt.inf + dt.d0/(P[0]-dt.c0)
    f = L(1)
    for root in dt.quartic_roots:
        f *= x-root**p
    assert f**(p**3) == f
    return f != 0 and f**((p**3-1)//2) == 1

rng = random.Random(receipt['seed'])
attempted = []
accepted = []
reject1 = reject2 = 0
for _ in range(receipt['attempt_cap']):
    a = Integer(rng.randrange(1, int(r)))
    attempted.append(a)
    P = a*G
    if not is_selected(P):
        reject1 += 1
    elif not is_selected(P+P):
        reject2 += 1
    else:
        accepted.append(a)
        if len(accepted) == receipt['accepted_queries']:
            break
sha = lambda xs: hashlib.sha256(','.join(map(str,xs)).encode()).hexdigest()
assert len(attempted) == receipt['attempted_scalars']
assert reject1 == receipt['first_character_rejects']
assert reject2 == receipt['second_character_rejects']
assert sha(attempted) == receipt['attempted_scalar_sha256_decimal_csv']
assert sha(accepted) == receipt['accepted_scalar_sha256_decimal_csv']
assert list(map(str,accepted)) == receipt['accepted_scalars_decimal']

triples = list(itertools.combinations_with_replacement(range(len(base)), 3))
assert len(triples) == 816
triple_values = [tuple(base[i][m]+base[j][m]+base[k][m] for m in range(3))
                 for i,j,k in triples]
N = (-2*G,-2*G,-2*G)
assert all(r*P == O and P != O for P in N)
hits = []
for query_index, a in enumerate(accepted):
    target = tuple(a*P for P in N)
    for indices, value in zip(triples,triple_values):
        if value == target:
            hits.append({'query_index':query_index, 'scalar':str(a), 'indices':list(indices)})
assert len(hits) == receipt['exact_three_sum_matches']
assert [(v['query_index'],v['scalar'],v['indices']) for v in receipt['relation_witnesses']] == [
    (v['query_index'],v['scalar'],v['indices']) for v in hits]
out = {'status':'verified_type_i_f3_character_query_probe',
       'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
       'producer_sha256':receipt['source_sha256'],
       'attempted_scalars_replayed':len(attempted),
       'accepted_queries_replayed':len(accepted),
       'unordered_triples_enumerated_per_query':len(triples),
       'exact_matches':len(hits)}
(HERE/'type_i_geom_character_verification.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps(out,indent=2,sort_keys=True))
