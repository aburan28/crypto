"""Independent point-coordinate replay of the P-192 Type I norm certificate."""
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, prod
import type_i_generator_map as gm

HERE = Path(__file__).resolve().parent
producer = json.loads((HERE / 'type_i_norm_quotients.json').read_text())
fibers = json.loads((HERE / 'type_i_divisor_transfer.json').read_text())
assert producer['source_sha256'] == hashlib.sha256(
    (HERE / 'type_i_norm_quotients.py').read_bytes()).hexdigest()
assert fibers['source_sha256'] == hashlib.sha256(
    (HERE / 'type_i_divisor_transfer.py').read_bytes()).hexdigest()
assert producer['divisor_fiber_sha256'] == fibers['fiber_sha256']

L, p = gm.L, gm.p
q = p**2
r = gm.Integer('FFFFFFFFFFFFFFFFFFFFFFFF99DEF836146BC9B1B4D22831', 16)
E = EllipticCurve(L, [L(gm.a), L(gm.b)])
G = E(L(gm.gx), L(gm.gy))
O = E(0)
assert r.is_prime() and r*G == O and G != O
T = -2*G
assert T != O and r*T == O

def decode(c):
    return sum((L(gm.Integer(v))*L.gen()**i for i,v in enumerate(c)), L(0))

def decode_point(row):
    return tuple(decode(c) for c in row)

positive = [decode_point(row) for row in fibers['fiber_G_Fp6']]
negative = [decode_point(row) for row in fibers['fiber_2G_Fp6']]
assert len(positive) == len(negative) == 2

def project(P, i):
    X,Y0,Y1,Y2 = P
    z = gm.z**(p**i)
    y = Y0 + z*Y1 + z**2*Y2
    rhs = prod(X-v**(p**i) for v in gm.rl)
    assert y**2 == rhs
    rho = gm.inf**(p**i)
    assert X != rho
    xE = gm.c0**(p**i) + gm.d0**(p**i)/(X-rho)
    yE = gm.kappa**(p**i)*y/(X-rho)**2
    return E(xE,yE)

def encode(P):
    if P == O:
        return None
    return [[str(v) for v in P[k].polynomial().list()] for k in (0,1)]

vectors = []
for j in range(3):
    pos = [tuple(c**(q**j) for c in P) for P in positive]
    neg = [tuple(c**(q**j) for c in P) for P in negative]
    row = [sum((project(P,i) for P in pos),O)
           - sum((project(P,i) for P in neg),O) for i in range(3)]
    assert row[(-j)%3] == T
    assert all(row[i] == O for i in range(3) if i != (-j)%3)
    assert [encode(P) for P in row] == producer['three_quotient_pushforward_vectors'][j]
    vectors.append(row)

norm = [sum((vectors[j][i] for j in range(3)),O) for i in range(3)]
assert norm == [T,T,T]
assert [encode(P) for P in norm] == producer['norm_pushforward_vector']
assert all(P != O and r*P == O for P in norm)
assert all(tuple(c**(q**3) for c in P) == P for P in positive+negative)

out = {
    'status':'independent_type_i_norm_quotient_replay_passed',
    'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'producer_receipt_sha256':hashlib.sha256(
        (HERE / 'type_i_norm_quotients.json').read_bytes()).hexdigest(),
    'fiber_receipt_sha256':hashlib.sha256(
        (HERE / 'type_i_divisor_transfer.json').read_bytes()).hexdigest(),
    'all_nine_component_pushforwards_replayed':True,
    'norm_vector_diagonal_nonzero_order_r':True,
    'three_step_frobenius_cycle_replayed':True,
}
(HERE / 'type_i_norm_verification.json').write_text(
    json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps(out,indent=2,sort_keys=True))
