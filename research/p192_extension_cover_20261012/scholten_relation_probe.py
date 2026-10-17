"""Bounded split-support probe on the exact descended P-192 Scholten class."""
import hashlib
import json
from pathlib import Path

from sage.all import GF, HyperellipticCurve, Integer, PolynomialRing
from scholten_twist import f, zeta as old_zeta, s as old_s, p, descended_K

HERE = Path(__file__).resolve().parent
transfer = json.loads((HERE/'scholten_transfer.json').read_text())
assert transfer['status'] == 'verified_scholten_pullback_and_norm_over_Fp6'
L = GF(p**6, 'k')
root = f.change_ring(L).roots(multiplicities=False)[0]
embed = lambda c: L(c.polynomial()(root))
zeta, s = embed(old_zeta), embed(old_s)
eta = L(-3).sqrt()
omega = (-1+eta)/2
P = PolynomialRing(L, 'T')
T = P.gen()
descended_L = P([embed(c) for c in descended_K.list()])
H = HyperellipticCurve(descended_L)
J = H.jacobian()
gx = Integer('188DA80EB03090F67CBF20EB43A18800F4FF0AFD82FF1012', 16)
gy = Integer('07192B95FFC8DA78631011ED6B24CDD573F977A11E794811', 16)
tG = (L(gx)-zeta)/s-1
U = P(omega*(eta+T)**2-tG*(eta-T)**2).monic()
V = P(L(gy)*(eta-T)**3/(2*eta)**3) % U
assert [str(c) for c in U.list()] == transfer['pullback_u_coefficients']
assert [str(c) for c in V.list()] == transfer['pullback_v_coefficients']
D = J([U,V])
sigma = lambda z: z**(p**3)
Ds = J([P([sigma(c) for c in U.list()]), P([sigma(c) for c in V.list()])])
N = D+Ds
assert N != J(0)

def is_square_in_subfield(z, size):
    if z == 0:
        return True
    return z**((size-1)//2) == 1

support = {'zero': 0, 'degree_one': 0, 'degree_two': 0,
           'split_over_Fp3': 0, 'split_over_Fp': 0}
Q = N
for a in range(1, 257):
    Uq = Q[0]
    degree = Uq.degree()
    if degree == 0:
        support['zero'] += 1
    elif degree == 1:
        support['degree_one'] += 1
        assert all(sigma(c) == c for c in Uq.list())
        support['split_over_Fp3'] += 1
        if all(c**p == c for c in Uq.list()):
            support['split_over_Fp'] += 1
    elif degree == 2:
        support['degree_two'] += 1
        assert all(sigma(c) == c for c in Uq.list())
        coeffs = Uq.list()
        disc = coeffs[1]**2-4*coeffs[0]*coeffs[2]
        if is_square_in_subfield(disc, p**3):
            support['split_over_Fp3'] += 1
            if all(c**p == c for c in coeffs) and is_square_in_subfield(disc, p):
                support['split_over_Fp'] += 1
    else:
        raise AssertionError(f'unexpected Mumford degree {degree}')
    Q += N

out = {
    'status': 'bounded_256_exact_scholten_split_support_probe',
    'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'transfer_source_sha256': transfer['source_sha256'],
    'query_law': 'successive a*N for a=1,...,256; N is verified nonzero r-torsion norm class',
    'ordinary_queries': 256,
    'support_counts': support,
    'Fp3_policy': 'reduced degree-2 Mumford polynomial splits over F_(p^3)',
    'Fp_policy': 'reduced Mumford support abscissae lie in F_p',
    'boundary': 'split-support diagnostic only; no indexed factor base, complete Jacobian relation equations, rank, LA or target descent',
}
(HERE/'scholten_relation_probe.json').write_text(json.dumps(out, indent=2, sort_keys=True)+'\n')
print(json.dumps(out, indent=2, sort_keys=True))
