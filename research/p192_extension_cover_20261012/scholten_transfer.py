"""Literal Jacobian pullback and norm gate for the descended Scholten cover."""
import hashlib
import json
from pathlib import Path

from sage.all import GF, HyperellipticCurve, Integer, PolynomialRing
from scholten_twist import f, K3, zeta as old_zeta, s as old_s, p, descended_K

HERE = Path(__file__).resolve().parent
r = Integer('FFFFFFFFFFFFFFFFFFFFFFFF99DEF836146BC9B1B4D22831', 16)
gx = Integer('188DA80EB03090F67CBF20EB43A18800F4FF0AFD82FF1012', 16)
gy = Integer('07192B95FFC8DA78631011ED6B24CDD573F977A11E794811', 16)
print('constructing native F_(p^6)', flush=True)
L = GF(p**6, 'k')
print('finding cubic embedding', flush=True)
roots = f.change_ring(L).roots(multiplicities=False)
assert len(roots) == 3
root = roots[0]
embed = lambda c: L(c.polynomial()(root))
zeta = embed(old_zeta)
s = embed(old_s)
eta = L(-3).sqrt()
omega = (-1+eta)/2
assert omega**2+omega+1 == 0 and omega**3 == 1
P = PolynomialRing(L, 'T')
T = P.gen()
descended_L = P([embed(c) for c in descended_K.list()])
print('constructing Jacobian', flush=True)
H = HyperellipticCurve(descended_L)
J = H.jacobian()
print('constructing pullback divisor', flush=True)
tG = (L(gx)-L(zeta))/L(s)-1
U = P(omega*(eta+T)**2-tG*(eta-T)**2).monic()
assert U.degree() == 2 and U.is_squarefree()
V = P(L(gy)*(eta-T)**3/(2*eta)**3) % U
assert (V**2-descended_L) % U == 0
D = J([U,V])
print('checking pullback order', flush=True)
assert D != J(0) and r*D == J(0)

sigma = lambda z: z**(p**3)
Us = P([sigma(c) for c in U.list()])
Vs = P([sigma(c) for c in V.list()])
assert (Vs**2-descended_L) % Us == 0
Ds = J([Us,Vs])
print('checking conjugate order', flush=True)
assert Ds != J(0) and r*Ds == J(0)
N = D+Ds
print('checking norm order', flush=True)
assert r*N == J(0)

out = {
    'status': 'verified_scholten_pullback_and_norm_over_Fp6',
    'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'r': str(r),
    'pullback_u_coefficients': [str(c) for c in U.list()],
    'pullback_v_coefficients': [str(c) for c in V.list()],
    'pullback_nonzero': bool(D != J(0)),
    'r_times_pullback_zero': bool(r*D == J(0)),
    'sigma_pullback_nonzero': bool(Ds != J(0)),
    'r_times_sigma_pullback_zero': bool(r*Ds == J(0)),
    'norm_nonzero': bool(N != J(0)),
    'r_times_norm_zero': bool(r*N == J(0)),
    'boundary': 'Jacobian class arithmetic over F_(p^6); descent of norm coefficients and complete relation pipeline need separate checks',
}
(HERE/'scholten_transfer.json').write_text(json.dumps(out, indent=2, sort_keys=True)+'\n')
print(json.dumps(out, indent=2, sort_keys=True))
