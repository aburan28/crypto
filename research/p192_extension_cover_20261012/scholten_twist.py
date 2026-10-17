"""Construct a nonfixed-coefficient Scholten genus-2 descent for P-192."""
import hashlib
import json
import random
from pathlib import Path

from sage.all import GF, HyperellipticCurve, Integer, PolynomialRing

HERE = Path(__file__).resolve().parent
p = Integer('FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF', 16)
a = p-3
b = Integer('64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1', 16)
F = GF(p)
R = PolynomialRing(F, 'x')
x = R.gen()
f = x**3+F(a)*x+F(b)
assert f.is_irreducible() and p % 3 == 2
K3 = F.extension(f, 'zeta')
zeta = K3.gen()
s = -(3*zeta**2+K3(a))/(3*zeta)
assert s != 0
A = s**3
B = 3*s**3+3*zeta*s**2

# Since p^3 = 2 mod 3, t^2+t+1 is irreducible over F_(p^3).
V = PolynomialRing(K3, 'v')
v = V.gen()
L = K3.extension(v**2+v+1, 'omega')
omega = L.gen()
assert omega**3 == 1 and omega != 1
sigma = lambda z: z**(p**3)
assert sigma(omega) == omega**2
eta = omega-omega**2
assert sigma(eta) == -eta and eta**2 == -3

P = PolynomialRing(L, 'T')
T = P.gen()
AL, BL = L(A), L(B)
twisted_cubic = AL*T**3+BL*omega**2*T**2+BL*omega*T+AL
sigma_cubic = P([sigma(c) for c in twisted_cubic.list()])
assert sigma_cubic == T**3*twisted_cubic(1/T)
assert BL*omega**2 != sigma(BL*omega**2)

# Twisted Galois action on v^2=f(u^2): (u,v)->(1/u,v/u^3).
# T=eta*(u-1)/(u+1) and W=v/(u+1)^3 are invariant coordinates.
numerator = (
    AL*(eta+T)**6 + BL*omega**2*(eta+T)**4*(eta-T)**2
    + BL*omega*(eta+T)**2*(eta-T)**4 + AL*(eta-T)**6
)
descended_L = numerator/(2*eta)**6
assert descended_L.degree() == 6
assert all(sigma(c) == c for c in descended_L.list())
descended_coeffs = []
for c in descended_L.list():
    cs = c.lift().list()
    assert len(cs) <= 1
    descended_coeffs.append(K3(cs[0]) if cs else K3(0))
RK = PolynomialRing(K3, 'T')
descended_K = RK(descended_coeffs)
assert descended_K.is_squarefree()
H = HyperellipticCurve(descended_K)
assert H.genus() == 2

pilot_rng = random.Random(0x19220261013)
xs = [Integer(pilot_rng.randrange(int(p))) for _ in range(256)]
pilot_digest = hashlib.sha256(','.join(str(v) for v in xs).encode()).hexdigest()
pilot_eligible = sum(bool(descended_K(K3(x0)).is_square()) for x0 in xs)

# Verify the rational map from the descended genus-2 curve to P-192 after
# base change to F_(p^6), by eliminating W^2 with its exact sextic equation.
FL = P.fraction_field()
Tf = FL(T)
u = (eta+Tf)/(eta-Tf)
elliptic_x = L(zeta)+L(s)*(omega*u**2+1)
elliptic_y_squared = FL(descended_L)*(2*eta/(eta-Tf))**6
assert elliptic_y_squared == elliptic_x**3+L(a)*elliptic_x+L(b)

out = {
    'status': 'verified_nonfixed_scholten_genus2_descent',
    'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'p': str(p), 'base_field': 'F_(p^3)', 'extension_field': 'F_(p^6)',
    'omega_equation': 'omega^2+omega+1=0',
    'eta_equation': 'eta=omega-omega^2, eta^2=-3, sigma(eta)=-eta',
    'twisted_scholten_a_coefficients_Fp3': [str(c) for c in A.polynomial().list()],
    'twisted_scholten_b_coefficients_Fp6_over_Fp3': [str(c) for c in (BL*omega**2).lift().list()],
    'b_not_fixed_by_sigma': True,
    'descended_sextic_coefficients_Fp3': [[str(v) for v in c.polynomial().list()] for c in descended_coeffs],
    'descended_sextic_squarefree': True, 'genus': 2,
    'descent_map': 'u=(eta+T)/(eta-T), v=W*(2*eta/(eta-T))^3',
    'elliptic_map_x': 'zeta+s*(omega*u^2+1)',
    'elliptic_map_identity_verified': True,
    'pilot_x_law': 'Python MT19937 randrange(p), seed 0x19220261013, 256 independent draws',
    'pilot_x_sha256_decimal_csv': pilot_digest,
    'pilot_Fp_abscissae_with_Fp3_point': pilot_eligible,
    'boundary': 'literal F_(p^3) genus-2 model and map checked; Jacobian conorm-norm image, relation pipeline and target cost unmeasured',
}
(HERE/'scholten_twist.json').write_text(json.dumps(out, indent=2, sort_keys=True)+'\n')
print(json.dumps(out, indent=2, sort_keys=True))
