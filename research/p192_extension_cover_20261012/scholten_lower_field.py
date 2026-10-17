"""Extract and independently replay the P-192 Scholten norm class over F_(p^3)."""
import hashlib
import json
from pathlib import Path

from sage.all import GF, HyperellipticCurve, Integer, PolynomialRing, matrix, vector
from scholten_twist import F, K3, descended_K, f, p, s as old_s, zeta as old_zeta

HERE = Path(__file__).resolve().parent
old = json.loads((HERE / 'scholten_transfer.json').read_text())
assert old['status'] == 'verified_scholten_pullback_and_norm_over_Fp6'
r = Integer('FFFFFFFFFFFFFFFFFFFFFFFF99DEF836146BC9B1B4D22831', 16)
gx = Integer('188DA80EB03090F67CBF20EB43A18800F4FF0AFD82FF1012', 16)
gy = Integer('07192B95FFC8DA78631011ED6B24CDD573F977A11E794811', 16)

print('constructing native extension and lower-field embedding', flush=True)
L = GF(p**6, 'k')
root = f.change_ring(L).roots(multiplicities=False)[0]
embed = lambda c: L(c.polynomial()(root))
zeta, s = embed(old_zeta), embed(old_s)
eta = L(-3).sqrt()
omega = (-1 + eta) / 2
P = PolynomialRing(L, 'T')
T = P.gen()
H6 = HyperellipticCurve(P([embed(c) for c in descended_K.list()]))
J6 = H6.jacobian()
tG = (L(gx) - zeta) / s - 1
U = P(omega * (eta + T)**2 - tG * (eta - T)**2).monic()
V = P(L(gy) * (eta - T)**3 / (2 * eta)**3) % U
assert [str(c) for c in U.list()] == old['pullback_u_coefficients']
assert [str(c) for c in V.list()] == old['pullback_v_coefficients']
D = J6([U, V])
sigma = lambda c: c**(p**3)
Ds = J6([P([sigma(c) for c in U.list()]), P([sigma(c) for c in V.list()])])
N = D + Ds
assert N != J6(0) and r * N == J6(0)

# The embedded prime-field coordinates of 1, zeta, zeta^2 form a 6-by-3
# matrix.  A sigma-fixed coefficient has a unique preimage in its span.
basis = [L(1), zeta, zeta**2]
def prime_coordinates(c):
    coefficients = c.polynomial().list()
    return vector(F, coefficients + [F(0)] * (6 - len(coefficients)))

embedding_matrix = matrix(F, [prime_coordinates(c) for c in basis]).transpose()
assert embedding_matrix.nrows() == 6 and embedding_matrix.rank() == 3

def descend(c):
    assert sigma(c) == c
    coords = embedding_matrix.solve_right(prime_coordinates(c))
    result = sum((K3(coords[i]) * old_zeta**i for i in range(3)), K3(0))
    assert embed(result) == c
    return result

print('extracting Mumford coefficients', flush=True)
R3 = PolynomialRing(K3, 'T')
U3 = R3([descend(c) for c in N[0].list()])
V3 = R3([descend(c) for c in N[1].list()])
assert (V3**2 - descended_K) % U3 == 0
J3 = HyperellipticCurve(descended_K).jacobian()
N3 = J3([U3, V3])
assert N3 != J3(0)
print('checking r-torsion in J(F_(p^3))', flush=True)
assert r * N3 == J3(0)
assert J6([P([embed(c) for c in U3.list()]),
           P([embed(c) for c in V3.list()])]) == N

coefficients = lambda poly: [
    [str(v) for v in c.polynomial().list()] for c in poly.list()
]
out = {
    'status': 'verified_scholten_norm_class_over_Fp3',
    'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'transfer_source_sha256': old['source_sha256'],
    'field': 'F_(p^3) with zeta^3+a*zeta+b=0',
    'u_coefficients_Fp3': coefficients(U3),
    'v_coefficients_Fp3': coefficients(V3),
    'u_degree': int(U3.degree()),
    'norm_nonzero': True,
    'r_times_norm_zero_over_Fp3': True,
    'base_change_matches_Fp6_norm': True,
    'boundary': 'explicit lower-field class; subgroup-usable factor base and relation pipeline are separate gates',
}
(HERE / 'scholten_lower_field.json').write_text(json.dumps(out, indent=2, sort_keys=True) + '\n')
print(json.dumps(out, indent=2, sort_keys=True))
