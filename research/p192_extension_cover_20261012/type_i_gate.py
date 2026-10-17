"""Exact Momose--Chao Type I discriminant gate for standard P-192.

Run via the checked cryptanalysis Sage launcher. This verifies the paper's
Legendre-parameter criterion; it does not construct the descended plane
quartic or a relation algorithm.
"""
import hashlib
import json
import random
from pathlib import Path
from sage.all import EllipticCurve, GF, Integer, PolynomialRing, matrix, prod

HERE = Path(__file__).resolve().parent
p = Integer('FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF', 16)
a = p - 3
b = Integer('64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1', 16)
F = GF(p)
R = PolynomialRing(F, 'x')
x = R.gen()
f = x**3 + F(a)*x + F(b)
assert f.is_irreducible()
K3 = F.extension(f, 'zeta')
zeta = K3.gen()
roots = [zeta, zeta**p, zeta**(p**2)]
assert len(set(roots)) == 3
lam = (roots[2]-roots[0])/(roots[1]-roots[0])
assert lam not in (K3(0), K3(1))
norm = lam * lam**p * lam**(p**2)
trace = lam + lam**p + lam**(p**2)
assert norm == -1 and trace**p == trace
delta = (trace-1)**2 - 4*norm
assert delta**p == delta and delta != 0

# The q=p^2 trace/norm from F_(p^6) to F_(p^2) equal these expressions:
# Frobenius p^2 and p^4 permute the same three F_(p^3) conjugates.
delta_F = F(delta.polynomial()[0])
assert K3(delta_F) == delta
assert p % 4 == 3
S = PolynomialRing(F, 'w')
w = S.gen()
F2 = F.extension(w**2+1, 'i')
assert F2(delta_F).is_square()

j_legendre = 256*(lam**2-lam+1)**3/(lam**2*(lam-1)**2)
E = EllipticCurve(F, [F(a), F(b)])
assert j_legendre == K3(E.j_invariant())

# Momose--Chao equations (20)--(25): solve the Frobenius--Möbius fixed-point
# equation for beta and verify it directly rather than relying on a square
# discriminant alone.
q = p**2
eps = zeta
eps_q = eps**q
mu = (lam*eps-eps_q)/(lam-1)
nu = (lam*eps_q-eps)/(lam-1)
assert mu+nu == eps+eps_q
A = matrix(K3, [[nu, -eps*eps_q], [1, -mu]])
frobenius_matrix = lambda M, power: matrix(K3, 2, 2, [v**power for v in M.list()])
B = frobenius_matrix(A, q**2)*frobenius_matrix(A, q)*A
assert B[1,0] != 0
fixed_discriminant = (B[1,1]-B[0,0])**2+4*B[1,0]*B[0,1]
assert fixed_discriminant.is_square()
disc_root = fixed_discriminant.sqrt()
beta_roots = [
    (B[0,0]-B[1,1]+sgn*disc_root)/(2*B[1,0])
    for sgn in (1, -1)
]
assert len(beta_roots) == 2
valid_beta = []
for beta in beta_roots:
    denominator = A[1,0]*beta+A[1,1]
    if denominator == 0:
        continue
    if beta**q != (A[0,0]*beta+A[0,1])/denominator:
        continue
    if len({eps, eps_q, beta, beta**q}) != 4:
        continue
    lambda_replay = ((beta-eps_q)/(beta-eps))*((beta**q-eps)/(beta**q-eps_q))
    if lambda_replay == lam:
        valid_beta.append(beta)
assert valid_beta
beta = valid_beta[0]
quartic_roots = [eps, eps_q, beta, beta**q]
assert len(set(quartic_roots)) == 4
Q = PolynomialRing(K3, 'X')
X = Q.gen()
quartic = prod(X-root for root in quartic_roots)
assert quartic.is_squarefree()

# An explicit F_p-descended affine model. Over an algebraic closure the
# K3-valued equation splits into three conjugate double covers. Their product
# is a square; the norm equation selects one connected degree-four component.
# Its smooth projective normalization has six branch points and genus three.
names = ('X', 'Y0', 'Y1', 'Y2')
PK = PolynomialRing(K3, names=names)
Xk, Y0k, Y1k, Y2k = PK.gens()
PF = PolynomialRing(F, names=names)
Xf, Y0f, Y1f, Y2f = PF.gens()
yK = Y0k+zeta*Y1k+zeta**2*Y2k
fk = prod(Xk-root for root in quartic_roots)
equation_K = yK**2-fk

def split_coefficients(g):
    parts = [PF(0), PF(0), PF(0)]
    for exponents, coeff in g.dict().items():
        monomial = prod(v**e for v,e in zip(PF.gens(), exponents))
        for i, c in enumerate(coeff.polynomial().list()):
            parts[i] += F(c)*monomial
    return parts

equations_F = split_coefficients(equation_K)
sigma_y = lambda power: Y0k+zeta**power*Y1k+(zeta**2)**power*Y2k
norm_y = sigma_y(1)*sigma_y(p)*sigma_y(p**2)
norm_roots = prod(Xk-root**(p**i) for root in (eps,beta) for i in range(3))
norm_parts = split_coefficients(norm_y-norm_roots)
assert norm_parts[1] == 0 and norm_parts[2] == 0
norm_equation_F = norm_parts[0]
assert len({root**(p**i) for root in (eps,beta) for i in range(3)}) == 6
fk_conjugates = [prod(Xk-root**(p**i) for root in quartic_roots) for i in range(3)]
assert prod(fk_conjugates) == norm_roots**2

def sparse_polynomial(poly):
    return [
        {'exponents': list(exponents), 'coefficient': str(coeff)}
        for exponents,coeff in sorted(poly.dict().items())
    ]

model = [sparse_polynomial(poly) for poly in equations_F+[norm_equation_F]]
model_digest = hashlib.sha256(json.dumps(model,sort_keys=True,separators=(',',':')).encode()).hexdigest()

# The x in F_p factor-base policy is screened on the actual descended cover.
# If the quartic RHS is a nonzero square in K3, exactly one of its two roots
# has norm N(x-alpha)N(x-beta), because [K3:F_p]=3 is odd.
pilot_rng = random.Random(0x19220261013)
xs = [Integer(pilot_rng.randrange(int(p))) for _ in range(256)]
pilot_digest = hashlib.sha256(','.join(str(v) for v in xs).encode()).hexdigest()
eligible = 0
for x0 in xs:
    rhs = quartic(K3(x0))
    U = (K3(x0)-eps)*(K3(x0)-beta)
    assert rhs == U**(q+1)
    y = U**((q+1)//2)
    assert y*y == rhs
    s_val = prod(K3(x0)-root**(p**i) for root in (eps,beta) for i in range(3))
    ny = y*y**p*y**(p**2)
    if ny != s_val:
        y = -y
    assert y*y == rhs and y*y**p*y**(p**2) == s_val
    eligible += 1

out = {
    'status': 'verified_type_i_branch_and_explicit_descended_affine_model',
    'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'p': str(p), 'q': str(q),
    'lambda_coefficients_Fp3': [str(v) for v in lam.polynomial().list()],
    'lambda_norm_Fp6_over_Fp2': str(norm),
    'lambda_trace_Fp6_over_Fp2': str(trace),
    'delta_reduced_Fp': str(delta_F),
    'delta_square_in_Fp': bool(delta_F.is_square()),
    'delta_nonzero': True,
    'delta_square_in_Fp2': True,
    'j_matches_P192': True,
    'paper_condition': 'Momose-Chao Lemma 7, Type I nonhyperelliptic cover criterion',
    'beta_valid_candidates': len(valid_beta),
    'beta_coefficients_Fp3': [str(v) for v in beta.polynomial().list()],
    'beta_frobenius_mobius_verified': True,
    'beta_cross_ratio_replay_verified': True,
    'quartic_roots_distinct': True,
    'quartic_coefficients_Fp3': [[str(v) for v in c.polynomial().list()] for c in quartic.list()],
    'descended_model_field': 'F_p (base-changes to F_(p^2))',
    'descended_model_coordinates': list(names),
    'descended_model_equations_sparse': model,
    'descended_model_sha256': model_digest,
    'descended_model_definition': 'three F_p coefficients of Y^2-f0(X)=0 with Y=Y0+zeta*Y1+zeta^2*Y2, plus Norm_Fp3/Fp(Y)=N(X-alpha)N(X-beta); take smooth projective normalization',
    'branch_union_size': 6,
    'genus_by_tame_degree_four_riemann_hurwitz': 3,
    'f0_f1_f2_equals_square_of_six_root_polynomial': True,
    'pilot_x_law': 'Python MT19937 randrange(p), seed 0x19220261013, 256 independent draws',
    'pilot_x_sha256_decimal_csv': pilot_digest,
    'pilot_descended_cover_eligible': eligible,
    'subfield_square_identity_replayed_all_samples': True,
    'boundary': 'exact descended affine model and 256-coordinate eligibility; literal Jacobian transfer, subgroup-usable base, relations, rank, LA and target descent unmeasured',
}
(HERE/'type_i_gate.json').write_text(json.dumps(out, indent=2, sort_keys=True)+'\n')
print(json.dumps(out, indent=2, sort_keys=True))
