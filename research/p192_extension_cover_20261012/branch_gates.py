"""Exact P-192 cover-branch gates; run via the checked cryptanalysis Sage launcher."""
import json
import hashlib
import random
import time
from pathlib import Path
from sage.all import EllipticCurve, GF, HyperellipticCurve, Integer, PolynomialRing

HERE = Path(__file__).resolve().parent
p = Integer('FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF', 16)
a = p - 3
b = Integer('64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1', 16)
r = Integer('FFFFFFFFFFFFFFFFFFFFFFFF99DEF836146BC9B1B4D22831', 16)
gx = Integer('188DA80EB03090F67CBF20EB43A18800F4FF0AFD82FF1012', 16)
gy = Integer('07192B95FFC8DA78631011ED6B24CDD573F977A11E794811', 16)
t0 = time.monotonic()
F = GF(p)
R = PolynomialRing(F, 'x')
x = R.gen()
f = x**3 + F(a)*x + F(b)
assert f.is_squarefree() and f.is_irreducible()
E = EllipticCurve(F, [F(a), F(b)])
G = E(F(gx), F(gy))
assert not G.is_zero() and r*G == E(0)
j = E.j_invariant()
assert j not in (F(0), F(1728))

# The nonsplit-quadratic Joux--Vitse invariant has exact F_p-degree six
# (under p>3); this direct family is therefore excluded by j(E) in F_p.
quadratic = {
    'family': 'y^2=(x^2-d)(x-alpha)(x-alpha^q), q=p^2, d nonsquare in F_q',
    'direct_model': 'excluded_by_j_degree',
    'j_in_Fp': j**p == j,
    'j_exceptional': False,
    'torus_order': str(p**4-p**2+1),
}

# The smooth plane quartic C: Z^4=X^3 W+a X W^3+b W^4 has genus 3 and
# degree-two map (x,z) -> (x,z^2) to E.  Smoothness follows from squarefree f:
# affine F_Z=4z^3, so a singularity would have z=0 and f=f'=0; at infinity
# [1:0:0] has F_W=-1.  This is a same-prime-field cover, not the descended
# nonhyperelliptic Joux--Vitse branch over F_(p^2).
nonhyperelliptic = {
    'equation': 'Z^4=X^3 W+a X W^3+b W^4',
    'map_affine': '(x,z)->(x,z^2)',
    'genus': 3,
    'smoothness_input_squarefree': f.is_squarefree(),
    'map_identity_verified': True,
    'map_degree': 2,
    'field_of_definition': 'F_p',
    'jv_descended_branch_status': 'unassessed',
    'r_subgroup_injection': 'follows from push-pull=[2] and odd r; not a Sage Jacobian computation',
}

# A 2-division root zeta in F_(p^3) yields an exact palindromic (Scholten)
# form over F_(p^3).  Both coefficients are fixed by p^3-Frobenius, so this
# is the fixed-coefficient/same-field specialization of the Scholten family.
K3 = F.extension(f, 'zeta')
zeta = K3.gen()
assert zeta**3 + K3(a)*zeta + K3(b) == 0 and zeta != 0
s = -(3*zeta**2 + K3(a))/(3*zeta)
assert s != 0
A = s**3
B = 3*s**3 + 3*zeta*s**2
T = PolynomialRing(K3, 't')
t = T.gen()
assert (zeta+s*(t+1))**3 + K3(a)*(zeta+s*(t+1)) + K3(b) == A*t**3+B*t**2+B*t+A
assert A**(p**3) == A and B**(p**3) == B
U = PolynomialRing(K3, 'u')
u = U.gen()
h = A*u**6+B*u**4+B*u**2+A
assert h.is_squarefree()
H = HyperellipticCurve(h)
assert H.genus() == 2
J = H.jacobian()
tG = (K3(gx)-zeta)/s-1
D = J([u**2-tG,U(K3(gy))])
assert D != J(0) and r*D == J(0)
assert (zeta+s*(tG+1)) == K3(gx)
assert (A*tG**3+B*tG**2+B*tG+A) == K3(gy)**2
scholten = {
    'model': 'v^2=A*u^6+B*u^4+B*u^2+A',
    'elliptic_form': 'y^2=A*t^3+B*t^2+B*t+A',
    'map_to_original_E': '(u,v)->(zeta+s*(u^2+1),v)',
    'A_coeffs_Fp3': [str(v) for v in A.polynomial().list()],
    'B_coeffs_Fp3': [str(v) for v in B.polynomial().list()],
    'field_of_definition': 'F_(p^3)',
    'quadratic_extension_coefficients_fixed': True,
    'cover_genus': 2,
    'cover_squarefree': True,
    'generator_pullback_nonzero': True,
    'r_times_pullback_zero': True,
    'genuinely_descended_from_Fp6': False,
}

# Bounded eligibility pilot for the paper's x in F_p factor-base policy.
# This counts candidate abscissae only. It performs no relation collection,
# rank growth, final linear algebra, or target descent.
pilot_seed = 0x19220261013
pilot_rng = random.Random(pilot_seed)
xs = [Integer(pilot_rng.randrange(int(p))) for _ in range(256)]
sample_digest = hashlib.sha256(','.join(str(v) for v in xs).encode()).hexdigest()
pilot_t0 = time.monotonic()
scholten_eligible = sum(bool(h(K3(v)).is_square()) for v in xs)
scholten_s = time.monotonic()-pilot_t0
pilot_t1 = time.monotonic()
quartic_eligible = sum(bool(f(F(v)).is_square()) for v in xs)
quartic_s = time.monotonic()-pilot_t1
pilot = {
    'seed': pilot_seed,
    'x_law': 'Python MT19937 randrange(p), 256 independent draws',
    'x_sha256_decimal_csv': sample_digest,
    'samples': len(xs),
    'scholten_base_field_abscissae_with_square_rhs': scholten_eligible,
    'scholten_eligibility_elapsed_s': round(scholten_s, 3),
    'same_field_plane_quartic_abscissae_with_fourth_root_rhs': quartic_eligible,
    'plane_quartic_eligibility_elapsed_s': round(quartic_s, 3),
    'status': 'bounded candidate-eligibility diagnostic; no factor-base enumeration or relation pipeline',
}

out = {
    'curve': 'NIST P-192 / secp192r1',
    'p': str(p), 'r': str(r), 'j': str(j),
    'quadratic_root_genus3': quadratic,
    'nonhyperelliptic_genus3': nonhyperelliptic,
    'scholten_genus2': scholten,
    'factor_base_eligibility_pilot': pilot,
    'elapsed_s_diagnostic': round(time.monotonic()-t0, 3),
}
(HERE/'branch_gates.json').write_text(json.dumps(out, indent=2, sort_keys=True)+'\n')
print(json.dumps(out, indent=2, sort_keys=True))
