"""Find an explicit quartic-to-P-192 isomorphism and map the standard generator.

This settles the elliptic-coordinate part of a Type I divisor transfer.  A
class in the genus-three Jacobian additionally needs a concrete conorm and
norm representation, which is not supplied by these coordinate identities.
"""
import hashlib
import itertools
import json
from pathlib import Path

from sage.all import GF, Integer, PolynomialRing, prod

HERE = Path(__file__).resolve().parent
gate = json.loads((HERE / 'type_i_gate.json').read_text())
assert gate['status'] == 'verified_type_i_branch_and_explicit_descended_affine_model'
p = Integer('FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF', 16)
a = p - 3
b = Integer('64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1', 16)
gx = Integer('188DA80EB03090F67CBF20EB43A18800F4FF0AFD82FF1012', 16)
gy = Integer('07192B95FFC8DA78631011ED6B24CDD573F977A11E794811', 16)
F = GF(p)
R = PolynomialRing(F, 'x')
x = R.gen()
f = x**3 + F(a)*x + F(b)
K3 = F.extension(f, 'zeta')
zeta = K3.gen()
beta = sum((K3(Integer(c))*zeta**i for i,c in enumerate(gate['beta_coefficients_Fp3'])), K3(0))
quartic_roots = [zeta, zeta**(p**2), beta, beta**(p**2)]
elliptic_roots = [zeta, zeta**p, zeta**(p**2)]
Q = PolynomialRing(K3, 'X')
X = Q.gen()
quartic = prod(X-v for v in quartic_roots)

models = []
for infinity in quartic_roots:
    finite = [v for v in quartic_roots if v != infinity]
    for targets in itertools.permutations(elliptic_roots):
        reciprocal = lambda v: 1/(v-infinity)
        D = (targets[0]-targets[1])/(reciprocal(finite[0])-reciprocal(finite[1]))
        C = targets[0] - D*reciprocal(finite[0])
        if C + D*reciprocal(finite[2]) == targets[2]:
            models.append((infinity, C, D, finite, targets))
assert models
infinity, C, D, finite, targets = models[0]
linear = C*(X-infinity)+D
cleared_cubic = (linear**3*(X-infinity)
                 + K3(a)*linear*(X-infinity)**3
                 + K3(b)*(X-infinity)**4)
kappa_squared_K3 = cleared_cubic.leading_coefficient()
assert cleared_cubic == kappa_squared_K3*quartic

print('constructing native F_(p^6) and quartic generator point', flush=True)
L = GF(p**6, 'k')
root = f.change_ring(L).roots(multiplicities=False)[0]
embed = lambda c: L(c.polynomial()(root))
z, be = embed(zeta), embed(beta)
rl = [embed(c) for c in quartic_roots]
inf, c0, d0 = embed(infinity), embed(C), embed(D)
assert L(gx) != c0
Xg = inf + d0/(L(gx)-c0)
Qg = prod(Xg-v for v in rl)
assert Xg != inf and Qg != 0

# The transformed cubic relation differs from the monic quartic by one
# constant square class.  Evaluate at the exact generator and verify the
# rational-function identity at a second independent coordinate.
kappa_squared = (L(gy)**2 * (Xg-inf)**4)/Qg
assert kappa_squared == embed(kappa_squared_K3)
assert kappa_squared.is_square()
kappa = kappa_squared.sqrt()
Yg = L(gy)*(Xg-inf)**2/kappa
assert Yg**2 == Qg
for offset in (1, 2, 3):
    test_x = L(gx+offset)
    Xt = inf + d0/(test_x-c0)
    assert (test_x**3+L(a)*test_x+L(b))*(Xt-inf)**4 == \
           kappa_squared*prod(Xt-v for v in rl)
assert c0+d0/(Xg-inf) == L(gx)
assert kappa*Yg/(Xg-inf)**2 == L(gy)

# The two further conjugate quartic equations govern rational points above
# (Xg,Yg) on the degree-two genus-three cover.
conjugate_square_classes = []
for i in (1, 2):
    values = [embed(v**(p**i)) for v in quartic_roots]
    rhs = prod(Xg-v for v in values)
    conjugate_square_classes.append(bool(rhs.is_square()))

encode = lambda c: [str(v) for v in c.polynomial().list()]
out = {
    'status': 'verified_type_i_quartic_generator_coordinates',
    'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'type_i_source_sha256': gate['source_sha256'],
    'mobius_models_found': len(models),
    'infinity_root_Fp3': encode(infinity),
    'mobius_C_Fp3': encode(C),
    'mobius_D_Fp3': encode(D),
    'elliptic_root_assignment_Fp3': [encode(v) for v in targets],
    'kappa_squared_Fp6': str(kappa_squared),
    'cleared_quartic_identity_verified': True,
    'kappa_Fp6': str(kappa),
    'generator_quartic_X_Fp6': str(Xg),
    'generator_quartic_Y_Fp6': str(Yg),
    'generator_quartic_equation_verified': True,
    'inverse_generator_replay_verified': True,
    'conjugate_quartic_rhs_square_over_Fp6': conjugate_square_classes,
    'boundary': 'elliptic quotient coordinate map only; genus-three Jacobian conorm and norm remain unimplemented',
}
(HERE / 'type_i_generator_map.json').write_text(json.dumps(out, indent=2, sort_keys=True) + '\n')
print(json.dumps(out, indent=2, sort_keys=True))
