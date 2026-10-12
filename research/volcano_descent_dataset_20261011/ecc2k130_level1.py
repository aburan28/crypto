#!/usr/bin/env python3
"""ECC2K-130 crater -> level 1 of the 263-volcano, over F_{2^131}.

The Koblitz curve y^2 + xy = x^3 + 1 has j = 1 (in F_2) and End = O_K,
K = Q(sqrt -7).  With pi_2^2 + pi_2 + 2 = 0 the Frobenius over F_{2^131} is
pi_2^131 = (t + f sqrt(-7))/2 with f = 263 * 146505763881528721.  263 | f, so
pi_q acts as a scalar on E[263]: all 264 subgroups of order 263 are rational.
Predicted: Phi_263(X, 1) splits into 264 linear factors over F_{2^131}; the
root j = 1 appears twice (263 splits in K: two horizontal isogenies) and the
other 262 roots are distinct and are the 262 curves with End = O_263
(h(O_263) = 262).

PARI's polmodular only takes an evaluation point in the prime subfield in
characteristic 2, so this script can take exactly one step (from j = 1).  The
second step from a level-1 curve needs Phi_263 as a bivariate polynomial over
F_2, which this PARI does not provide, and the P-volcano step needs a
degree-1.5e17 isogeny, which nothing provides.

Writes ecc2k130_level1.json and ecc2k130_level1_roots.txt next to this file.
"""
import json
import os
import time

from cypari import pari

HERE = os.path.dirname(os.path.abspath(__file__))
pari.allocatemem(6 * 10 ** 9)
M = 131
L = 263

# Frobenius of the Koblitz curve over F_{2^M} from pi_2 = (-1 + sqrt(-7))/2, exact integers.
def mul(x, y):                       # elements (a + b sqrt(-7))/2 as (a, b)
    a, b = x; c, d = y
    return ((a * c - 7 * b * d) // 2, (a * d + b * c) // 2)
P = (2, 0)
for _ in range(M):
    P = mul(P, (-1, 1))
t, f = P
q = 2 ** M
assert t * t - 4 * q == -7 * f * f
N = q + 1 - t
fac = pari.factor(f)
fac = {int(fac[i, 0]): int(fac[i, 1]) for i in range(int(fac.nrows()))}
assert fac.get(L) == 1
print("t =", t, " f =", f, "=", fac, flush=True)
print("N = 4 * r, r =", N // 4, " prime:", bool(pari.isprime(N // 4)), flush=True)
print("h(O_263) =", int(pari.qfbclassno(-7 * L * L)), flush=True)

pari("T = ffinit(2, %d, 't)" % M)
t0 = time.time()
Phi = pari.polmodular(L, 0, pari("Mod(1, 2)"))          # Phi_263(X, 1) over F_2
build = time.time() - t0
t1 = time.time()
F = pari.factorff(pari.lift(Phi), 2, pari("T"))
fac_s = time.time() - t1
degs = [int(pari.poldegree(F[0][i])) for i in range(len(F[0]))]
exps = [int(F[1][i]) for i in range(len(F[0]))]
roots = []
for i, d in enumerate(degs):
    if d == 1:
        roots.append(str(pari.liftall(pari.polcoef(F[0][i], 0))))   # char 2: X + r; liftall: Mod(1,2) -> 1
n_lin = sum(e for d, e in zip(degs, exps) if d == 1)
res = dict(
    curve="y^2 + xy = x^3 + 1 over F_(2^131)", j=1, trace_t=t, conductor_f=f, f_factorisation=fac,
    N=N, ell=L, class_number_O_263=int(pari.qfbclassno(-7 * L * L)),
    modular_polynomial_degree=int(pari.poldegree(Phi)), build_seconds=round(build, 1), factor_seconds=round(fac_s, 2),
    factors_degrees=sorted(set(degs)), linear_factors_with_multiplicity=n_lin,
    distinct_roots=len(set(roots)), horizontal_root_multiplicity=None,
)
res["horizontal_root_multiplicity"] = [e for i, (d, e) in enumerate(zip(degs, exps))
                                       if d == 1 and str(pari.liftall(pari.polcoef(F[0][i], 0))) == "1"]
res["descending_distinct_roots"] = len([r for r in set(roots) if r != "1"])
json.dump(res, open(os.path.join(HERE, "ecc2k130_level1.json"), "w"), indent=1)
open(os.path.join(HERE, "ecc2k130_level1_roots.txt"), "w").write("\n".join(sorted(set(roots))) + "\n")
print(json.dumps(res, indent=1))
assert n_lin == L + 1 and len(set(roots)) == L and res["descending_distinct_roots"] == res["class_number_O_263"]
print("OK: 264 rational kernels = 2 horizontal + 262 descending = h(O_263)")
