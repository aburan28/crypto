# E2 (conductor gap): build the floor curves one level down from K_0 at a landed rung by the CM route.
#
#   sage e2_floor_curves_polclass.sage N ELL [OUT.json]
#
# Ring class polynomial H_{-7 ell^2} (degree ell - chi(ell)), reduced mod 2, roots in F_{2^N} give the
# j-invariants of the curves with End = Z + ell*O_K in the isogeny class of K_0 (a = 0, y^2 + xy = x^3 + 1).
# For each j: b = 1/j (ordinary binary curve y^2 + xy = x^3 + a2 x^2 + b has j = 1/b), and a2 in {0, 1}
# selects the twist; keep the twist whose order is #K_0 = 2^N + 1 - t, checked on random points.
# Output: one JSON record with the modulus, every j as an integer bitmask, b, a2, and the Frobenius orbit id,
# matching the encoding of research/isogeny_conductor_gap_20261007/e1_floor_panel.py (--model-b).
import json, sys
from sage.all import *

N = int(sys.argv[1]); ELL = int(sys.argv[2])
OUT = sys.argv[3] if len(sys.argv) > 3 else f"e2_floor_n{N}_ell{ELL}.json"

# trace of K_0 over F_{2^N} by the Lucas recurrence tau + taubar = -1, tau*taubar = 2
V = [2, -1]
for k in range(1, N):
    V.append(-V[k] - 2 * V[k - 1])
t = V[N]; order = 2**N + 1 - t
f2 = (t * t - 4 * 2**N) // (-7); f = ZZ(f2).isqrt(); assert f * f == f2 and f % ELL == 0

D = -7 * ELL**2
H = hilbert_class_polynomial(D)
print(f"n={N} ell={ELL} deg H={H.degree()} expected {ELL - kronecker(-7, ELL)}", file=sys.stderr)
F = GF(2**N, 'z', modulus='minimal_weight')
mod = F.modulus()
Hbar = H.change_ring(GF(2))
roots = Hbar.change_ring(F).roots(multiplicities=False)
print(f"roots in F_{{2^{N}}}: {len(roots)}", file=sys.stderr)

def to_int(x):
    return int(sum(int(c) << i for i, c in enumerate(x.polynomial().list())))

def frob_orbit_rep(x):
    best = x; y = x
    for _ in range(N - 1):
        y = y**2
        if to_int(y) < to_int(best):
            best = y
    return to_int(best)

recs = []
for j in roots:
    if j == 0:
        continue  # supersingular, not in the class
    b = 1 / j
    chosen = None
    for a2 in (F(0), F(1)):
        E = EllipticCurve(F, [1, a2, 0, 0, b])
        ok = True
        for _ in range(3):
            P = E.random_point()
            if not (order * P).is_zero():
                ok = False; break
        if ok:
            chosen = a2; break
    recs.append({"j_int": to_int(j), "b_int": to_int(b), "a2": None if chosen is None else int(chosen),
                 "order_matches_K0": chosen is not None, "orbit_rep_j": frob_orbit_rep(j)})

orbits = sorted(set(r["orbit_rep_j"] for r in recs))
for r in recs:
    r["orbit"] = f"O{orbits.index(r['orbit_rep_j']):02d}"
doc = {"n": N, "ell": ELL, "discriminant": D, "trace_K0": int(t), "order_K0": int(order), "conductor": int(f),
       "modulus_int": to_int(F.gen()**0 * 0 + F(0)) if False else int(sum(int(c) << i for i, c in enumerate(mod.list()))),
       "modulus": str(mod), "class_polynomial_degree": H.degree(), "roots_in_Fq": len(roots),
       "frobenius_orbits": len(orbits), "curves": recs}
json.dump(doc, open(OUT, "w"), indent=1)
print(json.dumps({k: v for k, v in doc.items() if k != "curves"}, indent=1))
print(f"curves with #E = #K_0: {sum(r['order_matches_K0'] for r in recs)} of {len(recs)}")
