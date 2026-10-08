#!/usr/bin/env sage -python
"""E4(2)-(3): the 263-torsion of E0 over F_{2^262}, its Weil/Tate pairings, and the gamma-weight of x(T).

F = F_2(gamma) with gamma a primitive 263rd root of unity (modulus = a degree-131 factor of x^263-1),
K = F_{2^262} with an explicit embedding.  E0: y^2+xy=x^3+1.  pi = -1 on E0[263], so E0[263] is rational
over K and x(T) in F for every T in E0[263].
"""
import json, time, random
from sage.all import GF, PolynomialRing, EllipticCurve, ZZ, set_random_seed
set_random_seed(7); t0=time.time()
R = PolynomialRing(GF(2), "x"); x = R.gen()
facs = [g for g, e in (x**263 - 1).factor() if g.degree() == 131]
assert len(facs) == 2, [g.degree() for g,_ in (x**263-1).factor()]
F = GF(2**131, "g", modulus=facs[0]); g = F.gen()
assert g**263 == 1 and g != 1
K = GF(2**262, "w")
rt = F.modulus().change_ring(K).roots(multiplicities=False)[0]
emb = F.hom([rt], K)
E = EllipticCurve(K, [1, 0, 0, 0, 1])
q = ZZ(2)**131; t = ZZ(-22283658519494248867)
NK = q**2 + 1 - (t*t - 2*q)        # #E(F_{q^2})
v = 0; m = NK
while m % 263 == 0: m //= 263; v += 1
assert v >= 2, v
cof = NK // 263**v
def point_of_order_263():
    while True:
        S = cof * E.random_point()
        while S != E(0) and 263*S != E(0): S = 263*S
        if S != E(0): return S
T1 = point_of_order_263()
while True:
    T2 = point_of_order_263()
    w = T1.weil_pairing(T2, 263)
    if w != 1: break
# (2) pairings
w_in_F = (w**q == w); w_order_263 = (w**263 == 1 and w != 1)
r = ZZ(680564733841876926932320129493409985129)
Gq = EllipticCurve(F, [1,0,0,0,1]).random_point()
G = E(emb(Gq[0]), emb(Gq[1]))
G = (NK // r) * G          # project into the order-r part (over K the r-part is r x r? no: r || #E(K); use cofactor)
assert r*G == E(0) and G != E(0)
# reduced Tate pairing with k = 1 must be evaluated on a divisor (Q+S)-(S), not at the point Q:
# with k = 1 the factor from evaluating the Miller function at O is not killed by the final exponentiation.
expo = (ZZ(2)**262 - 1) // 263
def tate_div(P, Q):
    while True:
        S = E.random_point()
        if S in (E(0), P, -P, Q, -Q, Q+P, Q-P): continue
        try:
            f = P._miller_(Q + S, ZZ(263)) / P._miller_(S, ZZ(263))
            return f**expo
        except ZeroDivisionError:
            continue
tate_T1_G = tate_div(T1, G)
tate_T1_R = tate_div(T1, E.random_point())
tate_T1_T2 = tate_div(T1, T2)
tate_pointwise_G = T1.tate_pairing(G, 263, 1, q=2**262)   # the naive point evaluation, kept for the record
# (3) gamma-weight of x(T) for 40 torsion points, versus random F elements
def coords_weight(c): return sum(1 for a in c.polynomial().coefficients(sparse=False) if a)
def min_conj_weight(xK):
    mp = xK.minimal_polynomial()
    roots = mp.roots(F, multiplicities=False)
    return min(coords_weight(c) for c in roots), mp.degree()
ws=[]; xs_in_F=0
for i in range(40):
    a, b = random.randrange(1,263), random.randrange(0,263)
    T = a*T1 + b*T2
    if T == E(0): continue
    xT = T[0]; xs_in_F += int(xT**q == xT)
    wmin, deg = min_conj_weight(xT); ws.append(wmin)
# same statistic for random elements: minimum polynomial-basis weight over the 131 conjugates
def min_conj_weight_F(c):
    best=10**9; y=c
    for _ in range(131):
        best=min(best, coords_weight(y)); y=y**2
    return best
rand_ws=[min_conj_weight_F(F.random_element()) for _ in range(60)]
rec={"weil_pairing_in_Fq": bool(w_in_F), "weil_pairing_order_263": bool(w_order_263),
     "tate_divisor_T1_G_trivial": bool(tate_T1_G == 1), "tate_divisor_T1_randompoint_trivial": bool(tate_T1_R == 1),
     "tate_divisor_T1_T2_order_263": bool(tate_T1_T2 != 1 and tate_T1_T2**263 == 1),
     "tate_pointwise_T1_G_trivial(k=1 shortcut, invalid)": bool(tate_pointwise_G == 1),
     "x_of_torsion_in_Fq": f"{xs_in_F}/{len(ws)}",
     "gamma_weight_torsion_min_mean_max": [min(ws), sum(ws)/len(ws), max(ws)],
     "gamma_weight_random_min_mean_max": [min(rand_ws), sum(rand_ws)/len(rand_ws), max(rand_ws)],
     "wall_s": round(time.time()-t0,1), "status": "PASS" if (w_in_F and w_order_263 and tate_T1_G==1 and xs_in_F==len(ws)) else "CHECK"}
print(json.dumps(rec, indent=2))
