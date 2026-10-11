"""Isogeny-fiber probe: rational vs geometric fibers, Frobenius orbits on
fibers, the pullback factorisation of summation polynomials, and a verified
relation lifted across an isogeny.  Pure Python + sympy, toy prime fields.

Run:  python3 -I fiber_probe.py            (prints a report, writes JSON)
"""
import json, random, sys
from itertools import product

# ---------------------------------------------------------------- F_p^k
class Ext:
    """F_{p^k} = F_p[z]/(m(z)), elements as tuples of length k."""
    def __init__(self, p, k, modulus):
        self.p, self.k, self.m = p, k, modulus  # monic, m[k]==1
    def zero(self): return (0,)*self.k
    def one(self):  return (1,)+(0,)*(self.k-1)
    def emb(self, a): return (a % self.p,)+(0,)*(self.k-1)
    def add(self, a, b): return tuple((x+y) % self.p for x, y in zip(a, b))
    def sub(self, a, b): return tuple((x-y) % self.p for x, y in zip(a, b))
    def neg(self, a): return tuple((-x) % self.p for x in a)
    def mul(self, a, b):
        p, k = self.p, self.k
        c = [0]*(2*k-1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    c[i+j] = (c[i+j] + x*y) % p
        for d in range(2*k-2, k-1, -1):
            if c[d]:
                t = c[d]
                for j in range(k):
                    c[d-k+j] = (c[d-k+j] - t*self.m[j]) % p
        return tuple(c[:k])
    def pow(self, a, e):
        r, b = self.one(), a
        while e:
            if e & 1: r = self.mul(r, b)
            b = self.mul(b, b); e >>= 1
        return r
    def inv(self, a):
        return self.pow(a, self.p**self.k - 2)
    def frob(self, a): return self.pow(a, self.p)
    def is_zero(self, a): return all(x == 0 for x in a)
    def sqrt(self, a):
        q = self.p**self.k
        if self.is_zero(a): return a
        if self.pow(a, (q-1)//2) != self.one(): return None
        assert q % 4 == 3, "probe uses q = 3 mod 4 only"
        return self.pow(a, (q+1)//4)

def irreducible(p, k, rng):
    """Random monic irreducible of degree k over F_p (brute force, k<=3)."""
    while True:
        m = [rng.randrange(p) for _ in range(k)] + [1]
        if k == 1: return m
        if any(sum(c*pow(x, i, p) for i, c in enumerate(m)) % p == 0 for x in range(p)):
            continue
        if k == 2: return m
        # k == 3: no linear factor suffices
        return m

# ---------------------------------------------------------------- curves
class Curve:
    def __init__(self, F, a, b):
        self.F, self.a, self.b = F, a, b
    def on(self, P):
        if P is None: return True
        x, y = P; F = self.F
        lhs = F.mul(y, y)
        rhs = F.add(F.add(F.mul(F.mul(x, x), x), F.mul(self.a, x)), self.b)
        return lhs == rhs
    def neg(self, P): return None if P is None else (P[0], self.F.neg(P[1]))
    def add(self, P, Q):
        F = self.F
        if P is None: return Q
        if Q is None: return P
        if P[0] == Q[0]:
            if F.is_zero(F.add(P[1], Q[1])): return None
            num = F.add(F.mul(F.emb(3), F.mul(P[0], P[0])), self.a)
            lam = F.mul(num, F.inv(F.mul(F.emb(2), P[1])))
        else:
            lam = F.mul(F.sub(Q[1], P[1]), F.inv(F.sub(Q[0], P[0])))
        x3 = F.sub(F.sub(F.mul(lam, lam), P[0]), Q[0])
        y3 = F.sub(F.mul(lam, F.sub(P[0], x3)), P[1])
        return (x3, y3)
    def mul(self, n, P):
        R = None
        while n:
            if n & 1: R = self.add(R, P)
            P = self.add(P, P); n >>= 1
        return R
    def frob(self, P): return None if P is None else (self.F.frob(P[0]), self.F.frob(P[1]))
    def points_prime(self):
        """All F_p-rational points (F is the prime field, k=1)."""
        F, p = self.F, self.F.p
        pts = [None]
        for x in range(p):
            X = F.emb(x)
            rhs = F.add(F.add(F.mul(F.mul(X, X), X), F.mul(self.a, X)), self.b)
            s = F.sqrt(rhs)
            if s is None: continue
            pts.append((X, s))
            if not F.is_zero(s): pts.append((X, F.neg(s)))
        return pts

# ---------------------------------------------------------------- Velu
def velu(E, kernel_reps):
    """kernel_reps: representatives Q of (ker \\ {O}) / +-, as points over E.F.
    Returns (E2, phi) with phi evaluating the isogeny (Velu/Kohel form)."""
    F = E.F
    t = w = F.zero()
    data = []
    for Q in kernel_reps:
        xq, yq = Q
        gx = F.add(F.mul(F.emb(3), F.mul(xq, xq)), E.a)
        gy = F.mul(F.emb(-2), yq)
        two_tors = F.is_zero(yq)
        vq = gx if two_tors else F.mul(F.emb(2), gx)
        uq = F.mul(gy, gy)
        t = F.add(t, vq)
        w = F.add(w, F.add(uq, F.mul(xq, vq)))
        data.append((xq, yq, gx, gy, vq, uq))
    A = F.sub(E.a, F.mul(F.emb(5), t))
    B = F.sub(E.b, F.mul(F.emb(7), w))
    E2 = Curve(F, A, B)
    def phi(P):
        if P is None: return None
        x, y = P
        X, Y = x, y
        for (xq, yq, gx, gy, vq, uq) in data:
            d = F.sub(x, xq)
            if F.is_zero(d): return None           # P in kernel
            di = F.inv(d); di2 = F.mul(di, di); di3 = F.mul(di2, di)
            X = F.add(X, F.add(F.mul(vq, di), F.mul(uq, di2)))
            term = F.add(F.mul(uq, F.mul(F.mul(F.emb(2), y), di3)),
                         F.mul(vq, F.mul(F.sub(y, yq), di2)))
            term = F.sub(term, F.mul(F.mul(gx, gy), di2))
            Y = F.sub(Y, term)
        return (X, Y)
    return E2, phi

# ---------------------------------------------------------------- helpers
def tonelli(F, a):
    """Square root in F_{q} by Tonelli-Shanks (any q odd)."""
    q = F.p**F.k
    if F.is_zero(a): return a
    if F.pow(a, (q-1)//2) != F.one(): return None
    if q % 4 == 3: return F.pow(a, (q+1)//4)
    s, Q = 0, q-1
    while Q % 2 == 0: Q //= 2; s += 1
    rng = random.Random(1)
    while True:
        z = tuple(rng.randrange(F.p) for _ in range(F.k))
        if not F.is_zero(z) and F.pow(z, (q-1)//2) != F.one(): break
    M, c, t, R = s, F.pow(z, Q), F.pow(a, Q), F.pow(a, (Q+1)//2)
    while t != F.one():
        i, tt = 0, t
        while tt != F.one(): tt = F.mul(tt, tt); i += 1
        b = c
        for _ in range(M-i-1): b = F.mul(b, b)
        M, c, t, R = i, F.mul(b, b), F.mul(t, F.mul(b, b)), F.mul(R, b)
    return R

def factor_int(n):
    f, d = {}, 2
    while d*d <= n:
        while n % d == 0: f[d] = f.get(d, 0)+1; n //= d
        d += 1
    if n > 1: f[n] = f.get(n, 0)+1
    return f

def lift_curve(E, G):
    """Re-instantiate E (over F_p) over the extension G."""
    return Curve(G, G.emb(E.a[0]), G.emb(E.b[0]))

def embed_pt(P, G):
    return None if P is None else (G.emb(P[0][0]), G.emb(P[1][0]))

def poly_mulmod(F, A, B):
    """Polynomials over F_p as lists of ints, low degree first."""
    p = F.p
    c = [0]*(len(A)+len(B)-1)
    for i, x in enumerate(A):
        for j, y in enumerate(B):
            c[i+j] = (c[i+j] + x*y) % p
    return c

def fiber_report(E1, E2, phi, ker_pts_prime, pts1, pts2, label):
    """Rational fibers of every Q in E2(F_p); Lang/H^1 bookkeeping."""
    img = {}
    for P in pts1:
        Q = phi(P)
        key = 'O' if Q is None else Q
        img.setdefault(key, []).append(P)
    nk = len(ker_pts_prime)             # #ker(F_p), including O
    sizes = {}
    for Q in pts2:
        key = 'O' if Q is None else Q
        n = len(img.get(key, []))
        sizes[n] = sizes.get(n, 0)+1
    N1, N2 = len(pts1), len(pts2)
    assert N1 == N2, (N1, N2)
    n_img = len(img)
    out = dict(label=label, N=N1, ker_rational=nk, image_size=n_img,
               index=N2//n_img, rational_fiber_size_histogram=sizes)
    assert set(sizes) <= {0, nk}, sizes
    assert N2 == n_img*nk
    return out, img

# ---------------------------------------------------------------- main
def main(p=211, seed=3):
    rng = random.Random(seed)
    F = Ext(p, 1, [0, 1])
    res = {'p': p}
    # ---- find a curve with a rational 2-torsion point AND a 3-torsion point
    #      whose x is rational (so a rational 3-isogeny); classify kernel rationality.
    found = None
    for a in range(1, p):
        for b in range(1, p):
            E = Curve(F, F.emb(a), F.emb(b))
            d = (4*a**3 + 27*b**2) % p
            if d == 0: continue
            pts = E.points_prime()
            N = len(pts)
            fac = factor_int(N)
            r = max(fac)
            if r < 11 or N % 3: continue
            two = [P for P in pts if P is not None and F.is_zero(P[1])]
            # roots of psi_3 = 3x^4 + 6a x^2 + 12 b x - a^2
            psi3 = [x for x in range(p) if (3*x**4 + 6*a*x*x + 12*b*x - a*a) % p == 0]
            if not two or not psi3: continue
            if len(two) == 1 and len(psi3) >= 1:
                found = (a, b, E, pts, N, fac, r, two, psi3); break
        if found: break
    a, b, E1, pts1, N, fac, r, two, psi3 = found
    res['E1'] = dict(a=a, b=b, N=N, N_factored=fac, r=r, t=p+1-N)
    # ---- 2-isogeny with rational kernel
    T2 = two[0]
    E2a, phi2 = velu(E1, [T2])
    pts2a = E2a.points_prime()
    rep2, img2 = fiber_report(E1, E2a, phi2, [None, T2], pts1, pts2a, 'ell=2 rational kernel')
    res['iso2'] = rep2
    # homomorphism check
    for _ in range(50):
        P, Q = rng.choice(pts1), rng.choice(pts1)
        assert phi2(E1.add(P, Q)) == E2a.add(phi2(P), phi2(Q))
        assert E2a.on(phi2(P))
    # ---- 3-isogenies: one per psi_3 root; classify whether kernel points are rational
    res['iso3'] = []
    for x0 in psi3:
        X0 = F.emb(x0)
        rhs = F.add(F.add(F.mul(F.mul(X0, X0), X0), F.mul(E1.a, X0)), E1.b)
        y0 = tonelli(F, rhs)
        if y0 is not None:                      # rational kernel
            T3 = (X0, y0)
            E2b, phi3 = velu(E1, [T3])
            pts2b = E2b.points_prime()
            kp = [None, T3, E1.neg(T3)]
            rep, img = fiber_report(E1, E2b, phi3, kp, pts1, pts2b, 'ell=3 rational kernel')
            for _ in range(50):
                P, Q = rng.choice(pts1), rng.choice(pts1)
                assert phi3(E1.add(P, Q)) == E2b.add(phi3(P), phi3(Q))
            # geometric fiber of a Q outside the image: lives over F_{p^3}
            Qout = next(Q for Q in pts2b if Q is not None and Q not in img)
            # numerator of X(x) - x_Q  (X = x + v/(x-x0) + u/(x-x0)^2) : cubic in x
            gx = (3*x0*x0 + a) % p; v = 2*gx % p; u = 4*(x0**3 + a*x0 + b) % p
            xq = Qout[0][0]
            # (x - xq)(x-x0)^2 + v (x-x0) + u
            lin = [(-x0) % p, 1]
            sq = poly_mulmod(F, lin, lin)
            cub = poly_mulmod(F, [(-xq) % p, 1], sq)
            cub[0] = (cub[0] + (-v*x0) % p + u) % p; cub[1] = (cub[1] + v) % p
            has_root = any(sum(c*pow(x, i, p) for i, c in enumerate(cub)) % p == 0 for x in range(p))
            rep['cubic_over_Fp_has_root'] = has_root
            if not has_root:
                G = Ext(p, 3, cub)             # F_{p^3} = F_p[z]/(cubic): z is x(P)
                E1G, E2G = lift_curve(E1, G), lift_curve(E2b, G)
                _, phiG = velu(E1G, [embed_pt(T3, G)])
                z = (0, 1, 0)
                fz = G.add(G.add(G.mul(G.mul(z, z), z), G.mul(E1G.a, z)), E1G.b)
                yz = tonelli(G, fz)
                rep['y_in_Fp3'] = yz is not None
                if yz is not None:
                    P = (z, yz)
                    if phiG(P) != embed_pt(Qout, G): P = E1G.neg(P)   # other sign of y
                    assert E1G.on(P) and phiG(P) == embed_pt(Qout, G)
                    orbit = [P]
                    while True:
                        nxt = E1G.frob(orbit[-1])
                        if nxt == P: break
                        orbit.append(nxt)
                    TG = embed_pt(T3, G)
                    c = E1G.add(E1G.frob(P), E1G.neg(P))     # cocycle pi(P)-P in ker
                    rep['fiber_outside_image'] = dict(
                        frobenius_orbit_length=len(orbit),
                        cocycle_is_T=(c == TG), cocycle_is_minus_T=(c == E1G.neg(TG)),
                        fiber_equals_orbit=sorted(map(str, orbit)) == sorted(map(str, [P, E1G.add(P, TG), E1G.add(P, E1G.neg(TG))])))
            res['iso3'].append(rep)
            res['_E2b'] = (E2b, phi3, pts2b, T3)
        else:                                   # kernel {O, +-T}, T in E1(F_{p^2}), pi(T) = -T
            c = rhs[0]
            G = Ext(p, 2, [(-c) % p, 0, 1])     # z^2 = c
            E1G = lift_curve(E1, G)
            T3 = (G.emb(x0), (0, 1))
            assert E1G.on(T3)
            assert E1G.frob(T3) == E1G.neg(T3)
            E2G, phiG = velu(E1G, [T3])
            # the isogeny is F_p-rational: coefficients lie in F_p
            assert E2G.a[1] == 0 and E2G.b[1] == 0
            E2c = Curve(F, F.emb(E2G.a[0]), F.emb(E2G.b[0]))
            pts2c = E2c.points_prime()
            def phi3c(P):
                Q = phiG(embed_pt(P, G))
                if Q is None: return None
                assert Q[0][1] == 0 and Q[1][1] == 0
                return (F.emb(Q[0][0]), F.emb(Q[1][0]))
            rep, img = fiber_report(E1, E2c, phi3c, [None], pts1, pts2c, 'ell=3 non-rational kernel (pi T = -T)')
            # geometric fiber over F_{p^2}: {P, P+T, P-T}, Frobenius swaps P+-T
            P = rng.choice([Q for Q in pts1 if Q is not None])
            PG = embed_pt(P, G)
            fib = [PG, E1G.add(PG, T3), E1G.add(PG, E1G.neg(T3))]
            orbits = {str(E1G.frob(x)) for x in fib}
            rep['sample_fiber'] = dict(rational_points=sum(1 for x in fib if x[0][1] == 0 and x[1][1] == 0),
                                       frobenius_swaps_P_plus_minus_T=(E1G.frob(fib[1]) == fib[2]))
            res['iso3'].append(rep)
    # ---- verified relation lifted across the rational-kernel 3-isogeny
    E2b, phi3, pts2b, T3 = res.pop('_E2b')
    h = N // r
    B = 12
    FB = [P for P in pts1 if P is not None and P[0][0] < B]
    R = rng.choice([P for P in pts1 if P is not None])
    Qt = phi3(R)
    rel = None
    for i in range(len(FB)):
        for j in range(i, len(FB)):
            for k in range(j, len(FB)):
                S = E2b.add(E2b.add(phi3(FB[i]), phi3(FB[j])), phi3(FB[k]))
                if S == Qt: rel = (FB[i], FB[j], FB[k]); break
            if rel: break
        if rel: break
    assert rel, "no 3-term relation found; enlarge B"
    S1 = E1.add(E1.add(rel[0], rel[1]), rel[2])
    diff = E1.add(S1, E1.neg(R))
    kerpts = [None, T3, E1.neg(T3)]
    assert diff in kerpts
    lhs, rhs_ = E1.mul(h, S1), E1.mul(h, R)
    assert lhs == rhs_
    res['relation'] = dict(factor_base_bound_x=B, factor_base_size=len(FB),
                           R=str(R), phi_R=str(Qt), P=[str(x) for x in rel],
                           sum_on_E1_minus_R=('O' if diff is None else str(diff)),
                           kernel_index=kerpts.index(diff), cofactor=h,
                           cofactor_projection_equal=(lhs == rhs_))
    return res, (E1, E2a, E2b, a, b, T2, T3)

if __name__ == '__main__':
    res, objs = main()
    print(json.dumps({k: v for k, v in res.items()}, indent=1, default=str))
    with open(sys.argv[1] if len(sys.argv) > 1 else 'fiber_probe.json', 'w') as f:
        json.dump(res, f, indent=1, default=str)
