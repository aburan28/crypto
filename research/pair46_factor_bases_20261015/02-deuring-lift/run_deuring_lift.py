#!/usr/bin/env python3
"""02-deuring-lift: canonical (Deuring) lift of the pair-46 crater, rank obstruction
and Xedni-style relation search.  All arithmetic is done by PARI/GP (cypari2);
Python only orchestrates and writes results.json.  See PROTOCOL.md."""
import json, math, random, sys, time
import cypari2

pari = cypari2.Pari()
pari.allocatemem(600 * 10**6)
gp = pari
gp("default(realprecision,80)")  # cypari2 default is 15 digits, too coarse for the frozen 1e-20 Gram tolerance
T0 = time.time()
SEED = 20261015
BUDGET = 38 * 60
OUT = sys.argv[1] if len(sys.argv) > 1 else "results.json"
P46 = 5947698030651343
N46 = 5947697996161643
CRATER = (5757116383569871, 3560173328417442)

def log(*a):
    print(f"[{time.time()-T0:7.1f}s]", *a, flush=True)

gp("E=ellinit([0,0,1,-38,90]); A=-27*E.c4; B=-54*E.c6;")  # short model y^2=x^3+Ax+B
A = int(gp("A")); B = int(gp("B"))
res = {"protocol": "PROTOCOL.md", "seed": SEED, "curve_Q": "361a1: y^2+y=x^3-38x+90",
       "short_model": [A, B]}

# ---------------- exact: canonical lift over Q ----------------
ex = {}
ex["j"] = str(gp("E.j")); ex["conductor"] = int(gp("ellglobalred(E)[1]"))
ex["torsion"] = str(gp("elltors(E)"))
ex["ellrank"] = str(gp("ellrank(E)"))
ex["analytic_rank"] = str(gp("ellanalyticrank(E)"))
ex["generator_minimal_model"] = str(gp("ellrank(E)[4]"))
ex["regulator"] = str(gp("ellheight(E, ellrank(E)[4][1])"))
gp(f"p={P46}; ap=ellap(E,p);")
ex["a_p_pair46"] = int(gp("ap")); ex["E_mod_p_order"] = int(gp("p+1-ap"))
ex["order_equals_N"] = int(gp("p+1-ap")) == N46
ex["chi_-19(p)"] = int(gp("kronecker(-19,p)"))
# explicit isomorphism short model of E -> crater over F_p
gp(f"Ac=Mod({CRATER[0]},p); Bc=Mod({CRATER[1]},p); u2=(Bc/Mod(B,p))/(Ac/Mod(A,p));")
ex["u2_is_square"] = int(gp("issquare(u2)"))
gp("u=sqrt(u2); Es=ellinit([A,B]); Ecr=ellinit([Ac,Bc]);")
ex["iso_check_A"] = bool(int(gp("u^4*Mod(A,p)==Ac"))); ex["iso_check_B"] = bool(int(gp("u^6*Mod(B,p)==Bc")))
gp("Gm=ellrank(E)[4][1]; Gs=ellchangepoint(Gm, [1,0,-1/2,0]);")  # minimal -> short? verify below
# map minimal-model point to short model y^2=x^3+Ax+B: x_s=36x+3b2, y_s=108(2y+a1x+a3)
gp("Gs=[36*Gm[1]+3*E.b2, 108*(2*Gm[2]+E.a1*Gm[1]+E.a3)];")
ex["G_short"] = str(gp("Gs")); ex["G_short_on_curve"] = bool(int(gp("ellisoncurve(Es,Gs)")))
gp("Gc=[u^2*Mod(Gs[1],p), u^3*Mod(Gs[2],p)];")
ex["G_mod_p_on_crater"] = bool(int(gp("ellisoncurve(Ecr,Gc)")))
ex["G_mod_p_order"] = str(gp("ellorder(Ecr,Gc," + str(N46) + ")"))
ex["image_of_EQ_is_whole_group"] = ex["G_mod_p_order"] == str(N46)
hG = float(gp("ellheight(E,Gm)"))
ex["canonical_height_G"] = hG
# a uniformly random target T = kG has lift kG with height k^2 h(G); E[k^2] for k uniform on [0,N)
ex["expected_height_random_target_lift"] = f"{hG/3:.4g} * N^2 = {hG*N46**2/3:.4e} (nats, PARI normalization)"
ex["log2_expected_height_random_target_lift"] = math.log2(hG * N46**2 / 3)
res["exact_pair46"] = ex
log("exact pair46 done", ex["ellrank"], ex["order_equals_N"], ex["image_of_EQ_is_whole_group"])

# ---------------- twist ranks ----------------
def is_fund(d):
    return int(gp(f"isfundamental({d})")) == 1
D = []
seen = set()
for ad in range(2, 201):
    for d in (ad, -ad):
        if d == -19 or not is_fund(d):
            continue
        partner = int(gp(f"coredisc({-19*d})"))
        key = frozenset((d, partner))
        if key in seen:
            continue
        seen.add(key); D.append(d)
D.sort(key=lambda d: (abs(d), -d))
twists = {}
for d in D:
    if time.time() - T0 > 0.45 * BUDGET:
        log("twist budget reached at d", d); break
    gp(f"Ed=ellinit([A*({d})^2,B*({d})^3]);")
    t1 = time.time()
    r = gp("iferr(alarm(60, ellrank(Ed)), e, \"ERR\")")
    rec = {"d": d}
    if str(r) == "ERR":
        rec.update(rank="undetermined", note="ellrank timeout 60s")
    else:
        lo, hi = int(r[0]), int(r[1])
        rec["ellrank_bounds"] = [lo, hi]
        if lo != hi:
            ar = gp("iferr(alarm(60, ellanalyticrank(Ed)[1]), e, -1)")
            ar = int(ar)
            rec["analytic_rank"] = ar
            rec["rank"] = ar if ar >= 0 and ar <= hi and ar >= lo else "undetermined"
        else:
            rec["rank"] = lo
        pts = r[3]
        rec["points"] = [[str(P[0]), str(P[1])] for P in pts][: max(lo, 0)] if lo > 0 else []
    rec["j_check"] = str(gp("Ed.j")) == "-884736"
    rec["sec"] = round(time.time() - t1, 2)
    twists[d] = rec
res["twists"] = list(twists.values())
log("twists done", len(twists), "with rank>=1:", sum(1 for v in twists.values() if isinstance(v.get("rank"), int) and v["rank"] >= 1))
ranks_pos = [d for d in twists if isinstance(twists[d].get("rank"), int) and twists[d]["rank"] >= 1 and twists[d]["points"]]
res["exact_pair46"]["twists_chi_plus_at_pair46"] = [
    {"d": d, "rank": twists[d].get("rank")} for d in twists if int(gp(f"kronecker({d},{P46})")) == 1]

# ---------------- toy population ----------------
def toy_primes(b, k=8):
    out = []; q = 2**b
    while len(out) < k:
        q = int(gp(f"nextprime({q}+1)"))
        if int(gp(f"kronecker(-19,{q})")) != 1 or q == 19:
            continue
        n = q + 1 - int(gp(f"ellap(E,{q})"))
        if int(gp(f"isprime({n})")):
            out.append((q, n))
    return out

def red_twist_point(d, P, q):
    """Reduce point P on E_d: y^2=x^3+A d^2 x+B d^3 to E short mod q via (x/d, y/(d s)), s^2=d."""
    gp(f"s=sqrt(Mod({d},{q})); Pq=[Mod({P[0]},{q})/{d}, Mod({P[1]},{q})/({d}*s)];")
    return "Pq"

def kernel_rho(logs, n):
    R = len(logs)
    gp(f"l={list(logs)}; n={n};")
    gp("i0=0; for(i=1,#l, if(l[i]%n, i0=i; break));")
    # basis columns: n e_{i0}; e_i - (l_i/l_{i0}) e_{i0}
    gp("R=#l; M=matrix(R,R); M[i0,i0]=n; inv=lift(1/Mod(l[i0],n)); "
       "for(i=1,R, if(i!=i0, M[i,i]=1; M[i0,i]=-lift(Mod(l[i]*inv,n)))); "
       "Mr=M*qflll(M); Gm2=Mr~*Mr; nm=vecmin(vector(R,i,Gm2[i,i])); "
       "Q=qfminim(Gm2,nm,,0); mn=nm; V=Q[3]; for(j=1,#V, v=Mr*V[,j]; mn=min(mn, v~*v));")
    lam = math.sqrt(int(gp("mn")))
    gh = math.sqrt(R / (2 * math.pi * math.e)) * n ** (1.0 / R)
    # check: shortest vector found is in kernel
    return lam / gh, lam

def coverage(logs, n, B=2):
    s = {0}
    for l in logs:
        s = {(x + c * l) % n for x in s for c in range(-B, B + 1)}
    return len(s) / n

def gram_rank(Estr, pts):
    """rank of canonical-height Gram matrix of pts on curve Estr (PARI)."""
    gp(f"Ex={Estr}; PP={pts}; H=ellheightmatrix(Ex,PP); ev=qfjacobi(H)[1];")
    ev = [float(x) for x in gp("ev")]
    mx = max(abs(x) for x in ev)
    return sum(1 for x in ev if abs(x) > 1e-20 * mx and abs(x) > 1e-30), ev

XCOLS = ["k6", "k4", "k2", "k3", "k1"]
def xedni_trial(q, pts, extra=()):
    r = len(pts)
    rows, rhs = [], []
    for (x, y) in pts:
        m = (y * y - x**3 - A * x - B)
        assert m % q == 0
        m //= q
        full = {"k6": -1, "k4": -x, "k2": -x * x, "k3": y, "k1": x * y}
        rows.append([full[c] for c in XCOLS[:r]]); rhs.append(-m)
    gp(f"Mx=matrix({r},{r},i,j,{rows}[i][j]); bx={rhs}~;")
    if int(gp("matrank(Mx)")) < r:
        return {"r": r, "status": "singular_incidence"}
    gp("kx=matsolve(Mx,bx);")
    k = {c: gp(f"kx[{i+1}]") for i, c in enumerate(XCOLS[:r])}
    for c in XCOLS[r:]:
        k[c] = 0
    gp(f"Ex2=ellinit([{q}*({k['k1']}),{q}*({k['k2']}),{q}*({k['k3']}),{A}+{q}*({k['k4']}),{B}+{q}*({k['k6']})]);")
    allp = list(pts) + list(extra)
    ps = "[" + ",".join(f"[{x},{y}]" for x, y in allp) + "]"
    ok = all(int(gp(f"ellisoncurve(Ex2,[{x},{y}])")) for x, y in allp)
    pint = int(gp(f"denominator(kx) % {q} != 0"))
    rk, ev = gram_rank("Ex2", ps)
    return {"r": len(allp), "on_curve": ok, "p_integral": bool(pint), "gram_rank": rk,
            "lifted_relations": len(allp) - rk, "min_eig": min(ev),
            "log10_disc": float(gp("log(abs(Ex2.disc))/log(10)"))}

toy = []
cut = []
for b in (16, 20, 24):
    for (q, n) in toy_primes(b):
        if time.time() - T0 > BUDGET:
            cut.append(q); continue
        rng = random.Random(SEED ^ q)
        rec = {"p": q, "bits": b, "N": n}
        gp(f"Eq=ellinit([A,B],{q}); Ng={n};")
        gp(f"Gq=[Mod(Gs[1],{q}),Mod(Gs[2],{q})];")
        gp("setrand(" + str((SEED ^ q) % 2**31) + ");")
        gp("g0=random(Eq); while(g0==[0], g0=random(Eq));")
        def L(ptexpr):
            return int(gp(f"elllog(Eq,{ptexpr},g0,Ng)"))
        logG = L("Gq")
        rec["A1_G_mod_p_nonzero"] = logG % n != 0
        # A2: canonical twist family
        ds = [d for d in ranks_pos if int(gp(f"kronecker({d},{q})")) == 1 and q % abs(d) != 0][:6]
        rec["A2_twists"] = ds
        if len(ds) == 6:
            logs = []
            for d in ds:
                P = twists[d]["points"][0]
                red_twist_point(d, P, q)
                assert int(gp("ellisoncurve(Eq,Pq)")), (q, d, P, str(gp("Pq")), str(gp("Eq[4..5]")))
                logs.append(L("Pq"))
            rho, lam = kernel_rho(logs, n)
            rec["A2_rho"] = rho; rec["A2_lambda1"] = lam
            rec["A2_cov"] = coverage(logs, n)
            # A2 lifted relations: one generator per distinct twist -> block-diagonal Gram
            hs = []
            for d in ds:
                P = twists[d]["points"][0]
                hs.append(float(gp(f"ellheight(ellinit([A*({d})^2,B*({d})^3]),[{P[0]},{P[1]}])")))
            rec["A2_heights"] = hs
            rec["A2_lifted_relations"] = 6 - sum(1 for h in hs if h > 1e-12)
            # PC: G,2G,3G + first 3 twists
            gp("G2=elladd(Es,Gs,Gs); G3=elladd(Es,G2,Gs);")
            rk, _ = gram_rank("Es", "[Gs,G2,G3]")
            rec["PC_lifted_relations"] = 6 - (rk + 3)
            pl = [logG, L("[Mod(G2[1],%d),Mod(G2[2],%d)]" % (q, q)), L("[Mod(G3[1],%d),Mod(G3[2],%d)]" % (q, q))] + logs[:3]
            rec["PC_rho"], rec["PC_lambda1"] = kernel_rho(pl, n)
            # NC: 6 random points
            nl = [rng.randrange(1, n) for _ in range(6)]  # uniform random points = uniform logs (prime order)
            rec["NC_rho"], rec["NC_lambda1"] = kernel_rho(nl, n)
            rec["NC_cov"] = coverage(nl, n)
        # A3: generic Xedni lift
        xed = []
        for r in (2, 3, 4, 5):
            ptsF = []
            for _ in range(r):
                gp("z=random(Eq); while(z==[0], z=random(Eq));")
                ptsF.append((int(gp("lift(z[1])")), int(gp("lift(z[2])"))))
            xed.append(xedni_trial(q, ptsF))
        rec["A3_xedni"] = xed
        gp("z=random(Eq); while(z==[0], z=random(Eq)); z2=random(Eq); while(z2==[0], z2=random(Eq));")
        x1, y1 = int(gp("lift(z[1])")), int(gp("lift(z[2])"))
        x2, y2 = int(gp("lift(z2[1])")), int(gp("lift(z2[2])"))
        rec["PC_xedni"] = xedni_trial(q, [(x1, y1), (x2, y2)], extra=[(x1, -y1)])
        toy.append(rec)
        log("p", q, "A2_rho", rec.get("A2_rho"), "NC_rho", rec.get("NC_rho"), "PC_rho", rec.get("PC_rho"))
res["toy"] = toy
res["cut_primes"] = cut
json.dump(res, open(OUT, "w"), indent=1, default=str)
log("wrote", OUT)
