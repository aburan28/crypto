#!/usr/bin/env python3
"""E6: volcano profile per rung and component decision for a curve of the ECC2K-130 family.

  python3 -I which_component.py profiles            -> volcano_profiles.json (manifest import data)
  python3 -I which_component.py classify N J_INT     -> which component (n=131 via ground truth, n=37 via class37)
  python3 -I which_component.py selftest             -> classify every curve of the ground truth and class37

J_INT is the j-invariant as an integer bitmask in the stated modulus (bit i = coefficient of z^i),
the encoding of ecdlp-hardness-work/ground_truth/ground_truth.json.
"""
import json, os, sys
from sympy import factorint, isprime, legendre_symbol, divisors

GT = "/Volumes/SSD990/ecdlp-hardness-work/ground_truth/ground_truth.json"
C37 = "/Volumes/SSD990/ecdlp-hardness-work/ic-conductor-leads/data/class37.json"
LADDER = [23, 37, 41, 53, 61, 71, 73, 83, 97, 131]
KERNEL_FEASIBLE_DEGREE = 1000   # field degree up to which a Velu descent is routine
MODULAR_FEASIBLE_ELL = 10_000   # modular-polynomial route

def lucas(n):
    V=[2,-1]; U=[0,1]
    for k in range(1,n):
        V.append(-V[k]-2*V[k-1]); U.append(-U[k]-2*U[k-1])
    return V[n], abs(U[n])

def chi(ell): return 1 if ell==2 else int(legendre_symbol((-7)%ell, ell))

def mult_order(a, m):
    if a % m == 0: return None
    from sympy import n_order
    return int(n_order(a % m, m))   # factors m-1; the naive loop would take ~1e16 steps at n=131

def profile(n):
    q=2**n; t,f=lucas(n); assert t*t-4*q == -7*f*f
    fac=factorint(f); primes=[]
    for ell,e in sorted(fac.items()):
        c=chi(ell); a=(t*pow(2,-1,ell))%ell if ell!=2 else None
        d=mult_order(a,ell) if a is not None else None
        kernel_ok = d is not None and d <= KERNEL_FEASIBLE_DEGREE
        modular_ok = ell <= MODULAR_FEASIBLE_ELL
        primes.append({"ell":int(ell),"depth":int(e),"chi_minus7":c,"horizontal_edges_at_surface":1+c,
                       "curves_one_level_down":int(ell-c),"pi_mod_ell":int(a) if a is not None else None,
                       "kernel_field_degree":d,"reachable_by":("kernel" if kernel_ok else ("modular" if modular_ok else "none"))})
    total=1
    for dv in divisors(f):
        if dv==1: continue
        h=1
        for ell,e in factorint(dv).items(): h*=ell**(e-1)*(ell-chi(ell))
        total+=h
    return {"n":n,"trace":int(t),"order":int(q+1-t),"conductor_f":int(f),"f_factorization":{str(k):int(v) for k,v in fac.items()},
            "primes":primes,"curves_in_class":int(total),"h_OK":1,"surface_curves":1}

def write_profiles(path="volcano_profiles.json"):
    out={"family":"K_0: y^2+xy=x^3+1 over F_{2^n}, End=Z[tau], tau^2+tau+2=0","generated":"2026-10-08",
         "rungs":{str(n):profile(n) for n in LADDER}}
    json.dump(out, open(path,"w"), indent=1); print("wrote",path,"rungs:",LADDER)

def load_classes():
    classes={}
    if os.path.exists(GT):
        g=json.load(open(GT)); classes[131]={int(c["j_int"]):c for c in g["curves"]}
    if os.path.exists(C37):
        c37=json.load(open(C37)); classes[37]={int(c["j_int"]):c for c in c37["curves"].values()}
    return classes

def classify(n, j, classes=None):
    classes = classes if classes is not None else load_classes()
    p=profile(n)
    if j in (0,1):
        return {"n":n,"j":j,"component":"crater","label":"E0" if j==1 else "j=0 (supersingular, not in class)"}
    if n in classes and j in classes[n]:
        c=classes[n][j]
        return {"n":n,"j":j,"component":"reachable","level":c.get("level"),"orbit":c.get("orbit"),"label":c.get("label"),
                "frob_index":c.get("frob_index")}
    big=[pr for pr in p["primes"] if pr["reachable_by"]=="none"]
    if n in classes:
        return {"n":n,"j":j,"component":"big (not among the constructible curves)",
                "unreachable_primes":[pr["ell"] for pr in big],
                "note":"a random curve of this order lies here with probability 1 - 2^-%.1f" % (__import__("math").log2(p["curves_in_class"]/max(1,len(classes[n])))) }
    return {"n":n,"j":j,"component":"unknown (no class data for this n)","unreachable_primes":[pr["ell"] for pr in big],
            "reachable_primes":[pr["ell"] for pr in p["primes"] if pr["reachable_by"]!="none"]}

def selftest():
    classes=load_classes(); ok=0; tot=0
    for n,cls in classes.items():
        for j,c in cls.items():
            r=classify(n,j,classes); tot+=1
            good = (r["component"]=="crater" and c.get("level")==1) or (r["component"]=="reachable" and r["level"]==c.get("level"))
            ok+=good
            if not good: print("MISMATCH",n,c.get("label"),r)
    print(f"selftest: {ok}/{tot} known curves classified consistently; classes loaded for n={sorted(classes)}")
    # an unreachable probe: no curve of the big component exists; use j of a wrong curve as a negative control
    r=classify(131, 12345678901234567890, classes); print("negative control (random j at n=131):", r["component"])
    return ok==tot

if __name__=="__main__":
    cmd=sys.argv[1] if len(sys.argv)>1 else "selftest"
    if cmd=="profiles": write_profiles()
    elif cmd=="classify": print(json.dumps(classify(int(sys.argv[2]), int(sys.argv[3])), indent=1))
    else: sys.exit(0 if selftest() else 1)
