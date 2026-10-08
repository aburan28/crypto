from sympy import factorint, isprime, legendre_symbol, n_order, divisors
def lucas(n):
    V=[2,-1]; U=[0,1]
    for k in range(1,n):
        V.append(-V[k]-2*V[k-1]); U.append(-U[k]-2*U[k-1])
    return V[n], abs(U[n])
def chi(ell): return 1 if ell==2 else int(legendre_symbol((-7)%ell, ell))
def profile(n, a_curve=0, verbose=True):
    q=2**n; t,f=lucas(n)
    if a_curve==1: t=-t          # K_1 is the twist for odd n: trace flips
    NE=q+1-t; fac=factorint(f)
    rows=[]
    for ell,e in sorted(fac.items()):
        c=chi(ell); a=(t*pow(2,-1,ell))%ell if ell!=2 else None
        o=int(n_order(a,ell)) if a not in (None,0) else None
        rows.append((int(ell),int(e),c,int(ell-c),a,o))
    tot=1
    for d in divisors(f):
        h=1
        for ell,e in factorint(d).items(): h*=ell**(e-1)*(ell-chi(ell))
        tot+=h if d>1 else 0
    return t,f,fac,rows,NE,tot
t,f,fac,rows,NE,tot=profile(131)
print("n=131 a=0: f =",f,"=",dict(fac))
print("ell_big prime:",isprime(146505763881528721))
for ell,e,c,lvl1,a,o in rows:
    print(f"  ell={ell} depth={e} chi={c:+d} horiz_at_surface={1+c} curves_level1={lvl1} pi_mod_ell={a} kernel_field_degree=ord(a)={o}")
print("  total curves in class (incl. K_0):",tot," ~2^%.1f"%(tot.bit_length()-1))
print("  fraction with End=O_K: 1/%d"%tot)
print("\nToy ladder (a=0 unless noted): n | f | factorization | per-prime (depth, chi, kernel-degree)")
for n in [23,41,53,61,71,73,83,97]:
    t,f,fac,rows,NE,tot=profile(n)
    print(f"n={n}: f={f} fac={dict(fac)} rows={[(r[0],r[1],r[2],r[5]) for r in rows]} total_curves~2^{tot.bit_length()-1}")
