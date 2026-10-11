# ECC2K-130: K_0: y^2+xy=x^3+1 over F_{2^131}; tau^2 + tau + 2 = 0 (a=0 => trace of tau = -1).
import sys
from sympy import factorint, isprime, legendre_symbol, n_order
n=131; q=2**n
# Lucas sequences for tau,taubar: V_k = tau^k + taubar^k, U_k=(tau^k - taubar^k)/(tau-taubar)
# tau+taubar = -1, tau*taubar = 2  => V_{k+1} = -V_k - 2V_{k-1}; U_{k+1} = -U_k - 2U_{k-1}
V=[2,-1]; U=[0,1]
for k in range(1,n):
    V.append(-V[k]-2*V[k-1]); U.append(-U[k]-2*U[k-1])
t=V[n]; f=abs(U[n])
NE=q+1-t
r=680564733841876926932320129493409985129
print("trace t =",t); print("#E =",NE, " = 4*r:",NE==4*r, " r prime:",isprime(r))
print("t^2-4q = -7 f^2 check:", t*t-4*q == -7*f*f)
print("conductor f of Z[pi] in O_K =",f, " bits:",f.bit_length())
fac=factorint(f); print("factorization of f:",fac)
print("\nper-prime volcano profile (K_0 at the surface, End=O_K, h(O_K)=1):")
print("ell | v_ell(f)=depth | (-7/ell) | surface horiz edges | curves at level 1 | a=pi mod ell | ord(a) = kernel field degree")
for ell,e in sorted(fac.items()):
    if ell==2:
        chi = 1  # 2 splits in Q(sqrt-7): -7 = 1 mod 8
    else:
        chi = legendre_symbol(-7 % ell, ell)
    horiz = 1+chi
    lvl1 = ell - chi
    a = (t * pow(2, -1, ell)) % ell if ell!=2 else None
    o = n_order(a, ell) if (a not in (None,0)) else None
    print(f"{ell:>22} | {e} | {chi:+d} | {horiz} | {lvl1} | {a} | {o}")
# total number of curves in the isogeny class = sum over orders Z+f'O_K, f'|f of h(Z+f'O_K)
from sympy import divisors
def h_order(fp):
    # h(Z+f'O_K) = h(O_K) * f' * prod_{ell|f'} (1 - chi(ell)/ell) / [O_K^*:O^*], units index 1 for D=-7 (units +-1)
    h=1
    for ell,e in factorint(fp).items():
        chi = 1 if ell==2 else legendre_symbol(-7%ell,ell)
        h *= ell**(e-1) * (ell-chi)
    return h if fp>1 else 1
tot=sum(h_order(d) for d in divisors(f))
print("\ncurves (F_q-isomorphism classes) in the isogeny class:",tot, " bits:",tot.bit_length())
print("fraction at the surface (K_0 alone): 1/%d" % tot)
# 2-adic: rational 2-power torsion of K_0 is cyclic Z/4; depth at 2:
print("v_2(f) =", fac.get(2,0))
