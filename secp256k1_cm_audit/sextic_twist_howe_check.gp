\\ sextic_twist_howe_check.gp
\\ Thread 18: Howe (H1)+(H2)+(H3) for all 15 pairs of 6 sextic twists of secp256k1.
\\ Note: loop vars named 'ii','jj','kk' to avoid clash with polynomial indeterminates i,j,k.

default(parisize, 256000000);

p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
n = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;
t0 = p + 1 - n;

print("=================================================================");
print("Thread 18: Howe-gluing for all 15 pairs of j=0 sextic twists");
print("=================================================================");
print("p mod 6 = ", p % 6, "  (=1 => 6 distinct sextic twists exist)");
print("");

\\ --- CM decomposition 4p = t0^2 + 3*bval^2 ---
bsq = (4*p - t0^2) \ 3;
bval = sqrtint(bsq);
aval = (t0 + bval) \ 2;
print("CM decomposition: a=", aval, "  b=", bval);
print("  2a-b = t0? ", 2*aval - bval == t0);
print("  a^2-ab+b^2 = p? ", aval^2 - aval*bval + bval^2 == p);
print("");

\\ --- 6 traces algebraically ---
\\ Trace of twist by zeta_6^kk equals:
\\   kk=0: 2a-b, kk=1: a-2b, kk=2: -(a+b),
\\   kk=3: -(2a-b), kk=4: 2b-a, kk=5: a+b
a = aval; b = bval;
TR = [2*a-b, a-2*b, -(a+b), -(2*a-b), 2*b-a, a+b];
ORD = vector(6, ii, p + 1 - TR[ii]);

print("6 traces:");
for(ii=1,6, print("  t_",ii-1," = ", TR[ii]));
print("6 group orders:");
for(ii=1,6, print("  n_",ii-1," = ", ORD[ii]));
print("Sum checks: t[0]+t[2]+t[4] = ", TR[1]+TR[3]+TR[5]);
print("            t[1]+t[3]+t[5] = ", TR[2]+TR[4]+TR[6]);
print("Sum of squares = 6p? ", (2*a-b)^2+(a-2*b)^2+(a+b)^2 == 6*p);
print("");

\\ --- H1: all orders distinct ---
print("H1: all 6 orders distinct? ", #Set(ORD) == 6);
print("");

\\ --- Find h of order 6 in F_p* ---
\\ p ≡ 3 (mod 4), so sqrt(-3) = (-3)^{(p+1)/4} mod p.
sq3 = lift(Mod(-3, p)^((p+1)\4));
zinv = lift(Mod(2, p)^(-1));
z3 = lift(Mod((-1 + sq3) * zinv, p));       \\ zeta_3 = (-1+sqrt(-3))/2
h_int = lift(Mod(-z3, p));                   \\ h = -zeta_3, order 6
print("zeta3 order = ", znorder(Mod(z3, p)));
print("h = -zeta3, order = ", znorder(Mod(h_int, p)));
print("h^3 = -1 mod p? ", lift(Mod(h_int,p)^3) == p-1);
print("");

\\ --- 2-torsion polynomial factorizations ---
\\ Twist kk uses curve y^2 = x^3 + 7*h^kk.
\\ Compute h^kk once and factor x^3 + 7*h^kk over F_p.
print("2-torsion patterns (factoring x^3 + 7*h^kk mod p):");
DPAT = vector(6);
hpow = Mod(1, p);                  \\ h^0 = 1
{
for(kk=0, 5,
    b7hk = lift(Mod(7,p) * hpow);
    \\ factor x^3 + b7hk over F_p
    fac = factormod(x^3 + b7hk, p);
    nfac = matsize(fac)[1];
    degs = [];
    for(ii=1, nfac,
        dg = poldegree(fac[ii,1]);
        degs = concat(degs, [dg])
    );
    DPAT[kk+1] = degs;
    label = if(degs==[3], "irred", if(degs==[1,2], "1+irred2", "split3"));
    print("  kk=",kk,": pattern=",degs,"  (",label,")");
    hpow = hpow * Mod(h_int, p)     \\ update h^kk
)
}
print("");

\\ --- H2 for all 15 pairs ---
print("H2 check (same 2-torsion Galois structure):");
h2pass = 0; h2fail = 0;
{
for(ii=1, 6,
    for(jj=ii+1, 6,
        ok = (DPAT[ii] == DPAT[jj]);
        if(ok, h2pass++, h2fail++;
               print("  H2 FAIL pair (",ii-1,",",jj-1,"): ",DPAT[ii]," vs ",DPAT[jj]))
    )
)
}
print("H2: ", h2pass, " pass / 15   (", h2fail, " fail)");
print("");

\\ --- H3 for all 15 pairs ---
print("H3 check (gcd of group orders = 1):");
h3pass = 0; h3fail = 0;
{
for(ii=1, 6,
    for(jj=ii+1, 6,
        gc = gcd(ORD[ii], ORD[jj]);
        if(gc == 1, h3pass++,
                    h3fail++;
                    print("  H3 FAIL pair (",ii-1,",",jj-1,") gcd=",gc))
    )
)
}
print("H3: ", h3pass, " pass / 15   (", h3fail, " fail)");
print("");

\\ --- Full table ---
print("================================================================");
print("Full table (all 15 pairs):");
print("pair  H1    H2    H3    ALL");
allpass = 0;
{
for(ii=1, 6,
    for(jj=ii+1, 6,
        H1 = (ORD[ii] != ORD[jj]);
        H2 = (DPAT[ii] == DPAT[jj]);
        H3 = (gcd(ORD[ii], ORD[jj]) == 1);
        ALL = H1 && H2 && H3;
        if(ALL, allpass++);
        print("(",ii-1,",",jj-1,")  ", if(H1,"Y","N"), "     ",
              if(H2,"Y","N"), "     ", if(H3,"Y","N"), "     ", if(ALL,"YES","NO"),
              "  t_i=",TR[ii]," t_j=",TR[jj])
    )
)
}
print("");
print("Pairs satisfying ALL three Howe conditions: ", allpass, "/15");
print("================================================================");
