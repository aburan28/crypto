\\ Thread 18: Howe-gluing conditions for all 15 pairs of j=0 sextic twists
\\ secp256k1: y^2 = x^3 + 7 over F_p
\\ Run: gp -q thread18_sextic_twist_howe.gp

default(parisize, 512000000);

p     = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
n_ec  = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;
bval  = 7;

print("================================================================");
print("Thread 18: Howe-gluing for all 15 pairs of j=0 sextic twists");
print("secp256k1: y^2 = x^3 + 7");
print("================================================================");
print("");

{

\\ ---- 1. CM group orders (no SEA) ----
t0 = p + 1 - n_ec;
print("t0 = p+1-n = ", t0);
print("t0 mod 2   = ", t0 % 2);
print("t0 mod 3   = ", t0 % 3);
rhs = 4*p - t0^2;
if (rhs % 3 != 0, error("4p-t0^2 not divisible by 3"));
vv = sqrtint(rhs / 3);
if (vv^2 != rhs/3, error("v not integer"));
print("v (4p=t0^2+3v^2): v = ", vv);
print("Verify 4p=t0^2+3v^2: ", t0^2 + 3*vv^2 == 4*p);
print("");

\\ 6 traces in canonical quad-twist-pair order:
\\   (k=0,k=3)=(t0,-t0)  (k=1,k=4)=(c1,-c1)  (k=2,k=5)=(c2,-c2)
c1 = (t0 - 3*vv) / 2;
c2 = -(t0 + 3*vv) / 2;
if ((t0 - 3*vv) % 2 != 0, error("c1 not integer"));
if ((t0 + 3*vv) % 2 != 0, error("c2 not integer"));

TR  = [t0, c1, c2, -t0, -c1, -c2];
ORD = [p+1-TR[1], p+1-TR[2], p+1-TR[3], p+1-TR[4], p+1-TR[5], p+1-TR[6]];

print("6 CM traces and N=p+1-trace (quad-twist pairs: (0,3),(1,4),(2,5)):");
for (kk = 1, 6,
    print("  k=", kk-1, ": trace=", TR[kk], "  N mod 4=", ORD[kk] % 4)
);
print("");

\\ ---- 2. Sextic non-residue dnr ----
e6   = (p-1) / 6;
e3   = (p-1) / 3;
dnr  = 2;
while (1,
    gg = Mod(dnr, p)^e6;
    if (gg^3 != 1 && gg^2 != 1, break);
    dnr++
);
print("Sextic non-residue dnr = ", dnr);
print("");

\\ ---- 3. 2-torsion type for each sextic twist k=0..5 ----
\\ Twist k: y^2=x^3+7*dnr^k. 2-torsion poly: x^3+bk.
\\ For p=1 mod 3: no [1,2] factorisation. Only [3] or [1,1,1].
\\ Type B iff -bk is a cube in F_p iff Mod(-bk,p)^((p-1)/3)==1.

TYP = vector(6);
BKV = vector(6);

print("2-torsion types (A=irred [3], B=split [1,1,1]):");
for (kk = 0, 5,
    bk = lift(Mod(bval, p) * Mod(dnr, p)^kk);
    BKV[kk+1] = bk;
    cub_char = lift(Mod(-bk, p)^e3);
    is_B = (cub_char == 1);
    TYP[kk+1] = is_B;
    fac  = factor(Mod(x^3 + bk, p));
    nf   = matsize(fac)[1];
    d1   = poldegree(fac[1, 1]);
    tstr = if (is_B, "B [1,1,1]", "A [3]    ");
    print("  k=", kk, ": type ", tstr,
          "  cubic_char(-bk)=", cub_char,
          "  top_factor_deg=", d1)
);
print("");

\\ Cross-check: count type-B from factorisation vs from N%4==0
cnt_B_fact = 0;
cnt_B_cm   = 0;
for (kk = 0, 5, if (TYP[kk+1],   cnt_B_fact++));
for (kk = 0, 5, if (ORD[kk+1]%4 == 0, cnt_B_cm++));
print("Type-B count (factorisation): ", cnt_B_fact);
print("Type-B count (CM N%4==0):     ", cnt_B_cm);
print("Match: ", cnt_B_fact == cnt_B_cm);
print("");

\\ ---- 4. Howe H1/H2/H3 for all 15 pairs ----
\\ Use canonical CM ordering for group orders.
\\ H2 assignment: CM-index with N%4==0 -> B, else A.
\\ This gives correct H2 answer because both classifiers agree on count,
\\ and all 6 orders are distinct, so the within-type gcd analysis is exact.

print("================================================================");
print("H1/H2/H3 for all 15 pairs (canonical CM ordering)");
print("  pair  H1       H2          H3");
print("  ----  --       --          --");

npair_all3 = 0;
for (ii = 0, 4,
    for (jj = ii+1, 5,
        Ni = ORD[ii+1];
        Nj = ORD[jj+1];
        h1 = (Ni != Nj);
        Bi = (Ni % 4 == 0);
        Bj = (Nj % 4 == 0);
        h2 = (Bi == Bj);
        gg = gcd(Ni, Nj);
        h3 = (gg == 1);
        ok = h1 && h2 && h3;
        if (ok, npair_all3++);
        h2s = if (h2, if(Bi, "OK(B+B)", "OK(A+A)"), "FAIL(A+B)");
        h3s = if (h3, "gcd=1", concat("gcd=", Str(gg)));
        print("  (", ii, ",", jj, ")  ",
              if(h1,"H1=OK  ","H1=FAIL"),
              "  ", h2s,
              "  ", h3s,
              if(ok,"  <<< PASS",""))
    )
);
print("");
print("Pairs with ALL 3 Howe conditions: ", npair_all3, " / 15");
print("");

\\ ---- 5. The type-B gcd specifically ----
print("--- Type-B pair gcd ---");
NB0 = 0; NB1 = 0; cnt_B = 0;
for (kk = 0, 5,
    if (ORD[kk+1] % 4 == 0,
        cnt_B++;
        if (cnt_B == 1, NB0 = ORD[kk+1]);
        if (cnt_B == 2, NB1 = ORD[kk+1])
    )
);
if (cnt_B == 2,
    gB = gcd(NB0, NB1);
    print("NB0 = ", NB0);
    print("NB1 = ", NB1);
    print("gcd(NB0,NB1) = ", gB);
    print("H3 for B-B pair: ", if(gB==1,"PASS","FAIL"))
);
print("");

\\ ---- 6. H1 sanity: all 6 traces distinct ----
h1_ok = 1;
for (ii = 1, 5, for (jj = ii+1, 6, if (TR[ii]==TR[jj], h1_ok=0)));
print("All 6 CM traces distinct (H1 holds for all 15 pairs): ", h1_ok);
print("");

\\ ---- 7. Final count ----
nA = 0; nB = 0;
for (kk = 0, 5, if (TYP[kk+1], nB++, nA++));
pairsH2 = nA*(nA-1)/2 + nB*(nB-1)/2;
print("================================================================");
print("Summary");
print("================================================================");
print("Type-A twists (irred 2-tor): ", nA, " of 6");
print("Type-B twists (split 2-tor): ", nB, " of 6");
print("Pairs satisfying H2: C(", nA, ",2)+C(", nB, ",2) = ", pairsH2, " of 15");
print("H1 satisfied for all 15 pairs: ", h1_ok);
print("Pairs satisfying all 3 conditions: ", npair_all3, " of 15");
print("");
print("ECDLP implication: each passing pair (E_i, E_j) admits a Howe (2,2)-gluing");
print("to a smooth genus-2 Jacobian J/F_p with |J(F_p)| ~ p^2. Pollard-rho");
print("cost on J is O(p), strictly worse than ECDLP on E (cost O(sqrt(p))).");
print("No ECDLP speedup from any Howe gluing of j=0 sextic twist pairs.");
print("================================================================");
print("Thread 18 complete.");
print("================================================================");

}
