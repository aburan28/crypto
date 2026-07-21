\\ ============================================================
\\ Howe (2,2)-gluing: all 15 pairs of secp256k1 j=0 sextic twists
\\ ============================================================
\\ secp256k1 has j=0; twists classified by F_p*/(F_p*)^6 (|class|=6,
\\ since p ≡ 1 mod 6).  Check all C(6,2)=15 pairs for H1+H2+H3.
\\
\\ (H1)  #E_i ≠ #E_j
\\ (H2)  E_i[2] ≃ E_j[2] as F_p-Galois modules
\\          ⟺  d_j/d_i is a cube in F_p*  (for j=0 curves)
\\ (H3)  gcd(#E_i, #E_j) = 1
\\
\\ Run: gp -q howe_sextic_twists.gp
\\ ============================================================

default(parisize, 256000000);
default(parisizemax, 2000000000);

p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
n = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;

print("=================================================================");
print("Howe gluing: all 15 pairs of secp256k1 j=0 sextic twists");
print("=================================================================");

\\ ---- Step 1: CM decomposition ----
\\ Find a,b in Z with a^2+ab+b^2=p and 2a+b = t (Frobenius trace).
\\ From 3a^2-3at+t^2=p, discriminant/3 = (4p-t^2)/3 must be a square.

t_secp = p + 1 - n;

disc3 = (4*p - t_secp^2) / 3;
sq = sqrtint(disc3);
if (sq^2 != disc3, error("disc3 not a perfect square"));

\\ Two candidate pairs; pick the one satisfying 2a+b=t
a1 = (t_secp + sq) / 2;
b1 = t_secp - 2*a1;
a2 = (t_secp - sq) / 2;
b2 = t_secp - 2*a2;

a_cm = 0; b_cm = 0;
if (a1^2 + a1*b1 + b1^2 == p, a_cm = a1; b_cm = b1,
                               a_cm = a2; b_cm = b2);

if (a_cm^2 + a_cm*b_cm + b_cm^2 != p, error("CM decomp failed"));

print("\nCM decomp: a=", a_cm, ", b=", b_cm);
print("Check: a^2+ab+b^2=p? ", a_cm^2+a_cm*b_cm+b_cm^2==p);
print("Check: 2a+b=t?       ", 2*a_cm+b_cm == t_secp);

\\ ---- Step 2: The 6 Frobenius traces ----
\\ Twisting by omega^k (primitive 6th root in Z[omega]) gives:
\\  k=0: tr = 2a+b       k=3: tr = -(2a+b)
\\  k=1: tr = a-b        k=4: tr = -(a-b) = b-a
\\  k=2: tr = -a-2b      k=5: tr = -(-a-2b) = a+2b

a = a_cm; b = b_cm;
T = [2*a+b, a-b, -a-2*b, -(2*a+b), -(a-b), a+2*b];

print("\n6 Frobenius traces (indexed 0..5):");
for (k=1, 6, print("  T[", k-1, "] = ", T[k], "   #E = ", p+1-T[k]));

if (T[1] != t_secp, error("T[0] should equal secp256k1 trace"));
print("Sanity T[0]=t: OK");

\\ ---- Step 3: Sextic-twist b-coefficients ----
\\ Primitive 6th root of unity omega_p in F_p:
\\ omega satisfies x^2-x+1=0 in Z[omega], so 2*omega-1 = sqrt(-3).
\\ Since p ≡ 3 mod 4, sqrt(-3) mod p = (-3)^((p+1)/4) mod p.

sq3 = lift(Mod(-3,p)^((p+1)/4));
if (sq3^2 % p != p-3, error("sqrt(-3) failed"));
\\ omega_p = (1 + sqrt(-3)) / 2 mod p
omega_p = lift(Mod(1 + sq3, p) * Mod(2,p)^(p-2));
if ((omega_p^2 - omega_p + 1) % p != 0, error("omega_p not a 6th root"));
print("\nomega_p = g^((p-1)/6) mod p = ", omega_p);
print("Check (omega_p)^2 - omega_p + 1 ≡ 0 mod p: ", (omega_p^2 - omega_p + 1) % p == 0);

\\ d_k = 7 * omega_p^k mod p  (representative b-coefficient of k-th sextic twist)
d = vector(6, k, lift(Mod(7,p) * Mod(omega_p,p)^(k-1)));
print("\nd_k = 7 * omega_p^k mod p:");
for (k=1, 6, print("  d[", k-1, "] = ", d[k]));

\\ Verify: d[0]=7 (original secp256k1), d[3] is quadratic twist
\\ (omega_p^3 = -1 mod p since omega_p has order 6)
print("omega_p^3 mod p = ", lift(Mod(omega_p,p)^3), "  (should be p-1 = ", p-1, ")");

\\ ---- Step 4: 2-torsion splitting types ----
\\ E_k: y^2 = x^3 + d[k].  2-torsion poly = x^3 + d[k].
\\ Irreducible iff -d[k] is not a cube mod p (no rational 2-torsion).

print("\n2-torsion poly x^3+d[k] splitting types:");
e3 = (p-1)/3;
split_irr = vector(6);
for (k=1, 6, {
    f = factormod(Pol([1,0,0,d[k]]), p);
    nr = matsize(f)[1];
    degs = [];
    for (i=1, nr, degs = concat(degs, [poldegree(f[i,1])]));
    degs = vecsort(degs);
    split_irr[k] = (degs == [3]);
    print("  d[", k-1, "]: degree pattern ", degs,
          if(split_irr[k]," (irreducible)"," (reducible)"))
});

\\ ---- Step 5: Elliptic curve group orders via ellcard ----
\\ Direct computation (PARI SEA algorithm for 256-bit primes).

print("\nComputing ellcard for 5 non-secp twists (SEA, ~10-60s each)...");
orders = vector(6);
orders[1] = n;   \\ d[0]=7 is secp256k1

\\ Validate: confirm secp256k1 order matches
E0 = ellinit([0,0,0,0,Mod(7,p)]);
ord0 = ellcard(E0);
if (ord0 != n, error("ellcard(secp256k1) != n"));
print("  d[0]=7: ellcard = n ✓");

for (k=2, 6, {
    Ek = ellinit([0,0,0,0,Mod(d[k],p)]);
    orders[k] = ellcard(Ek);
    print("  d[", k-1, "]: #E = ", orders[k], "   trace = ", p+1-orders[k])
});

print("\nGroup orders and traces:");
for (k=1, 6, print("  d[", k-1, "]: #E=", orders[k], "  t=", p+1-orders[k]));

\\ Verify all traces belong to the CM set {T[1]..T[6]}
trace_set = Set(T);
for (k=1, 6, {
    tr_k = p + 1 - orders[k];
    if (!setsearch(trace_set, tr_k), error("trace d[", k-1, "] not in CM set"))
});
print("All 6 traces lie in {T[0]..T[5]}: OK");

\\ ---- Step 6: All 15 pairs ----
print("\n=================================================================");
print("All 15 pairs — Howe (H1)+(H2)+(H3) check");
print("=================================================================\n");

howe_count = 0;
h2_ok_count = 0;
for (i=0, 4, {
    for (j=i+1, 5, {
        ni = orders[i+1];
        nj = orders[j+1];

        \\ H1: different group orders
        h1 = (ni != nj);

        \\ H2: d[j]/d[i] is a cube in F_p
        ratio = lift(Mod(d[j+1],p) * Mod(d[i+1],p)^(p-2));
        h2 = (lift(Mod(ratio,p)^e3) == 1);

        \\ H3: gcd(#E_i, #E_j) = 1
        g_ij = gcd(ni, nj);
        h3 = (g_ij == 1);

        ok = h1 && h2 && h3;
        if (h2, h2_ok_count++);
        if (ok, howe_count++);

        print("  Pair (", i, ",", j, "):  H1=", if(h1,"✓","✗"),
              "  H2=", if(h2,"✓ [ratio is cube]","✗ [ratio not cube]"),
              "  H3=", if(h3,"✓","✗ [gcd="),
              if(!h3, g_ij, ""),
              if(!h3, "]", ""),
              "  =>  ", if(ok,"HOWE-GLUABLE","not gluable"))
    })
});

print("\n=================================================================");
print("Summary");
print("=================================================================");
print("  Pairs with (H2)=same Galois module:   ", h2_ok_count, " / 15");
print("  Pairs with (H1)+(H2)+(H3) all pass:   ", howe_count, " / 15");

if (howe_count == 3,
    print("\nResult: exactly 3 Howe-gluable pairs (one per cubic-twist class)."),
if (howe_count > 0,
    print("\nResult: ", howe_count, " Howe-gluable pair(s)."),
    print("\nResult: no Howe-gluable pairs among 6 sextic twists.")));

print("\nStructural note:");
print("  (H2) holds iff d[j]/d[i] = omega_p^(j-i) is a cube,");
print("  iff (j-i) ≡ 0 mod 3.  Exactly 3 pairs: {0,3},{1,4},{2,5}.");
print("  Each is the pair (cubic-twist-class representative,");
print("       its quadratic twist within that class).");
print("  Howe gluing cannot cross cubic-twist classes.");
print("\n=================================================================");
print("howe_sextic_twists.gp complete.");
print("=================================================================");
