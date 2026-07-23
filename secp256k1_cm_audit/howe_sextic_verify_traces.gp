\\ ============================================================
\\ Empirically verify trace assignment for all 6 sextic twists
\\ of secp256k1, then re-run all 15 Howe checks.
\\
\\ The original script assumed b_k <-> T_k by CM ordering without
\\ empirical verification.  This script uses ellcard() to get the
\\ true order of each twist, then checks H1/H2/H3 correctly.
\\
\\ Run: gp -q howe_sextic_verify_traces.gp

default(parisize, 512000000);
default(timer, 1);

p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
n_known = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;
t_known = p + 1 - n_known;

print("================================================================");
print("Empirical trace + Howe check for 6 sextic twists of secp256k1");
print("================================================================");
print("p = ", p);
print("t (secp256k1) = ", t_known);
print("");

\\ s from 4p = t^2 + 3s^2
s = sqrtint((4*p - t_known^2) / 3);
print("s = ", s, "  (verify: 4p=t^2+3s^2? ", 4*p == t_known^2 + 3*s^2, ")");
print("");

\\ Primitive 6th root of unity u in F_p*
gg = znprimroot(p);
u = lift(gg ^ ((p - 1) / 6));
print("u = g^{(p-1)/6} = ", u);
print("u^6 = 1? ", Mod(u,p)^6 == 1, "  u^3 = -1? ", Mod(u,p)^3 == p-1);
print("");

\\ The 6 b-values: b[k] = 7 * u^(k-1) mod p  (k = 1..6, i.e., k-1 = 0..5)
bv = vector(6, k, lift(Mod(7,p) * Mod(u,p)^(k-1)));
print("b-values (1-indexed, k-1 = twist index):");
for (k=1,6, print("  twist ", k-1, ": b = ", bv[k]));
print("");

\\ ---- True orders via ellcard ----
print("----------------------------------------------------------------");
print("Computing true orders via ellcard (SEA algorithm, ~10s each)...");
print("----------------------------------------------------------------");
Nv = vector(6);
Tv = vector(6);
{
for (k = 1, 6,
    local(E_k, nk);
    print("  Computing order for twist k=", k-1, " ...");
    E_k = ellinit([0, bv[k]], p);
    nk = ellcard(E_k);
    Nv[k] = nk;
    Tv[k] = p + 1 - nk;
    print("    order = ", nk, "   trace = ", Tv[k])
);
}
print("");
print("Verified orders:");
for (k=1,6, print("  twist ", k-1, ": N = ", Nv[k]));
print("");

\\ Sanity checks
print("Sanity checks:");
print("  N[0] = secp256k1 order?  ", Nv[1] == n_known);
print("  N[3] = quad-twist order? ", Nv[4] == p + 1 + t_known);
print("");

\\ ---- H2: 2-torsion factorisation ----
print("----------------------------------------------------------------");
print("H2: 2-torsion polynomial x^3 + b_k factorisation over F_p");
print("----------------------------------------------------------------");
dpat = vector(6);
{
for (k = 1, 6,
    local(fac, nrows, pat, note);
    fac = factor(Mod(1,p) * Pol([1, 0, 0, bv[k]]));
    nrows = matsize(fac)[1];
    pat = [];
    for (i=1,nrows, pat = concat(pat, [poldegree(fac[i,1])]));
    dpat[k] = pat;
    note = if(pat==[3], "no F_p-rat 2-tor",
           if(pat==[1,1,1], "E[2] in E(F_p)",
           if(pat==[1,2], "1 F_p-rat 2-tor pt", "other")));
    printf("  twist %d: deg-pattern=%s  (%s)  [4|N? %d]\n",
           k-1, Str(pat), note, Nv[k]%4==0)
);
}
print("");

\\ ---- 15 pairwise Howe checks ----
print("================================================================");
print("All 15 pairwise Howe (H1)+(H2)+(H3) checks (verified traces)");
print("================================================================");
print("");
print("  (i,j) | H1:n_i!=n_j | H2:same2-tor | H3:gcd=1 | Glueable?");
print("  ------|-------------|--------------|----------|----------");
cnt = 0;
glue_pairs = List([]);
{
for (i=1,6,
    for (j=i+1,6,
        local(ni, nj, h1, h2, h3, g, gl);
        ni = Nv[i]; nj = Nv[j];
        h1 = (ni != nj);
        h2 = (dpat[i] == dpat[j]);
        g  = gcd(ni, nj);
        h3 = (g == 1);
        gl = h1 && h2 && h3;
        if (gl, cnt++; listput(glue_pairs, [i-1, j-1]));
        printf("  (%d,%d)  | %s | %s | %s (gcd=%d) | %s\n",
               i-1, j-1,
               if(h1,"YES","NO "),
               if(h2,"YES","NO "),
               if(h3,"YES","NO "), g,
               if(gl,"YES <- Howe-glueable","no"))
    )
);
}
print("");
print("Howe-glueable pairs (empirically verified): ", cnt, " / 15");
if (#glue_pairs > 0,
    for (r=1,#glue_pairs,
         print("  (", glue_pairs[r][1], ",", glue_pairs[r][2], ")")
    )
);
print("");
print("================================================================");
print("Done.");
print("================================================================");
