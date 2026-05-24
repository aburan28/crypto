\\ ============================================================
\\ Howe conditions for all 15 pairs of j=0 sextic twists
\\ secp256k1: y² = x³ + 7  over F_p  (j=0, CM disc = -3)
\\ ============================================================
\\
\\ For p ≡ 1 (mod 6), there are exactly 6 sextic twist isomorphism
\\ classes of j=0 curves over F_p.  They are y²=x³+b_k where
\\ b_k = 7 * u^k mod p, k=0,...,5, and u = g^{(p-1)/6} is a
\\ primitive 6th root of unity in F_p*.
\\
\\ For each of the C(6,2)=15 pairs (i,j), we check:
\\   (H1)  #E_i ≠ #E_j            (E_i ≁ E_j over F_p)
\\   (H2)  x³+b_i, x³+b_j have   (same F_p-Galois module structure
\\          the same degree pattern  of 2-torsion)
\\   (H3)  gcd(#E_i, #E_j) = 1   (no shared F_p-rational order)
\\
\\ Order-identification strategy (avoids slow ellcard for large p):
\\   - We know T[1]=t, T[4]=-t for k=0,3 (secp256k1 and its quad twist).
\\   - CM formula: 4p = t²+3s² gives s; then all 6 traces are
\\     ±t, ±(t+3s)/2, ±(t-3s)/2 (see derive_traces()).
\\   - For each unknown twist (k=1,2,4,5), we find a random affine
\\     point on the curve and check which candidate trace T satisfies
\\     (p+1-T)*P = O (point at infinity).  This is O(log p) per check.
\\
\\ Run: gp -q howe_sextic_twists_all15.gp

default(parisize, 1024000000);
default(timer, 0);

p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
n_known = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;
t_known = p + 1 - n_known;

print("================================================================");
print("Howe (2,2)-conditions for all 15 pairs of secp256k1 sextic twists");
print("================================================================");
print("");
print("p   = ", p);
print("t   = ", t_known, "   (trace of secp256k1)");
print("");

\\ ---- 1. Verify p ≡ 1 (mod 6) → exactly 6 sextic twists ----
print("---- Step 1: Verify p ≡ 1 (mod 6) ----");
r6 = p % 6;
print("p mod 6 = ", r6);
if (r6 != 1, error("p not ≡ 1 mod 6; sextic-twist count is not 6"));
print("p ≡ 1 (mod 6) → there are exactly 6 sextic twist classes ✓");
print("");

\\ ---- 2. CM structure: 4p = t² + 3s² ----
print("---- Step 2: CM structure 4p = t² + 3s² ----");
rem = 4*p - t_known^2;
if (rem % 3 != 0, error("4p - t^2 not divisible by 3"));
s_sq = rem / 3;
s = sqrtint(s_sq);
if (s^2 != s_sq, error("s is not an integer; CM formula failed"));
print("4p - t² = 3 * ", s_sq);
print("s = sqrt((4p-t²)/3) = ", s);
print("Verify 4p = t² + 3s²: ", 4*p == t_known^2 + 3*s^2);
print("");

\\ ---- 3. Six traces from CM ----
\\ Traces of the 6 sextic twists (order corresponds to the ζ₆^k twist):
\\ T_0 =  t (secp256k1 itself)
\\ T_1 = (t - 3s)/2
\\ T_2 = -(t + 3s)/2
\\ T_3 = -t  (quadratic twist)
\\ T_4 = -(t - 3s)/2 = (3s - t)/2
\\ T_5 =  (t + 3s)/2
{
t_plus  = (t_known + 3*s) / 2;
t_minus = (t_known - 3*s) / 2;
}
T = [t_known, t_minus, -(t_known + 3*s)/2, -t_known, (3*s - t_known)/2, t_plus];
N = vector(6, k, p + 1 - T[k]);

print("---- Step 3: Six CM traces and curve orders ----");
print("");
print("  k | trace T_k                       | order #E_k");
print("  --|----------------------------------|-----------------------------------");
for (k = 1, 6,
    printf("  %d | %30d | %d\n", k-1, T[k], N[k])
);
print("");

\\ Sanity: k=0 matches secp256k1, k=3 matches quadratic twist
print("Sanity: N[1] = n_secp256k1?       ", N[1] == n_known);
print("Sanity: N[4] = p+1+t (quad twist)?", N[4] == p + 1 + t_known);
print("");

\\ ---- 4. Find primitive 6th root u in F_p* ----
print("---- Step 4: Find primitive 6th root u mod p ----");
\\  Find a generator g of F_p* (primitive root), then u = g^{(p-1)/6}.
\\  PARI's znprimroot gives the smallest primitive root, but it operates
\\  on Mod objects; use znprimroot(p) directly.
g_gen = znprimroot(p);
u = lift(g_gen ^ ((p - 1) / 6));
print("Primitive root g = ", g_gen);
print("u = g^{(p-1)/6} = ", u);
print("u^6 ≡ 1 (mod p)?  ", Mod(u,p)^6 == 1);
print("u^3 ≡ -1 (mod p)? ", Mod(u,p)^3 == p-1);
print("u^2 ≢ 1 (mod p)?  ", Mod(u,p)^2 != 1);
print("");

\\ ---- 5. List 6 b-values (sextic twist coefficients) ----
print("---- Step 5: b-values of 6 sextic twists y²=x³+b_k ----");
b = vector(6, k, lift(Mod(7, p) * Mod(u, p)^(k-1)));
for (k = 1, 6,
    printf("  k=%d: b_%d = %d\n", k-1, k-1, b[k])
);
print("");
print("b_3 = -7 mod p (i.e., p-7)?  ", b[4] == p - 7);
print("");

\\ ---- 6. Identify trace of each twist via scalar-multiplication test ----
\\
\\ Strategy: for each b_k (k=1,2,4,5; k=0 and k=3 already known),
\\ find a random affine point P = (x, y) with y² ≡ x³+b_k (mod p),
\\ then test which candidate trace T satisfies (p+1-T)*P = O.
\\
\\ For k=0: T = T[1] = t_known (known)
\\ For k=3: T = T[4] = -t_known (known)
\\
\\ For the rest, we iterate over x = 2, 3, ... until y²=x³+b_k has a solution.

print("---- Step 6: Match b_k to trace via point-order test ----");
print("");

trace_of = vector(6, k, 0);  \\ trace_of[k] = trace of E_k
trace_of[1] = T[1];  \\ k=0: secp256k1
trace_of[4] = T[4];  \\ k=3: quadratic twist

{
for (k = 1, 6,
    local(bk, found_x, ysq, y, P, n_cand, test_result);
    if (k == 1 || k == 4, next);  \\ already known
    bk = b[k];
    \\ Find a random affine point on y² = x³ + bk mod p
    found_x = 0;
    for (x_try = 1, 1000,
        ysq = Mod(x_try^3 + bk, p);
        if (issquare(ysq),
            y = lift(sqrt(ysq));
            P = [Mod(x_try, p), Mod(y, p)];
            found_x = x_try;
            break
        )
    );
    if (found_x == 0,
        print("  k=", k-1, ": could not find affine point in x=1..1000 (unlikely)");
        next
    );
    printf("  k=%d: found affine point (%d, %d)\n", k-1, found_x, lift(P[2]));
    \\ Test which trace candidate gives (p+1-T)*P = O
    local(E_k, match_found);
    E_k = ellinit([Mod(0,p), Mod(bk,p)]);
    match_found = 0;
    for (ti = 1, 6,
        n_cand = p + 1 - T[ti];
        test_result = ellmul(E_k, P, n_cand);
        \\ test_result is O iff point has order dividing n_cand
        if (test_result == ellinf(),
            printf("  k=%d: matches T[%d] = %d (order %d)\n",
                   k-1, ti, T[ti], n_cand);
            trace_of[k] = T[ti];
            match_found = 1;
            break
        )
    );
    if (!match_found,
        printf("  k=%d: WARNING — no trace match found!\n", k-1)
    )
);
}

print("");
print("Trace assignments:");
for (k = 1, 6,
    printf("  b_%d = %-70d  →  trace %d\n", k-1, b[k], trace_of[k])
);
print("");

\\ ---- 7. H2: factorisation of x³+b_k for each k ----
print("---- Step 7: (H2) 2-torsion factorisation of x³+b_k ----");
print("");
print("  k | deg-pattern | irreducible? | note");
print("  --|-------------|--------------|------------------------------");
deg_pat = vector(6);
{
for (k = 1, 6,
    local(bk, poly_mod, fac, nrows, pat, irreducible_flag, note);
    bk = b[k];
    poly_mod = Mod(1, p) * (x^3 + bk);
    fac = factor(poly_mod);
    nrows = matsize(fac)[1];
    pat = [];
    for (i = 1, nrows, pat = concat(pat, [poldegree(fac[i,1])]));
    deg_pat[k] = pat;
    irreducible_flag = (pat == [3]);
    note = "";
    if (pat == [1,1,1], note = "splits completely; E[2] all F_p-rational");
    if (pat == [1,2],   note = "1 F_p-rational 2-torsion pt; Gal swaps other two");
    if (pat == [3],     note = "no F_p-rational 2-torsion; Gal acts cyclically (Z/3)");
    printf("  %d | %-11s | %-12s | %s\n", k-1, Str(pat), Str(irreducible_flag), note)
);
}
print("");

\\ ---- 8. All 15 pairwise Howe checks ----
print("================================================================");
print("Pairwise Howe conditions for all C(6,2)=15 pairs");
print("================================================================");
print("");
print("  (i,j) | H1: n_i≠n_j | H2: same 2-tor | H3: gcd=1 | Glueable?");
print("  ------|-------------|----------------|-----------|----------");

count_yes = 0;
count_no  = 0;
glueable_pairs = [];

{
for (i = 1, 6,
    for (j = i+1, 6,
        local(ni, nj, ti, tj, h1, h2, h3, glueable);
        ti = trace_of[i];
        tj = trace_of[j];
        ni = p + 1 - ti;
        nj = p + 1 - tj;
        h1 = (ni != nj);
        h2 = (deg_pat[i] == deg_pat[j]);
        h3 = (gcd(ni, nj) == 1);
        glueable = h1 && h2 && h3;
        printf("  (%d,%d)  | %-11s | %-14s | %-9s | %s\n",
               i-1, j-1,
               if(h1,"✓","✗"),
               if(h2,"✓","✗"),
               if(h3,"✓","✗"),
               if(glueable,"YES ← Howe-glueable","no"));
        if (glueable,
            count_yes++;
            glueable_pairs = concat(glueable_pairs, [[i-1, j-1]])
        ,
            count_no++
        )
    )
);
}
print("");
printf("Howe-glueable pairs: %d / 15\n", count_yes);
if (#glueable_pairs > 0,
    printf("Glueable pair indices: ");
    for (r = 1, #glueable_pairs, printf("(%d,%d) ", glueable_pairs[r][1], glueable_pairs[r][2]));
    print("")
);
print("");

\\ ---- 9. Structural summary ----
print("================================================================");
print("Structural summary");
print("================================================================");
print("");
print("The 6 sextic twists of secp256k1 (j=0, p ≡ 1 mod 6):");
print("");
print("  2-torsion factorisation patterns among the 6 twists:");
{
local(unique_pats);
unique_pats = [];
for (k = 1, 6,
    local(found);
    found = 0;
    for (r = 1, #unique_pats,
        if (unique_pats[r] == deg_pat[k], found = 1; break)
    );
    if (!found, unique_pats = concat(unique_pats, [deg_pat[k]]))
);
printf("  %d distinct patterns: ", #unique_pats);
for (r = 1, #unique_pats, printf("%s ", Str(unique_pats[r])));
print("")
}
print("");
print("  (H2) holds between two twists iff they have the same");
print("  2-torsion pattern, i.e., their b-values are in the same");
print("  cubic residue class mod p (since -1 is a cube when p ≡ 1 mod 6).");
print("");
print("  (H3): gcd(n_i, n_j) for all 15 pairs — see table above.");
print("");
print("  Implications for the isogeny-graph ECDLP analysis:");
print("  If Howe-glueable pairs exist, there exist genus-2 curves C/F_p");
print("  with Jac(C) → E_i × E_j via a (2,2)-isogeny, for each such pair.");
print("  Per the structural-completeness paper, no such cover provides");
print("  a sub-√p ECDLP algorithm.");
print("");
print("================================================================");
print("Done.");
print("================================================================");
