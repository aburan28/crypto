\\ ============================================================
\\ (3,3)-Howe gluing check for j=0 sextic twists
\\ ============================================================
\\
\\ Thread 3 (2026-07-08) confirmed that 5 of 15 pairs of secp256k1's
\\ sextic twists are Howe-(2,2)-glueable.  This script asks whether
\\ any pair is also (3,3)-glueable — i.e., whether there exists a
\\ smooth genus-2 Jacobian (3,3)-isogenous to E_i × E_j for some
\\ sextic-twist pair.
\\
\\ Motivation: Decru-Kunzweiler (eprint.iacr.org/2026/334) provide
\\ efficient tripling formulas on Hessian (j=0) Kummer lines via
\\ (3,3)-isogeny decomposition.  The analogous Howe-gluing condition
\\ for (3,3) is:
\\
\\   (A1) Hom_{F_p}(E_i, E_j) = 0  (not F_p-isogenous; equiv #E_i ≠ #E_j)
\\   (A2) E_i[3] ≃ E_j[3] as F_p-Galois modules (3-torsion compatibility)
\\   (A3) (N,N)-kernel compatibility: existence of F_p-equivariant iso
\\        φ: E_i[3] → E_j[3] preserving the Weil e_3 pairing
\\
\\ For j=0 curves, the 3-torsion polynomial of y²=x³+b is:
\\   ψ_3(x) = 3x(x^3 + 2b)   [classical 3-division polynomial]
\\ The roots of x^3 + 2b give the non-identity 3-torsion x-coordinates.
\\
\\ Key observation for j=0 sextic twists:
\\ All 6 sextic twists have the same j-invariant j=0, and the CM field
\\ Q(√-3) acts on their 3-torsion via the same Galois representation
\\ (since the 3-torsion field depends only on j for CM curves).
\\ Hence (A2) and (A3) may be automatic for all pairs — in which case
\\ gluability reduces to (A1) alone (non-isogenous, i.e., different
\\ group orders).
\\
\\ This script verifies these conditions for a SMALL prime p with
\\ j=0 curve, then checks the secp256k1 trace analytically.
\\
\\ Run: gp -q howe_33_gluing_check.gp
\\
\\ PARI/GP not installed in CI container; expected output is NOT
\\ committed.  This is a research prototype.

default(parisize, 256000000);
default(timer, 0);

print("================================================================");
print("(3,3)-Howe gluing check on j=0 sextic twists");
print("================================================================");
print("");
print("Background: Decru-Kunzweiler (ePrint 2026/334) give efficient");
print("(3,3)-isogeny formulas on Hessian (j=0) curves.  We check the");
print("Howe-gluing conditions for (3,3)-covers of secp256k1-type curves.");
print("");

\\ ---- 1. Small prime setup ----
\\ Find a prime p ≡ 1 (mod 6) with j=0 curve E: y^2 = x^3 + b
\\ having prime group order (to match secp256k1 regime).
print("---- 1. Small prime p ≡ 1 (mod 6) with prime-order j=0 curve ----");
found_p = 0; found_b = 0; found_t = 0;
{
forprime(p0 = 7, 3000,
    if (Mod(p0, 6) != 1, next);
    for (b0 = 1, 20,
        local(E0, n0);
        E0 = ellinit([0, b0], p0);
        n0 = ellcard(E0);
        if (isprime(n0),
            found_p = p0; found_b = b0;
            found_t = p0 + 1 - n0;
            break(2)
        )
    )
);
}
p = found_p;
print("p = ", p, "   b = ", found_b, "   t (Frobenius trace) = ", found_t);
print("n = p+1-t = ", p + 1 - found_t, " (prime order)");
print("");

\\ ---- 2. Enumerate the 6 sextic twists of E: y^2 = x^3 + b ----
\\ Sextic twists of j=0 curves over F_p (p ≡ 1 mod 6):
\\ Twist by d ∈ F_p* / (F_p*)^6.  The subgroup (F_p*)^6 has index 6 in F_p*.
\\ Representatives: pick d_i with different sextic residue characters.
print("---- 2. Sextic twists of E ----");

\\ Find a generator of F_p*
gen = 0;
{
for (g0 = 2, p - 1,
    if (znorder(Mod(g0, p)) == p - 1,
        gen = g0;
        break
    )
);
}
print("Generator of F_p* = ", gen);

\\ The 6 sextic twist representatives are gen^k for k = 0..5
\\ Twist d: E_d: y^2 = x^3 + b * d^5   (standard sextic twist formula for j=0)
\\ (Note: y^2=x^3+b twisted by d gives y^2=x^3+b*d^(-3)… depends on convention.
\\  We use: E_d: y^2 = x^3 + b*d  where d ∈ F_p*/(F_p*)^6.)
\\  Under substitution (x,y) -> (d^2 x, d^3 y) this becomes y^2 = x^3 + b*d^(-5).
\\  Standard: twists of y^2 = x^3 + c are y^2 = x^3 + c*d where d runs over
\\  F_p*/(F_p*)^6.  We take representatives d = gen^k, k=0..5.)

twist_reps = vector(6);
twist_traces = vector(6);
twist_orders = vector(6);
{
for (k = 0, 5,
    local(d_rep, E_k, n_k);
    d_rep = Mod(gen, p)^k;
    E_k = ellinit([0, found_b * lift(d_rep)], p);
    n_k = ellcard(E_k);
    twist_reps[k+1] = lift(d_rep);
    twist_orders[k+1] = n_k;
    twist_traces[k+1] = p + 1 - n_k;
    print("  Twist k=", k, ": d=", lift(d_rep), "   #E = ", n_k, "   t = ", p+1-n_k)
);
}
print("");

\\ ---- 3. Check condition (A1) for all 15 pairs ----
print("---- 3. Condition (A1): non-isogenous pairs ----");
a1_pass = 0; a1_fail = 0;
{
for (i = 1, 6,
    for (j = i+1, 6,
        local(same_order);
        same_order = (twist_orders[i] == twist_orders[j]);
        if (!same_order,
            a1_pass++;
            \\print("  (", i, ",", j, "): #E_i≠#E_j  (A1) ✓")
        ,
            a1_fail++;
            print("  (", i, ",", j, "): #E_i=#E_j  (A1) FAIL: ", twist_orders[i])
        )
    )
);
}
print("(A1) passes: ", a1_pass, " / 15 pairs");
print("(A1) fails:  ", a1_fail, " / 15 pairs  (isogenous pairs share order)");
print("");

\\ ---- 4. Check condition (A2): 3-torsion Galois module structure ----
\\ E_k: y^2 = x^3 + b*d_k.  3-torsion polynomial: x^3 + 2*b*d_k.
\\ Galois module E_k[3] determined by splitting of x^3 + 2*b*d_k over F_p.
print("---- 4. Condition (A2): 3-torsion Galois module structure ----");
three_tors_patterns = vector(6);
{
for (k = 1, 6,
    local(b_k, poly_3, fact_3, pat);
    b_k = found_b * twist_reps[k];
    poly_3 = Mod(x^3 + 2*b_k, p);
    fact_3 = factor(poly_3);
    pat = [];
    for (r = 1, matsize(fact_3)[1],
        local(deg_r, mult_r);
        deg_r = poldegree(fact_3[r, 1]);
        mult_r = fact_3[r, 2];
        for (s = 1, mult_r, pat = concat(pat, [deg_r]))
    );
    three_tors_patterns[k] = pat;
    print("  Twist k=", k-1, " (d=", twist_reps[k], "): x^3+2b*d factors as deg-pattern ", pat)
);
}
print("");

\\ Check pairwise compatibility
print("3-torsion pattern pairwise matches:");
a2_pass = 0; a2_fail = 0;
{
for (i = 1, 6,
    for (j = i+1, 6,
        local(match);
        match = (three_tors_patterns[i] == three_tors_patterns[j]);
        if (match, a2_pass++, a2_fail++);
        print("  (", i-1, ",", j-1, "): ", if(match, "(A2) ✓ same pattern", "(A2) ✗ DIFFERENT"))
    )
);
}
print("(A2) compatible pairs: ", a2_pass, " / 15");
print("(A2) incompatible:     ", a2_fail, " / 15");
print("");

\\ ---- 5. Pairs satisfying both (A1) and (A2) ----
print("---- 5. Pairs satisfying (A1) AND (A2) ----");
both_pass = 0;
{
for (i = 1, 6,
    for (j = i+1, 6,
        local(a1_ok, a2_ok);
        a1_ok = (twist_orders[i] != twist_orders[j]);
        a2_ok = (three_tors_patterns[i] == three_tors_patterns[j]);
        if (a1_ok && a2_ok,
            both_pass++;
            print("  (", i-1, ",", j-1, "): (A1)+(A2) satisfied — candidate for (3,3)-Howe-gluing")
        )
    )
);
}
print("Pairs satisfying (A1)+(A2): ", both_pass, " / 15");
print("");
print("Note: (A3) — Weil pairing compatibility of the isomorphism");
print("E_i[3] → E_j[3] — requires explicit 3-torsion point computation");
print("and is not checked here.  For j=0 curves, the CM action of Z[ζ_3]");
print("on the 3-torsion is the same for all sextic twists, so (A3) may");
print("follow automatically from (A2).");
print("");

\\ ---- 6. Structural observation for secp256k1 ----
print("---- 6. Implication for secp256k1 ----");
print("");
print("For secp256k1: p = 2^256 - 2^32 - 977, p ≡ 1 (mod 6).");
print("The 6 sextic twists have traces t_k = ζ_6^k * t_0 where ζ_6 is");
print("a primitive 6th root of unity in Z (i.e., ±1, ±t_0/2 ± ... ).");
print("More precisely, for CM by Z[ζ_3], the traces cycle as:");
print("  {t_0, -t_0, t_0/2+..., ...} — the 6 Weil numbers.");
print("");
print("The 3-torsion polynomial for E_k: y^2=x^3+b_k is x^3 + 2*b_k.");
print("For secp256k1: b_0 = 7, so x^3 + 14 mod p.");
print("All sextic twists have 3-torsion poly x^3 + 2*b*d_k where d_k = gen^k.");
print("Since d_k = gen^k, we have 2*b*d_k = 2*7*gen^k.");
print("The splitting of x^3 + c over F_p depends on whether c is a cube:");
print("  x^3 + c irreducible over F_p iff c is not a cube mod p.");
print("  For k=0: c = 14.  For k=1: c = 14*gen.  The cube character of c");
print("  changes with k, so different twists may have different 3-torsion patterns.");
print("");
print("Prediction: among 6 sextic twists, 2 will have c = cube (splitting),");
print("4 will have c = non-cube (irreducible) OR split differently.");
print("Exact distribution depends on whether 14 is a cubic residue mod p_secp256k1.");
print("This script computes the exact pattern for the small prime p above.");
print("");

print("================================================================");
print("(3,3)-Howe gluing check complete.");
print("See RESEARCH_AUTOLAB_LOG.md for context and next-step proposal.");
print("================================================================");
