\\ ============================================================
\\ Correct Howe cover via F_{p^3} 2-torsion computation
\\ ============================================================
\\
\\ chlrs_igusa_formula.gp Rosenhain formula gives match=0 for p=1009.
\\
\\ ROOT CAUSE: the Mobius transform coefficients involve alpha (cube root
\\ of -b in F_{p^3}), which is NOT in F_p. But when computing cross-ratios
\\ T(z) for the remaining branch points, alpha DOES cancel — giving F_p values.
\\
\\ So the lambda_i ARE in F_p, but the Rosenhain curve built from them
\\ is not (2,2)-isogenous to E1 x E2 over F_p. This script diagnoses why.
\\
\\ Run: gp -q chlrs_fp3_howe_correct.gp

default(parisize, 256000000);
default(timer, 0);

print("================================================================");
print("Howe cover diagnosis: why Rosenhain formula gives wrong Jacobian");
print("================================================================");
print("");

p_toy = 1009;
b1 = 11;
d_ns = 2;
b2 = (d_ns^3 * b1) % p_toy;

E1 = ellinit([0, b1], p_toy);
E2 = ellinit([0, b2], p_toy);
n1 = ellcard(E1); t1 = p_toy + 1 - n1;
n2 = ellcard(E2); t2 = p_toy + 1 - n2;

print("E1: y^2 = x^3 + ", b1, "   #E1 = ", n1, "   trace t1 = ", t1);
print("E2: y^2 = x^3 + ", b2, "  #E2 = ", n2, "   trace t2 = ", t2);
print("t2 = -t1? ", t2 == -t1, "   n1*n2 = ", n1*n2, "  (target #Jac)");
print("");

\\ Cube root of unity omega in F_p
om_roots = polrootsmod(x^2 + x + 1, p_toy);
om = Mod(lift(om_roots[1]), p_toy);
d_mod = Mod(d_ns, p_toy);
print("omega (primitive cube root of unity) = ", lift(om));
print("");

\\ ============================================================
\\ Compute the 3 Rosenhain lambda values (same as chlrs_igusa_formula.gp)
\\ ============================================================
\\ Mobius transform T: alpha->0, d*alpha->1, om*alpha->inf
\\ For z = om^2*alpha: T(z) = (om^2-1)(d-om) / (om*(om-1)*(d-1))
\\ For z = d*om*alpha: T(z) = (d*om-1)(d-om) / (om*(d-1)^2)
\\ For z = d*om^2*alpha: T(z) = (d*om^2-1)(d-om) / (om*(d*om-1)*(d-1))
lam1 = (om^2 - 1)*(d_mod - om) / (om*(om - 1)*(d_mod - 1));
lam2 = (d_mod*om - 1)*(d_mod - om) / (om*(d_mod - 1)^2);
lam3 = (d_mod*om^2 - 1)*(d_mod - om) / (om*(d_mod*om - 1)*(d_mod - 1));

l1 = lift(lam1); l2 = lift(lam2); l3 = lift(lam3);
print("Rosenhain lambdas: lambda1=", l1, "  lambda2=", l2, "  lambda3=", l3);
print("(Same as chlrs_igusa_formula.gp variant 1: ", [l1,l2,l3] == [661, 637, 954], ")");
print("");

\\ Frobenius of Rosenhain curve C_R: y^2 = x(x-1)(x-l1)(x-l2)(x-l3)
h_R = Mod(x,p_toy) * Mod(x-1,p_toy) * Mod(x-l1,p_toy) * Mod(x-l2,p_toy) * Mod(x-l3,p_toy);
cp_R = hyperellcharpoly(h_R);
n_R = subst(cp_R, variable(cp_R), 1);
print("Rosenhain C_R: y^2 = x(x-1)(x-", l1, ")(x-", l2, ")(x-", l3, ")");
print("  Frobenius char poly: ", cp_R);
print("  #Jac(C_R) = ", n_R, "   vs n1*n2 = ", n1*n2, "   MATCH: ", n_R == n1*n2);
print("");

\\ Frobenius of NAIVE cover C_N: y^2 = (x^3+b1)(x^3+b2)
h_N = Mod((x^3+b1)*(x^3+b2), p_toy);
cp_N = hyperellcharpoly(h_N);
n_N = subst(cp_N, variable(cp_N), 1);
print("Naive C_N: y^2 = (x^3+", b1, ")(x^3+", b2, ")");
print("  Frobenius char poly: ", cp_N);
print("  #Jac(C_N) = ", n_N, "   vs n1*n2 = ", n1*n2, "   MATCH: ", n_N == n1*n2);
print("");

\\ ============================================================
\\ Try ALL F_p-rational choices of 3 basepoints from the 6 branch locus.
\\ The 6 branch pts {alpha, om*alpha, om^2*alpha, d*alpha, d*om*alpha, d*om^2*alpha}
\\ are in F_{p^3}. Over F_p, no individual branch point is rational.
\\ So any F_p-rational Mobius T cannot fix any branch point to 0, 1, or inf.
\\ A Rosenhain form over F_p requires choosing 3 OF THE 6 branch pts
\\ as {0, 1, inf} — but this forces those 3 to be F_p-rational.
\\ CONCLUSION: no F_p-rational Rosenhain form exists for this curve.
\\
\\ The genus-2 curve Jac^{-1}(E1 x E2) exists over F_p (by descent),
\\ but its model is NOT of Rosenhain form over F_p.
\\ ============================================================

print("---- Can the correct Howe Jac be found via brute force over F_p? ----");
print("");
print("Strategy: scan degree-6 polynomials h(x) over F_p with:");
print("  (1) h smooth (J10 != 0)");
print("  (2) Frobenius char poly = (X^2 - t1*X + p)(X^2 + t1*X + p)");
print("  (3) h has the right structural form (6 roots forming 2 conjugate triples)");
print("");

\\ Expected char poly if Jac ~ E1 x E2:
expected_cp = (x^2 - t1*x + p_toy) * (x^2 + t1*x + p_toy);
print("Expected char poly: ", expected_cp);
expected_n = subst(expected_cp, x, 1);
print("Expected #Jac: ", expected_n);
print("");

\\ Scan h(x) = x^6 + a*x^3 + c for a, c in F_p (the most natural j=0 form)
\\ This parameterizes sextic polynomials with the j=0 symmetry x -> om*x.
\\ Symmetry: h(om*x) = om^6*h(x)/lead = h(x) iff h(x) = x^6 + a*x^3 + c.
\\ We expect the correct cover to have this form (CM structure).

print("Scanning h(x) = x^6 + a*x^3 + c over F_", p_toy, " ...");
found_a = -1; found_c = -1; found_cp = 0;
scan_count = 0;
{
forstep(a = 0, p_toy - 1, 1,
    forstep(c = 1, p_toy - 1, 1,
        scan_count = scan_count + 1;
        h_cand = Mod(x^6 + a*x^3 + c, p_toy);
        \\ Quick check: discriminant / J10 != 0
        disc_check = poldisc(lift(h_cand));
        if (disc_check % p_toy == 0, next);
        cp_cand = hyperellcharpoly(h_cand);
        cp_cand_X = subst(cp_cand, variable(cp_cand), x);
        if (cp_cand_X == expected_cp,
            found_a = a; found_c = c; found_cp = cp_cand_X;
            break(2)
        )
    )
)
}
print("Scanned ", scan_count, " (a,c) pairs.");
if (found_a >= 0,
    print("FOUND: h(x) = x^6 + ", found_a, "*x^3 + ", found_c);
    print("  Frobenius char poly: ", found_cp);
    print("  #Jac = ", subst(found_cp, x, 1));
    ,
    print("NOT FOUND in x^6 + a*x^3 + c form.");
    print("The correct Howe Jac is not of this symmetric form over F_p.");
    print("Requires Mestre reconstruction (needs Sage/Magma).")
);
print("");

\\ ============================================================
\\ Summary
\\ ============================================================
print("================================================================");
print("SUMMARY");
print("================================================================");
print("");
print("1. Rosenhain formula (chlrs_igusa_formula.gp): WRONG for F_p-model.");
print("   Lambdas are in F_p but Jac(C_R) is NOT isogenous to E1 x E2.");
print("   Reason: Mobius T is defined over F_{p^3}; the Rosenhain model");
print("   it produces is F_{p^3}-rational, not an F_p-model of the cover.");
print("");
print("2. Naive cover (x^3+b1)(x^3+b2): also NOT isogenous to E1 x E2 over F_p.");
print("   (Already established in howe_explicit_cover.gp.)");
print("");
print("3. Brute-force scan over x^6 + a*x^3 + c: see result above.");
print("   If not found: the correct F_p-model is NOT of this symmetric form.");
print("");
print("4. Correct Howe cover: needs Mestre algorithm (Sage/Magma).");
print("   THREAD 2 STATUS: BLOCKED (no Sage/Magma available).");
print("");
print("5. ECDLP security (B5): UNAFFECTED.");
print("   Any genus-2 cover C of secp256k1 over F_p has");
print("   HCDLP cost >= sqrt(#Jac(C)(F_p)) >= sqrt(p^2) >> sqrt(n_secp).");
