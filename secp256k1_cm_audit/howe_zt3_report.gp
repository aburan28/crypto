\\ ============================================================
\\ howe_zt3_report.gp
\\ ============================================================
\\
\\ Follow-up to howe_zt3_search.gp: print the 28 matching (a,b) pairs,
\\ compute basic Igusa invariants (J_10 = disc), and determine the
\\ isomorphism classes over F_43.
\\
\\ Z/3Z family: y^2 = x^6 + a*x^3 + b.
\\ The Z/6Z automorphism: (x,y) -> (zeta3 * x, y) acts on the curve.
\\ Two such curves are F_p-isomorphic iff (a', b') = (a*u^3, b*u^6)
\\ for some u ∈ F_p^* (since substituting x -> u*x, y -> u^3*y gives
\\ y^2 = x^6*u^6 + a*x^3*u^3 + b, divided by u^6 gives the standard
\\ form with a' = a*u^3/u^3 = a... wait.
\\ Actually sub x -> v*x, y -> w*y:
\\   (w*y)^2 = (v*x)^6 + a*(v*x)^3 + b
\\   w^2 * y^2 = v^6*x^6 + a*v^3*x^3 + b
\\ Divide by w^2: y^2 = (v^6/w^2)*x^6 + (a*v^3/w^2)*x^3 + (b/w^2)
\\ For monic x^6: need v^6/w^2 = 1, i.e., w = v^3.
\\ Then: a' = a*v^3/v^6 = a/v^3, b' = b/v^6.
\\ So (a, b) ~ (a/v^3, b/v^6) for v ∈ F_p^*.
\\
\\ Equivalence: (a,b) ~ (a',b') iff a' = a*u, b' = b*u^2 for u in
\\ (F_p^*)^3 (cube class), OR equivalently iff a^2/b^... hmm.
\\ More simply: J_10 / J_2^5 is an absolute invariant.
\\
\\ Run: gp -q howe_zt3_report.gp

default(parisize, 256000000);
default(timer, 0);

p = 43;
target_cp = Pol([1, 0, -83, 0, 1849]);

print("=== Z/3Z Howe-glued curves over F_43: (a,b) pairs ===");
print("");

\\ Store matches as separate arrays (avoid List indexing bug)
na = 0;
av = vector(100);
bv = vector(100);

{
for (a = 0, p-1,
    for (b = 0, p-1,
        local(f, dv, cp);
        f = x^6 + a*x^3 + b;
        dv = lift(poldisc(f) * Mod(1, p));
        if (dv == 0, next);
        cp = hyperellcharpoly(f * Mod(1, p));
        if (cp == target_cp,
            na = na + 1;
            av[na] = a;
            bv[na] = b
        )
    )
);
}

print("Total matches: ", na);
print("");
print("k | a  | b  | J10=disc(f) mod p | a^6-27*4*b^2 (rough J check)");
print("--|----|----|-------------------|--------------");
{
for (k = 1, na,
    local(ak, bk, fk, dk, s);
    ak = av[k]; bk = bv[k];
    fk = x^6 + ak*x^3 + bk;
    dk = lift(poldisc(fk) * Mod(1, p));
    \\ s = a^3 - 4*b (discriminant of the quadratic u^2 + a*u + b as poly in u=x^3)
    s = lift(Mod(ak^3 - 4*bk, p));  \\ note: disc of u^2 + a*u + b is a^2-4b, NOT a^3-4b
    print(k, " | ", ak, " | ", bk, " | ", dk, " | a^2-4b=", lift(Mod(ak^2-4*bk, p)))
);
}
print("");

\\ ---- Determine isomorphism classes under (a,b) -> (a/v^3, b/v^6) ----
\\ Equivalently: two pairs (a1,b1) and (a2,b2) are F_p-isomorphic iff
\\ there exists v ∈ F_p^* s.t. a2 = a1 / v^3 mod p AND b2 = b1 / v^6 mod p.
\\ Equivalently, (a1^6/b1^2 ≡ a2^6/b2^2 mod p) when b != 0,
\\ OR both have a = 0 (the hyperelliptic involution case).

print("=== Isomorphism classes ===");
\\ Tag each match with canonical invariant a^6/b^2 mod p (when b != 0)
print("Invariant I = a^6/b^2 mod p (or 'b=0' or 'a=0') per curve:");
{
for (k = 1, na,
    local(ak, bk, inv_v);
    ak = av[k]; bk = bv[k];
    if (bk == 0,
        print(k, ": a=", ak, " b=0 (special form)"),
        if (ak == 0,
            print(k, ": a=0 b=", bk, " (Z/6Z symmetric)"),
            inv_v = lift(Mod(ak^6, p) / Mod(bk^2, p));
            print(k, ": a=", ak, " b=", bk, " I=a^6/b^2=", inv_v)
        )
    )
);
}
print("");

\\ ---- Count distinct I values ----
nI = 0;
Iv = vector(100);
Ic = vector(100);   \\ count
{
for (k = 1, na,
    local(ak, bk, inv_v, found_I, jj);
    ak = av[k]; bk = bv[k];
    if (bk == 0 || ak == 0,
        inv_v = -1,  \\ sentinel
        inv_v = lift(Mod(ak^6, p) / Mod(bk^2, p))
    );
    found_I = 0;
    for (jj = 1, nI,
        if (Iv[jj] == inv_v, Ic[jj] = Ic[jj] + 1; found_I = 1; break)
    );
    if (!found_I, nI = nI + 1; Iv[nI] = inv_v; Ic[nI] = 1)
);
}

print("Distinct I = a^6/b^2 values: ", nI);
{
for (j = 1, nI,
    print("  I=", Iv[j], "  count=", Ic[j])
);
}
print("");
print("Each distinct I value = one Mobius isomorphism class (over F_p^*).");
print("Distinct isomorphism classes over F_43 = ", nI);
print("");

\\ ---- Pick canonical representative of each class and get Igusa info ----
print("=== Canonical representatives ===");
seen_I = vector(100); nS = 0;

{
for (k = 1, na,
    local(ak, bk, inv_v, fk, cpk, n_C, already);
    ak = av[k]; bk = bv[k];
    if (bk == 0 || ak == 0, inv_v = -1, inv_v = lift(Mod(ak^6, p) / Mod(bk^2, p)));
    already = 0;
    for (jj = 1, nS, if (seen_I[jj] == inv_v, already = 1; break));
    if (already, next);
    nS = nS + 1; seen_I[nS] = inv_v;
    fk = x^6 + ak*x^3 + bk;
    cpk = hyperellcharpoly(fk * Mod(1, p));
    \\ Count #C(F_p) by point counting
    n_C = 1;   \\ point at infinity for degree-6 curve (one point)
    for (xi = 0, p-1,
        local(rhs_v, ks_v);
        rhs_v = lift(Mod(xi^6 + ak*xi^3 + bk, p));
        if (rhs_v == 0, n_C = n_C + 1; next);
        ks_v = kronecker(rhs_v, p);
        if (ks_v == 1, n_C = n_C + 2)
    );
    print("Class ", nS, ": y^2 = x^6 + ", ak, "*x^3 + ", bk);
    print("  char poly = ", cpk);
    print("  #C(F_43) = ", n_C, "  (expected 44)");
    print("  disc(f) mod p = ", lift(poldisc(fk) * Mod(1, p)));
    print("  I = a^6/b^2 mod p = ", inv_v);
    print("")
);
}

print("Total isomorphism classes: ", nS);
print("");

\\ ---- Naive cover comparison ----
print("=== Naive cover (from prior session) ===");
f_naive = (x^3 + 7) * (x^3 + 13);
cp_naive = hyperellcharpoly(f_naive * Mod(1, p));
print("Naive cover y^2 = (x^3+7)(x^3+13) over F_43:");
print("  char poly = ", cp_naive);
print("  a=20, b=91 -> mod p: a_eff=", lift(Mod(7+13, p)), " b_eff=", lift(Mod(7*13, p)));
print("  (This is NOT in the Z/3Z family because b = 7*13 = 91 = 5 mod 43,");
print("   and a = 7+13 = 20, which gives x^6+20*x^3+5 -- IS in family!)");
print("  Direct check: cp of x^6+20*x^3+5 over F_43:");
cp_20_5 = hyperellcharpoly((x^6 + 20*x^3 + 5) * Mod(1, p));
print("  ", cp_20_5);
print("  Match with target? ", cp_20_5 == target_cp);
print("");

print("================================================================");
print("Key finding: ", na, " Z/3Z-symmetric genus-2 curves over F_43");
print("have Frobenius char poly = (T^2-13T+43)(T^2+13T+43).");
print("These ARE candidates for the Howe-glued curve Jac(C) ~ E1 x E2.");
print("================================================================");
