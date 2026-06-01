\\ ============================================================
\\ howe_zt3_search.gp
\\ ============================================================
\\
\\ Search for the Howe-glued curve over the proxy prime p=43
\\ in the Z/3Z-symmetric family y^2 = x^6 + a*x^3 + b.
\\
\\ Motivation: when x^3+7 is irreducible over F_p (as for secp256k1
\\ and the proxy p=43), the 2-torsion of E1 sits in F_{p^3}.
\\ The Frobenius acts as a 3-cycle on the 2-torsion, so the
\\ Howe-glued curve C inherits the Z/3Z automorphism (x,y)->(ζ_3 x, y).
\\ All such curves over F_p have the form y^2 = x^6 + a*x^3 + b.
\\
\\ The 2026-05-31 exhaustive search over 120 Rosenhain orderings
\\ (individual branch-point cross-ratios) found NO match.
\\ That search is not equivalent to this one: Rosenhain form requires
\\ the sextic to split into 6 distinct F_{p^3} linear factors arranged
\\ as individual cross-ratios; the Z/3Z form has pairs of conjugate
\\ cubics as the natural unit.
\\
\\ Target (proxy prime p=43):
\\   E1: y^2 = x^3 + 7,  #E1 = 31, trace = 13
\\   E2: quadratic twist of E1, #E2 = 57, trace = -13
\\   Howe char poly: (T^2-13T+43)(T^2+13T+43) = T^4-83T^2+1849
\\
\\ Run: gp -q howe_zt3_search.gp

default(parisize, 256000000);
default(timer, 0);

print("================================================================");
print("Z/3Z-symmetric Howe-glued curve search over proxy prime p=43");
print("================================================================");
print("");

p = 43;

\\ ---- Verify E1 and its quadratic twist (from prior audit) ----
print("---- Verify E1/E2 traces over F_p ----");
E1 = ellinit([0, 7], p);
n1 = ellcard(E1);
t1 = p + 1 - n1;
print("E1: y^2 = x^3 + 7  over F_43:  #E1 = ", n1, "   trace = ", t1);

\\ Find quadratic non-residue to construct twist
qnr = 2;
{
while (kronecker(qnr, p) != -1, qnr = qnr + 1);
}
print("Quadratic non-residue mod p = ", qnr);

\\ For y^2 = x^3 + c, quadratic twist is y^2 = x^3 + c * qnr^3 mod p
c2 = lift(Mod(7 * qnr^3, p));
E2 = ellinit([0, c2], p);
n2 = ellcard(E2);
t2 = p + 1 - n2;
print("E2 (quad twist): y^2 = x^3 + ", c2, "  over F_43:  #E2 = ", n2, "   trace = ", t2);
print("H1: #E1 ≠ #E2? ", n1 != n2, " | H3: gcd(#E1,#E2) = ", gcd(n1, n2));
print("");

if (t2 != -t1, print("WARNING: twist trace not negated; adjust qnr"));
t_target = t1;  \\ positive trace, should be 13
if (t1 < 0, t_target = -t1);

\\ ---- Target Frobenius char poly ----
\\ (T^2 - t_target*T + p)(T^2 + t_target*T + p) = T^4 - (t^2 - 2p)T^2 + p^2
a2_coeff = -(t_target^2 - 2*p);  \\ = -(169 - 86) = -83
target_cp_Z = Pol([1, 0, a2_coeff, 0, p^2]);
print("Target char poly (over Z): ", target_cp_Z);
print("Expected a2 = ", a2_coeff, "   p^2 = ", p^2);
print("");

\\ ---- Search over Z/3Z family: y^2 = x^6 + a*x^3 + b ----
print("---- Exhaustive search: y^2 = x^6 + a*x^3 + b, a,b in F_", p, " ----");
print("");

found_list = List();
total_checked = 0;
total_smooth = 0;

{
for (a = 0, p-1,
    for (b = 0, p-1,
        local(f, disc_mod, cp);
        f = x^6 + a*x^3 + b;
        disc_mod = lift(poldisc(f) * Mod(1, p));
        total_checked = total_checked + 1;
        if (disc_mod == 0, next);
        total_smooth = total_smooth + 1;
        cp = hyperellcharpoly(f * Mod(1, p));
        if (cp == target_cp_Z,
            listput(found_list, [a, b, cp])
        )
    )
);
}

print("Checked: ", total_checked, " pairs (a,b)");
print("Smooth curves: ", total_smooth);
print("Matches (char poly = target): ", #found_list);
print("");

if (#found_list == 0,
    print("NEGATIVE RESULT: No Z/3Z-symmetric curve over F_", p,
          " has the target Howe char poly.");
    print("Implication: the Howe-glued curve for E1/E2 at p=43 does NOT");
    print("lie in the family y^2 = x^6 + a*x^3 + b.");
    print("Next: extend search to general sextic or use Mestre algorithm."),

    print("POSITIVE RESULT: Found ", #found_list, " curve(s):");
    for (k = 1, #found_list,
        local(entry, ak, bk, cpk, disc_val, J2v, J4v, J6v, J10v, fk);
        entry = found_list[k];
        ak = entry[1]; bk = entry[2]; cpk = entry[3];
        fk = x^6 + ak*x^3 + bk;
        print("  Curve #", k, ": y^2 = x^6 + ", ak, "*x^3 + ", bk);
        print("    char poly = ", cpk);
        disc_val = lift(poldisc(fk) * Mod(1, p));
        print("    disc mod p = ", disc_val);
        \\ Compute #C(F_p) for sanity
        local(cnt_pts);
        cnt_pts = 0;
        for (xi = 0, p-1,
            local(rhs, yi);
            rhs = lift(Mod(xi^6 + ak*xi^3 + bk, p));
            if (rhs == 0, cnt_pts = cnt_pts + 1; next);
            \\ Count y s.t. y^2 = rhs mod p
            if (kronecker(rhs, p) == 1, cnt_pts = cnt_pts + 2)
        );
        cnt_pts = cnt_pts + 1;  \\ point at infinity
        print("    #C(F_p) = ", cnt_pts, "  (expected: ", 1 + p - 0, " = 44)");
        print("")
    )
);
print("");

\\ ---- Extended search: include all a's giving E1-trace on sub-cubics ----
print("================================================================");
print("Extended: for which (a,b) do sub-cubics give E1/E2 traces?");
print("(cross-check with naive cover search from prior session)");
print("================================================================");
print("");

\\ Find all a s.t. ellcard(y^2 = x^3 + a) gives trace +t_target or -t_target
e1_values = List();
e2_values = List();
{
for (c = 1, p-1,
    local(Ec, nc, tc);
    Ec = ellinit([0, c], p);
    nc = ellcard(Ec);
    tc = p + 1 - nc;
    if (tc == t_target, listput(e1_values, c));
    if (tc == -t_target, listput(e2_values, c))
);
}
print("Values c with ellcharpoly E_c trace = +", t_target, " (E1-type): ", Vec(e1_values));
print("Values c with ellcharpoly E_c trace = -", t_target, " (E2-type): ", Vec(e2_values));
print("");

print("Naive-cover pairs (c1,c2) with y^2=(x^3+c1)(x^3+c2), c1 in E1-type, c2 in E2-type:");
naive_found = 0;
{
for (i1 = 1, #e1_values,
    for (i2 = 1, #e2_values,
        local(c1, c2, fv, dv, cpv);
        c1 = e1_values[i1]; c2 = e2_values[i2];
        fv = (x^3 + c1) * (x^3 + c2);
        dv = lift(poldisc(fv) * Mod(1, p));
        if (dv == 0, next);
        cpv = hyperellcharpoly(fv * Mod(1, p));
        naive_found = naive_found + 1;
        print("  (c1=", c1, ", c2=", c2, "): charpoly = ", cpv,
              if(cpv == target_cp_Z, "  <-- MATCH", ""))
    )
);
}
print("Total naive-cover pairs checked: ", naive_found);
print("");

\\ ---- Broader search: y^2 = x^5 + a*x^3 + b*x^2 + c*x + d (4-param family) ----
print("================================================================");
print("Targeted degree-5 search (Mobius: fix leading term, 4 free coeffs)");
print("Restrict to a,b in small range, d fixed = 0..42, c fixed:");
print("================================================================");
print("");
print("Running a,b in 0..42, c in {0,1,2}, d in {0,1,2} to sample...");

deg5_matches = 0;
deg5_checked = 0;
{
for (a = 0, p-1,
    for (b = 0, p-1,
        for (c = 0, 2,
            for (d = 0, 2,
                local(f5, d5, cp5);
                f5 = x^5 + a*x^3 + b*x^2 + c*x + d;
                d5 = lift(poldisc(f5) * Mod(1, p));
                deg5_checked = deg5_checked + 1;
                if (d5 == 0, next);
                cp5 = hyperellcharpoly(f5 * Mod(1, p));
                if (cp5 == target_cp_Z,
                    deg5_matches = deg5_matches + 1;
                    print("MATCH deg5: y^2 = x^5 + ", a, "*x^3 + ", b,
                          "*x^2 + ", c, "*x + ", d, "  cp=", cp5)
                )
            )
        )
    )
);
}
print("Degree-5 sample: checked ", deg5_checked, " curves, found ", deg5_matches, " matches.");
print("");

print("================================================================");
print("Summary");
print("================================================================");
print("Proxy prime p = 43");
print("E1 trace = +", t_target, "   E2 trace = -", t_target);
print("Target char poly = ", target_cp_Z);
print("Z/3Z family (x^6+ax^3+b): matches = ", #found_list);
print("Naive-cover family (E1xE2): all pairs checked above");
print("Degree-5 sample: matches = ", deg5_matches);
print("");
print("If #found_list = 0 AND deg5_matches = 0:");
print("  The Howe-glued curve likely requires a full sextic equation");
print("  y^2 = x^6 + a5 x^5 + a4 x^4 + a3 x^3 + a2 x^2 + a1 x + a0");
print("  (5 free params mod Mobius, search space 43^5 = 147M).");
print("  Recommended: use Mestre algorithm (port conic+sextic resolution).");
print("================================================================");
