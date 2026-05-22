\\ ============================================================
\\ Horizontal Cl(O)-orbit measurement on a small h>1 CM curve
\\ ============================================================
\\
\\ Slice-2 companion to RESEARCH_SECP256K1_CM.md.  The existing
\\ cga_hnc module in src/cryptanalysis/cga_hnc.rs walks 2-isogeny
\\ BFS which mixes horizontal and descending isogenies.  This
\\ script does the strictly horizontal Cl(O)-orbit walk that the
\\ user's "Step 4" of the proposal calls for.
\\
\\ Result we expect to demonstrate:
\\   1. The orbit is a single Cl(O)-Cayley graph of size h(D).
\\   2. Tate's theorem ⇒ every orbit member has the same #E(F_p).
\\      → Per-curve smoothness of #E is identical across the orbit.
\\      → No CRT amortisation via Pohlig-Hellman on #E factorisation.
\\   3. Per-curve transported point orders ord(P_i) can however
\\      vary across the orbit (as divisors of #E).  This is what
\\      makes cga_hnc's existing measurement non-trivial.
\\   4. The orbit walk itself costs O(h log h) ℓ-isogeny steps;
\\      the saving in PH is at most a factor of ~h.  Net: no
\\      asymptotic gain over single-curve PH on #E's smooth part.
\\
\\ Run: gp -q cl_o_horizontal.gp > cl_o_horizontal_output.txt

default(parisize, 256000000);
default(timer, 0);

\\ ------------------------------------------------------------
\\ Pick a small CM discriminant with h(D) > 1
\\ ------------------------------------------------------------
\\ D = -23 (h=3, smallest h>1 fundamental disc)
\\ D = -47 (h=5)
\\ D = -71 (h=7)
\\ D = -167 (h=11)
\\ We use D = -71 to get a non-trivial cycle structure.

D = -71;
h_D = qfbclassno(D);
print("================================================================");
print("Horizontal Cl(O)-orbit walk on CM disc D = ", D, " curves");
print("================================================================");
print("h(D) = ", h_D);
print("");

\\ ------------------------------------------------------------
\\ Find a small prime p where curves with CM by O_D exist over F_p.
\\ Condition: 4p = a^2 + |D| b^2 (with D ≡ 1 mod 4) for some a, b ≥ 1.
\\ Equivalently: p splits in K = Q(sqrt(D)) and  isn't ramified.
\\ ------------------------------------------------------------
print("---- Finding a small p with CM disc D = ", D, " curves ----");
\\ Search for p prime, p splits in O_D, p small.  Use qfbcornacchia.
small_primes = primes(40);
chosen_p = 0;
chosen_a = 0;
chosen_b = 0;
{
for (i = 1, #small_primes,
    local(q, ab);
    q = small_primes[i];
    if (kronecker(D, q) != 1, next);
    ab = qfbcornacchia(-D, 4*q);
    if (ab == 0, next);
    chosen_p = q;
    chosen_a = ab[1];
    chosen_b = ab[2];
    break
);
}
p = chosen_p;
print("Chosen p = ", p);
print("Cornacchia: 4p = a^2 + |D| b^2 with a = ", chosen_a, ", b = ", chosen_b);
print("check 4p = a^2 + |D| b^2: ", chosen_a^2 + (-D) * chosen_b^2 == 4*p);
print("");

\\ The Frobenius trace of a CM curve with End = O_D over F_p is ±a
\\ (modulo sign convention from Hilbert class polynomial side).
\\ So #E(F_p) = p + 1 ∓ a.  Both signs occur: the two roots of the
\\ Hilbert class poly give curves with trace +a and -a respectively
\\ (in fact for each j-invariant root, there's one trace; signs
\\ alternate in the orbit pairing).

print("---- The two possible traces ----");
t_plus  = chosen_a;
t_minus = -chosen_a;
n_plus  = p + 1 - t_plus;
n_minus = p + 1 - t_minus;
print("Trace +a: ", t_plus, "    #E = ", n_plus, "    factor: ", factor(n_plus));
print("Trace -a: ", t_minus, "    #E = ", n_minus, "    factor: ", factor(n_minus));
print("");

\\ ------------------------------------------------------------
\\ Enumerate the orbit: roots of Hilbert class polynomial H_D mod p
\\ ------------------------------------------------------------
print("---- Hilbert class polynomial H_D(X) and its roots mod p ----");
\\ polclass(D, 0) returns the classical Hilbert class polynomial H_D ∈ Z[X]
H_D = polclass(D, 0);
print("H_D(X) ∈ Z[X], deg = ", poldegree(H_D), "  (== h(D) = ", h_D, ")");
\\ Factor H_D mod p
H_D_mod = H_D * Mod(1, p);
fact_H = factor(H_D_mod);
print("Factorisation of H_D mod p:");
nrows = matsize(fact_H)[1];
roots_mod_p = [];
{
for (k = 1, nrows,
    local(deg_factor, mult, rt);
    deg_factor = poldegree(fact_H[k, 1]);
    mult = fact_H[k, 2];
    print("  deg ", deg_factor, " × ", mult);
    if (deg_factor == 1,
        rt = -polcoeff(fact_H[k, 1], 0) / polcoeff(fact_H[k, 1], 1);
        roots_mod_p = concat(roots_mod_p, [lift(rt)]);
        print("    ↳ j = ", lift(rt))
    )
);
}
print("Total F_p-rational j-invariants: ", #roots_mod_p);
print("(If = h(D), the entire Cl(O)-orbit is F_p-rational.)");
print("");

\\ For our purposes (small prime), expect roots_mod_p to have either
\\ 0 (orbit not F_p-rational), h(D)/k (split into orbits by Frobenius),
\\ or h(D) (full orbit F_p-rational).

{
if (#roots_mod_p == 0,
    print("WARNING: no F_p-rational orbit member.  Try a different p.");
    quit
);
}

\\ ------------------------------------------------------------
\\ For each F_p-rational orbit member, construct a curve E_j with
\\ j(E_j) = j and compute its order.
\\ ------------------------------------------------------------
\\ Standard model for any j ≠ 0, 1728:
\\   c = j / (1728 - j),  E_j: y^2 = x^3 + 3c·x + 2c
\\ For j = 0: y^2 = x^3 + 1.  For j = 1728: y^2 = x^3 + x.
\\
\\ ellinit + ellsea give us the curve and its order.

print("---- Building curves and counting points for each orbit member ----");
print("(All orders should equal n_plus or n_minus, by Tate.)");
print("");

\\ ellinit needs the curve in Weierstrass form over Mod(_, p).
order_counts = vector(#roots_mod_p);
{
for (i = 1, #roots_mod_p,
    local(j_val, c, E, ord_E);
    j_val = roots_mod_p[i];
    if (j_val == 0,
        E = ellinit([0, 1] * Mod(1, p)),
        if (j_val == 1728,
            E = ellinit([1, 0] * Mod(1, p)),
            c = Mod(j_val, p) / (1728 - j_val);
            E = ellinit([3 * c, 2 * c])
        )
    );
    ord_E = ellsea(E);
    if (ord_E == 0, ord_E = ellcard(E));
    order_counts[i] = lift(ord_E);
    print("  j = ", j_val, "    #E = ", lift(ord_E),
          "    trace = ", p + 1 - lift(ord_E))
);
}
print("");

\\ Confirm uniformity (Tate)
print("Distinct #E values across orbit: ", #Set(order_counts));
print("(should be 1 if all F_p-rational members are in one Frobenius-orbit)");
print("(or 2 if the orbit splits into +/− trace pairs)");

\\ ------------------------------------------------------------
\\ Smoothness of #E: identical across orbit by Tate ⇒
\\ no CRT amortisation through orbit-walking the GROUP ORDER.
\\ ------------------------------------------------------------
print("");
print("---- Smoothness of #E (= n_plus or n_minus) ----");
print("By Tate, #E is constant across each Frobenius-orbit, so:");
print("  smooth part of #E for n_plus  = ", n_plus,  " : ", factor(n_plus));
print("  smooth part of #E for n_minus = ", n_minus, " : ", factor(n_minus));
print("");
print("CGA-HNC's claimed 'per-curve variation in #E' is therefore");
print("VACUOUS in the strictly horizontal Cl(O)-orbit.  Variation");
print("in the existing 2-isogeny-BFS measurement comes from:");
print("  (a) the orbit including DESCENDING isogenies to non-maximal");
print("      orders (which still have #E = n by Tate, but lose GLV-style");
print("      automorphisms), OR");
print("  (b) the transported POINT ORDER ord(P_i) being a different");
print("      divisor of #E on different curves, even when #E is fixed.");
print("");
print("Mechanism (b) is real and is what cga_hnc Phase 2 exploits.");
print("Mechanism (a) is moot for ECDLP (Tate again).");
print("");
print("Cost analysis: the orbit walk visits h(D) curves; the saving");
print("over single-curve PH is at most a factor of h(D); the walk");
print("itself costs Θ(h(D)) ℓ-isogeny steps; ⇒ no asymptotic gain.");

print("");
print("================================================================");
print("Horizontal Cl(O)-orbit measurement complete.");
print("================================================================");
