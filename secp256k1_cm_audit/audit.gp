\\ ============================================================
\\ secp256k1 CM / isogeny-graph structural audit (PARI/GP)
\\ ============================================================
\\ Run: gp -q audit.gp > audit_output.txt
\\
\\ NB: PARI/GP terminates statements at EOL outside { }.  All
\\ multi-line for/if bodies below are wrapped in { }.

default(parisize, 256000000);
default(timer, 0);

\\ secp256k1 parameters
p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
n = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;
t_actual = p + 1 - n;

print("================================================================");
print("secp256k1 CM / isogeny-graph structural audit");
print("================================================================");
print("p = ", p);
print("n = ", n);
print("t (Frobenius trace) = ", t_actual);
print("|t| bit-length      = ", #binary(abs(t_actual)));
print("t mod 2  = ", Mod(t_actual, 2));
print("t mod 3  = ", Mod(t_actual, 3));
print("");

print("---- Frobenius discriminant ----");
D_pi = t_actual^2 - 4*p;
print("t^2 - 4p = ", D_pi);
print("sign = ", sign(D_pi));
print("(t^2 - 4p) / -3 = ", D_pi / -3);
print("is square?  ", issquare(D_pi / -3));
print("sqrt        = ", sqrtint(D_pi / -3));
print("");

print("---- Cornacchia: 4p = L^2 + 27 M^2 ----");
LM = qfbcornacchia(27, 4*p);
L = LM[1];
M = LM[2];
print("L (raw) = ", L);
print("M       = ", M);
print("check 4p = L^2 + 27 M^2: ", L^2 + 27*M^2 == 4*p);
print("L mod 3 = ", Mod(L, 3));
L = if(Mod(L, 3) == 2, -L, L);
print("L (normalised, L ≡ 1 mod 3) = ", L);
print("L odd? ", L % 2 != 0, "    M odd? ", M % 2 != 0);
print("");

print("---- Six possible Frobenius traces ----");
traces = [L, -L, (L+3*M)/2, -((L+3*M)/2), (L-3*M)/2, -((L-3*M)/2)];
{
for (i = 1, 6,
    print("  t_", i, " = ", traces[i])
);
}
print("Sum of six traces (must be 0): ", vecsum(traces));
print("");

print("---- Which twist is secp256k1 ----");
matched = 0;
{
for (i = 1, 6,
    if (traces[i] == t_actual, matched = i)
);
}
if (matched == 0, print("WARNING: secp256k1's trace not found in predicted set"));
if (matched != 0, print("secp256k1 = twist #", matched, "  (trace = ", traces[matched], ")"));
print("");

print("---- Six twist orders ----");
orders = vector(6, i, p + 1 - traces[i]);
{
for (i = 1, 6,
    print("  #E_", i, " = ", orders[i], if(i == matched, " <-- secp256k1", ""))
);
}
print("Sum of six orders = 6(p+1)? ", vecsum(orders) == 6*(p+1));
print("");

\\ ------------------------------------------------------------
\\ Partial factorisation of twist orders.
\\ Use factor(N, 10^7) for trial-division to 10^7, then ECM the
\\ residue with B1 = 10^4 (factorint flag 4).  Report cofactor.
\\ ------------------------------------------------------------
print("---- Twist-order partial factorisation ----");
print("(trial division to 10^7, then short ECM; primality probabilistic)");
print("");
{
for (i = 1, 6,
    local(ord, fac, smooth_part, max_pf, cof, k, q, e, nrows);
    ord = orders[i];
    print("Twist #", i, if(i == matched, " (secp256k1)", ""), ":  N = ", ord);
    fac = factor(ord, 10^7);
    smooth_part = 1;
    max_pf = 1;
    nrows = matsize(fac)[1];
    for (k = 1, nrows,
        q = fac[k, 1];
        e = fac[k, 2];
        if (q < 10^7,
            print("    factor ", q, "^", e);
            smooth_part = smooth_part * q^e;
            if (q > max_pf, max_pf = q)
        )
    );
    cof = ord / smooth_part;
    print("  smooth-part (primes < 10^7): ", smooth_part);
    print("  cofactor: ", cof);
    print("  cofactor bit-len: ", #binary(cof));
    print("  cofactor pseudoprime? ", ispseudoprime(cof));
    if (max_pf > 1,
        print("  largest small prime factor: ", max_pf,
              "    (Pohlig-Hellman on it: ~", sqrtint(max_pf), " ops)")
    );
    print("")
);
}

\\ ------------------------------------------------------------
\\ Frobenius in O_K = Z[w], w^2 + w + 1 = 0
\\ ------------------------------------------------------------
print("---- Frobenius pi in O_K = Z[w] ----");
K = nfinit(y^2 + y + 1);
print("K = Q(w), disc(K) = ", K.disc);
fact_p = idealprimedec(K, p);
print("primes above p in K: ", #fact_p);
{
for (i = 1, #fact_p,
    local(pp);
    pp = fact_p[i];
    print("  prime ", i, "  norm=", idealnorm(K, pp), " e=", pp.e, " f=", pp.f)
);
}
print("");

A_inner = (4*p - t_actual^2) / 3;
print("(4p - t^2)/3 = ", A_inner, "   is square? ", issquare(A_inner));
sqrt_val = sqrtint(A_inner);
print("sqrt        = ", sqrt_val);
a_plus  = (t_actual + sqrt_val)/2;
a_minus = (t_actual - sqrt_val)/2;
b_plus  = 2*a_plus  - t_actual;
b_minus = 2*a_minus - t_actual;
print("Two Frobenius candidates pi = a + b*w:");
print("  pi_+ = ", a_plus,  " + (", b_plus,  ") w");
print("    norm = ", a_plus^2  - a_plus*b_plus   + b_plus^2);
print("    trace = ", 2*a_plus  - b_plus);
print("  pi_- = ", a_minus, " + (", b_minus, ") w   [complex conjugate]");
print("    norm = ", a_minus^2 - a_minus*b_minus + b_minus^2);
print("    trace = ", 2*a_minus - b_minus);
norm_pi_minus_1 = (a_plus - 1)^2 - (a_plus - 1)*b_plus + b_plus^2;
print("Norm(pi - 1)     = ", norm_pi_minus_1, "    == n?  ", norm_pi_minus_1 == n);
norm_pibar_minus_1 = (a_minus - 1)^2 - (a_minus - 1)*b_minus + b_minus^2;
print("Norm(pi-bar - 1) = ", norm_pibar_minus_1);
{
mt = 0;
for (i = 1, 6, if (orders[i] == norm_pibar_minus_1, mt = i));
print("    matches twist #", mt);
}
print("");

\\ ------------------------------------------------------------
\\ Small-ℓ isogeny graph from j = 0
\\ ------------------------------------------------------------
print("---- Small-ℓ isogeny neighbours of j = 0 (mod p) ----");
print("");

isogeny_audit(ell) = {
    my(phi, phi0, fact, kron, split_type, num_roots, self_loops, nrows,
       deg_factor, mult, rt, rt_lift, k, distinct_neighbours);
    print("---- ℓ = ", ell, " ----");
    kron = if(ell == 3, 0, kronecker(-3, ell));
    split_type = if(ell == 3, "ramified", if(kron == 1, "split", "inert"));
    print("  (-3/ℓ) = ", kron, "  ⇒  ℓ is ", split_type, " in Z[ω]");

    phi = polmodular(ell);
    phi0 = subst(phi, x, 0) * Mod(1, p);
    print("  deg Φ_ℓ(0, Y) = ", poldegree(phi0), "  (expected ", ell+1, ")");

    fact = factor(phi0);
    nrows = matsize(fact)[1];
    num_roots = 0;
    self_loops = 0;
    distinct_neighbours = 0;
    print("  factorisation of Φ_ℓ(0, Y) over F_p:");
    for (k = 1, nrows,
        deg_factor = poldegree(fact[k, 1]);
        mult = fact[k, 2];
        if (deg_factor == 1,
            rt = -polcoeff(fact[k, 1], 0) / polcoeff(fact[k, 1], 1);
            rt_lift = lift(rt);
            num_roots = num_roots + mult;
            distinct_neighbours = distinct_neighbours + 1;
            if (rt_lift == 0,
                self_loops = self_loops + mult;
                print("    deg 1 × ", mult, "    ↳ j' = 0 (self-loop)"),
                print("    deg 1 × ", mult, "    ↳ j' = ", rt_lift)
            ),
            print("    deg ", deg_factor, " × ", mult, "    [not F_p-rational]")
        )
    );
    print("  F_p-rational neighbour j-invariants (w/mult): ", num_roots);
    print("  distinct F_p-rational j-neighbours:           ", distinct_neighbours);
    print("  self-loops to j = 0:                          ", self_loops);
};

for_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47];
{
for (i = 1, #for_primes,
    isogeny_audit(for_primes[i]);
    print("")
);
}

\\ ------------------------------------------------------------
\\ Disc-(-3 ℓ²) class numbers for the descending neighbours
\\ ------------------------------------------------------------
print("---- Class numbers h(-3·ℓ²) of descending orders ----");
print("(non-maximal order Z + ℓZ[ω] has discriminant -3·ℓ²)");
{
for (i = 1, #for_primes,
    local(ell, disc, qf, h);
    ell = for_primes[i];
    disc = -3 * ell^2;
    qf = qfbclassno(disc);
    h = qf;
    print("  ℓ = ", ell, "    disc = ", disc, "    h(disc) = ", h)
);
}

print("");
print("================================================================");
print("Audit done.");
print("================================================================");
