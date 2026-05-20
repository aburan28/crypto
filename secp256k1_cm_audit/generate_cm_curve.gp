\\ ============================================================
\\ Generate a fresh CM curve for the cargo-test audit.
\\ ============================================================
\\
\\ Expects: D_input (negative fundamental discriminant), bit_target
\\ (approximate desired bit length of p).
\\
\\ Output (KEY=VALUE format consumed by Rust harness):
\\   FRESH_P=<hex>
\\   FRESH_N=<hex>
\\   FRESH_A=<hex>
\\   FRESH_B=<hex>
\\   FRESH_DISC=<int>
\\   FRESH_CLASS_NUMBER=<int>
\\   GENERATED=<0|1>
\\
\\ Algorithm:
\\   1. Compute Hilbert class polynomial H_D(X) ∈ Z[X].
\\   2. Search for prime p with p ≈ 2^bit_target, kronecker(D, p) = 1
\\      (so D-CM curves exist over F_p), and H_D mod p has an F_p-root.
\\   3. Take j_0 = first F_p-root of H_D mod p; build E with j(E) = j_0.
\\   4. Compute #E(F_p) = ellcard; report parameters.

default(parisize, 256000000);
default(timer, 0);

D = D_input;
print_h = qfbclassno(D);

\\ Build Hilbert class polynomial
H_D = polclass(D, 0);

p_search_low  = 2^bit_target;
p_search_high = p_search_low * 4;

found_p = 0;
found_j = 0;
{
forprime (p = p_search_low, p_search_high,
    local(H_mod, factored, n_rows, root_found, k, q, deg);
    if (kronecker(D, p) != 1, next);
    H_mod = H_D * Mod(1, p);
    factored = factor(H_mod);
    n_rows = matsize(factored)[1];
    root_found = 0;
    for (k = 1, n_rows,
        if (poldegree(factored[k, 1]) == 1,
            found_p = p;
            found_j = lift(-polcoeff(factored[k, 1], 0) / polcoeff(factored[k, 1], 1));
            root_found = 1;
            break
        )
    );
    if (root_found, break)
);
}

if (found_p == 0, print("GENERATED=0"); quit);

\\ Build E with j(E) = found_j.  Standard formula:
\\   j = 0    ⇒ E: y² = x³ + 1
\\   j = 1728 ⇒ E: y² = x³ + x
\\   otherwise: c = j/(1728 - j) mod p; E: y² = x³ + 3c x + 2c
p_val = found_p;
j_val = found_j;
{
if (j_val == 0, a_val = 0; b_val = 1, if (j_val == 1728, a_val = 1; b_val = 0, c_val = lift(Mod(j_val, p_val) / Mod(1728 - j_val, p_val)); a_val = lift(Mod(3 * c_val, p_val)); b_val = lift(Mod(2 * c_val, p_val))));
}

E_fresh = ellinit([a_val, b_val], p_val);
n_val = ellcard(E_fresh);

print("FRESH_P=", Strprintf("%x", p_val));
print("FRESH_N=", Strprintf("%x", n_val));
print("FRESH_A=", Strprintf("%x", a_val));
print("FRESH_B=", Strprintf("%x", b_val));
print("FRESH_DISC=", D);
print("FRESH_CLASS_NUMBER=", print_h);
print("FRESH_J=", j_val);
print("GENERATED=1");
