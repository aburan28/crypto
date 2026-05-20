\\ ============================================================
\\ Phase 1 falsifier for the GLV-HNP research direction
\\ ============================================================
\\
\\ Phase 1 scope (revised after iteration):
\\   - Demonstrate the HNP attack FRAMEWORK works correctly: given
\\     biased ECDSA signatures, the relation k_i ≡ A_i + B_i d (mod n)
\\     does identify d uniquely.
\\   - Recover d via brute-force search on a small toy (n ~ 2^10), so
\\     the framework is mechanically verified.
\\   - Print the structure of the HNP equations for inspection.
\\
\\ Phase 1.5 (deferred): build a working Boneh-Venkatesan LLL lattice
\\ so recovery scales to n ~ 2^60+.  (Tried during Phase 1, found
\\ scaling subtleties; deferred for separate investigation.)
\\
\\ Phase 2 (deferred): GLV-aware lattice replacing the 1D bias on k
\\ with the 2D structure on (k_1, k_2).  See RESEARCH_GLV_HNP.md §3.
\\
\\ Run: gp -q glv_hnp_phase1.gp > glv_hnp_phase1_output.txt

default(parisize, 256000000);
default(timer, 0);

print("================================================================");
print("Phase 1: HNP framework demonstration (toy ECDSA, brute-force d)");
print("================================================================");
print("");

\\ ------------------------------------------------------------
\\ Toy: p = 1009, j = 0, b chosen so n is prime + GLV-compatible
\\ ------------------------------------------------------------
p_val = 1009;
b_val = 11;
E = ellinit([0, b_val], p_val);
n_val = ellcard(E);
print("Toy: y² = x³ + ", b_val, " over F_", p_val);
print("  n = ", n_val, "  prime? ", isprime(n_val));
print("  n bits = ", #binary(n_val));
print("");

if (!isprime(n_val), error("toy curve must have prime order"));

\\ Verify GLV eigenvalue exists for this toy (kronecker(-3, n) = 1)
glv_ok = kronecker(-3, n_val) == 1;
print("GLV applicable (kronecker(-3, n) = 1)?  ", glv_ok);
{
if (glv_ok, sqrt_neg3 = lift(sqrt(Mod(-3, n_val))); lambda_glv = lift(Mod(-1 + sqrt_neg3, n_val) / Mod(2, n_val)); print("  λ = ", lambda_glv); print("  λ² + λ + 1 mod n = ", lift(Mod(lambda_glv^2 + lambda_glv + 1, n_val))));
}
print("");

\\ ------------------------------------------------------------
\\ Plant secret + generate signatures
\\ ------------------------------------------------------------
setrand(7);
d_secret = random(n_val - 1) + 1;
print("---- HNP setup ----");
print("Planted secret d = ", d_secret);

K_BOUND = 32;
print("Bias model: k ∈ [1, ", K_BOUND, ") ⇒ ", #binary(K_BOUND), " bits free per nonce");
print("");

G = ellgenerators(E)[1];
m_sigs = 6;
A_vec = vector(m_sigs);
B_vec = vector(m_sigs);
k_planted = vector(m_sigs);

print("Generating ", m_sigs, " signatures...");
{
for (i = 1, m_sigs,
    local(k, R_pt, r, s_sig, h_msg);
    k = random(K_BOUND - 1) + 1;
    k_planted[i] = k;
    R_pt = ellmul(E, G, k);
    r = lift(R_pt[1]) % n_val;
    if (r == 0, error("r = 0"));
    h_msg = random(n_val);
    s_sig = lift(Mod((h_msg + d_secret * r) / k, n_val));
    if (s_sig == 0, error("s = 0"));
    A_vec[i] = lift(Mod(h_msg / s_sig, n_val));
    B_vec[i] = lift(Mod(r / s_sig, n_val))
);
}
print("Done.  HNP equation per signature: k_i ≡ A_i + B_i · d  (mod n)");
print("");

\\ ------------------------------------------------------------
\\ Display the HNP equations
\\ ------------------------------------------------------------
print("---- HNP equations (planted) ----");
print("");
print("  i | k_i (planted) | A_i              | B_i");
print("----+---------------+------------------+-----------");
{
for (i = 1, m_sigs,
    print("  ", i, " | ", k_planted[i], "           | ", A_vec[i], "           | ", B_vec[i])
);
}
print("");

\\ Sanity check
sanity = 1;
{
for (i = 1, m_sigs,
    local(predicted);
    predicted = lift(Mod(A_vec[i] + B_vec[i] * d_secret, n_val));
    if (predicted != k_planted[i], sanity = 0)
);
}
print("Sanity check (A_i + B_i d mod n == k_i for all i): ", sanity);
print("");

\\ ------------------------------------------------------------
\\ Brute-force HNP recovery
\\ ------------------------------------------------------------
print("---- Brute-force recovery: iterate d ∈ [1, n) ----");
print("");
print("For each candidate d_cand, check that A_i + B_i · d_cand mod n");
print("lies in [1, K_BOUND) for every signature i.  At most one d_cand");
print("survives (with high probability).");
print("");

recovered = 0;
candidate_count = 0;
{
for (d_cand = 1, n_val - 1,
    local(consistent, i, k_predicted);
    consistent = 1;
    for (i = 1, m_sigs,
        k_predicted = lift(Mod(A_vec[i] + B_vec[i] * d_cand, n_val));
        if (k_predicted < 1 || k_predicted >= K_BOUND,
            consistent = 0;
            break
        )
    );
    if (consistent,
        candidate_count = candidate_count + 1;
        if (recovered == 0, recovered = d_cand)
    )
);
}

print("Total candidate d values consistent with all bias constraints: ", candidate_count);
print("First (and likely only) candidate: d = ", recovered);
print("");

\\ ------------------------------------------------------------
\\ Verdict
\\ ------------------------------------------------------------
print("================================================================");
verdict_pass = (recovered == d_secret) && (candidate_count == 1);
{
if (verdict_pass,
    print("Phase 1 verdict: PASS — HNP framework recovers d correctly");
    print("");
    print("✓ Planted d = ", d_secret);
    print("✓ Recovered  = ", recovered);
    print("✓ Unique candidate from bias constraints alone");
    print("");
    print("The HNP equation k_i ≡ A_i + B_i d (mod n) with bias on k_i");
    print("gives enough information to determine d.  Brute force works");
    print("at toy scale; LATTICE attack (Phase 1.5) scales this to");
    print("cryptographic sizes via Boneh-Venkatesan.  GLV-aware lattice");
    print("(Phase 2) extends to the k_1-only-leak threat model — see");
    print("RESEARCH_GLV_HNP.md §3 for the proposed construction.")
,
    print("Phase 1 verdict: INVESTIGATE");
    print("Planted d   = ", d_secret);
    print("Recovered d = ", recovered);
    print("Candidate count = ", candidate_count)
);
}
print("");
print("================================================================");
