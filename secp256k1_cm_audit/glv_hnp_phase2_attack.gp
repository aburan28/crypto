\\ ============================================================
\\ Phase 2 GLV-aware HNP toy ATTACK demonstration
\\ Extends glv_hnp_phase2_toy.gp with three concrete findings:
\\
\\ Finding 1: brute-force d recovery using public key Q (r-comparison).
\\
\\ Finding 2: k_1-only-leak feasibility check is UNDERCONSTRAINED:
\\   for ANY d_cand, some k_{i,2} makes k_{i,1}<K1_BOUND. Zero power.
\\
\\ Finding 3: both-bounded GLV model (k_1 AND k_2 << sqrt(n))
\\   admits unique recovery WITHOUT public key via (K_BOUND^2/n)^m << 1.
\\
\\ Run: gp -q glv_hnp_phase2_attack.gp

default(parisize, 256000000);
default(timer, 0);

print("================================================================");
print("Phase 2 GLV-HNP toy ATTACK demonstration");
print("================================================================\n");

\\ ---- Setup: same toy j=0 curve as glv_hnp_phase2_toy.gp ----
print("---- Setup: toy j=0 curve with GLV ----");
found_p = 0; found_b = 0; found_n = 0;
{
forprime(p = 200, 2000,
    local(E0, n0);
    if (Mod(p, 6) != 1, next);
    for (b = 1, 50,
        E0 = ellinit([0, b], p);
        n0 = ellcard(E0);
        if (isprime(n0) && kronecker(-3, n0) == 1,
            found_p = p; found_b = b; found_n = n0;
            break(2)
        )
    )
);
}
p_val = found_p; b_val = found_b; n_val = found_n;
E = ellinit([0, b_val], p_val);
G = ellgenerators(E)[1];
sqrt_neg3 = lift(sqrt(Mod(-3, n_val)));
lambda = lift(Mod(-1 + sqrt_neg3, n_val) / Mod(2, n_val));
lambda_inv = lift(Mod(1, n_val) / Mod(lambda, n_val));
print("E: y^2 = x^3 + ", b_val, " over F_", p_val);
print("n (prime) = ", n_val);
print("lambda = ", lambda, "   lambda^{-1} mod n = ", lambda_inv);
print("lambda^2 + lambda + 1 mod n = ", lift(Mod(lambda^2 + lambda + 1, n_val)));
print("");

\\ Helper: compute r-value from k (x-coord of k*G mod n)
compute_r(k_full) = lift(ellmul(E, G, k_full)[1]) % n_val;

\\ ================================================================
\\ PART 1: k_1-only-bias model (k_2 full-range, as in Phase 2 toy)
\\ ================================================================
print("================================================================");
print("PART 1: k_1-only-bias model (k_2 full range)");
print("================================================================\n");

setrand(11);
d_secret = random(n_val - 1) + 1;
Q_pub = ellmul(E, G, d_secret);
K1_BOUND = max(2, n_val \ 100);
print("Planted d = ", d_secret);
print("Public key Q = ", Q_pub);
print("K1_BOUND = ", K1_BOUND, "   (k_1 in [0, K1_BOUND), k_2 full-range [0,n))");
print("");

m_sigs = 5;
A_vec = vector(m_sigs);
B_vec = vector(m_sigs);
r_vec = vector(m_sigs);
k1_planted = vector(m_sigs);

{
i = 1;
while(i <= m_sigs,
    local(k1, k2, k_full, R_pt, r_sig, s_sig, h_msg);
    k1 = random(K1_BOUND - 1) + 1;
    k2 = random(n_val - 1) + 1;
    k_full = lift(Mod(k1 + lambda * k2, n_val));
    if (k_full == 0, next);
    R_pt = ellmul(E, G, k_full);
    r_sig = lift(R_pt[1]) % n_val;
    if (r_sig == 0, next);
    h_msg = random(n_val);
    s_sig = lift(Mod((h_msg + d_secret * r_sig) / k_full, n_val));
    if (s_sig == 0, next);
    k1_planted[i] = k1;
    r_vec[i] = r_sig;
    A_vec[i] = lift(Mod(h_msg / s_sig, n_val));
    B_vec[i] = lift(Mod(r_sig / s_sig, n_val));
    i = i + 1
);
}
print("m = ", m_sigs, " signatures.  Stored r values: ", r_vec);
print("Planted k_1 values: ", k1_planted);
print("");

\\ ---- FINDING 1: Brute-force with public key (r-comparison) ----
print("---- FINDING 1: Brute-force via public-key r-comparison ----");
found_d = 0;
n_ec_mults = 0;
{
for(d_cand = 1, n_val - 1,
    local(match, k_full_cand, r_cand);
    match = 1;
    for(i = 1, m_sigs,
        k_full_cand = lift(Mod(A_vec[i] + B_vec[i] * d_cand, n_val));
        r_cand = compute_r(k_full_cand);
        n_ec_mults = n_ec_mults + 1;
        if (r_cand != r_vec[i], match = 0; break)
    );
    if (match && found_d == 0, found_d = d_cand)
);
}
print("Brute-force over ", n_val - 1, " candidates, ", n_ec_mults, " EC multiplications.");
print("Recovered d = ", found_d, "   Planted d = ", d_secret);
print("Recovery correct? ", if(found_d == d_secret, "YES", "NO"));
print("");
print("Note: r-comparison correctly discriminates via the ECDSA commitment.");
print("The check r_cand == r_vec[i] fails for wrong d because only k_full=k_i");
print("satisfies k_i*G x-coord == r_i (the other solution n-k_i gives neg. G pt");
print("with same x-coord, corresponding to a different d by CRT).\n");

\\ ---- FINDING 2: k_1-only feasibility = vacuous (all d pass) ----
print("---- FINDING 2: k_1-only feasibility check (no public key) ----");
print("For each d_cand: does any k_{i,2} in [0,n) give k_{i,1}<K1_BOUND?");
print("Since lambda is invertible mod n, as k_{i,2} ranges over [0,n),");
print("(A_i+B_i*d - lambda*k_{i,2}) mod n takes ALL n residues.");
print("Exactly K1_BOUND of them are in [0, K1_BOUND).  Hence EVERY d passes.\n");

feasible_all = 0;
{
for(d_cand = 1, n_val - 1,
    local(all_ok, k1_check);
    all_ok = 1;
    for(i = 1, m_sigs,
        \\ Try k_{i,1}=0: k_{i,2} = lambda_inv*(A_i+B_i*d) mod n (in [0,n), always valid)
        k1_check = lift(Mod(A_vec[i] + B_vec[i]*d_cand, n_val));
        k1_check = k1_check % K1_BOUND;   \\ this k_{i,1} is in [0,K1_BOUND) by construction
        \\ k_{i,2} = lambda_inv*(A_i+B_i*d - k1_check) mod n, always in [0,n)
        \\ All d pass trivially.
        if (0, all_ok = 0; break)   \\ never fail — the check is vacuous
    );
    if (all_ok, feasible_all = feasible_all + 1)
);
}
print("d_cand values passing k_1-only feasibility: ", feasible_all, " / ", n_val - 1);
print("(All pass — the check gives zero discriminating power.)\n");

\\ Explicit demonstration: for d_cand = (d_secret+1)%n (a wrong candidate),
\\ show we can find k_{i,1}<K1_BOUND and k_{i,2} in [0,n) for EVERY signature.
d_wrong = (d_secret % (n_val-1)) + 1;
print("Explicit wrong candidate d_wrong = ", d_wrong, ":");
{
for(i = 1, m_sigs,
    local(c, k1_ex, k2_ex, check);
    c = lift(Mod(A_vec[i] + B_vec[i]*d_wrong, n_val));
    k1_ex = c % K1_BOUND;   \\ k_{i,1} in [0, K1_BOUND)
    k2_ex = lift(Mod((c - k1_ex) * lambda_inv, n_val));  \\ k_{i,2} in [0,n)
    check = lift(Mod(k1_ex + lambda * k2_ex, n_val));
    print("  sig ", i, ":  k_1=", k1_ex, " (< ", K1_BOUND, ")  k_2=", k2_ex,
          "  verify A+Bd=k1+lambda*k2: ", if(check==c, "OK", "FAIL"))
);
}
print("\nConclusion: k_1-only bias gives ZERO lattice attack power without Q.\n");

\\ ================================================================
\\ PART 2: Both-bounded GLV model (k_1 AND k_2 small)
\\ ================================================================
print("================================================================");
print("PART 2: Both-bounded GLV model (k_1 AND k_2 < K_BB << sqrt(n))");
print("================================================================\n");

\\ For discrimination: need (K_BB^2/n)^m << 1/(n-1)
\\ K_BB = 4: (16/199)^5 * 198 ≈ 198 * 8.2e-7 ≈ 1.6e-4 (< 1, unique)
\\ K_BB = 5: (25/199)^5 * 198 ≈ 198 * 1.0e-5 ≈ 0.002 (< 1, unique)
K_BB = 4;
print("K_BB (both-bounded bound) = ", K_BB);
expected_fp_bb = (K_BB^2.0 / n_val)^m_sigs * (n_val - 1);
print("Expected false positives (theory): (K_BB^2/n)^m * (n-1) = ", expected_fp_bb);
print("Compare: sqrt(n) = ", sqrt(n_val*1.0), " — K_BB << sqrt(n) ensures discrimination.\n");

setrand(77);
d_bb = random(n_val - 1) + 1;
print("Planted d (both-bounded) = ", d_bb);

A_bb = vector(m_sigs);
B_bb = vector(m_sigs);
k1_bb = vector(m_sigs);
k2_bb = vector(m_sigs);
{
i = 1;
while(i <= m_sigs,
    local(k1, k2, k_full, R_pt, r_sig, s_sig, h_msg);
    k1 = random(K_BB - 1) + 1;
    k2 = random(K_BB - 1) + 1;
    k_full = lift(Mod(k1 + lambda * k2, n_val));
    if (k_full == 0, next);
    R_pt = ellmul(E, G, k_full);
    r_sig = lift(R_pt[1]) % n_val;
    if (r_sig == 0, next);
    h_msg = random(n_val);
    s_sig = lift(Mod((h_msg + d_bb * r_sig) / k_full, n_val));
    if (s_sig == 0, next);
    k1_bb[i] = k1; k2_bb[i] = k2;
    A_bb[i] = lift(Mod(h_msg / s_sig, n_val));
    B_bb[i] = lift(Mod(r_sig / s_sig, n_val));
    i = i + 1
);
}
print("k_1 values: ", k1_bb, "  k_2 values: ", k2_bb);
print("Combined nonces k_full: ", vector(m_sigs, i, lift(Mod(k1_bb[i]+lambda*k2_bb[i],n_val))));
print("");

\\ Brute-force WITHOUT public key: check if any (k1,k2) in [0,K_BB)^2
\\ decomposes c_i = A_i + B_i*d_cand for ALL i.
bb_feasible_count = 0;
bb_found_d = 0;
{
for(d_cand = 1, n_val - 1,
    local(all_ok, c_i, k1_try, k2_try);
    all_ok = 1;
    for(i = 1, m_sigs,
        local(found_pair);
        found_pair = 0;
        c_i = lift(Mod(A_bb[i] + B_bb[i] * d_cand, n_val));
        for(k1_try = 1, K_BB - 1,
            k2_try = lift(Mod((c_i - k1_try) * lambda_inv, n_val));
            if (k2_try >= 1 && k2_try < K_BB,
                found_pair = 1; break
            )
        );
        if (!found_pair, all_ok = 0; break)
    );
    if (all_ok,
        bb_feasible_count = bb_feasible_count + 1;
        if (bb_found_d == 0, bb_found_d = d_cand)
    )
);
}
print("Both-bounded feasibility (K_BB=", K_BB, ", m=", m_sigs, ", NO public key):");
print("  Feasible candidates: ", bb_feasible_count, " / ", n_val - 1);
print("  First feasible: ", bb_found_d, "   Planted d: ", d_bb);
print("  Recovery correct? ", if(bb_found_d == d_bb, "YES", "NO"));
print("  Actual false positives: ", bb_feasible_count - if(bb_found_d == d_bb, 1, 0));
print("");

\\ ================================================================
\\ SUMMARY
\\ ================================================================
print("================================================================");
print("SUMMARY of Findings");
print("================================================================\n");
print("Curve: y^2=x^3+", b_val, " over F_", p_val, ", n=", n_val, ", lambda=", lambda);
print("");
f1_ok = if(found_d == d_secret, "YES", "NO");
print("Finding 1 (k_1-only, WITH public key Q):");
print("  Brute-force O(n) EC mults recovers d uniquely.  Needs Q.");
print("  Correct? ", f1_ok, ".  d=", found_d, " (planted ", d_secret, ").");
print("");
print("Finding 2 (k_1-only, NO public key):");
print("  ALL ", n_val-1, " candidates pass k_1-only feasibility check.");
print("  k_1-only lattice attack is VACUOUS without public key.");
print("  The Phase 2 lattice design (RESEARCH_GLV_HNP_PHASE2.md §2) is");
print("  UNDERCONSTRAINED: needs Q for any discrimination.");
print("");
bb_ok = if(bb_found_d == d_bb, "YES", "NO");
bb_fp = bb_feasible_count - if(bb_found_d == d_bb, 1, 0);
print("Finding 3 (both-bounded, K_BB=", K_BB, "<<sqrt(n)=", floor(sqrt(n_val*1.0)), ", NO Q):");
print("  Only ", bb_feasible_count, " candidates feasible.  Recovery correct? ", bb_ok, ".");
print("  Expected FP: ", expected_fp_bb, ".  Actual FP: ", bb_fp, ".");
print("");
print("Conclusion:");
print("  The Phase 2 threat model (k_1-only-leak, k_2 full-range) needs Q.");
print("  For a lattice-only attack (no Q), BOTH k_1 AND k_2 must be small");
print("  (K_BB < n^{(m-1)/(2m)} for m sigs, here n^0.4=", n_val^0.4, ").");
print("  The realistic GLV model (k_1,k_2 ~ sqrt(n)) is NOT directly exploitable");
print("  via a standard lattice attack — contradiction with Phase 2 design note.");
f3_ok = if(bb_found_d == d_bb, "YES", "NO");
print("  Finding 3 recovery confirmed: ", f3_ok, " (d=", bb_found_d, " planted=", d_bb, ").");
print("");
print("Next step: implement (2m+1)-dim lattice for both-bounded model at");
print("64-bit scale; determine minimum bias (K_BB vs n ratio) for LLL success.");
