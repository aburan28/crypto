\\ ============================================================
\\ Thread 20: λ/n threshold bisection for GLV-HNP Phase 2 attack
\\
\\ Known from 2026-07-26 run:
\\   λ/n = 0.07 (p=2677, n=2647): LLL FAILS (BKZ-40 also fails)
\\   λ/n = 0.34 (p=524347, n=523969): LLL SUCCEEDS at m=9
\\
\\ Goal: find the threshold between 0.07 and 0.34.
\\ Method: find j=0 prime-order curves at λ/n ≈ 0.10, 0.15, 0.20, 0.25, 0.30
\\         and test the column-balanced Phase 2 lattice attack at each.
\\
\\ Lattice: glv_hnp_phase2_attack.py balanced-column construction (S_K1=N//K1, etc.)
\\ PARI port: uses qflll(M~) for LLL reduction.
\\
\\ Run: gp -q glv_hnp_phase2_lambda_threshold.gp

default(parisize, 512000000);
default(timer, 0);

print("================================================================");
print("Thread 20: λ/n threshold bisection for GLV-HNP Phase 2 lattice");
print("================================================================");
print("");

\\ ============================================================
\\ Build column-balanced lattice (mirrors glv_hnp_phase2_attack.py)
\\
\\ dim = 2m+2.  Column layout (1-indexed):
\\   cols 1..m      : k1_i slots, scale S_K1 = N // K1_BOUND
\\   col  m+1       : d slot, scale S_D = 1
\\   cols m+2..2m+1 : k2_i slots, scale S_K2 = N // K2_BOUND
\\   col  2m+2      : Kannan slot, scale S_KANNAN = N
\\
\\ Rows (1-indexed):
\\   rows 1..m      : mod-n constraints: M[i,i] = N * S_K1
\\   row  m+1       : d-row:  M[m+1,i] = B_i * S_K1;  M[m+1,m+1] = S_D
\\   rows m+2..2m+1 : k2-rows: M[m+1+i,i] = -λ * S_K1; M[m+1+i,m+1+i] = S_K2
\\   row  2m+2      : Kannan: M[2m+2,i] = A_i * S_K1;  M[2m+2,2m+2] = S_KANNAN
\\ ============================================================
build_glv_lat(n_val, lam, K1, K2, A_vec, B_vec, m) = {
    local(dim, M, S_K1, S_K2, S_KANNAN);
    dim = 2*m + 2;
    S_K1     = n_val \ K1;
    S_K2     = n_val \ K2;
    S_KANNAN = n_val;
    M = matrix(dim, dim);
    for (i = 1, m, M[i, i] = n_val * S_K1);
    for (i = 1, m, M[m+1, i] = B_vec[i] * S_K1);
    M[m+1, m+1] = 1;       \\ S_D = 1
    for (i = 1, m,
        M[m+1+i, i]       = -lam * S_K1;
        M[m+1+i, m+1+i]   = S_K2
    );
    for (i = 1, m, M[2*m+2, i] = A_vec[i] * S_K1);
    M[2*m+2, dim] = S_KANNAN;
    M
}

\\ ============================================================
\\ Run one attack trial.  Returns 1 if d recovered, 0 otherwise.
\\ seed_val: setrand seed.  m_sigs: number of signatures.
\\ ============================================================
glv_attack(E, G, n_val, lam, K1, K2, d_secret, m_sigs, seed_val) = {
    local(A_vec, B_vec, i, k1, k2, k_full, R, r_sig, h_msg, s_sig, s_inv);
    local(M, U, Red, dim, last, sign_v, d_cand);
    setrand(seed_val);
    A_vec = vector(m_sigs);
    B_vec = vector(m_sigs);

    i = 1;
    while (i <= m_sigs,
        k1     = random(K1);
        k2     = random(K2);
        k_full = lift(Mod(k1 + lam*k2, n_val));
        if (k_full == 0, next);
        R      = ellmul(E, G, k_full);
        r_sig  = lift(R[1]) % n_val;
        if (r_sig == 0, next);
        h_msg  = random(n_val);
        s_sig  = lift(Mod((h_msg + d_secret * r_sig) / k_full, n_val));
        if (s_sig == 0, next);
        s_inv  = lift(Mod(1/s_sig, n_val));
        A_vec[i] = lift(Mod(h_msg  * s_inv, n_val));
        B_vec[i] = lift(Mod(r_sig  * s_inv, n_val));
        i += 1
    );

    M   = build_glv_lat(n_val, lam, K1, K2, A_vec, B_vec, m_sigs);
    dim = 2*m_sigs + 2;

    \\ LLL: qflll(M~) treats columns of M~ (= rows of M) as basis vectors.
    \\ Returns unimodular U; reduced rows = rows of U~ * M.
    U   = qflll(M~);
    Red = U~ * M;

    \\ Scan for d
    for (r = 1, dim,
        last = Red[r, dim];
        if (abs(last) != n_val, next);      \\ Kannan marker = ±S_KANNAN = ±n
        sign_v = last / n_val;
        d_cand = lift(Mod(sign_v * Red[r, m_sigs+1], n_val));
        if (d_cand == d_secret, return(1))
    );
    0
}

\\ ============================================================
\\ Compute effective λ/n ratio (canonical in [0, 0.5])
\\ lam and n-1-lam are the two roots of x^2+x+1=0 mod n.
\\ We use the smaller one.
\\ ============================================================
eff_ratio(lam, n_val) = {
    local(lam2);
    lam2 = n_val - 1 - lam;
    if (lam <= lam2, lam * 1.0 / n_val, lam2 * 1.0 / n_val)
}

\\ canonical lam (the smaller of lam and n-1-lam)
canon_lam(lam, n_val) = {
    local(lam2);
    lam2 = n_val - 1 - lam;
    if (lam <= lam2, lam, lam2)
}

\\ ============================================================
\\ Curve discovery: find j=0 prime-order curves at target λ/n
\\ ratios.  Search p from p_min to p_max with p ≡ 1 mod 6.
\\ ============================================================
\\ Target ratios (inclusive brackets of width 0.035):
\\   slot 1: [0.055, 0.090) — validates known fail (n=2647, r=0.070)
\\   slot 2: [0.090, 0.130) — new test
\\   slot 3: [0.130, 0.175) — new test
\\   slot 4: [0.175, 0.225) — new test
\\   slot 5: [0.225, 0.275) — new test
\\   slot 6: [0.275, 0.320) — new test
\\   slot 7: [0.320, 0.360) — validates known success (n=523969, r=0.340)
\\   (all 7 slots in [0.05, 0.36))

N_SLOTS = 7;
slot_lo  = [0.055, 0.090, 0.130, 0.175, 0.225, 0.275, 0.320];
slot_hi  = [0.090, 0.130, 0.175, 0.225, 0.275, 0.320, 0.360];
slot_label = ["0.07 (fail?)", "0.10", "0.15", "0.20", "0.25", "0.30", "0.34 (ok?)"];

\\ Slots already filled (p, b, n, lam)
slot_p   = vector(N_SLOTS, i, 0);
slot_b   = vector(N_SLOTS, i, 0);
slot_n   = vector(N_SLOTS, i, 0);
slot_lam = vector(N_SLOTS, i, 0);
slot_r   = vector(N_SLOTS, i, 0.0);

\\ Pre-fill known data points to save search time
\\ Known fail: p=2677, b=2, n=2647, λ=185  → r≈0.070
slot_p[1] = 2677; slot_b[1] = 2; slot_n[1] = 2647;
E_tmp = ellinit([0, 2], 2677);
sq_tmp = lift(sqrt(Mod(-3, 2647)));
lam_tmp = lift(Mod((-1+sq_tmp)/2, 2647));
if (lift(Mod(lam_tmp^2+lam_tmp+1,2647)) != 0, lam_tmp = lift(Mod((-1-sq_tmp)/2,2647)));
slot_lam[1] = canon_lam(lam_tmp, 2647);
slot_r[1] = eff_ratio(lam_tmp, 2647);

\\ Known success: p=524347, b=2, n=523969, λ=177902 → r≈0.340
\\ (too slow for ellmul at toy scale; instead use a smaller curve in same range)
\\ We'll search for a curve in slot 7 range from small primes first.
\\ (The 524347 curve is correct but searching small primes is much faster.)

print("---- Searching for curves at each λ/n slot ----");
print("");

P_MIN = 1000; P_MAX = 300000;
slots_done = 1;   \\ slot 1 pre-filled; remaining: 2..7

{
forprime (p = P_MIN, P_MAX,
    if (slots_done == N_SLOTS, break);
    if (Mod(p, 6) != 1, next);
    for (b = 1, 30,
        local(E_try, n_try, sq_try, lam_try, r_try, sl);
        E_try = ellinit([0, b], p);
        n_try = ellcard(E_try);
        if (!isprime(n_try) || n_try < 1000, next);
        if (kronecker(-3, n_try) != 1, next);
        sq_try  = lift(sqrt(Mod(-3, n_try)));
        lam_try = lift(Mod((-1+sq_try)/2, n_try));
        if (lift(Mod(lam_try^2+lam_try+1,n_try)) != 0,
            lam_try = lift(Mod((-1-sq_try)/2,n_try)));
        lam_try = canon_lam(lam_try, n_try);
        r_try = lam_try * 1.0 / n_try;
        \\ Find a free slot that matches this ratio
        for (sl = 2, N_SLOTS,
            if (slot_p[sl] == 0 && r_try >= slot_lo[sl] && r_try < slot_hi[sl],
                slot_p[sl] = p; slot_b[sl] = b;
                slot_n[sl] = n_try; slot_lam[sl] = lam_try; slot_r[sl] = r_try;
                slots_done += 1;
                break
            )
        )
    )
);
}

print("Curves found:");
print(Str("  Slot | target  |    p      |   n     | λ      | λ/n"));
print(Str("  -----|---------|-----------|---------|--------|--------"));
{
for (sl = 1, N_SLOTS,
    if (slot_p[sl] == 0,
        print("  ", sl, "    | ", slot_label[sl], " | NOT FOUND")
    ,
        printf("  %d    | %-11s | %7d   | %6d  | %6d | %.4f\n",
            sl, slot_label[sl], slot_p[sl], slot_n[sl], slot_lam[sl], slot_r[sl])
    )
);
}
print("");

\\ ============================================================
\\ Attack parameters: K1, K2, m sweep, seeds
\\ ============================================================
\\ K1_BOUND = sqrtint(n) \ 3  (about n^0.5 / 3: small but non-trivial)
\\ K2_BOUND = sqrtint(n) + 1   (GLV domain)
\\ Seeds: 3 per (m, slot) pair
SEEDS = [11, 42, 12345];

\\ ============================================================
\\ Run attack for each slot
\\ ============================================================
print("---- Attack results ----");
print("");
print("For each slot: test m = m_thresh-1 to m_thresh+4 with 3 seeds each.");
print("Report first m where 3/3 seeds succeed, or 'FAIL' if never.");
print("");

{
for (sl = 1, N_SLOTS,
    if (slot_p[sl] == 0,
        print("Slot ", sl, " (", slot_label[sl], "): CURVE NOT FOUND — skip");
        next
    );

    local(p_c, b_c, n_c, lam_c, r_c, E_c, G_c, K1_c, K2_c, eff_c, m_thresh_c);
    local(d_c, m_val, wins, all_m_results);

    p_c   = slot_p[sl];
    b_c   = slot_b[sl];
    n_c   = slot_n[sl];
    lam_c = slot_lam[sl];
    r_c   = slot_r[sl];

    E_c = ellinit([0, b_c], p_c);

    \\ Find generator: prime-order curve, so any non-O point works.
    \\ Iterate x=0,1,... and take first point on curve.
    G_c = 0;
    for (x_try = 0, p_c-1,
        local(rhs_try, y2_try, y_try);
        rhs_try = lift(Mod(x_try^3 + b_c, p_c));
        if (!issquare(Mod(rhs_try, p_c)), next);
        y_try = lift(sqrt(Mod(rhs_try, p_c)));
        G_c = [x_try, y_try];
        break
    );
    if (G_c == 0, print("No generator for slot ", sl, "!"); next);

    \\ Verify order
    if (ellmul(E_c, G_c, n_c) != [0],
        print("Generator order check failed for slot ", sl, "!"); next);

    \\ K1 ~ sqrt(n)/15 gives eff = K1*K2/n ~ (sqrt(n)/15)*sqrt(n)/n = 1/15 ~ 0.067
    \\ This matches the Python attack.py regime (eff=0.05-0.10, m_thresh=3-4)
    K1_c = max(2, sqrtint(n_c) \ 15);
    K2_c = sqrtint(n_c) + 1;
    eff_c = K1_c * K2_c * 1.0 / n_c;
    m_thresh_c = ceil(log(n_c * 1.0) / log(1.0 / eff_c));

    setrand(999);
    d_c = random(n_c - 2) + 1;   \\ d in [1, n-1]

    printf("Slot %d: λ/n=%.4f  n=%d  K1=%d  K2=%d  eff=%.4f  m_thresh=%d\n",
        sl, r_c, n_c, K1_c, K2_c, eff_c, m_thresh_c);

    all_m_results = "";
    \\ Test up to m_thresh+10 to clearly distinguish "fails even with many sigs"
    for (m_val = max(2, m_thresh_c - 1), m_thresh_c + 10,
        wins = 0;
        for (si = 1, #SEEDS,
            wins += glv_attack(E_c, G_c, n_c, lam_c, K1_c, K2_c, d_c, m_val, SEEDS[si])
        );
        all_m_results = Str(all_m_results, " m=", m_val, ":", wins, "/", #SEEDS);
        if (wins == #SEEDS,
            printf("  FIRST FULL SUCCESS at m=%d (3/3)\n", m_val);
            break
        ,
            if (m_val == m_thresh_c + 10,
                printf("  FAILED at all tested m up to %d\n", m_thresh_c + 10)
            )
        )
    );
    printf("  Results:%s\n", all_m_results);
    print("")
);
}

print("================================================================");
print("Thread 20 complete.");
print("================================================================");
