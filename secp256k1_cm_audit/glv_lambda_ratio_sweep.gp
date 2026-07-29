\\ Thread 20 Part 2: λ/n ratio sweep with controlled n size
\\
\\ Fixes n in [2000,4000] (12-bit) to control K1, K2 across curves.
\\ Tests j=0 prime-order curves with varying λ/n and records m*.
\\
\\ Key question: when n is fixed, does larger λ/n give smaller m*?
\\ Answer determines if λ/n is the driver or just a correlate of n-size.
\\
\\ Run: gp -q glv_lambda_ratio_sweep.gp

default(parisize, 256000000);
default(timer, 0);

print("================================================================");
print("Thread 20 Part 2: λ/n ratio sweep, fixed n in [2000,4000]");
print("================================================================");
print("");

\\ ---- Lattice builder (no embedded braces) ----
build_glv_lat(n_val, lam, K1, K2, A_vec, B_vec, m) = {
    local(dim, M, S_K1, S_K2);
    dim = 2*m + 2;
    S_K1 = n_val \ K1;
    S_K2 = n_val \ K2;
    M = matrix(dim, dim);
    for (i = 1, m, M[i, i] = n_val * S_K1);
    for (i = 1, m, M[m+1, i] = B_vec[i] * S_K1);
    M[m+1, m+1] = 1;
    for (i = 1, m, M[m+1+i, i] = -lam * S_K1; M[m+1+i, m+1+i] = S_K2);
    for (i = 1, m, M[2*m+2, i] = A_vec[i] * S_K1);
    M[2*m+2, 2*m+2] = n_val;
    M
}

\\ ---- Attack runner: 3 seeds, returns win count ----
glv_3wins(E_c, G_c, n_val, lam, K1, K2, d_secret, m_sigs) = {
    local(wins, i, k1, k2, k_full, R, r_sig, h_msg, s_sig, s_inv);
    local(A_vec, B_vec, M, U, Red, dim, last, sign_v, d_cand, ok);
    wins = 0;
    for (sd = 1, 3,
        setrand(sd * 1000 + m_sigs * 7 + n_val);
        A_vec = vector(m_sigs); B_vec = vector(m_sigs);
        i = 1;
        while (i <= m_sigs,
            k1 = random(K1); k2 = random(K2);
            k_full = lift(Mod(k1 + lam*k2, n_val));
            if (k_full == 0, next);
            R = ellmul(E_c, G_c, k_full);
            r_sig = lift(R[1]) % n_val;
            if (r_sig == 0, next);
            h_msg = random(n_val);
            s_sig = lift(Mod((h_msg + d_secret*r_sig) / k_full, n_val));
            if (s_sig == 0, next);
            s_inv = lift(Mod(1/s_sig, n_val));
            A_vec[i] = lift(Mod(h_msg * s_inv, n_val));
            B_vec[i] = lift(Mod(r_sig * s_inv, n_val));
            i += 1
        );
        M = build_glv_lat(n_val, lam, K1, K2, A_vec, B_vec, m_sigs);
        dim = 2*m_sigs + 2;
        U = qflll(M~);
        Red = U~ * M;
        ok = 0;
        for (r = 1, dim,
            last = Red[r, dim];
            if (abs(last) != n_val, next);
            sign_v = last / n_val;
            d_cand = lift(Mod(sign_v * Red[r, m_sigs+1], n_val));
            if (d_cand == d_secret, ok = 1; break)
        );
        wins += ok
    );
    wins
}

\\ ---- Scan j=0 prime-order curves with n in [2000,4000] ----
N_MIN = 2000; N_MAX = 4000;
print("Scanning j=0 prime-order curves n in [", N_MIN, ",", N_MAX, "] ...");
print("");

curves_found = 0;
MAX_CURVES = 50;
res_n     = vector(MAX_CURVES);
res_ratio = vector(MAX_CURVES);
res_mstar = vector(MAX_CURVES);

{
forprime (p = N_MIN, 30000,
    if (curves_found >= MAX_CURVES, break);
    if (Mod(p, 6) != 1, next);
    for (b = 1, 30,
        local(E2, n2, sq2, lam2, lam_c, r2, G2, K1_c, K2_c, d_c, m_val, wins, ms);
        if (curves_found >= MAX_CURVES, break(2));
        E2 = ellinit([0, b], p);
        n2 = ellcard(E2);
        if (!isprime(n2) || n2 < N_MIN || n2 > N_MAX, next);
        if (kronecker(-3, n2) != 1, next);
        sq2 = lift(sqrt(Mod(-3, n2)));
        lam2 = lift(Mod((-1+sq2)/2, n2));
        if (lift(Mod(lam2^2+lam2+1,n2)) != 0, lam2 = lift(Mod((-1-sq2)/2,n2)));
        lam_c = if (lam2 <= n2-1-lam2, lam2, n2-1-lam2);
        r2 = lam_c * 1.0 / n2;
        \\ Find generator
        G2 = 0;
        for (x = 0, p,
            local(rhs2, y2v);
            rhs2 = lift(Mod(x^3+b, p));
            if (!issquare(Mod(rhs2,p)), next);
            G2 = [x, lift(sqrt(Mod(rhs2, p)))]; break
        );
        if (G2 == 0, next);
        K1_c = max(2, sqrtint(n2) \ 15);
        K2_c = sqrtint(n2) + 1;
        setrand(7777 + n2 + b);
        d_c = random(n2 - 2) + 1;
        ms = -1;
        for (m_val = 2, 14,
            wins = glv_3wins(E2, G2, n2, lam_c, K1_c, K2_c, d_c, m_val);
            if (wins == 3, ms = m_val; break)
        );
        curves_found += 1;
        res_n[curves_found]     = n2;
        res_ratio[curves_found] = r2;
        res_mstar[curves_found] = ms;
        printf("  [%2d] n=%4d  λ/n=%.4f  K1=%d  K2=%d  m*=%2d\n",
            curves_found, n2, r2, K1_c, K2_c, ms)
    )
);
}

print("");
print("Total curves: ", curves_found);
print("");

\\ Pearson correlation between λ/n and m*
n_ok = 0;
sum_r = 0.0; sum_m = 0.0; sum_r2 = 0.0; sum_m2 = 0.0; sum_rm = 0.0;
for (i = 1, curves_found,
    if (res_mstar[i] < 0, next);
    n_ok += 1;
    sum_r  += res_ratio[i];
    sum_m  += res_mstar[i];
    sum_r2 += res_ratio[i]^2;
    sum_m2 += res_mstar[i]^2;
    sum_rm += res_ratio[i] * res_mstar[i]
);
if (n_ok > 2,
    local(mr, mm, cov, vr, vm, pearson);
    mr = sum_r / n_ok; mm = sum_m / n_ok;
    cov = sum_rm/n_ok - mr*mm;
    vr  = sum_r2/n_ok - mr^2;
    vm  = sum_m2/n_ok - mm^2;
    pearson = if (vr*vm > 0.0, cov / sqrt(vr*vm), 0);
    printf("Pearson r(λ/n, m*) = %.4f  (n_curves=%d)\n", pearson, n_ok);
    printf("Mean λ/n=%.4f  Mean m*=%.2f\n", mr, mm);
    if (abs(pearson) < 0.3,
        print("FINDING: weak correlation — λ/n is NOT the main driver of m*");
        print("         n-size and eff dominate. 2026-07-26 threshold was n-size artifact."),
        if (pearson < -0.3,
            print("FINDING: negative correlation — larger λ/n → smaller m* (easier attack)");
            print("         λ/n IS a significant driver of lattice difficulty."),
            print("FINDING: positive correlation — larger λ/n → larger m* (harder attack)");
            print("         Unexpected result, needs investigation.")
        )
    )
);

\\ Summary sorted by λ/n
print("");
print("Sorted by λ/n:");
print("  λ/n    | m*");
print("  -------|---");
\\ Insertion sort
for (i = 2, curves_found,
    local(j, ri, ni, mi);
    ri = res_ratio[i]; ni = res_n[i]; mi = res_mstar[i];
    j = i - 1;
    while (j >= 1 && res_ratio[j] > ri,
        res_ratio[j+1] = res_ratio[j];
        res_n[j+1] = res_n[j];
        res_mstar[j+1] = res_mstar[j];
        j -= 1
    );
    res_ratio[j+1] = ri; res_n[j+1] = ni; res_mstar[j+1] = mi
);
for (i = 1, curves_found,
    printf("  %.4f  | %2d  (n=%d)\n", res_ratio[i], res_mstar[i], res_n[i])
);

print("");
print("================================================================");
print("Thread 20 Part 2 complete.");
print("================================================================");
