\\ ============================================================
\\ Thread 18: Verify trace assignments for b_1 and b_4
\\ (sextic twists of secp256k1) via scalar multiplication.
\\
\\ Background: howe_sextic_twists_all15.gp assigns traces to b_k
\\ using CM formula ordering without scalar-mult verification.
\\ The 2026-07-08 Python script found different trace assignments
\\ for k=1,2,4,5, giving glueable pair set {(0,2),(0,3),(0,5),(2,3),(2,5)}
\\ vs. PARI's {(0,2),(0,3),(0,5),(1,4),(2,3)}.
\\
\\ Here we verify via scalar multiplication which candidate trace
\\ is actually correct for b_1 and b_4.
\\
\\ Run: gp -q thread18_trace_verify.gp

default(parisize, 512000000);
default(timer, 0);

p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
n_known = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;
t_known = p + 1 - n_known;
s = sqrtint((4*p - t_known^2) / 3);

\\ The 6 candidate traces
t_plus  = (t_known + 3*s) / 2;
t_minus = (t_known - 3*s) / 2;
T_cands = [t_known, t_minus, -(t_known + 3*s)/2, -t_known, (3*s - t_known)/2, t_plus];
N_cands = vector(6, jj, p + 1 - T_cands[jj]);

\\ Find primitive 6th root of unity
g_gen = znprimroot(p);
u = lift(g_gen ^ ((p - 1) / 6));

\\ b values
b_vals = vector(6, jj, lift(Mod(7, p) * Mod(u, p)^(jj-1)));

print("================================================================");
print("Thread 18: Trace verification for sextic twists via scalar mult");
print("================================================================");
print("");
print("6 candidate N values (p+1-T_k):");
for (k = 1, 6, print("  N_", k-1, " = ", N_cands[k]));
print("");

\\ --- Verify trace of b_1 (k=1) ---
\\ Find a random affine point on E_{b_1}: y^2 = x^3 + b_1 over F_p.
\\ Then find which N_cand satisfies N_cand * P = O.

verify_trace = func(bval, label) = {
    local(E, P, correct_idx, correct_N, correct_T);
    E = ellinit([0, bval], p);
    \\ Find a point (try x = 2, 3, ... until y^2 = x^3 + bval is a QR)
    P = 0;
    for (x = 2, 1000,
        local(rhs, y2, sq);
        rhs = (x^3 + bval) % p;
        sq = Mod(rhs, p)^((p-1)/2);
        if (lift(sq) == 1,
            y2 = lift(Mod(rhs, p)^((p+1)/4));
            P = [x, y2];
            break
        )
    );
    if (P == 0, error("No point found for b = ", bval));
    print("Curve ", label, " (b = ", bval % 10^10, "...): point P = (", P[1], ", ", P[2] % 10^10, "...)");
    correct_idx = 0;
    for (k = 1, 6,
        local(Pk);
        Pk = ellmul(E, P, N_cands[k]);
        if (Pk == [0], correct_idx = k; break)
    );
    if (correct_idx == 0,
        print("  ERROR: no candidate N annihilates P!")
    ,
        correct_N = N_cands[correct_idx];
        correct_T = T_cands[correct_idx];
        print("  Correct trace index: T[", correct_idx, "] = ", correct_T);
        print("  Correct order: N[", correct_idx, "] = ", correct_N);
        print("  N mod 4  = ", correct_N % 4, "   N mod 3 = ", correct_N % 3);
        print("  gcd(N,4) = ", gcd(correct_N, 4), "  gcd(N,3) = ", gcd(correct_N, 3));
    );
    print("");
    return(correct_idx);
};

print("---- Verifying b_1 (k=1) ----");
idx1 = verify_trace(b_vals[2], "b_1");

print("---- Verifying b_4 (k=4) ----");
idx4 = verify_trace(b_vals[5], "b_4");

print("---- Verifying b_2 (k=2) and b_5 (k=5) for completeness ----");
idx2 = verify_trace(b_vals[3], "b_2");
idx5 = verify_trace(b_vals[6], "b_5");

\\ Summary of correct assignments
print("================================================================");
print("SUMMARY: Correct CM-trace assignment (scalar-mult verified)");
print("================================================================");
print("");
print("k=0: N[1] =", N_cands[1], "(secp256k1)");
print("k=1: N[", idx1, "] =", N_cands[idx1]);
print("k=2: N[", idx2, "] =", N_cands[idx2]);
print("k=3: N[4] =", N_cands[4], "(quad twist)");
print("k=4: N[", idx4, "] =", N_cands[idx4]);
print("k=5: N[", idx5, "] =", N_cands[idx5]);
print("");

\\ Recompute gcd(N_1, N_4) and gcd(N_2, N_5) with correct assignments
N_b1 = N_cands[idx1]; N_b4 = N_cands[idx4];
N_b2 = N_cands[idx2]; N_b5 = N_cands[idx5];
print("gcd(N_{b1}, N_{b4}) = ", gcd(N_b1, N_b4), " → (1,4) pair H3: ", if(gcd(N_b1,N_b4)==1,"PASS","FAIL"));
print("gcd(N_{b2}, N_{b5}) = ", gcd(N_b2, N_b5), " → (2,5) pair H3: ", if(gcd(N_b2,N_b5)==1,"PASS","FAIL"));
print("");
print("Correct glueable pair set (combining [3]-class non-gcd pairs with [1,1,1] class):");
\\ We already know from howe_sextic_twists_all15.gp that within [3]-class:
\\ (0,2),(0,3),(0,5),(2,3) all pass H3.  (3,5) fails because gcd>1 regardless.
\\ The only uncertainty was between (1,4) and (2,5).
if (gcd(N_b2,N_b5)==1,
    print("  Glueable: {(0,2),(0,3),(0,5),(2,3),(2,5)}  [as Python found]")
);
if (gcd(N_b1,N_b4)==1,
    print("  Glueable: {(0,2),(0,3),(0,5),(1,4),(2,3)}  [as PARI CM-formula found]")
);
print("");
print("================================================================");
print("Done: Thread 18 trace verification complete.");
print("================================================================");
