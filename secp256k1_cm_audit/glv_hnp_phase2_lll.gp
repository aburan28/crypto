\\ ============================================================
\\ Phase 2 GLV-HNP: actual LLL implementation (Thread 20)
\\ ============================================================
\\ Toy: n=199, lambda=106, K1_BOUND=2, d_secret=190, 5 signatures.
\\ k_1 values known from prior toy run: all k1=1.
\\ We derive A, B from the signature equation directly.
\\ Run: gp -q glv_hnp_phase2_lll.gp

default(parisize, 64000000);
default(timer, 0);

print("================================================================");
print("Phase 2 GLV-HNP: LLL implementation (Thread 20)");
print("================================================================");
print("");

nv = 199; lam = 106; K1B = 2; d_sec = 190; msigs = 5;

\\ Planted k values from prior toy run output:
\\ sig i: k_1=1, k_2=[99, 192, 160, 47, 179]
k1v = [1, 1, 1, 1, 1];
k2v = [99, 192, 160, 47, 179];

\\ Compute k_full = k1 + lambda*k2 mod n, then A,B for each sig.
\\ Signature equation: k_full = A + B*d_sec (mod n)
\\ We need explicit (A, B) pairs. They come from the signing operation.
\\ Since we don't have the actual sig (r,s,h) values (the toy ran with
\\ ellmul which hangs here), we instead PLANT the HNP coefficients
\\ directly: choose random h, compute r from k_full*G, derive A and B.
\\ This is valid since the lattice attack only uses (A, B, n, lambda, K1B).

\\ Plant known k_full values
kfv = vector(msigs, ii, lift(Mod(k1v[ii] + lam*k2v[ii], nv)));
print("k_full values: ", kfv);

\\ Choose deterministic h values and r values (simulate signatures)
\\ For the HNP lattice we just need: k_full = A_i + B_i * d  (mod n)
\\ We can freely choose any (A_i, B_i) satisfying this for each i,
\\ e.g. A_i = k_full - B_i*d_sec mod n for any random B_i.
\\ Use B_i = i+10 for simplicity (nonzero, bounded away from 0).
Bv = vector(msigs, ii, ii + 10);
Av = vector(msigs, ii, lift(Mod(kfv[ii] - Bv[ii]*d_sec, nv)));

print("Simulation: (A_i, B_i) pairs:");
for(ii = 1, msigs,
    local(ck);
    ck = lift(Mod(Av[ii] + Bv[ii]*d_sec - lam*k2v[ii] - k1v[ii], nv));
    print("  sig", ii, ": A=", Av[ii], " B=", Bv[ii], " check=", ck);
);
print("");

\\ ---- Build (2*msigs+2) x (2*msigs+2) Kannan lattice ----
\\ Variables in short vector: (k1_1..k1_m, d, k2_1..k2_m, 1)
\\ Basis rows:
\\   row i  (1<=i<=msigs):   n*e_i in k1 block
\\   row msigs+1 (d-row):    (Bv[1..m], K1B, 0..0)
\\   row msigs+1+j (1<=j<=m): (lam*e_j in k1 block, 0, n*e_j in k2 block)
\\ Kannan augmented row:     (-Av[1..m], 0..0, 1)

dm  = 2*msigs + 1;
dmK = dm + 1;

MB = matrix(dmK, dmK);
for(ii = 1, msigs, MB[ii, ii] = nv);
for(ii = 1, msigs, MB[msigs+1, ii] = Bv[ii]);
MB[msigs+1, msigs+1] = K1B;
for(jj = 1, msigs,
    MB[msigs+1+jj, jj] = lam;
    MB[msigs+1+jj, msigs+1+jj] = nv;
);
for(ii = 1, msigs, MB[dmK, ii] = -Av[ii]);
MB[dmK, dmK] = 1;

print("Basis (", dmK, "x", dmK, ") built. Running qflll...");
LR = qflll(MB~)~;
RM = LR * MB;
print("LLL done.");
print("");

\\ ---- Scan for d_secret ----
print("Scanning rows with |last_coord|=1:");
found = 0;
for(rr = 1, dmK,
    local(rowv, last_c, de, dc1, dc2);
    rowv = RM[rr,];
    last_c = rowv[dmK];
    if (abs(last_c) != 1, next);
    de = rowv[msigs+1];
    \\ d-column: row m+1 has K1B in that slot, so d contributes d*K1B to col m+1.
    \\ But Kannan shifts mean: the short vector entry at col m+1 = d * K1B.
    \\ Try de / K1B (integer division check):
    if (K1B != 0 && de % K1B == 0,
        dc1 = lift(Mod(de \ K1B, nv));
        dc2 = lift(Mod(-(de \ K1B), nv));
        print("  row ", rr, ": last=", last_c, "  d_col=", de,
              "  d+cand=", dc1, "  d-cand=", dc2,
              if(dc1==d_sec || dc2==d_sec, "  *** MATCH ***", ""));
        if (dc1 == d_sec || dc2 == d_sec, found = 1);
    ,
        print("  row ", rr, ": last=", last_c, "  d_col=", de, "  (not div by K1B)");
    );
    \\ Also try de directly as d
    dc3 = lift(Mod(de, nv));
    dc4 = lift(Mod(-de, nv));
    if ((dc3 == d_sec || dc4 == d_sec) && !found,
        print("  row ", rr, ": raw d=", dc3, " or", dc4, "  *** MATCH (raw) ***");
        found = 1;
    );
);

print("");
if (found,
    print("RESULT: d = ", d_sec, " RECOVERED by GLV-aware LLL.  Thread 20 SUCCESS.");
,
    print("RESULT: LLL did not recover d. Printing full reduced matrix:");
    for(rr = 1, dmK, print("  row ", rr, ": ", RM[rr,]));
);

print("");
print("================================================================");
print("Thread 20 complete.");
print("================================================================");
