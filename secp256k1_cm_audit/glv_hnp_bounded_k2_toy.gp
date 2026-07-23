\\ ============================================================
\\ GLV-aware HNP: bounded-k_2 threat model demonstration
\\ Thread 19 (2026-07-23 autolab run)
\\ ============================================================
\\
\\ glv_hnp_phase2_toy.gp used k_2 in [0, n) (unbounded):
\\   for ANY target T, exists k_1<K1_BOUND, k_2 in [0,n) with
\\   k_1 + lambda*k_2 == T mod n  =>  attack impossible.
\\
\\ This script uses k_2 in [0, sqrt(n)] (Babai GLV bound)
\\ and shows brute-force d recovery works in the correct model.
\\
\\ Run: gp -q glv_hnp_bounded_k2_toy.gp

default(parisize, 256000000);
default(timer, 0);

print("================================================================");
print("GLV-HNP Thread 19: bounded-k_2 threat model");
print("================================================================");
print("");

\\ Same toy curve as glv_hnp_phase2_toy.gp
\\ y^2 = x^3 + 2 over F_211, #E = 199 (prime), lambda = 106
E = ellinit([0, 2], 211);
p_toy = 211;
n_toy = 199;
lambda = 106;
print("Curve: y^2 = x^3 + 2  over  F_", p_toy, ",  n = ", n_toy, " (prime)");
print("lambda = ", lambda, "  (lambda/n = ", lambda * 1.0 / n_toy, ")");
print("Verify lambda^2+lambda+1 mod n = ", lift(Mod(lambda^2 + lambda + 1, n_toy)));
print("");

K1_BOUND = 5;
K2_BOUND = sqrtint(n_toy);
print("K1_BOUND = ", K1_BOUND);
print("K2_BOUND = floor(sqrt(n)) = ", K2_BOUND, " (Babai GLV bound)");
print("");

\\ ---- Step 1: Bias set S = {k1 + lambda*k2 mod n} ----
print("---- Step 1: Bias set S ----");
{
local(s_tmp, v);
s_tmp = [];
for (k1 = 0, K1_BOUND - 1,
    for (k2 = 0, K2_BOUND - 1,
        v = lift(Mod(k1 + lambda * k2, n_toy));
        s_tmp = concat(s_tmp, [v])
    )
);
S_vec = Set(s_tmp);
}
S_size = #S_vec;
density = S_size * 1.0 / n_toy;
print("|S| = ", S_size, " out of n = ", n_toy, " (K1*K2 = ", K1_BOUND * K2_BOUND, ")");
print("density = ", density);
print("");

\\ ---- Step 2: Uniqueness analysis ----
print("---- Step 2: Expected false positives per run ----");
{
local(m_try, fp);
m_try = 5;  fp = (n_toy-1) * density^m_try;
print("  m=5  : density^5  = ", density^m_try, "  FP = ", fp);
m_try = 8;  fp = (n_toy-1) * density^m_try;
print("  m=8  : density^8  = ", density^m_try, "  FP = ", fp);
m_try = 10; fp = (n_toy-1) * density^m_try;
print("  m=10 : density^10 = ", density^m_try, "  FP = ", fp);
m_try = 12; fp = (n_toy-1) * density^m_try;
print("  m=12 : density^12 = ", density^m_try, "  FP = ", fp);
}
m_sigs = 10;
print("Chosen: m = ", m_sigs, " signatures");
print("");

\\ ---- Step 3: Plant d, generate m signatures with bounded k1, k2 ----
print("---- Step 3: Generate ", m_sigs, " signatures ----");
setrand(42);
d_secret = random(n_toy - 2) + 1;
print("Planted d = ", d_secret);

A_vec = vector(m_sigs);
B_vec = vector(m_sigs);
k1_arr = vector(m_sigs);
k2_arr = vector(m_sigs);
kfull_arr = vector(m_sigs);
G = ellgenerators(E)[1];

{
local(i, att, k1, k2, kf, Rpt, r_s, h_m, s_s);
i = 1;
att = 0;
while (i <= m_sigs && att < 20000,
    att++;
    k1 = random(K1_BOUND);
    k2 = random(K2_BOUND);
    kf = lift(Mod(k1 + lambda * k2, n_toy));
    if (kf == 0, next);
    Rpt = ellmul(E, G, kf);
    r_s = lift(Rpt[1]) % n_toy;
    if (r_s == 0, next);
    h_m = random(n_toy);
    s_s = lift(Mod((h_m + d_secret * r_s) / kf, n_toy));
    if (s_s == 0, next);
    k1_arr[i] = k1;
    k2_arr[i] = k2;
    kfull_arr[i] = kf;
    A_vec[i] = lift(Mod(h_m / s_s, n_toy));
    B_vec[i] = lift(Mod(r_s / s_s, n_toy));
    i++
);
if (i <= m_sigs, error("not enough valid signatures"))
}
print("Generated ", m_sigs, " signatures.");
print("Sample (k1, k2, k_full):");
{
for (j = 1, m_sigs,
    print("  sig ", j, ": k1=", k1_arr[j], " k2=", k2_arr[j], " k_full=", kfull_arr[j])
)
}
print("");

\\ Sanity check
sanity = 1;
{
for (j = 1, m_sigs,
    if (lift(Mod(A_vec[j] + B_vec[j] * d_secret, n_toy)) != lift(Mod(k1_arr[j] + lambda * k2_arr[j], n_toy)), sanity = 0)
)
}
print("Sanity (A_i+B_i*d == k1_i+lambda*k2_i mod n for all i): ", sanity);
kfull_max = 0;
for (j = 1, m_sigs, if (kfull_arr[j] > kfull_max, kfull_max = kfull_arr[j]));
print("max(k_full) = ", kfull_max, "  K1_BOUND = ", K1_BOUND, "  -> std HNP fails (k_full >> K1_BOUND)");
print("");

\\ ---- Step 4: Brute-force d recovery ----
print("================================================================");
print("Step 4: Brute-force d recovery");
print("================================================================");
print("For each d' in [1,n): check A_i+B_i*d' mod n in S for all i.");
print("");

{
local(cand_list, d_cand, all_ok, j, target);
cand_list = [];
for (d_cand = 1, n_toy - 1,
    all_ok = 1;
    for (j = 1, m_sigs,
        target = lift(Mod(A_vec[j] + B_vec[j] * d_cand, n_toy));
        if (!setsearch(S_vec, target), all_ok = 0; break)
    );
    if (all_ok, cand_list = concat(cand_list, [d_cand]))
);
cands = cand_list;
}

print("Candidates found: ", #cands);
print("Candidate d values: ", cands);
if (#cands == 1 && cands[1] == d_secret, print("RESULT: UNIQUE RECOVERY CONFIRMED  d = ", d_secret, "  OK"));
if (#cands == 1 && cands[1] != d_secret, print("RESULT: WRONG d found -- bug!"));
if (#cands == 0, print("RESULT: NO candidate found -- bug!"));
if (#cands > 1, print("RESULT: AMBIGUOUS -- ", #cands, " candidates"));
print("");

\\ ---- Step 5: Unbounded k_2 comparison ----
print("================================================================");
print("Step 5: Unbounded k_2 -- attack impossible");
print("================================================================");
\\ With k_2 in [0,n), set k_1=0, k_2=target*lambda^{-1} mod n.
\\ Then k_1 + lambda*k_2 = target for ANY target. So every d passes.
lambda_inv = lift(Mod(lambda, n_toy)^(-1));
d_wrong = d_secret + 1;
if (d_wrong >= n_toy, d_wrong = 1);

pass_bnd = 0;
pass_unb = 0;
{
local(tgt, k2u);
for (j = 1, m_sigs,
    tgt = lift(Mod(A_vec[j] + B_vec[j] * d_wrong, n_toy));
    if (setsearch(S_vec, tgt), pass_bnd++);
    k2u = lift(Mod(tgt * lambda_inv, n_toy));
    if (lift(Mod(lambda * k2u, n_toy)) == tgt, pass_unb++)
)
}
print("Wrong d = ", d_wrong, " (true d = ", d_secret, "):");
print("  Bounded (k2 in [0,", K2_BOUND, ")): ", pass_bnd, " / ", m_sigs, " sigs pass");
print("  Unbounded (k2 in [0,n)): ", pass_unb, " / ", m_sigs, " sigs pass");
print("Bounded: wrong d rejected. Unbounded: any d accepted (0 info gain).");
print("");

\\ ---- Step 6: Lattice sketch ----
print("================================================================");
print("Step 6: Lattice formulation (LLL = future work)");
print("================================================================");
dim_lat = 2 * m_sigs + 2;
print("Lattice dimension: 2m+2 = ", dim_lat);
print("Target short vector: (k_{i,1},..., d*K1, k_{i,2}*K2,..., K1)");
print("  k_{i,1} entries : ~", K1_BOUND);
print("  d * K1_BOUND    : ~", n_toy * K1_BOUND, " (DOMINATES vs K2^2 = ", K2_BOUND^2, ")");
print("  k_{i,2}*K2_BOUND: ~", K2_BOUND^2);
{
local(tn, sv, det);
tn = sqrt(m_sigs * K1_BOUND^2 + (n_toy * K1_BOUND / 2)^2 * 1.0 + m_sigs * (K2_BOUND^2/2.0)^2);
det = n_toy^m_sigs * K1_BOUND * K2_BOUND^m_sigs * 1.0;
sv = det^(1.0 / dim_lat);
print("  target norm ~ ", tn, "  Minkowski ~ ", sv);
print("  Balancing: replace d by d mod (n/K2) reduces d-entry to ~ ", (n_toy \ K2_BOUND) * K1_BOUND, " ~ K2^2 = ", K2_BOUND^2);
}
print("Full balanced lattice: Thread 20 (next session).");
print("");

\\ ---- Summary ----
print("================================================================");
print("THREAD 19 SUMMARY");
print("================================================================");
print("1. BASELINE (glv_hnp_phase2_toy.gp): setup verified, lattice deferred.");
print("   k_2 in [0,n) -> attack impossible (trivial consistency for any d).");
print("");
print("2. BOUNDED-k_2 BRUTE FORCE:");
print("   Planted d = ", d_secret);
print("   |S| = ", S_size, " / ", n_toy, "  density = ", density);
print("   Expected FP at m=", m_sigs, " : ", (n_toy-1)*density^m_sigs);
print("   Candidates found: ", cands);
if (#cands == 1 && cands[1] == d_secret, print("   RESULT: UNIQUE RECOVERY CONFIRMED"));
if (#cands != 1 || cands[1] != d_secret, print("   RESULT: see above"));
print("");
print("3. THREAT MODEL CORRECTION:");
print("   GLV-HNP requires k_2 in [0,sqrt(n)] (Babai bound), NOT [0,n).");
print("   Unbounded: any d passes all m checks (S=Z_n, zero discrimination).");
print("   Bounded:   |S|/n = ", S_size, "/", n_toy, " -> unique d recovery.");
print("");
print("4. NEXT (Thread 20): balanced 22-dim PARI LLL lattice attack.");
print("   Key insight: reduce d-entry via d mod (n/K2_BOUND) transformation.");
print("================================================================");
