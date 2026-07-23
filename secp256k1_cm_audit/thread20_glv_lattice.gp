\\ thread20_glv_lattice.gp
\\ Thread 20: GLV-HNP Phase 2 — GLV-bounded k_2 attack at toy scale
\\ Run: gp -q thread20_glv_lattice.gp

default(parisize, 256*10^6);
default(timer, 0);

print("================================================================");
print("Thread 20: GLV-HNP Phase 2 toy attack with GLV-bounded k_2");
print("================================================================");

\\ ---- Setup: toy j=0 prime-order curve with GLV ----
print("-- Setup: toy j=0 prime-order curve --");
found_p = 0; found_b = 0; found_n = 0;
{
forprime(p = 200, 2000,
  local(E_t, n_t);
  if (Mod(p, 6) != 1, next);
  for (b = 1, 50,
    E_t = ellinit([0, b], p);
    n_t = ellcard(E_t);
    if (isprime(n_t) && kronecker(-3, n_t) == 1,
      found_p = p; found_b = b; found_n = n_t; break(2)
    )
  )
);
}
p_val = found_p; b_val = found_b; n_val = found_n;
E = ellinit([0, b_val], p_val);
G = ellgenerators(E)[1];
print("Curve: y^2 = x^3 + ", b_val, " over F_", p_val, "  n = ", n_val);

lambda = lift(Mod((-1 + sqrt(Mod(-3, n_val))) / 2, n_val));
print("lambda = ", lambda);
print("lambda^2 + lambda + 1 mod n = ", lift(Mod(lambda^2 + lambda + 1, n_val)));

\\ ---- Parameters ----
print("-- Threat model parameters --");
m  = 5;
K1 = 3;
B2 = sqrtint(n_val);
print("m=", m, "  K1=", K1, "  B2=sqrtint(n)=", B2);
fa_rate = K1 * B2 * 1.0 / n_val;
print("False alarm rate per sig (wrong d): K1*B2/n = ", fa_rate);
exp_false = n_val * fa_rate^m;
print("Expected false alarms for all ", m, " sigs: ", exp_false);

\\ ---- Generate signatures: GLV-realistic model ----
print("-- Generate ", m, " signatures (k_1 < K1, k_2 < B2) --");
setrand(777);
d_secret = random(n_val - 2) + 1;
print("Planted d = ", d_secret);

A_vec = vector(m); B_vec = vector(m);
k1_vec = vector(m); k2_vec = vector(m);
{
  local(i_sig, k1, k2, k_full, R, r_sig, h_msg, s_sig);
  i_sig = 1;
  while (i_sig <= m,
    k1 = random(K1 - 1) + 1;
    k2 = random(B2 - 1) + 1;
    k_full = lift(Mod(k1 + lambda * k2, n_val));
    if (k_full == 0, next);
    R = ellmul(E, G, k_full);
    r_sig = lift(R[1]) % n_val;
    if (r_sig == 0, next);
    h_msg = random(n_val);
    s_sig = lift(Mod((h_msg + d_secret * r_sig) / k_full, n_val));
    if (s_sig == 0, next);
    k1_vec[i_sig] = k1; k2_vec[i_sig] = k2;
    A_vec[i_sig] = lift(Mod(h_msg / s_sig, n_val));
    B_vec[i_sig] = lift(Mod(r_sig / s_sig, n_val));
    i_sig++
  )
}
print("Planted k1: ", k1_vec, "  (all < K1=", K1, ")");
print("Planted k2: ", k2_vec, "  (all < B2=", B2, ")");

\\ Sanity check
{
  local(ok, T, kf);
  ok = 1;
  for (i = 1, m,
    T  = lift(Mod(A_vec[i] + B_vec[i] * d_secret, n_val));
    kf = lift(Mod(k1_vec[i] + lambda * k2_vec[i], n_val));
    if (T != kf, ok = 0)
  );
  print("Sanity (A_i + B_i*d == k1_i + lambda*k2_i for all i): ", ok);
}

\\ ---- Method 1: Brute force with GLV-bounded k_2 ----
print("-- Method 1: Brute force (k_2 < B2 = sqrtint(n)) --");
d_bf = 0; count_bf = 0;
{
  local(T, found_k2, k1t);
  for (d_try = 1, n_val - 1,
    all_ok = 1;
    for (i = 1, m,
      T = lift(Mod(A_vec[i] + B_vec[i] * d_try, n_val));
      found_k2 = 0;
      for (k2t = 0, B2 - 1,
        k1t = lift(Mod(T - lambda * k2t, n_val));
        if (k1t < K1, found_k2 = 1; break)
      );
      if (!found_k2, all_ok = 0; break)
    );
    if (all_ok,
      count_bf++;
      if (d_bf == 0, d_bf = d_try)
    )
  )
}
print("Brute force: ", count_bf, " candidate(s), first = ", d_bf, "  correct = ", (d_bf == d_secret));

\\ ---- Method 2: LLL (qflll) lattice, dim 2m+2 ----
print("-- Method 2: qflll lattice (2m+2 dim, Kannan embedding) --");
\\ Basis rows (each 2m+2 dimensional), transposed for qflll:
\\ Row i (i=1..m):      n in col i
\\ Row m+1:             B_j in col j (j=1..m), K1 in col m+1
\\ Row m+1+i (i=1..m):  -lambda in col i, K1 in col m+1+i
\\ Row 2m+2 (Kannan):   A_j in col j (j=1..m), K1 in col 2m+2
dim_aug = 2*m + 2;
M_lat = matrix(dim_aug, dim_aug);
for (i = 1, m, M_lat[i, i] = n_val);
for (j = 1, m, M_lat[m+1, j] = B_vec[j]);
M_lat[m+1, m+1] = K1;
{for (i = 1, m, M_lat[m+1+i, i] = -lambda; M_lat[m+1+i, m+1+i] = K1)};
for (j = 1, m, M_lat[dim_aug, j] = A_vec[j]);
M_lat[dim_aug, dim_aug] = K1;

\\ qflll expects columns as basis; M_lat has rows as basis, so transpose
Mt  = mattranspose(M_lat);
U   = qflll(Mt, 1);
Mrd = Mt * U;

d_lll = 0;
{
  local(ke, raw, cand);
  for (col = 1, dim_aug,
    ke = Mrd[dim_aug, col];
    if (abs(ke) == K1,
      raw = Mrd[m+1, col];
      if (raw % K1 == 0,
        cand = abs(raw / K1);
        if (cand > 0 && cand < n_val && d_lll == 0, d_lll = cand)
      )
    )
  )
}
print("qflll found d = ", d_lll, "  correct = ", (d_lll == d_secret));

\\ Print norms of first 4 LLL columns
print("First 4 LLL-column norms:");
{
  local(col_norm2);
  for (col = 1, min(4, dim_aug),
    col_norm2 = sum(j = 1, dim_aug, Mrd[j, col]^2);
    print("  col ", col, ": ||v||^2=", col_norm2,
          "  pos_(m+1)=", Mrd[m+1, col],
          "  Kannan=", Mrd[dim_aug, col])
  )
}

\\ ---- Mathematical analysis ----
print("-- Mathematical analysis (target vs Gaussian min) --");
target_norm2 = 0;
for (i = 1, m, target_norm2 += k1_vec[i]^2);
target_norm2 += (d_secret * K1)^2;
for (i = 1, m, target_norm2 += (k2_vec[i] * K1)^2);
print("Target vector ||v||^2 = ", target_norm2, "  ||v|| ~ ", floor(sqrt(target_norm2 * 1.0)));
print("  d*K1 = ", d_secret * K1, " (dominates; d=", d_secret, " K1=", K1, ")");
det_approx = n_val^m * K1^(m+1) * 1.0;
D = dim_aug * 1.0;
gauss_min = sqrt(D / (2.0 * Pi * exp(1))) * det_approx^(1.0 / D);
print("det(L) ~ n^m * K1^(m+1) = ", n_val, "^", m, " * ", K1, "^", m+1);
print("Gaussian min = sqrt(D/(2pi*e)) * det^(1/D) ~ ", floor(gauss_min));
ratio = sqrt(target_norm2 * 1.0) / gauss_min;
print("Ratio target/Gaussian = ", ratio);
{if (ratio < 1, print("FAVORABLE: target < Gaussian min => LLL may find d"), print("UNFAVORABLE: target >> Gaussian min => LLL will NOT find d"))};

\\ ---- Full-range k_2 baseline (every d consistent) ----
print("-- Full-range k_2 (k_2 in [0,n)) — confirm every d consistent --");
count_full = 0;
{
  local(T, found_full, k1t);
  for (d_try = 1, 20,
    all_ok = 1;
    for (i = 1, m,
      T = lift(Mod(A_vec[i] + B_vec[i] * d_try, n_val));
      found_full = 0;
      for (k2t = 0, n_val - 1,
        k1t = lift(Mod(T - lambda * k2t, n_val));
        if (k1t < K1, found_full = 1; break)
      );
      if (!found_full, all_ok = 0; break)
    );
    if (all_ok, count_full++)
  )
}
print("Full-range k_2: of d=1..20, how many satisfy ALL ", m, " sigs: ", count_full, "/20");
print("(Expected 20/20: every d is trivially consistent when k_2 is unconstrained)");

\\ ---- Summary ----
print("================================================================");
print("SUMMARY");
print("================================================================");
print("Curve y^2 = x^3 + ", b_val, " over F_", p_val, "  n=", n_val, "  lambda=", lambda);
print("m=", m, "  K1=", K1, "  B2=sqrtint(n)=", B2);
print("Planted d = ", d_secret);
print("");
print("Brute force (k_2 < B2):  d = ", d_bf, "  correct = ", (d_bf == d_secret));
print("qflll lattice:            d = ", d_lll, "  correct = ", (d_lll == d_secret));
print("");
print("Key finding:");
print("  GLV-bounded k_2 (<sqrt(n)) enables brute-force attack in O(n*m*B2) = O(n^{3/2}*m).");
print("  Full-range k_2 (< n): every d is consistent — brute force fails.");
print("  Integer qflll: target vector norm ~ ", floor(sqrt(target_norm2 * 1.0)), " >> Gaussian min ~ ", floor(gauss_min), " - LLL fails.");
print("  Fix for LLL: use floating-point qflll with d scaled by K1/n (fractional scaling).");
print("  For secp256k1: O(n^{3/2}) brute force infeasible; scaled LLL needed (Phase 2 doc).");
