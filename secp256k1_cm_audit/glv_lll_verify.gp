\\ Minimal test: verify PARI qflll recovers d for n=199, lam=106, K1=2, K2=15
\\ This replicates the Python attack.py success case in PARI.
\\ If this works, the PARI LLL implementation is sound.
default(parisize, 256000000);
default(timer, 0);

p_c = 211; b_c = 2; n_c = 199; lam_c = 106;
K1 = 2; K2 = 15;
S_K1 = n_c \ K1;
S_K2 = n_c \ K2;
S_KANNAN = n_c;

print("n=199, lam=106, lam/n=0.5327, K1=2, K2=15");
print("S_K1=", S_K1, " S_K2=", S_K2, " S_KANNAN=", S_KANNAN);
E = ellinit([0, b_c], p_c);

\\ Find generator (first x where x^3+b is QR mod p)
G_c = 0;
for (x_try = 0, p_c,
  rhs = lift(Mod(x_try^3 + b_c, p_c));
  if (!issquare(Mod(rhs, p_c)), next);
  y_try = lift(sqrt(Mod(rhs, p_c)));
  G_c = [x_try, y_try];
  break
);
print("Generator G=", G_c);
if (ellmul(E, G_c, n_c) != [0], error("bad generator"));

success_total = 0;
fail_total = 0;

for (m_sigs = 2, 8,
  local(wins);
  wins = 0;
  for (seed_val = 1, 5,
    local(A_vec, B_vec, k1_vec, k2_vec, i, k1, k2, k_full);
    local(R, r_sig, h_msg, s_sig, s_inv, M, dim, U, Red, d_cand, found);
    setrand(seed_val * 100 + m_sigs);
    d_secret = random(n_c - 2) + 1;
    A_vec = vector(m_sigs); B_vec = vector(m_sigs);
    k1_vec = vector(m_sigs); k2_vec = vector(m_sigs);
    i = 1;
    while (i <= m_sigs,
      k1 = random(K1); k2 = random(K2);
      k_full = lift(Mod(k1 + lam_c*k2, n_c));
      if (k_full == 0, next);
      R = ellmul(E, G_c, k_full);
      r_sig = lift(R[1]) % n_c;
      if (r_sig == 0, next);
      h_msg = random(n_c);
      s_sig = lift(Mod((h_msg + d_secret * r_sig) / k_full, n_c));
      if (s_sig == 0, next);
      s_inv = lift(Mod(1/s_sig, n_c));
      A_vec[i] = lift(Mod(h_msg * s_inv, n_c));
      B_vec[i] = lift(Mod(r_sig * s_inv, n_c));
      k1_vec[i] = k1; k2_vec[i] = k2;
      i += 1
    );
    dim = 2*m_sigs + 2;
    M = matrix(dim, dim);
    for (i = 1, m_sigs, M[i, i] = n_c * S_K1);
    for (i = 1, m_sigs, M[m_sigs+1, i] = B_vec[i] * S_K1);
    M[m_sigs+1, m_sigs+1] = 1;
    for (i = 1, m_sigs,
      M[m_sigs+1+i, i] = -lam_c * S_K1;
      M[m_sigs+1+i, m_sigs+1+i] = S_K2
    );
    for (i = 1, m_sigs, M[2*m_sigs+2, i] = A_vec[i] * S_K1);
    M[2*m_sigs+2, dim] = S_KANNAN;
    U = qflll(M~);
    Red = U~ * M;
    found = 0;
    for (r = 1, dim,
      last = Red[r, dim];
      if (abs(last) != S_KANNAN, next);
      sign_v = last / S_KANNAN;
      d_cand = lift(Mod(sign_v * Red[r, m_sigs+1], n_c));
      if (d_cand == d_secret, found = 1; break)
    );
    wins += found
  );
  print("m=", m_sigs, ": ", wins, "/5 recovered")
);
