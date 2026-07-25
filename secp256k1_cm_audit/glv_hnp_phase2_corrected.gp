\\ Phase 2 GLV-aware HNP — corrected lattice with GLV-domain k2 bound.
\\ Key fix vs. glv_hnp_phase2_lattice.gp:
\\   Old: k2 full-range [0, n)  → lattice unbalanced, LLL fails
\\   New: k2 ∈ [0, K2_BOUND=ceil(sqrt(n))) (GLV domain) + matching column scales
\\
\\ LLL finding: planted vector (norm~451) is NOT the shortest in the lattice —
\\ the spurious vector [0,..,n,..,0] (norm=n=199) is shorter.  Direct LLL recovery
\\ therefore fails for this toy.  Brute-force EC-verification recovers d in O(n*m).
\\ For secp256k1 scale, the Python/fpylll BKZ version (glv_hnp_phase2_attack.py)
\\ succeeds at m=6 (5/5 seeds).
\\
\\ Run: gp -q glv_hnp_phase2_corrected.gp

default(parisize, 512000000);

print("================================================================");
print("Phase 2 GLV-aware HNP — corrected k2 bound + column scaling");
print("================================================================");
print("");

\\ ---- Find toy j=0 curve ----
found_p = 0; found_b = 0; found_n = 0;
{
forprime(p = 200, 2000,
  if(Mod(p,6) != 1, next);
  for(b = 1, 50,
    local(E, n);
    E = ellinit([0,b], p);
    n = ellcard(E);
    if(isprime(n) && kronecker(-3,n) == 1,
      found_p = p; found_b = b; found_n = n; break(2)
    )
  )
);
}
p_val  = found_p; b_val = found_b; n_val = found_n;
E_main = ellinit([0, b_val], p_val);
G_main = ellgenerators(E_main)[1];
sqrt_neg3 = lift(sqrt(Mod(-3, n_val)));
lam = lift(Mod(-1 + sqrt_neg3, n_val) / Mod(2, n_val));
print("Curve: y^2 = x^3 + ", b_val, " over F_", p_val);
print("n = ", n_val, "   lambda = ", lam);
print("Check: lam^2+lam+1 mod n = ", lift(Mod(lam^2+lam+1, n_val)));
print("");

\\ ---- Bounds and scales ----
K1_BOUND = 2;
K2_BOUND = ceil(sqrt(n_val * 1.0));
S_K1 = n_val \ K1_BOUND;
S_D  = 1;
S_K2 = n_val \ K2_BOUND;
S_KANNAN = n_val;
ratio = K1_BOUND * K2_BOUND * 1.0 / n_val;
m_thresh = ceil(log(n_val * 1.0) / log(1.0 / ratio));
m_sigs = m_thresh + 1;
print("K1_BOUND=", K1_BOUND, "  K2_BOUND=", K2_BOUND, " (ceil(sqrt(n)))");
print("Scales: S_K1=", S_K1, "  S_K2=", S_K2, "  S_D=", S_D, "  S_KANNAN=", S_KANNAN);
print("Info-thresh m >= ", m_thresh, "  using m=", m_sigs);
print("");

\\ ---- Single trial function (returns: 1=LLL, 2=brute-force, 0=fail) ----
trial(seed_in) = {
  local(d_s, A_s, B_s, r_s, i_s, k1t, k2t, kft, Rp, rs, hs, ss);
  local(dim, M, U, L, r, v, last, sv, dc, Ti, Rc, aok);
  setrand(seed_in);
  d_s = random(n_val - 1) + 1;
  A_s = vector(m_sigs); B_s = vector(m_sigs); r_s = vector(m_sigs);
  i_s = 1;
  while(i_s <= m_sigs,
    k1t = random(K1_BOUND - 1) + 1;
    k2t = random(K2_BOUND - 1) + 1;
    kft = lift(Mod(k1t + lam * k2t, n_val));
    if(kft == 0, next);
    Rp  = ellmul(E_main, G_main, kft);
    rs  = lift(Rp[1]) % n_val;
    if(rs == 0, next);
    hs  = random(n_val);
    ss  = lift(Mod((hs + d_s * rs) / kft, n_val));
    if(ss == 0, next);
    r_s[i_s] = rs;
    A_s[i_s] = lift(Mod(hs / ss, n_val));
    B_s[i_s] = lift(Mod(rs / ss, n_val));
    i_s++
  );

  dim = 2 * m_sigs + 2;
  M = matrix(dim, dim);
  for(i = 1, m_sigs, M[i, i] = n_val * S_K1);
  for(i = 1, m_sigs, M[m_sigs+1, i] = B_s[i] * S_K1);
  M[m_sigs+1, m_sigs+1] = S_D;
  for(i = 1, m_sigs,
    M[m_sigs+1+i, i] = -lam * S_K1;
    M[m_sigs+1+i, m_sigs+1+i] = S_K2
  );
  for(i = 1, m_sigs, M[2*m_sigs+2, i] = A_s[i] * S_K1);
  M[2*m_sigs+2, dim] = S_KANNAN;

  \\ LLL via qflll(M~): reduces columns of M~ (= rows of M)
  U = qflll(M~);
  L = U~ * M;

  \\ Check all 10 reduced rows for d recovery
  for(r = 1, dim,
    v    = L[r,];
    last = v[dim];
    if(abs(last) != S_KANNAN, next);
    sv = if(last > 0, 1, -1);
    dc = lift(Mod(sv * v[m_sigs+1] / S_D, n_val));
    if(dc == 0, next);
    aok = 1;
    for(i = 1, m_sigs,
      Ti = lift(Mod(A_s[i] + B_s[i] * dc, n_val));
      if(Ti == 0, aok = 0; break);
      Rc = ellmul(E_main, G_main, Ti);
      if(lift(Rc[1]) % n_val != r_s[i], aok = 0; break)
    );
    if(aok && dc == d_s, return(1))  \\ LLL success
  );

  \\ Brute-force fallback
  for(d_try = 1, n_val - 1,
    aok = 1;
    for(i = 1, m_sigs,
      Ti = lift(Mod(A_s[i] + B_s[i] * d_try, n_val));
      if(Ti == 0, aok = 0; break);
      Rc = ellmul(E_main, G_main, Ti);
      if(lift(Rc[1]) % n_val != r_s[i], aok = 0; break)
    );
    if(aok, return(if(d_try == d_s, 2, 0)))
  );
  return(0);
}

\\ ---- Verbose trial: seed=11 ----
setrand(11);
d_verbose = random(n_val - 1) + 1;
print("---- Verbose trial: seed=11, d=", d_verbose, " ----");
print("");
setrand(11);
A_v = vector(m_sigs); B_v = vector(m_sigs); r_v = vector(m_sigs);
k1_v = vector(m_sigs); k2_v = vector(m_sigs);
i_v = 1;
{
while(i_v <= m_sigs,
  local(k1t, k2t, kft, Rp, rs, hs, ss);
  k1t = random(K1_BOUND - 1) + 1;
  k2t = random(K2_BOUND - 1) + 1;
  kft = lift(Mod(k1t + lam * k2t, n_val));
  if(kft == 0, next);
  Rp = ellmul(E_main, G_main, kft);
  rs = lift(Rp[1]) % n_val;
  if(rs == 0, next);
  hs = random(n_val);
  ss = lift(Mod((hs + d_verbose * rs) / kft, n_val));
  if(ss == 0, next);
  k1_v[i_v] = k1t; k2_v[i_v] = k2t; r_v[i_v] = rs;
  A_v[i_v] = lift(Mod(hs / ss, n_val));
  B_v[i_v] = lift(Mod(rs / ss, n_val));
  i_v++
);
}
print("Generated signatures (k2 in GLV domain [1,", K2_BOUND, ")):");
for(j = 1, m_sigs,
  print("  sig ", j, ": k1=", k1_v[j], "  k2=", k2_v[j])
);
print("");

\\ Build lattice
dim_v = 2 * m_sigs + 2;
M_v = matrix(dim_v, dim_v);
for(i = 1, m_sigs, M_v[i, i] = n_val * S_K1);
for(i = 1, m_sigs, M_v[m_sigs+1, i] = B_v[i] * S_K1);
M_v[m_sigs+1, m_sigs+1] = S_D;
for(i = 1, m_sigs,
  M_v[m_sigs+1+i, i] = -lam * S_K1;
  M_v[m_sigs+1+i, m_sigs+1+i] = S_K2
);
for(i = 1, m_sigs, M_v[2*m_sigs+2, i] = A_v[i] * S_K1);
M_v[2*m_sigs+2, dim_v] = S_KANNAN;

\\ Planted vector
q_v = vector(m_sigs);
for(i = 1, m_sigs,
  local(unred);
  unred = A_v[i] + B_v[i] * d_verbose - lam * k2_v[i];
  q_v[i] = (unred - k1_v[i]) / n_val
);
pv = vector(dim_v, j, 0);
for(i = 1, m_sigs, for(j = 1, dim_v, pv[j] -= q_v[i] * M_v[i, j]));
for(j = 1, dim_v, pv[j] += d_verbose * M_v[m_sigs+1, j]);
for(i = 1, m_sigs, for(j = 1, dim_v, pv[j] += k2_v[i] * M_v[m_sigs+1+i, j]));
for(j = 1, dim_v, pv[j] += M_v[2*m_sigs+2, j]);
print("Planted short vector:  ", pv);
print("  Norm = ", sqrt(sum(j = 1, dim_v, pv[j]^2) * 1.0));
print("  d-slot (should be d=", d_verbose, "): ", pv[m_sigs+1]);
print("  Kannan slot (should be ", S_KANNAN, "): ", pv[dim_v]);
print("");

\\ Spurious short vector: take c_{m+1} = n_val, all other c_i chosen to zero out k1 columns
print("Spurious shorter vector exists: [0,...,0, n, 0,...,0] with norm = n = ", n_val);
print("  This is the combination c_{d-row} = n => d-slot = n, Kannan-slot = 0.");
print("  Its norm (", n_val, ") < planted norm (", sqrt(sum(j=1,dim_v,pv[j]^2)*1.0), ").");
print("  LLL finds this spurious vector first → direct LLL recovery fails.");
print("");

\\ LLL
print("LLL rows (all ", dim_v, "):");
U_v = qflll(M_v~);
L_v = U_v~ * M_v;
for(r = 1, dim_v,
  print("  row ", r, ": ", L_v[r, ], "  ||·||=",
    sqrt(sum(j = 1, dim_v, L_v[r, j]^2) * 1.0))
);
print("");

\\ Multi-seed sweep (seeds 11,22,33,44,55) — unrolled
print("---- Multi-seed sweep (5 seeds, m=", m_sigs, ") ----");
res11 = trial(11);
res22 = trial(22);
res33 = trial(33);
res44 = trial(44);
res55 = trial(55);
n_lll = (res11 == 1) + (res22 == 1) + (res33 == 1) + (res44 == 1) + (res55 == 1);
n_brf = (res11 == 2) + (res22 == 2) + (res33 == 2) + (res44 == 2) + (res55 == 2);
n_ttl = (res11 > 0) + (res22 > 0) + (res33 > 0) + (res44 > 0) + (res55 > 0);
print("  seed=11: ", if(res11==1,"LLL ✓",if(res11==2,"brute-force ✓","FAILED")));
print("  seed=22: ", if(res22==1,"LLL ✓",if(res22==2,"brute-force ✓","FAILED")));
print("  seed=33: ", if(res33==1,"LLL ✓",if(res33==2,"brute-force ✓","FAILED")));
print("  seed=44: ", if(res44==1,"LLL ✓",if(res44==2,"brute-force ✓","FAILED")));
print("  seed=55: ", if(res55==1,"LLL ✓",if(res55==2,"brute-force ✓","FAILED")));
print("  LLL direct: ", n_lll, "/5   Brute-force: ", n_brf, "/5   Total: ", n_ttl, "/5");
print("");
print("================================================================");
print("SUMMARY");
print("  K1_BOUND=", K1_BOUND, "  K2_BOUND=", K2_BOUND, " (GLV domain, ~sqrt(n))");
print("  Column scales: S_K1=", S_K1, "  S_K2=", S_K2, "  S_KANNAN=", S_KANNAN);
print("  LLL direct recovery: ", n_lll, "/5");
print("  Brute-force recovery: ", n_brf, "/5");
print("  Total (LLL+BF): ", n_ttl, "/5");
print("");
print("  Root cause of LLL failure:");
print("    Spurious lattice vector [0,...,n,...,0] has norm=n=", n_val);
print("    Planted vector has norm ~sqrt(m*K1^2 + n^2 + m*(K2/2)^2) > n");
print("    Fix: BKZ at larger block size (fpylll), or restate as CVP with Babai.");
print("    Python/fpylll BKZ version (glv_hnp_phase2_attack.py): 5/5 at m=6.");
print("================================================================");
