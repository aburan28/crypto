\\ Features H7-H11 of PREREGISTRATION.md for the 21 accepted NIST seed candidates.
\\ Input: the recorded census orders (../nist_seed_provenance_20261016/census).
\\ Output: one JSON object per candidate on stdout.
default(parisizemax, 2^31);
CENSUS = "research/nist_seed_provenance_20261016/census/";
{RUNS = [
  [863711343, 0x7fffffffffffffffffffffff7fffffffffff8000000000007fffffffffff, [36,1562,2185,2320,2491,2682,2863,3048], [2185,2320,3048]],
  [863966734, 0xfffffffffffffffffffffffffffffffeffffffffffffffff, [139,2198,5228,5405,6013], [139,5228,6013]],
  [871002724, 2^521-1, [1136], [1136]],
  [871696354, 0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff, [876,1842], [1842]],
  [872952608, 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffeffffffff0000000000000000ffffffff, [2688], [2688]],
  [891607927, 0xffffffffffffffffffffffffffffffff000000000000000000000001, [777,1485,5736,7026], [7026]]];}

\\ Largest prime factor bit length as [lo, hi]: lo = bits of the largest proven
\\ prime factor, hi = max(lo, bits of an unfactored composite cofactor).
lpf(N) = {
  my(F, done = 0, lo = 0, hi = 0);
  F = alarm(600, factor(N));
  if (type(F) == "t_ERROR", F = factor(N, 2^24); done = 1);
  for (k = 1, #F~,
    my(q = F[k,1], b = #binary(q));
    if (ispseudoprime(q), lo = max(lo, b); hi = max(hi, b), hi = max(hi, b)));
  [lo, hi, !done];  \\ third entry: 1 = fully factored
}
\\ odd primes l < 10^6 with l^2 | D
sqdiv(D) = { my(c = List()); forprime(l = 3, 10^6, if (D % l^2 == 0, listput(c, l))); Vec(c); }

kr(x, m) = kronecker(x, m);
J(v) = if (type(v) == "t_VEC", Str("[", strjoin(apply(x -> Str(x), v), ","), "]"), Str(v));

{
if (type(RSEL) == "t_POL", RSEL = [1..#RUNS]);
for (r = 1, #RUNS, if (!setsearch(Set(RSEL), r), next);
  my(run = RUNS[r][1], P = RUNS[r][2], L = readstr(Str(CENSUS, "orders_", run, ".txt")), ord = Map());
  for (k = 1, #L, my(v = apply(eval, strsplit(L[k], " "))); mapput(ord, v[1], v[2]));
  foreach (RUNS[r][3], i,
    my(n = mapget(ord, i), t = P + 1 - n, D = t^2 - 4*P, nt = 2*P + 2 - n, pb = #binary(P));
    my(sq = sqdiv(D), Ft = lpf(nt), Fm = lpf(n - 1), Fp = lpf(n + 1));
    my(fn1 = factor(n - 1), k = znorder(Mod(P, n), [n - 1, fn1]));
    print("{\"run\":", run, ",\"i\":", i, ",\"label\":\"", if (setsearch(Set(RUNS[r][4]), i), "P", "S"), "\"",
      ",\"p_bits\":", pb, ",\"n_low20\":", n % 2^20,
      ",\"D_sq_primes\":", J(sq),
      ",\"lpf_twist\":", J(Ft), ",\"lpf_nm1\":", J(Fm), ",\"lpf_np1\":", J(Fp),
      ",\"twist_cofactor_bits\":", if (Ft[3], #binary(nt) - Ft[1], -1),
      ",\"kr_p_n\":", kr(P, n), ",\"kr_mp_n\":", kr(-P, n),
      ",\"emb_ratio\":", J((n - 1) / k),
      ",\"kr_m1_n\":", kr(-1, n), ",\"kr_2_n\":", kr(2, n), ",\"kr_3_n\":", kr(3, n),
      ",\"kr_t_n\":", kr(t, n), ",\"kr_nt_n\":", kr(nt, n),
      ",\"kr_n_p\":", kr(n, P), ",\"kr_t_p\":", kr(t, P), ",\"kr_D_p\":", kr(D, P), ",\"kr_nt_p\":", kr(nt, P),
      "}")));
}
quit;
