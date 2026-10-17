\\ Decision rule of PREREGISTRATION.md on out/features.jsonl.
\\ Prints a table and writes out/analysis.json.
OUT = "research/nist_seed_skip_filter_20261017/out/";
\\ value of "key" in a flat JSON line (numbers, "strings" or [a,b,...])
jget(s, key) = {
  my(k = Str("\"", key, "\":"), c = Vec(s), kc = Vec(k), i = 0, j, d = 0, v = "");
  for (a = 1, #c - #kc + 1, my(ok = 1); for (b = 1, #kc, if (c[a + b - 1] != kc[b], ok = 0; break)); if (ok, i = a + #kc; break));
  if (!i, error("missing ", key));
  j = i; while (j <= #c, if (c[j] == "[", d++); if (c[j] == "]", d--); if (d == 0 && (c[j] == "," || c[j] == "}"), break); j++);
  for (a = i, j - 1, v = Str(v, c[a]));
  if (Vec(v)[1] == "\"", Vec(v)[2], eval(v));
}
{
L = readstr(Str(OUT, "features.jsonl"));
N = #L;
lab = vector(N, r, jget(L[r], "label") == "P");
pb = vector(N, r, jget(L[r], "p_bits"));
\\ interval threshold: 1 / 0 decided, -1 undecided
thr(iv, x) = if (iv[1] >= x, 1, if (iv[2] < x, 0, -1));
lt = vector(N, r, jget(L[r], "lpf_twist")); lm = vector(N, r, jget(L[r], "lpf_nm1")); lq = vector(N, r, jget(L[r], "lpf_np1"));
B = List();  \\ [name, vector of 0/1/-1]
addb(nm, v) = listput(B, [nm, v]);
addb("H7 D has no l^2, l<1e6", vector(N, r, #jget(L[r], "D_sq_primes") == 0));
addb("H8 lpf(twist) >= 160", vector(N, r, thr(lt[r], 160)));
addb("H8 lpf(twist) >= pbits-32", vector(N, r, thr(lt[r], pb[r] - 32)));
addb("H8 twist cofactor <= 32 bits", vector(N, r, my(c = jget(L[r], "twist_cofactor_bits")); if (c < 0, -1, c <= 32)));
addb("H9 lpf(n-1) >= 160", vector(N, r, thr(lm[r], 160)));
addb("H9 lpf(n-1) >= pbits-32", vector(N, r, thr(lm[r], pb[r] - 32)));
addb("H9 lpf(n+1) >= 160", vector(N, r, thr(lq[r], 160)));
addb("H9 lpf(n+1) >= pbits-32", vector(N, r, thr(lq[r], pb[r] - 32)));
addb("H10 (p/n) = 1", vector(N, r, jget(L[r], "kr_p_n") == 1));
addb("H10 (-p/n) = 1", vector(N, r, jget(L[r], "kr_mp_n") == 1));
addb("H10 p primitive root mod n", vector(N, r, jget(L[r], "emb_ratio") == 1));
addb("H10 (n-1)/k <= 2", vector(N, r, jget(L[r], "emb_ratio") <= 2));
foreach (["kr_m1_n", "kr_2_n", "kr_3_n", "kr_t_n", "kr_nt_n", "kr_n_p", "kr_t_p", "kr_D_p", "kr_nt_p"], f,
  addb(Str("H11 ", f, " = 1"), vector(N, r, jget(L[r], f) == 1)));
B = Vec(B);
perfect(v, y) = { my(a = 1, b = 1); for (r = 1, N, if (v[r] < 0, return(0)); if (v[r] != y[r], a = 0); if (v[r] == y[r], b = 0)); a || b; };
cnt(v, y, s) = sum(r = 1, N, v[r] == 1 && y[r] == s);
print("binary feature                      P:true/10  S:true/11  undecided  perfect");
res = List();
for (k = 1, #B, my(v = B[k][2], u = #select(x -> x < 0, v), pf = perfect(v, lab));
  printf("%-36s %5d      %5d      %5d      %s\n", B[k][1], cnt(v, lab, 1), cnt(v, lab, 0), u, if (pf, "YES", "no"));
  listput(res, Str("{\"feature\":\"", B[k][1], "\",\"P_true\":", cnt(v, lab, 1), ",\"S_true\":", cnt(v, lab, 0), ",\"undecided\":", u, ",\"perfect\":", pf, "}")));
\\ continuous features (undecided intervals: midpoint, flagged)
C = [["H8 lpf(twist)/pbits", vector(N, r, (lt[r][1] + lt[r][2]) / 2 / pb[r])],
     ["H9 lpf(n-1)/pbits", vector(N, r, (lm[r][1] + lm[r][2]) / 2 / pb[r])],
     ["H9 lpf(n+1)/pbits", vector(N, r, (lq[r][1] + lq[r][2]) / 2 / pb[r])],
     ["H7 #{l<1e6: l^2|D}", vector(N, r, #jget(L[r], "D_sq_primes"))],
     ["POSCTRL n mod 2^20", vector(N, r, jget(L[r], "n_low20"))]];
\\ Mann-Whitney U (ties 1/2) of group y=1 vs y=0
U(x, y) = sum(i = 1, N, if (y[i], sum(j = 1, N, if (!y[j], (x[i] > x[j]) + (x[i] == x[j]) / 2))));
splits(x, y) = { my(a = vecmax(vector(N, r, if (y[r], x[r], -oo))), b = vecmin(vector(N, r, if (y[r], x[r], oo))),
  c = vecmax(vector(N, r, if (!y[r], x[r], -oo))), d = vecmin(vector(N, r, if (!y[r], x[r], oo)))); a < d || c < b; };
\\ exact two-sided null of U over all C(21,10) labellings
Unull(x) = { my(m = sum(r = 1, N, lab[r]), h = List()); forsubset([N, m], s, my(y = vector(N)); for (k = 1, #s, y[s[k]] = 1); listput(h, U(x, y))); Vec(h); };
print("\ncontinuous feature        median P   median S   U     exact p (2-sided)   perfect threshold");
cres = List();
for (k = 1, #C - 1, my(x = C[k][2], u = U(x, lab), m = sum(r = 1, N, lab[r]) * (N - sum(r = 1, N, lab[r])) / 2, h = Unull(x),
    pv = #select(z -> abs(z - m) >= abs(u - m), h) / #h);
  my(xp = vecsort(vector(N, r, if (lab[r], x[r], oo))), xq = vecsort(vector(N, r, if (!lab[r], x[r], oo))));
  printf("%-24s %9.3f  %9.3f   %5.1f  %.4f              %s\n", C[k][1], xp[5], xq[6], u, pv, if (splits(x, lab), "YES", "no"));
  listput(cres, Str("{\"feature\":\"", C[k][1], "\",\"U\":", u, ",\"p_exact\":", pv * 1., ",\"p_bonferroni\":", min(1, pv * (#C - 1)) * 1., ",\"perfect\":", splits(x, lab), "}")));
\\ positive control: plant P := 10 smallest n mod 2^20
x = C[#C][2]; o = vecsort(x, , 1); yp = vector(N); for (k = 1, 10, yp[o[k]] = 1);
pc = splits(x, yp) && !splits(x, lab);
print("\npositive control (planted label on n mod 2^20): detected = ", splits(x, yp), "; same feature on real labels splits = ", splits(x, lab));
\\ negative control: random relabellings
setrand(20261017); T = 100000; hit = 0;
for (it = 1, T, my(s = numtoperm(N, random(N!)), y = vector(N)); for (k = 1, 10, y[s[k]] = 1);
  if (#select(b -> perfect(b[2], y), B) > 0, hit++));
print("negative control: ", hit, " of ", T, " random labellings split perfectly on some binary feature (FWER ", hit * 1. / T, ")");
write(Str(OUT, "analysis.json"), Str("{\"binary\":[", strjoin(Vec(res), ","), "],\"continuous\":[", strjoin(Vec(cres), ","),
  "],\"positive_control_detected\":", splits(x, yp), ",\"negative_control\":{\"labellings\":", T, ",\"any_perfect\":", hit, "}}"));
}
quit;
