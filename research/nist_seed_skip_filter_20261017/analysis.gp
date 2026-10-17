\\ Decision rule of PREREGISTRATION.md on out/features.jsonl.
\\ Prints a table and writes out/analysis.json.
OUT = "research/nist_seed_skip_filter_20261017/out/";

\\ value of "key" in a flat JSON line (numbers, "strings" or [a,b,...])
jget(s, key) =
{
  my(k = Str("\"", key, "\":"), c = Vec(s), kc = Vec(k), i = 0, j, d = 0, v = "", ok);
  for (a = 1, #c - #kc + 1, ok = 1; for (b = 1, #kc, if (c[a + b - 1] != kc[b], ok = 0; break)); if (ok, i = a + #kc; break));
  if (!i, error("missing ", key));
  j = i;
  while (j <= #c, if (c[j] == "[", d++); if (c[j] == "]", d--); if (d == 0 && (c[j] == "," || c[j] == "}"), break); j++);
  for (a = i, j - 1, v = Str(v, c[a]));
  if (Vec(v)[1] == "\"", Vec(v)[2], eval(v));
}
\\ interval threshold: 1 / 0 decided, -1 undecided
thr(iv, x) = if (iv[1] >= x, 1, if (iv[2] < x, 0, -1));
\\ v splits labels y perfectly (an undecided entry blocks a split)
perfect(v, y) =
{
  my(a = 1, b = 1);
  for (r = 1, #y, if (v[r] < 0, return(0)); if (v[r] != y[r], a = 0); if (v[r] == y[r], b = 0));
  a || b;
}
cnt(v, y, s) = sum(r = 1, #y, v[r] == 1 && y[r] == s);
\\ Mann-Whitney U (ties 1/2) of group y=1 against y=0
U(x, y) = sum(i = 1, #y, if (y[i], sum(j = 1, #y, if (!y[j], (x[i] > x[j]) + (x[i] == x[j]) / 2))));
\\ some threshold on x separates y perfectly
splits(x, y) =
{
  my(N = #y, a, b, c, d);
  a = vecmax(vector(N, r, if (y[r], x[r], -oo))); b = vecmin(vector(N, r, if (y[r], x[r], oo)));
  c = vecmax(vector(N, r, if (!y[r], x[r], -oo))); d = vecmin(vector(N, r, if (!y[r], x[r], oo)));
  a < d || c < b;
}
\\ exact two-sided p of U over all labellings with the same group size
Up(x, y) =
{
  my(N = #y, m = sum(r = 1, N, y[r]), mid = m * (N - m) / 2, u = U(x, y), tot = 0, ge = 0, z);
  forsubset([N, m], s, z = vector(N); for (k = 1, #s, z[s[k]] = 1); tot++; if (abs(U(x, z) - mid) >= abs(u - mid), ge++));
  ge / tot;
}
med(v) = my(w = vecsort(v), n = #w); if (n % 2, w[(n + 1) / 2], (w[n / 2] + w[n / 2 + 1]) / 2);

L = readstr(Str(OUT, "features.jsonl"));
N = #L;
lab = vector(N, r, jget(L[r], "label") == "P");
pb = vector(N, r, jget(L[r], "p_bits"));
lt = vector(N, r, jget(L[r], "lpf_twist"));
lm = vector(N, r, jget(L[r], "lpf_nm1"));
lq = vector(N, r, jget(L[r], "lpf_np1"));
B = List();
listput(B, ["H7 D has no l^2, l<1e6", vector(N, r, #jget(L[r], "D_sq_primes") == 0)]);
listput(B, ["H8 lpf(twist) >= 160", vector(N, r, thr(lt[r], 160))]);
listput(B, ["H8 lpf(twist) >= pbits-32", vector(N, r, thr(lt[r], pb[r] - 32))]);
listput(B, ["H8 twist cofactor <= 32 bits", vector(N, r, my(c = jget(L[r], "twist_cofactor_bits")); if (c < 0, -1, c <= 32))]);
listput(B, ["H9 lpf(n-1) >= 160", vector(N, r, thr(lm[r], 160))]);
listput(B, ["H9 lpf(n-1) >= pbits-32", vector(N, r, thr(lm[r], pb[r] - 32))]);
listput(B, ["H9 lpf(n+1) >= 160", vector(N, r, thr(lq[r], 160))]);
listput(B, ["H9 lpf(n+1) >= pbits-32", vector(N, r, thr(lq[r], pb[r] - 32))]);
listput(B, ["H10 (p/n) = 1", vector(N, r, jget(L[r], "kr_p_n") == 1)]);
listput(B, ["H10 (-p/n) = 1", vector(N, r, jget(L[r], "kr_mp_n") == 1)]);
listput(B, ["H10 p primitive root mod n", vector(N, r, jget(L[r], "emb_ratio") == 1)]);
listput(B, ["H10 (n-1)/k <= 2", vector(N, r, jget(L[r], "emb_ratio") <= 2)]);
F11 = ["kr_m1_n", "kr_2_n", "kr_3_n", "kr_t_n", "kr_nt_n", "kr_n_p", "kr_t_p", "kr_D_p", "kr_nt_p"];
for (k = 1, #F11, listput(B, [Str("H11 ", F11[k], " = 1"), vector(N, r, jget(L[r], F11[k]) == 1)]));
B = Vec(B);

print("binary feature                       P true/10  S true/11  undecided  perfect");
res = vector(#B);
for (k = 1, #B, v = B[k][2]; u = #select(x -> x < 0, v); pf = perfect(v, lab); printf("%-36s %6d     %6d     %6d      %s\n", B[k][1], cnt(v, lab, 1), cnt(v, lab, 0), u, if (pf, "YES", "no")); res[k] = Str("{\"feature\":\"", B[k][1], "\",\"P_true\":", cnt(v, lab, 1), ",\"S_true\":", cnt(v, lab, 0), ",\"undecided\":", u, ",\"perfect\":", pf, "}"));

\\ continuous features (an undecided interval enters at its midpoint)
C = [["H8 lpf(twist)/pbits", vector(N, r, (lt[r][1] + lt[r][2]) / 2 / pb[r])], ["H9 lpf(n-1)/pbits", vector(N, r, (lm[r][1] + lm[r][2]) / 2 / pb[r])], ["H9 lpf(n+1)/pbits", vector(N, r, (lq[r][1] + lq[r][2]) / 2 / pb[r])], ["H7 #{l<1e6: l^2|D}", vector(N, r, #jget(L[r], "D_sq_primes"))]];
print("\ncontinuous feature        median P   median S   U      exact p   Bonferroni   perfect threshold");
cres = vector(#C);
for (k = 1, #C, x = C[k][2]; u = U(x, lab); pv = Up(x, lab); printf("%-24s %9.3f  %9.3f   %5.1f   %.4f    %.4f       %s\n", C[k][1], med(vecextract(x, select(z -> z, lab, 1))), med(vecextract(x, select(z -> !z, lab, 1))), u, pv, min(1, pv * #C), if (splits(x, lab), "YES", "no")); cres[k] = Str("{\"feature\":\"", C[k][1], "\",\"U\":", u * 1., ",\"p_exact\":", pv * 1., ",\"p_bonferroni\":", min(1, pv * #C) * 1., ",\"perfect\":", splits(x, lab), "}"));

\\ positive control: plant P := the 10 candidates with the smallest n mod 2^20
x = vector(N, r, jget(L[r], "n_low20")); o = vecsort(x, , 1); yp = vector(N); for (k = 1, 10, yp[o[k]] = 1);
pcd = splits(x, yp);
print("\npositive control (planted label on n mod 2^20): detected = ", pcd, "; that feature splits the real labels = ", splits(x, lab));

\\ negative control: random relabellings, 10 P of 21
setrand(20261017); T = 100000; hit = 0;
for (it = 1, T, s = numtoperm(N, random(N!)); y = vector(N); for (k = 1, 10, y[s[k]] = 1); if (#select(b -> perfect(b[2], y), B) > 0, hit++));
print("negative control: ", hit, " of ", T, " random labellings split perfectly on some binary feature (family-wise rate ", hit * 1. / T, ")");
write(Str(OUT, "analysis.json"), Str("{\"binary\":[", strjoin(res, ","), "],\"continuous\":[", strjoin(cres, ","), "],\"positive_control_detected\":", pcd, ",\"negative_control\":{\"labellings\":", T, ",\"any_perfect\":", hit, "}}"));
quit;
