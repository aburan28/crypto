\\ Class-size predictors for ISO-1 trace rows (2026-10-16).
\\ Input: env ISO1_P and ISO1_TRACES (file, one positive trace per line).
\\ Output CSV: t,D_K,f_pi,h_K,H_quarter,H_sixteenth where
\\   H_quarter   = Hurwitz class number of |t^2-4Q|/4  (curves with full 2-torsion in the class),
\\   H_sixteenth = Hurwitz class number of |t^2-4Q|/16 (curves with Z[(pi-1)/4] in End),
\\ both with the usual 2/|O^*| weights, from H(f'^2 D_K) = h_K/(w_K/2) * sum_{d|f'} d prod_{l|d}(1-(D_K/l)/l).
\\ h_K is quadclassunit(D_K).no (GRH-conditional for large |D_K|).
p = eval(getenv("ISO1_P")); Q = p^6;
hur(DK, hK, g) = {
  my(wk = if(DK == -3, 3, DK == -4, 2, 1), s = 0);
  fordiv(g, d, s += d * prod(i = 1, #factor(d)[,1], my(l = factor(d)[i,1]); 1 - kronecker(DK, l)/l));
  hK / wk * s
};
ts = readvec(getenv("ISO1_TRACES"));
print("t,D_K,f_pi,h_K,H_quarter,H_sixteenth");
{
for (i = 1, #ts,
  my(t = ts[i], D = t^2 - 4*Q, DK = coredisc(D), f = sqrtint(D / DK), hK = quadclassunit(DK).no);
  if (f^2 * DK != D, error("bad conductor at t=", t));
  my(Hq = if(f % 2 == 0, hur(DK, hK, f/2), -1), Hs = if(f % 4 == 0, hur(DK, hK, f/4), -1));
  printf("%d,%d,%d,%d,%.6f,%.6f\n", t, DK, f, hK, Hq, Hs));
}
quit
