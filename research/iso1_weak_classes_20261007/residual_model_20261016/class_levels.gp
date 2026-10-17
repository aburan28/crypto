\\ Hurwitz class numbers H(|t^2-4Q|/4^j), j = 1..10, for ISO-1 trace rows (2026-10-16).
\\ Input: env ISO1_P, ISO1_TRACES (one positive trace per line). -1 marks j with 2^j not dividing f_pi.
\\ H_odd = H for the odd part of f_pi (the 2-adic part removed): h_K/(w_K/2) * sum_{d | f_odd} ...
\\ Same weighted formula and GRH-conditional h_K as class_size.gp.
p = eval(getenv("ISO1_P")); Q = p^6;
hur(DK, hK, g) = {
  my(wk = if(DK == -3, 3, DK == -4, 2, 1), s = 0);
  fordiv(g, d, s += d * prod(i = 1, #factor(d)[,1], my(l = factor(d)[i,1]); 1 - kronecker(DK, l)/l));
  hK / wk * s
};
ts = readvec(getenv("ISO1_TRACES"));
print("t,D_K,f_pi,h_K,H_odd,H1,H2,H3,H4,H5,H6,H7,H8,H9,H10");
{
for (i = 1, #ts,
  my(t = ts[i], D = t^2 - 4*Q, DK = coredisc(D), f = sqrtint(D / DK), hK = quadclassunit(DK).no, s);
  if (f^2 * DK != D, error("bad conductor at t=", t));
  my(fo = f >> valuation(f, 2));
  s = Str(t, ",", DK, ",", f, ",", hK, ",", hur(DK, hK, fo));
  for (j = 1, 10, s = Str(s, ",", if(f % 2^j == 0, hur(DK, hK, f / 2^j), -1)));
  print(s));
}
quit
