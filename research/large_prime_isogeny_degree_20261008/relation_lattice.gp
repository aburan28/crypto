\\ E2a instrument: class group of the Frobenius order and relation norms.
\\ See PROTOCOL.md, section E2a.  Run from a directory holding an empty
\\ `results/`:
\\     gp -q relation_lattice.gp              (registered P-256; 6 h cap outside)
\\     E2_TOY=1 gp -q relation_lattice.gp     (toy discriminant self-test)
\\
\\ Computes Cl(Z[pi]) for D = t^2 - 4p under GRH, the discrete logarithms
\\ of the factor-base primes q <= 200 with (D/q) = +1, the relation
\\ lattice, and for each target ell a Babai-reduced relation vector e with
\\ prod q_i^{e_i} ~ l (the prime over ell), its l1 norm and its charge.
\\ The second solution for lbar is checked to be -e mod L.
\\
\\ Discrete logs are taken with bnfisprincipal on the maximal order, which
\\ is Z[pi] exactly when D is fundamental (settled for P-256 in PROTOCOL.md
\\ step 0).  A non-fundamental D is refused rather than silently computed
\\ in the wrong order.

default(parisizemax, 2^32);

toy = (getenv("E2_TOY") == "1");

\\ Toy class: p = 1000003 and the first odd t >= 1235 with t^2 - 4p fundamental.
{
if (toy,
  p = 1000003; t = 1235; while (!isfundamental(t^2 - 4*p), t += 2); n = p + 1 - t
,
  p = 115792089210356248762697446949407573530086143415290314195533631308867097853951;
  n = 115792089210356248762697446949407573529996955224135760342422259061068512044369
);
}
t = p + 1 - n;
D = t^2 - 4*p;
if (D % 4 == 2 || D % 4 == 3, error("D must be 0 or 1 mod 4"));

\\ Charge sheet for one small-degree step (PROTOCOL.md).
charge_step(q) = 2*q^2 + 4*q*#binary(p) + 40*q^2;

\\ ---- step 1: class group of the order of discriminant D (GRH) ----
T0 = getabstime();
fund = isfundamental(D);
if (!fund, error("non-fundamental D: refused; see PROTOCOL.md E2 step 0"));
K = bnfinit(y^2 - D, 1);
h = K.no; cyc = K.cyc; k = #cyc;
printf("D bits %d, h = %d (log2 %.2f), cyc = %s, %.1f s\n", #binary(abs(D)), h, log(h)/log(2), cyc, (getabstime()-T0)/1000.0);

\\ Coordinates of the prime over q in the cyclic decomposition, as a row.
dl(q) = bnfisprincipal(K, idealprimedec(K, q)[1], 0)~;

\\ ---- step 2: factor base of Elkies primes q <= 200 ----
B = select(q -> kronecker(D, q) == 1, primes([3, 200]));
m = #B;
V = matrix(m, k, i, j, dl(B[i])[j]);   \\ row i = coordinates of q_i
printf("factor base |B| = %d: %s\n", m, B);

\\ ---- step 3: relation lattice L = { e : e*V = 0 mod cyc } ----
\\ Integer kernel of [V~ | diag(cyc)], projected to the e-part, LLL-reduced.
\\ `cols` selects the factor-base columns in use: a target that is itself in
\\ the base is solved with its own column removed, so the relation is never
\\ the trivial one-step vector.
{
lattice(cols) =
  my(Vc = vecextract(V, cols, Str("..")), Mc = matconcat([Vc~, matdiagonal(cyc)]),
     Lc = matkerint(Mc)[1..#cols, ]);
  Lc * qflll(Lc);
}
Lfull = lattice([1..m]);
printf("relation lattice rank %d, min basis l1 norm %d\n", #Lfull, vecmin(vector(#Lfull, j, vecsum(abs(Lfull[, j])))));

\\ Babai rounding of a target against an LLL-reduced basis (square, full rank).
reduce(Lc, e) = e - Lc * round(matsolve(Lc, e));

\\ One Babai-reduced solution of e*V = target mod cyc over the columns `cols`,
\\ returned as a column of length m with zeros outside `cols`.
{
solve_target(cols, target) =
  my(Vc = vecextract(V, cols, Str("..")), sol = matsolvemod(Vc~, cyc~, target~, 1), e, full);
  if (sol == 0, error("no solution: target not in the span of the factor base"));
  e = reduce(lattice(cols), sol[1][1..#cols]);
  full = vectorv(m); for (i = 1, #cols, full[cols[i]] = e[i]);
  full;
}

\\ ---- step 4: targets ----
{
first_elkies_above(x) =
  my(q = nextprime(x));
  while (kronecker(D, q) != 1, q = nextprime(q + 1));
  q;
}
{
T = if (toy,
  select(q -> kronecker(D, q) == 1, primes([3, 60])),
  concat([11, 13, 17, 23, 29, 37, 41, 43, 47, 59],
         vector(6, i, first_elkies_above([10^3, 10^4, 10^5, 10^6, 2^32, 2^64][i]))));
}

out = fileopen("results/e2a_relations.jsonl", "w");
{
for (i = 1, #T,
  my(ell = T[i], cols = select(j -> B[j] != ell, [1..m]), v = dl(ell),
     e = solve_target(cols, v), ebar = solve_target(cols, -v),
     sum_check = reduce(lattice(cols), vecextract(e + ebar, cols)),
     l1 = vecsum(abs(e)),
     charge = sum(j = 1, m, abs(e[j]) * charge_step(B[j])));
  printf("ell = %d: |e|_1 = %d, charge = %.3g, conj consistent: %d\n", ell, l1, charge, sum_check == 0);
  filewrite(out, strprintf("{\"ell\": %d, \"in_base\": %d, \"e\": %s, \"l1\": %d, \"charge\": %.6g, \"conjugate_consistent\": %d}", ell, #cols < m, Vec(e), l1, charge, sum_check == 0));
);
}
fileclose(out);

out = fileopen("results/e2a_classgroup.json", "w");
filewrite(out, strprintf("{\"p_bits\": %d, \"D_bits\": %d, \"fundamental\": %d, \"h\": %d, \"log2_h\": %.4f, \"cyc\": %s, \"factor_base\": %s, \"seconds\": %.1f}", #binary(p), #binary(abs(D)), fund, h, log(h)/log(2), cyc, B, (getabstime()-T0)/1000.0));
fileclose(out);
