\\ Curve shortlist, 2026-10-16: 46-59-bit ECDLP targets and deep isogeny-volcano curves.
\\ Reproduce:  gp -q -s 2000000000 research/curve_shortlist_20261016/shortlist.gp
\\ Writes research/curve_shortlist_20261016/curves.json (one JSON object per curve record)
\\ and research/curve_shortlist_20261016/searches.jsonl (the negative / scan results).
\\ Everything here is exact (PARI/GP); nothing is sampled.

OUT = "research/curve_shortlist_20261016/curves.json";
SRCH = "research/curve_shortlist_20261016/searches.jsonl";
RECS = List();
default(parisizemax, 4000000000);

\\ ---------- helpers ----------
hexs(n) = { my(d = Vec("0123456789abcdef"), s = "", v = n); if (v == 0, return("0x0"));
  while (v > 0, s = concat(d[v % 16 + 1], s); v \= 16); concat("0x", s); }
q(x) = Str("\"", x, "\"");
\\ FF element <-> integer (bit i = coefficient of w^i)
ff2int(e) = subst(lift(e.pol), 'w, 2);
int2ff(n, g) = subst(Pol(binary(n), 'z), 'z, g);
sha_mod(s, bits) = eval(concat("0x", externstr(Str("printf %s '", s, "' | sha256sum | cut -c1-64"))[1])) % 2^bits;
\\ roots with multiplicity of a univariate polynomial over F_p or F_2^m
rts(f) = { my(F = factor(f), R = List());
  for (i = 1, #F~, if (poldegree(F[i,1]) == 1,
    my(r = -polcoef(F[i,1], 0) / polcoef(F[i,1], 1)); for (k = 1, F[i,2], listput(R, r))));
  Vec(R); }
nbr(l, j) = rts(subst(subst(polmodular(l), 'x, j), 'y, 'Y));
nbstep(l, prev, cur) = { my(N = nbr(l, cur), used = 0);
  for (i = 1, #N, if (N[i] == prev && !used, used = 1, return(N[i])));
  if (#N == 1 && N[1] != prev, return(N[1])); 0; }
walkh(l, j, first, cap) = { my(prev = j, cur = first, h = 1);
  while (h < cap, my(N = nbr(l, cur)); if (#N <= 1, return(h));
    my(nx = nbstep(l, prev, cur)); if (nx == 0, return(h)); prev = cur; cur = nx; h++); cap; }
\\ height above the floor of the l-volcano (0 = floor): shortest non-backtracking walk to a degree-1 vertex
height(l, j, cap) = { my(N = nbr(l, j), best = cap, seen = List());
  if (#N <= 1, return(0));
  for (i = 1, #N, if (setsearch(Set(seen), N[i]), next); listput(seen, N[i]);
    best = min(best, walkh(l, j, N[i], cap)); if (#seen >= 3, break));
  best; }
\\ walk from j to the floor; returns [floor j, path]
descend(l, j, cap) = { my(h = height(l, j, cap), path = List([j]), nx = 0);
  if (h == 0, return([j, Vec(path)]));
  my(N = nbr(l, j)); for (i = 1, #N, if (height(l, N[i], cap) == h - 1, nx = N[i]; break));
  my(prev = j, cur = nx); listput(path, cur);
  for (s = 2, h, my(t = nbstep(l, prev, cur)); prev = cur; cur = t; listput(path, cur));
  [cur, Vec(path)]; }
frob(qq, t) = { my(D = t^2 - 4*qq, DK = coredisc(D), v = sqrtint(D / DK)); [D, DK, v]; }
fstr(f) = { my(s = "{"); for (i = 1, #f~, s = concat(s, Str(if (i > 1, ",", ""), q(f[i,1]), ":", f[i,2]))); concat(s, "}"); }
twopow(N) = 2^valuation(N, 2);
\\ deterministic generator: cofactor times the point with the least x (integer encoding), first ordinate
gen_prime(E, p, n) = { my(h = ellcard(E) / n);
  for (x = 0, 10^6, my(Y = ellordinate(E, Mod(x, p))); if (#Y,
    my(G = ellmul(E, [Mod(x, p), Y[1]], h)); if (G != [0] && ellmul(E, G, n) == [0], return(G)))); }
gen_bin(E, g, n) = { my(h = ellcard(E) / n);
  for (x = 1, 10^6, my(Y = ellordinate(E, int2ff(x, g))); if (#Y,
    my(G = ellmul(E, [int2ff(x, g), Y[1]], h)); if (G != [0] && ellmul(E, G, n) == [0], return(G)))); }
emb(qq, n) = znorder(Mod(qq, n));

rec_bin(name, role, m, low, g, E, n, extra) = { my(N = ellcard(E), G = gen_bin(E, g, n), t = 2^m + 1 - N, fr = frob(2^m, t));
  listput(RECS, Str("{\"name\":", q(name), ",\"role\":", q(role),
    ",\"field\":{\"kind\":\"binary\",\"degree\":", m, ",\"polynomial_low_terms\":", low, "}",
    ",\"form\":\"y^2+xy=x^3+a*x^2+b\",\"a\":", q(hexs(ff2int(E.a2))), ",\"b\":", q(hexs(ff2int(E.a6))),
    ",\"order\":", q(N), ",\"trace\":", q(t), ",\"subgroup_order\":", q(n), ",\"cofactor\":", N / n,
    ",\"order_factorization\":", fstr(factor(N)),
    ",\"generator\":[", q(hexs(ff2int(G[1]))), ",", q(hexs(ff2int(G[2]))), "]",
    ",\"frobenius_disc\":", q(fr[1]), ",\"D_K\":", q(fr[2]), ",\"conductor_gap_v\":", q(fr[3]),
    ",\"v_factorization\":", fstr(factor(fr[3])), ",\"embedding_degree\":", q(emb(2^m, n)), extra, "}")); }
rec_prime(name, role, p, E, n, extra) = { my(N = ellcard(E), G = gen_prime(E, p, n), t = p + 1 - N, fr = frob(p, t));
  listput(RECS, Str("{\"name\":", q(name), ",\"role\":", q(role), ",\"p\":", q(p),
    ",\"form\":\"y^2=x^3+a*x+b\",\"a\":", q(lift(E.a4)), ",\"b\":", q(lift(E.a6)),
    ",\"order\":", q(N), ",\"trace\":", q(t), ",\"subgroup_order\":", q(n), ",\"cofactor\":", N / n,
    ",\"order_factorization\":", fstr(factor(N)),
    ",\"generator\":[", q(lift(G[1])), ",", q(lift(G[2])), "]",
    ",\"frobenius_disc\":", q(fr[1]), ",\"D_K\":", q(fr[2]), ",\"conductor_gap_v\":", q(fr[3]),
    ",\"v_factorization\":", fstr(factor(fr[3])), ",\"embedding_degree\":", q(emb(p, n)), extra, "}")); }
heights_json(l_list, j, label) = { my(s = Str(",\"", label, "\":{"));
  for (i = 1, #l_list, s = concat(s, Str(if (i > 1, ",", ""), q(l_list[i][1]), ":", height(l_list[i][1], j, l_list[i][2] + 3))));
  concat(s, "}"); }

system(Str("rm -f ", OUT, " ", SRCH));

\\ ---------- 1. Koblitz orders, m = 44..89 (negative search for prime m in 46..70) ----------
kobN(m, a) = { my(mu = if (a == 1, 1, -1), t0 = 2, t1 = mu); for (k = 2, m, [t0, t1] = [t1, mu*t1 - 2*t0]); [2^m + 1 - t1, t1]; }
{
for (m = 44, 89, for (a = 0, 1, my(r = kobN(m, a), f = factor(r[1]), n = f[#f~, 1]);
  write(SRCH, Str("{\"search\":\"koblitz\",\"m\":", m, ",\"a\":", a, ",\"m_prime\":", isprime(m),
    ",\"order\":", q(r[1]), ",\"largest_prime_bits\":", #binary(n), ",\"v\":", q(sqrtint((2^(m+2) - r[2]^2) / 7)), "}"))));
}

\\ ---------- 2. The three general targets ----------
\\ Koblitz a = 0 over GF(2^73), modulus x^73+x^4+x^3+x^2+1 (find_irreducible_sparse(73), docs/curves/ICV1.md)
g73 = ffgen(Mod(1, 2) * ('w^73 + 'w^4 + 'w^3 + 'w^2 + 1), 'w);
E = ellinit([1, 0, 0, 0, g73^0]);
rec_bin("shortlist-koblitz-73", "target", 73, "[4,3,2,0]", g73, E, 86020738150056119, ",\"koblitz\":true");

\\ random binary over GF(2^59), modulus x^59+x^7+x^4+x^2+1; b = SHA-256("binary59-<c>") mod 2^59, least c with #E = 2*prime
g59 = ffgen(Mod(1, 2) * ('w^59 + 'w^7 + 'w^4 + 'w^2 + 1), 'w);
{
for (c = 0, 10^4, my(bi = sha_mod(Str("binary59-", c), 59)); if (bi == 0, next);
  E = ellinit([1, 1, 0, 0, int2ff(bi, g59)]); my(N = ellcard(E));
  if (N % 2 == 0 && isprime(N / 2),
    rec_bin("shortlist-binary-59", "target", 59, "[7,4,2,0]", g59, E, N / 2, Str(",\"seed\":\"binary59-", c, "\""));
    break));
}

\\ prime field p = 2^56 - 5, y^2 = x^3 - 3x + b, b = SHA-256("prime56-<c>") mod 2^56 mod p, least c with prime #E
p56 = precprime(2^56);
{
for (c = 0, 10^4, my(b = sha_mod(Str("prime56-", c), 56) % p56);
  E = ellinit([-3, b], p56); my(N = ellcard(E));
  if (isprime(N),
    rec_prime("shortlist-prime-56", "target", p56, E, N, Str(",\"seed\":\"prime56-", c, "\",\"twist_order_factorization\":", fstr(factor(2*p56 + 2 - N))));
    break));
}

\\ ---------- 3. Deep prime-field volcanoes by CM ----------
\\ find p = (t^2 + |DK| v^2)/4 prime in [2^lo, 2^hi), #E = 2^k * prime with 2^k <= 16
findp(DK, v, lo, hi, start) = { my(a = abs(DK), t);
  start = max(start, sqrtint(max(0, 4*2^lo - a*v^2)) + 1) + start % 1000;
  t = if ((start^2 + a*v^2) % 4 == 0, start, start + 1);
  while (1, my(num = t^2 + a*v^2);
    if (num % 4 == 0, my(p = num / 4); if (p >= 2^hi, return(0));
      if (p >= 2^lo && isprime(p), foreach([t, -t], tt, my(N = p + 1 - tt, c = twopow(N));
        if (c <= 16 && isprime(N / c), return([p, tt, N]))))); t++); }
twist_prime(a4, a6, p) = { my(d = 2); while (kronecker(d, p) != -1, d++); [a4*d^2, a6*d^3]; }
curve_with_order(j, p, N) = { my(E = ellinit(ellfromj(j)));
  if (ellcard(E) != N, my(tw = twist_prime(E.a4, E.a6, p)); E = ellinit([tw[1], tw[2]], p));
  if (ellcard(E) != N, error("no twist with order ", N)); E; }
deep_prime(name, DK, v, lo, hi, st) = { my(r = findp(DK, v, lo, hi, st), p = r[1], t = r[2], N = r[3], j0, jj, L, steps = List(), cyc = "");
  if (DK == -3, j0 = Mod(0, p), my(cr = rts(polclass(DK) * Mod(1, p))); j0 = cr[1];
    \\ crater cycles: follow crater-to-crater edges back to the start
    cyc = ",\"crater_cycles\":{";
    my(fl = factor(v)[,1]~, first = 1);
    foreach(fl, l, my(nb = nbr(l, j0), crn = [x | x <- nb, setsearch(Set(cr), x)], len = 0);
      if (#crn, my(prev = j0, cur = crn[1], k = 1);
        while (cur != j0 && k < 50, my(nn = [x | x <- nbr(l, cur), setsearch(Set(cr), x) && x != prev]); prev = cur; cur = nn[1]; k++);
        len = k);
      cyc = concat(cyc, Str(if (first, "", ","), q(l), ":{\"crater_neighbours\":", #crn, ",\"cycle_length\":", len, "}")); first = 0);
    cyc = concat(cyc, "}"));
  L = [[f[1], f[2]] | f <- Vec(factor(v)~)];
  my(hc = heights_json(L, j0, "height_at_crater"));
  jj = j0;
  foreach(L, le, my(d = descend(le[1], jj, le[2] + 3)); listput(steps, [le[1], #d[2] - 1]); jj = d[1]);
  my(E = curve_with_order(jj, p, N), n = N / twopow(N), st_json = "{");
  for (i = 1, #steps, st_json = concat(st_json, Str(if (i > 1, ",", ""), q(steps[i][1]), ":", steps[i][2])));
  st_json = concat(st_json, "}");
  rec_prime(name, "deep-volcano floor", p, E, n,
    Str(",\"construction\":{\"D_K\":", DK, ",\"h_K\":", if (DK < -4, qfbclassno(DK), 1), ",\"crater_j\":", q(lift(j0)),
      ",\"descent_steps\":", st_json, "}", hc, heights_json(L, jj, "height_at_floor"), cyc)); }

deep_prime("shortlist-deep-tall", -3, 2^27, 55, 57, 2^27 + 12345);
deep_prime("shortlist-deep-broad", -3, 2^4*3^2*5*7*11*13*17*19, 56, 59, 2^28 + 777);
deep_prime("shortlist-deep-cycle", -71, 2^6*3^4*5^2*7, 55, 57, 2^27 + 4242);

\\ ---------- 4. Binary CM scan: m | h(D_K), t^2 + |D_K| v^2 = 2^(m+2), t odd ----------
best = 0; bestscore = [-1, 0, 0];
{
foreach([53, 59, 61], m,
  forstep (a = 7, 400000, 8, my(DK = -a);
    if (!isfundamental(DK), next); my(h = qfbclassno(DK)); if (h % m, next);
    my(S = qfbsolve(Qfb(1, 0, a), 2^(m+2), 1));
    foreach(S, s, my(t = abs(s[1]), v = abs(s[2])); if (t % 2 == 0 || v == 0, next);
      foreach([t, -t], tt, my(N = 2^m + 1 - tt, f = factor(N), n = f[#f~, 1], fv = factor(v),
          small = sum(i = 1, #fv~, if (fv[i,1] <= 200, fv[i,2], 0)), nsmall = sum(i = 1, #fv~, fv[i,1] <= 200));
        write(SRCH, Str("{\"search\":\"binary-cm\",\"m\":", m, ",\"D_K\":", DK, ",\"h\":", h, ",\"t\":", q(tt),
          ",\"v\":", q(v), ",\"v_factorization\":", fstr(fv), ",\"largest_prime_bits\":", #binary(n), ",\"small_l_depth\":", small, "}"));
        if (#binary(n) >= 50 && lex([small, nsmall, #binary(n)], bestscore) > 0, bestscore = [small, nsmall, #binary(n)]; best = [m, DK, h, tt, N, n, v])))));
}

\\ ---------- 5. Deep binary: floor of the best scan entry ----------
[m, DK, h, t, N, n, v] = best;
if (m != 61, error("scan winner changed: ", best));
g61 = ffgen(Mod(1, 2) * ('w^61 + 'w^5 + 'w^2 + 'w + 1), 'w);
cr = rts(polclass(DK) * g61^0); j0 = cr[1];
curve_bin(j) = { my(E = ellinit([1, 0, 0, 0, 1/j])); if (ellcard(E) != N, E = ellinit([1, 1, 0, 0, 1/j])); if (ellcard(E) != N, error("twist")); E; }
d7 = descend(7, j0, 6);
E = curve_bin(d7[1]);
{
rec_bin("shortlist-deep-binary", "deep-volcano floor", 61, "[5,2,1,0]", g61, E, n,
  Str(",\"construction\":{\"D_K\":", DK, ",\"h_K\":", h, ",\"crater_j\":", q(hexs(ff2int(j0))), ",\"descent_steps\":{\"7\":", #d7[2] - 1, "}}",
      ",\"height_at_crater\":{\"7\":", height(7, j0, 6), "},\"height_at_floor\":{\"7\":", height(7, d7[1], 6), "}",
      ",\"unwalked\":\"l=31543 volcano (depth 1) not walked: Phi_31543 is out of reach; this curve is on the l=7 floor but at the l=31543 crater, so End(E) = Z + 7^3 O_K, which strictly contains Z[pi]\""));
}
{
  my(body = "");
  for (i = 1, #RECS, body = concat(body, Str(if (i > 1, ",\n", ""), RECS[i])));
  write(OUT, Str("{\"schema\":\"curve-shortlist/v1\",\"generated_by\":\"research/curve_shortlist_20261016/shortlist.gp\",\"curves\":[\n", body, "\n]}"));
}
print("done");
quit;
