\\ PARI/GP side of the external comparison: the workloads, inputs and digests of
\\ examples/ext_compare.rs (see that file for the definitions), one JSON line per case.
\\
\\   gp -q -D parisizemax=4G pari_bench.gp < /dev/null  with OP and TAG set, e.g.
\\   echo 'OP="poly"; TAG="p127"; read("pari_bench.gp")' | gp -q -D parisizemax=4G
\\
\\ OP: "poly" (products and squares), "powx" (x^p mod g), "phi" (Phi_l mod p for all primes
\\ l <= MAXELL, polmodular), "sea" (ellcard on CURVES curves). Set NOSEADATA=1 to hide the
\\ seadata package (PARI then computes its modular polynomials itself). Times are CPU time
\\ (getabstime, ms resolution), so short calls are repeated ITERS times inside each timed block.

PR = Map(["p61", 2^61-1; "p127", 2^127-1; "p255", 2^255-19; "p511", 5326738796327623094747867617954605554069371494832722337612446642054009560026576537626892113026381253624626941643949444792662881241621373288942880288065659]);
p = mapget(PR, TAG);
if (type(NOSEADATA) == "t_INT" && NOSEADATA, default(datadir, "/nonexistent"));
\\ coefficient i (from 0) is G^(i + s) mod p (full-size residues), G the 64-bit golden-ratio
\\ constant of examples/ext_compare.rs; Polrev puts v[1] at degree 0
seqpol(n, s) = { my(g = Mod(11400714819323198485, p), c = g^s); Polrev(vector(n, i, my(v = c); c *= g; v), 'x) };
dig(P) = lift(subst(P, 'x, Mod(123456789, p)));
\\ enough calls per timed block that the 1 ms timer resolution stays under ~2 % of the block
itersfor(n) = if (n <= 64, 2000, if (n <= 256, 200, if (n <= 2048, 40, if (n <= 4096, 20, 10))));
\\ best and median, in microseconds per call, of 5 timed blocks of `it` calls
timeit(f, it) = {
  my(v = vector(5), r);
  for (k = 1, 5, my(t0 = getabstime()); for (i = 1, it, r = f()); v[k] = (getabstime() - t0) * 1000.0 / it);
  v = vecsort(v);
  [v[1], v[3], r];
}
jline(op, n, it, T, d) = printf("{\"tool\":\"pari-%s\",\"op\":\"%s\",\"p\":\"%s\",\"n\":%d,\"reps\":%d,\"best_us\":%.2f,\"median_us\":%.2f,\"digest\":\"%s\"}\n", PARIVERSION, op, TAG, n, 5 * it, T[1], T[2], Str(d));
PARIVERSION = Str(version()[1], ".", version()[2], ".", version()[3]);

{
if (OP == "poly",
  foreach ([16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384], n,
    if (TAG == "p511" && n > 8192, next);
    my(a = seqpol(n, 16), b = seqpol(n, 1000003), it = itersfor(n), T);
    T = timeit(() -> a * b, it); jline("mul", n, it, T, dig(T[3]));
    T = timeit(() -> sqr(a), it); jline("sqr", n, it, T, dig(T[3]))));
}

{
if (OP == "powx",
  foreach ([32, 64, 90, 128, 256, 512, 1024], d,
    my(g = seqpol(d, 2000003) + 'x^d, it = if (d <= 128, 3, 1), T);
    T = timeit(() -> lift(Mod(Mod(1, p) * 'x, g)^p), it); jline("powx", d, it, T, dig(T[3]))));
}

{
if (OP == "phi",
  \\ phi_spec: Phi_l(j0, Y) mod p directly (what SEA needs, one j); phi_full: Phi_l over Z
  \\ (polmodular(l)), reduced mod p. Digest: Phi_l(j0, r) mod p for both.
  my(ts = 0, tf = 0, j0 = Mod(987654321, p), r = Mod(123456789, p));
  forprime (l = 3, MAXELL,
    my(t0 = getabstime(), S = polmodular(l, 0, j0), F);
    ts += getabstime() - t0;
    t0 = getabstime(); F = polmodular(l) * Mod(1, p); tf += getabstime() - t0;
    printf("{\"tool\":\"pari-%s\",\"op\":\"phi_digest\",\"p\":\"%s\",\"ell\":%d,\"digest\":\"%s\",\"digest_full\":\"%s\"}\n", PARIVERSION, TAG, l, Str(lift(subst(S, 'y, r))), Str(lift(subst(subst(F, 'x, j0), 'y, r)))));
  printf("{\"tool\":\"pari-%s\",\"op\":\"phi_spec_all\",\"p\":\"%s\",\"max_ell\":%d,\"total_us\":%d}\n", PARIVERSION, TAG, MAXELL, ts * 1000);
  printf("{\"tool\":\"pari-%s\",\"op\":\"phi_full_all\",\"p\":\"%s\",\"max_ell\":%d,\"total_us\":%d}\n", PARIVERSION, TAG, MAXELL, tf * 1000));
}

{
if (OP == "sea",
  my(src = if (type(NOSEADATA) == "t_INT" && NOSEADATA, "nosea", "seadata"));
  for (k = 0, CURVES - 1,
    my(E = ellinit([k + 3, 7 * k + 11], p), t0, N);
    if (E == [], next);
    t0 = getabstime(); N = ellcard(E);
    printf("{\"tool\":\"pari-%s-%s\",\"op\":\"sea\",\"p\":\"%s\",\"k\":%d,\"us\":%d,\"order\":\"%s\"}\n", PARIVERSION, src, TAG, k, (getabstime() - t0) * 1000, Str(N))));
}
