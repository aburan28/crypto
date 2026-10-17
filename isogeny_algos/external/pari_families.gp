\\ PARI/GP side of the per-family comparison (examples/ext_families.rs defines the inputs;
\\ `ext_families gen` writes them as KCASES / CCASES): one JSON line per case.
\\
\\   echo 'OP="codomain"; read("cases.gp"); read("external/pari_families.gp")' | gp -q
\\
\\ OP "codomain": ellisogeny(E, P, 1) (Vélu, codomain only); "maps": ellisogeny(E, P) (Vélu
\\ with the rational maps, l <= 1009); "chain": the l^n-isogeny as n l-isogenies: the kernel of
\\ step i is [l^(n-1-i)] P_i, and P_(i+1) = ellisogenyapply(phi_i, P_i).
\\ Times are CPU time (getabstime, ms resolution); short calls run ITERS times per block.

PARIVERSION = Str(version()[1], ".", version()[2], ".", version()[3]);
\\ best and median, in microseconds per call, of `blocks` timed blocks of `it` calls
timeit(f, blocks, it) = {
  my(v = vector(blocks), r);
  for (k = 1, blocks, my(t0 = getabstime()); for (i = 1, it, r = f()); v[k] = (getabstime() - t0) * 1000.0 / it);
  v = vecsort(v);
  [v[1], v[(blocks + 1) \ 2], r];
}
jline(op, c, reps, T, out) = printf("{\"tool\":\"pari-%s\",\"op\":\"%s\",\"bits\":%d,\"ell\":%d,\"n\":%d,\"reps\":%d,\"best_us\":%.2f,\"median_us\":%.2f,%s}\n", PARIVERSION, op, c[1], c[2], c[3], reps, T[1], T[2], out);
cod(E2) = Str("\"a4\":\"", lift(E2[4]), "\",\"a6\":\"", lift(E2[5]), "\"");
\\ blocks and calls per block: enough calls that 1 ms resolution stays small
plan(l, bits) = if (l <= 101, [5, 400], if (l <= 1009 && bits <= 255, [5, 40], [3, 4]));

{
if (OP == "codomain" || OP == "maps",
  foreach (KCASES, c,
    my(p = c[4], E = ellinit([1, 0], Mod(1, p)), P = [Mod(c[5], p), Mod(c[6], p)], pl = plan(c[2], c[1]), T);
    if (OP == "codomain",
      T = timeit(() -> ellisogeny(E, P, 1), pl[1], pl[2]);
      jline("codomain_velu", c, pl[1] * pl[2], T, cod(T[3])));
    \\ the full maps cost O(l^2): far fewer calls
    if (OP == "maps" && c[2] <= 1009,
      pl = if (c[2] <= 101, [5, 10], [3, 1]);
      T = timeit(() -> ellisogeny(E, P), pl[1], pl[2]);
      my(m = T[3][2], r = Mod(123456789, p));
      jline("maps_velu", c, pl[1] * pl[2], T, Str(cod(T[3][1]), ",\"x_map_at_r\":\"", lift(subst(m[1], 'x, r) / subst(m[3], 'x, r)^2), "\"")))));
}

{
if (OP == "chain",
  foreach (CCASES, c,
    my(p = c[4], l = c[2], n = c[3], E0 = ellinit([1, 0], Mod(1, p)), P0 = [Mod(c[5], p), Mod(c[6], p)], it = if (c[1] <= 127, 20, if (c[1] <= 255, 4, 1)), T);
    T = timeit(() -> my(E = E0, P = P0, r); for (i = 0, n - 1, my(K = ellmul(E, P, l^(n - 1 - i)), s = ellisogeny(E, K)); if (i < n - 1, P = ellisogenyapply(s[2], P)); E = ellinit(s[1])); E, 5, it);
    jline("chain_naive", c, 5 * it, T, cod(T[3]))));
}
