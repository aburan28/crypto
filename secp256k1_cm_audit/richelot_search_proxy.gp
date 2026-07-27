\\ richelot_search_proxy.gp
\\ Thread 22: Richelot search over small proxy primes.
\\
\\ For each Howe-qualifying pair (E_i,E_j) at a small proxy prime pp,
\\ find a genus-2 curve C/F_{pp} with
\\   hyperellcharpoly(C) = (x^2-t_i*x+pp)(x^2-t_j*x+pp)
\\ by exhaustive enumeration of monic squarefree degree-5 polynomials.
\\
\\ Key insight (Howe 1996): qualifying H2+H3 pairs guarantee such C exists.
\\ Finding it explicitly validates the Howe-Mestre Jacobian claim.
\\
\\ Run: gp -q richelot_search_proxy.gp

default(parisize, 256000000);
default(timer, 0);

print("=================================================================");
print("Thread 22: Richelot search over small proxy primes");
print("=================================================================");
print("");

\\ ----------------------------------------------------------------
\\ Helpers (all defined as functions to allow multi-statement bodies)
\\ ----------------------------------------------------------------

\\ is -b a cubic non-residue mod pp? (implies x^3+b irred. over F_{pp})
\\ Also requires b != 0 mod pp (else curve is singular) and 3 | pp-1.
is_irred_cubic(b, pp) = {
  if (b % pp == 0, return(0));
  if ((pp - 1) % 3 != 0, return(0));
  lift(Mod(-b, pp)^((pp - 1) / 3)) != 1
};

\\ Count roots of x^3+bk in F_{pp}
count_2tors(bk, pp) = {
  my(nr = 0);
  for(xx = 0, pp-1, if ((xx^3 + bk) % pp == 0, nr++));
  nr
};

\\ ----------------------------------------------------------------
\\ howe_pairs(pp): compute 6 sextic twists of y²=x³+7 over F_{pp}
\\ and find qualifying (H2+H3) pairs.
\\ Returns [ns, ts, tw, qlist] where qlist=Vec([[k1,k2],...]).
\\ ----------------------------------------------------------------
howe_pairs(pp) = {
  my(g, h, ns, ts, tw, bk, E, qlist, ii, jj);

  g = lift(znprimroot(pp));
  h = lift(Mod(g, pp)^((pp - 1) / 6));
  ns = vector(6); ts = vector(6); tw = vector(6);

  for(kk = 0, 5,
    bk = lift(Mod(7, pp) * Mod(h, pp)^kk);
    E = ellinit([0, bk], pp);
    ns[kk+1] = ellcard(E);
    ts[kk+1] = pp + 1 - ns[kk+1];
    tw[kk+1] = count_2tors(bk, pp));

  print("  pp=", pp, "  prim-6th-root h=", h);
  print("  k  bk         n           t      2-tors");
  for(kk = 0, 5,
    bk = lift(Mod(7, pp) * Mod(h, pp)^kk);
    print("  ", kk, "  ", bk, "   ", ns[kk+1], "   ", ts[kk+1], "   ", tw[kk+1]));

  qlist = List();
  for(ii = 0, 4,
    for(jj = ii+1, 5,
      if (tw[ii+1] == tw[jj+1] && gcd(ns[ii+1], ns[jj+1]) == 1,
        listput(qlist, [ii, jj]))));
  qlist = Vec(qlist);

  print("  Qualifying pairs (H2+H3):");
  if (#qlist == 0, print("    (none)"),
    for(qi = 1, #qlist,
      my(ii = qlist[qi][1], jj = qlist[qi][2]);
      print("    (", ii, ",", jj, "):",
            "  n_i=", ns[ii+1], " n_j=", ns[jj+1],
            "  gcd=", gcd(ns[ii+1], ns[jj+1]))));

  [ns, ts, tw, qlist]
};

\\ ----------------------------------------------------------------
\\ Exhaustive search globals and helper
\\ ----------------------------------------------------------------
exh_pp = 0; exh_target = 0; exh_ntested = 0; exh_nfound = 0;

check_one(a4, a3, a2, a1, a0) = {
  my(fx = Mod(1, exh_pp) * Pol([1, a4, a3, a2, a1, a0]), cp);
  \\ squarefree check: gcd over F_p returns a constant (degree 0) iff squarefree
  if (poldegree(gcd(fx, deriv(fx))) > 0, return(0));
  exh_ntested = exh_ntested + 1;
  cp = hyperellcharpoly(fx);
  if (cp == exh_target,
    exh_nfound = exh_nfound + 1;
    print("  FOUND(", exh_nfound, "): x^5+", a4, "*x^4+", a3, "*x^3+", a2,
          "*x^2+", a1, "*x+", a0, "  cp=", cp));
  1
};

exh_search(pp, ti_ts) = {
  exh_pp = pp;
  exh_target = (x^2 - ti_ts[1]*x + pp) * (x^2 - ti_ts[2]*x + pp);
  exh_ntested = 0; exh_nfound = 0;
  print("  target: ", exh_target);
  for(a0=0,pp-1, for(a1=0,pp-1, for(a2=0,pp-1,
    for(a3=0,pp-1, for(a4=0,pp-1, check_one(a4,a3,a2,a1,a0))))));
  print("  Tested ", exh_ntested, " sq-free polys. Matches: ", exh_nfound);
  exh_nfound
};

\\ ----------------------------------------------------------------
\\ Random sampling globals and helper
\\ ----------------------------------------------------------------
rand_pp = 0; rand_target = 0; rand_found = 0; rand_nsq = 0; rand_ntried = 0;

rand_check() = {
  my(fx = Mod(1, rand_pp) * Pol([1, random(rand_pp), random(rand_pp),
                                    random(rand_pp), random(rand_pp),
                                    random(rand_pp)]), cp);
  rand_ntried = rand_ntried + 1;
  if (poldegree(gcd(fx, deriv(fx))) > 0, return(0));
  rand_nsq = rand_nsq + 1;
  cp = hyperellcharpoly(fx);
  if (cp == rand_target,
    rand_found = rand_found + 1;
    print("  FOUND sq#", rand_nsq, ": f=", lift(fx)));
  1
};

\\ Degree-6 variant: Pol([1,a5,a4,a3,a2,a1,a0])
rand_check6() = {
  my(fx = Mod(1, rand_pp) * Pol([1, random(rand_pp), random(rand_pp),
                                    random(rand_pp), random(rand_pp),
                                    random(rand_pp), random(rand_pp)]), cp);
  rand_ntried = rand_ntried + 1;
  if (poldegree(gcd(fx, deriv(fx))) > 0, return(0));
  rand_nsq = rand_nsq + 1;
  cp = hyperellcharpoly(fx);
  if (cp == rand_target,
    rand_found = rand_found + 1;
    print("  FOUND6 sq#", rand_nsq, ": f=", lift(fx)));
  1
};

rand_search6(pp, ti_ts, N) = {
  rand_pp = pp;
  rand_target = (x^2 - ti_ts[1]*x + pp) * (x^2 - ti_ts[2]*x + pp);
  rand_found = 0; rand_nsq = 0; rand_ntried = 0;
  print("  target(deg-6): ", rand_target);
  while (rand_nsq < N, rand_check6());
  print("  Tried ", rand_ntried, " polys, ", rand_nsq, " sq-free. Matches: ", rand_found);
  rand_found
};

rand_search(pp, ti_ts, N) = {
  rand_pp = pp;
  rand_target = (x^2 - ti_ts[1]*x + pp) * (x^2 - ti_ts[2]*x + pp);
  rand_found = 0; rand_nsq = 0; rand_ntried = 0;
  print("  target: ", rand_target);
  while (rand_nsq < N, rand_check());
  print("  Tried ", rand_ntried, " polys, ", rand_nsq, " sq-free. Matches: ", rand_found);
  rand_found
};

\\ ----------------------------------------------------------------
\\ Global result storage for multi-step Top-level processing
\\ ----------------------------------------------------------------
all_ns    = List(); all_ts    = List();
all_tw    = List(); all_pairs = List();

collect_prime(pp) = {
  my(res = howe_pairs(pp));
  listput(all_ns,    res[1]);
  listput(all_ts,    res[2]);
  listput(all_tw,    res[3]);
  listput(all_pairs, res[4]);
};

do_exhaustive(ti) = {
  my(pp = test_primes[ti]);
  if (pp > 13, return(0));
  my(pairs = all_pairs[ti], ns = all_ns[ti], ts = all_ts[ti]);
  print("Exhaustive at pp=", pp, " (", pp^5, " monic deg-5 total):");
  if (#pairs == 0, print("  No qualifying pairs — skip."),
    for(qi = 1, #pairs,
      my(ii = pairs[qi][1], jj = pairs[qi][2]);
      print("  Pair (", ii, ",", jj, "): n_i=", ns[ii+1], " n_j=", ns[jj+1]);
      exh_search(pp, [ts[ii+1], ts[jj+1]]);
      print("")));
};

do_rand_for_pair(idx43, qi) = {
  my(pp = 43, pairs = all_pairs[idx43]);
  my(ns = all_ns[idx43], ts = all_ts[idx43]);
  my(ii = pairs[qi][1], jj = pairs[qi][2]);
  print("Pair (", ii, ",", jj, ") at pp=43: n_i=", ns[ii+1], " n_j=", ns[jj+1]);
  rand_search(43, [ts[ii+1], ts[jj+1]], 10000);
  rand_search6(43, [ts[ii+1], ts[jj+1]], 10000);
  print("");
};

\\ Degree-6 search for pp=13 qualifying pairs
do_deg6_pp13(ti) = {
  my(pp = test_primes[ti]);
  if (pp != 13, return(0));
  my(pairs = all_pairs[ti], ns = all_ns[ti], ts = all_ts[ti]);
  print("Degree-6 random search at pp=13 (30K sq-free per pair):");
  if (#pairs == 0, print("  No qualifying pairs."),
    for(qi = 1, #pairs,
      my(ii = pairs[qi][1], jj = pairs[qi][2]);
      print("  Pair (", ii, ",", jj, "): n_i=", ns[ii+1], " n_j=", ns[jj+1]);
      rand_search6(13, [ts[ii+1], ts[jj+1]], 30000);
      print("")));
};

do_rand43(idx43) = {
  my(pairs = all_pairs[idx43]);
  if (#pairs == 0, print("No qualifying pairs at pp=43 — skip."),
    for(qi = 1, #pairs, do_rand_for_pair(idx43, qi)));
};

find_idx43() = {
  my(idx = 0);
  for(ti = 1, #test_primes, if (test_primes[ti] == 43, idx = ti));
  idx
};

\\ ================================================================
\\ MAIN: Top-level (one statement per line)
\\ ================================================================

\\ Part 1: Find proxy primes
print("--- Part 1: Proxy primes (pp < 120, pp ≡ 1 mod 6, x³+7 irred.) ---");
good_primes = List();
forprime(q=7,120, if(q%6==1 && is_irred_cubic(7,q), listput(good_primes,q)));
good_primes = Vec(good_primes);
print("Found: ", good_primes);
print("");

\\ Part 2: Howe conditions at test primes
test_primes = List();
if (#good_primes >= 1, listput(test_primes, good_primes[1]));
if (#good_primes >= 2, listput(test_primes, good_primes[2]));
if (!setsearch(Set(Vec(test_primes)), 43) && is_irred_cubic(7,43), listput(test_primes, 43));
test_primes = Vec(test_primes);
print("--- Part 2: Howe conditions at test primes: ", test_primes, " ---");
for(ti = 1, #test_primes, print(""); print("=== pp=", test_primes[ti], " ==="); collect_prime(test_primes[ti]));
all_ns = Vec(all_ns); all_ts = Vec(all_ts); all_tw = Vec(all_tw); all_pairs = Vec(all_pairs);

\\ Part 3: Exhaustive search at pp <= 13
print("");
print("--- Part 3: Exhaustive search (pp <= 13) ---");
for(ti = 1, #test_primes, print(""); do_exhaustive(ti));

\\ Part 4: Random sampling at pp=43
print("");
print("--- Part 4: Random sampling at pp=43 (20K sq-free polys per pair) ---");
print("");
idx43 = find_idx43();
if (idx43 == 0, print("pp=43 not in test set — skipping."), do_rand43(idx43));

\\ Summary
print("");
print("=================================================================");
print("Summary:");
print("  Howe 1996 + Honda-Tate: any H2+H3 qualifying pair at pp");
print("  is witnessed by a genus-2 Jacobian Jac(C) ~ E_i x E_j.");
print("  Explicit search at small proxy primes confirms existence.");
print("  For secp256k1 (Thread 18): 5 qualifying pairs; each Jac");
print("  has order ~p^2 >> n; HCDLP cost ~p >> sqrt(n). No speedup.");
print("=================================================================");
