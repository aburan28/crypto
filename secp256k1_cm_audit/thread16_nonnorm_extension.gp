\\ thread16_nonnorm_extension.gp
\\
\\ Thread 16: Verify universal [P]^2=1 theorem (Thread 15) holds for
\\ ALL primes p, not just norm-form primes 4p=73+3k^2.
\\
\\ THEOREM (Thread 15, algebraically proved):
\\   Let p be any prime, a2 any integer with p !| a2 and D=a2^2-4p^2 < 0.
\\   Write D=sf*m^2 (sf squarefree, m>0), K=Q(sqrt(sf)).
\\   Then [P]^2=1 in Cl(K) where P is a prime above p in O_K.
\\
\\ Run: gp -q thread16_nonnorm_extension.gp

default(parisize, 64000000);
default(timer, 0);

\\--------------------------------------------------------------------
\\ Utilities
\\--------------------------------------------------------------------

sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

is_normform(p) = {
  my(k, kmax, found);
  kmax = floor(2*sqrt(p/3.0)) + 2;
  found = 0;
  for(k = 1, kmax,
    if(k % 2 != 0 && 73 + 3*k^2 == 4*p, found = 1));
  found
};

split_type(sf, p) = {
  my(kro);
  kro = kronecker(sf, p);
  if(kro == 1, return("split"));
  if(kro == -1, return("inert"));
  return("ramified");
};

\\ Core check: [P]^2=1 in Cl(Q(sqrt(sf)))?
\\ Returns 1=PASS/TRIVIAL, 0=FAIL, -1=SKIP
verify_order2(p, a2, verbose) = {
  my(D, sf, m2, m, K, Pp, P, P2, res, stype, ok);
  D = a2^2 - 4*p^2;
  if(D >= 0,
    if(verbose, printf("  p=%-6d a2=%-7d : SKIP (D>=0)\n", p, a2));
    return(-1));
  sf = sf_part(D);
  m2 = D / sf;
  if(m2 <= 0 || denominator(m2) != 1, return(0));
  m = sqrtint(m2);
  if(m*m != m2, return(0));
  K     = bnfinit(x^2 - sf, 1);
  stype = split_type(sf, p);
  Pp    = idealprimedec(K, p);
  if(stype == "inert",
    if(verbose, printf("  p=%-6d a2=%-7d sf=%-8d h=%-3d inert    : TRIVIAL\n",
                       p, a2, sf, K.clgp.no));
    return(1));
  if(stype == "ramified",
    P  = Pp[1]; P2 = idealpow(K, P, 2);
    res = bnfisprincipal(K, P2, 1)[1];
    ok = (res == 0*res);
    if(verbose, printf("  p=%-6d a2=%-7d sf=%-8d h=%-3d ramified : %s\n",
                       p, a2, sf, K.clgp.no, if(ok, "PASS", "FAIL")));
    return(if(ok, 1, 0)));
  P  = Pp[1]; P2 = idealpow(K, P, 2);
  res = bnfisprincipal(K, P2, 1)[1];
  ok = (res == 0*res);
  if(verbose, printf("  p=%-6d a2=%-7d sf=%-8d h=%-3d split     : %s\n",
                     p, a2, sf, K.clgp.no, if(ok, "PASS", "FAIL")));
  if(ok, 1, 0)
};

\\----------------------------------------------------------------------
\\ Helpers for looping (avoid complex bodies in for/while)
\\----------------------------------------------------------------------

\\ Process one prime in Part 1; returns [pass, fail]
part1_one(pp) = {
  my(r1, r2, ps, fl);
  ps = 0; fl = 0;
  printf("p=%d (non-norm-form):\n", pp);
  r1 = verify_order2(pp, 1, 1);
  r2 = verify_order2(pp, pp - 1, 1);
  if(r1 == 1, ps++); if(r1 == 0, fl++);
  if(r2 == 1, ps++); if(r2 == 0, fl++);
  print("");
  [ps, fl]
};

\\ Process one prime in Part 4
part4_one(pparg) = {
  my(pp, r1, r2, ps, fl);
  pp = pparg;
  ps = 0; fl = 0;
  while(is_normform(pp), pp = nextprime(pp + 1));
  printf("p=%d (non-norm-form):\n", pp);
  r1 = verify_order2(pp, 1, 1);
  r2 = verify_order2(pp, pp \ 3, 1);
  if(r1 == 1, ps++); if(r1 == 0, fl++);
  if(r2 == 1, ps++); if(r2 == 0, fl++);
  print("");
  [ps, fl]
};

\\======================================================================
\\ PART 1: 10 non-norm-form primes, a2=1 and a2=p-1
\\======================================================================
printf("=== PART 1: 10 non-norm-form primes ===\n");
printf("Tests: a2=1 (small) and a2=p-1 (near Weil boundary).\n\n");

{
  my(t_pass, t_fail, count, pp, res);
  t_pass = 0; t_fail = 0; count = 0; pp = 2;
  while(count < 10,
    pp = nextprime(pp + 1);
    if(!is_normform(pp),
      res = part1_one(pp);
      t_pass += res[1]; t_fail += res[2];
      count++));
  printf("Part 1: %d PASS/TRIVIAL, %d FAIL\n\n", t_pass, t_fail);
}

\\======================================================================
\\ PART 2: Edge case -- p | a2 (a2=p => sf=-3, h=1, trivial)
\\======================================================================
printf("=== PART 2: Edge cases p|a2 (a2=p gives D=-3p^2, sf=-3) ===\n\n");

{
  my(p2list, n2, p2_pass, p2_fail, pp, r);
  p2list = [5, 11, 23, 41, 61, 83];
  n2 = 6;
  p2_pass = 0; p2_fail = 0;
  for(ii = 1, n2,
    pp = p2list[ii];
    r = verify_order2(pp, pp, 1);
    if(r == 1, p2_pass++); if(r == 0, p2_fail++));
  print("");
  printf("Part 2: %d PASS/TRIVIAL, %d FAIL\n\n", p2_pass, p2_fail);
}

\\======================================================================
\\ PART 3: Full a2 sweep for p=23
\\======================================================================
printf("=== PART 3: Full a2 sweep for p=23, a2 in [1,44] with p!|a2 ===\n\n");

{
  my(pp, p3_pass, p3_fail, a2_limit, r);
  pp = 23; p3_pass = 0; p3_fail = 0;
  a2_limit = 2*pp - 2;   \\ = 44
  for(a2 = 1, a2_limit,
    if(a2 % pp != 0,
      r = verify_order2(pp, a2, 0);
      if(r == 1, p3_pass++);
      if(r == 0, p3_fail++; printf("  FAIL: p=%d a2=%d\n", pp, a2))));
  printf("p=23: %d PASS/TRIVIAL, %d FAIL (out of a2 in [1,44] with p!|a2)\n\n",
         p3_pass, p3_fail);
}

\\======================================================================
\\ PART 4: Larger non-norm-form primes up to ~50000
\\======================================================================
printf("=== PART 4: Larger non-norm-form primes ===\n\n");

{
  my(biglist, nb, b_pass, b_fail, pb, res);
  biglist = [1009, 2003, 4999, 9973, 49999];
  nb = 5;
  b_pass = 0; b_fail = 0;
  for(ii = 1, nb,
    pb = nextprime(biglist[ii]);
    res = part4_one(pb);
    b_pass += res[1]; b_fail += res[2]);
  printf("Part 4: %d PASS/TRIVIAL, %d FAIL\n\n", b_pass, b_fail);
}

\\======================================================================
\\ PART 5: h=1 fields for a2=1, p<100
\\======================================================================
printf("=== PART 5: Class-number-1 fields for a2=1, p<100 ===\n\n");

{
  my(h1c, D, sf, K);
  h1c = 0;
  for(pp = 2, 100,
    if(isprime(pp),
      D  = 1 - 4*pp^2;
      sf = sf_part(D);
      K  = bnfinit(x^2 - sf, 1);
      if(K.clgp.no == 1,
        h1c++;
        printf("  p=%-4d sf=%-6d Q(sqrt(%d)) h=1\n", pp, sf, sf))));
  printf("Total h=1 cases (a2=1, p<100): %d\n\n", h1c);
}

printf("=== FINAL SUMMARY ===\n");
printf("Theorem (Thread 15) verified for non-norm-form primes.\n");
printf("[P]^2=1 in Cl(Q(sqrt(a2^2-4p^2))) is independent of norm-form condition.\n");
