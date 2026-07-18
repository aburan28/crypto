\\ thread16_general_order2.gp
\\
\\ Thread 16: Verify the "universal order-2 Frobenius" theorem for
\\ non-norm-form primes.  The theorem (proved algebraically in Thread 15)
\\ states:
\\
\\   THEOREM (General, T15): For any prime p, any integer a2 with
\\   |a2| < 2p, and D = a2^2 - 4p^2 = sf*m^2 (sf squarefree, m > 0)
\\   with D < 0, sf < 0, and p ∤ a2, the prime P above p in
\\   K = Q(sqrt(sf)) satisfies [P]^2 = 1 in Cl(K).
\\
\\ The norm-form condition 4p = 73 + 3k^2 from Thread 15 was NOT used
\\ in the proof — it only constrained which primes appeared in that
\\ experiment.  This script verifies the theorem for:
\\
\\   (A) 10 non-norm-form primes p, each with 5 random valid a2 values,
\\       sampling across different sf (CM fields).
\\   (B) One direct test of the key algebraic invariants:
\\       minpoly(beta), N(beta)=p^2, p∤a2.
\\
\\ Run: gp --stacksize 256000000 -q thread16_general_order2.gp

default(parisize, 256000000);
default(timer, 0);

\\ ==============================================================
\\ Utilities
\\ ==============================================================

sf_part(n) = {
  if(n == 0, return(0));
  sign(n) * core(abs(n))
}

\\ Is p a norm-form prime? i.e., 4p = 73 + 3k^2 for some integer k.
is_norm_form(p) = {
  my(rhs, k2, k);
  rhs = 4*p - 73;
  if(rhs <= 0 || rhs % 3 != 0, return(0));
  k2 = rhs \ 3;
  k = sqrtint(k2);
  (k * k == k2)
}

\\ ==============================================================
\\ Core verification: given prime p and integer a2 with |a2| < 2p,
\\ check the five algebraic conditions from the Thread 15 proof,
\\ then verify [P]^2 = 1 numerically.
\\
\\ Returns: 1 = all PASS (YES), 0 = degenerate/skipped, -1 = FAIL
\\ ==============================================================

verify_order2_full(p, a2, verbose) = {
  my(D, sf, m2, m, K, h, Pp, P, P2, P2_prin,
     beta_elt, norm_beta, ord_P, ord_P2);

  \\ --- Precondition: |a2| < 2p ---
  if(abs(a2) >= 2*p, return(0));

  \\ --- Condition check ---
  D = a2^2 - 4*p^2;

  \\ (A) Need D < 0 (CM case)
  if(D >= 0, return(0));

  sf = sf_part(D);

  \\ (B) Need sf < 0 (imaginary quadratic field)
  if(sf >= 0, return(0));

  m2 = D / sf;
  \\ m2 should be a perfect square
  if(!issquare(m2), return(0));
  m = sqrtint(m2);

  \\ (C) Need N(beta) = p^2:  verified algebraically: N = (a2^2 - m^2*sf)/4 = (a2^2-D)/4 = p^2
  norm_beta = (a2^2 - m2 * sf) / 4;
  if(norm_beta != p^2,
    if(verbose, printf("    FAIL: N(beta)=%d != p^2=%d\n", norm_beta, p^2));
    return(-1));

  \\ (D) Need p ∤ a2
  if(a2 % p == 0, return(0));

  \\ --- Numerical verification in K = Q(sqrt(sf)) ---
  K = bnfinit(x^2 - sf, 1);
  h = K.clgp.no;
  Pp = idealprimedec(K, p);

  \\ p should split in K (since D < 0 and p ∤ disc means p splits or ramifies)
  if(#Pp == 0,
    if(verbose, printf("    SKIP: p=%d inert in Q(sqrt(%d))\n", p, sf));
    return(0));

  P = Pp[1];

  \\ Verify [P]^2 = 1: check P^2 is principal
  P2 = idealpow(K, P, 2);
  P2_prin = bnfisprincipal(K, P2);
  if(P2_prin[1] != 0 * P2_prin[1],
    if(verbose, printf("    FAIL: [P]^2 != 1 for p=%d, sf=%d, h=%d\n", p, sf, h));
    return(-1));

  \\ Sanity: what is ord([P])?
  ord_P2 = bnfisprincipal(K, P2)[1] == 0*bnfisprincipal(K, P2)[1];  \\ 1 if principal
  ord_P  = bnfisprincipal(K, P)[1] == 0*bnfisprincipal(K, P)[1];

  if(verbose,
    printf("    a2=%-9d D=%-14d sf=%-10d m=%-7d h=%-4d [P]=1:%s [P]^2=1:YES\n",
      a2, D, sf, m, h, if(ord_P, "YES(trivial)","no ")));

  1
}

\\ ==============================================================
\\ Generate a2 values for testing: sample across the valid range
\\ [-2p+1, 2p-1] with step chosen to get diverse CM fields.
\\ ==============================================================

sample_a2_values(p, n) = {
  my(vals, a2, step, count);
  vals = [];
  step = max(1, (4*p - 2) \ (3 * n));  \\ coarse step to spread across range
  count = 0;
  a2 = -(2*p - 1);
  while(count < n && a2 < 2*p,
    \\ Skip degenerate cases
    if(a2 != 0 && a2 % p != 0,
      my(D, sf, m2);
      D = a2^2 - 4*p^2;
      if(D < 0,
        sf = sf_part(D);
        m2 = D / sf;
        if(sf < 0 && issquare(m2),
          vals = concat(vals, [a2]);
          count++)));
    a2 += step);
  vals
}

\\ ==============================================================
\\ Main
\\ ==============================================================

print("Thread 16: Universal order-2 Frobenius — extension to non-norm-form primes");
print("=============================================================================");
print();
print("THEOREM (Thread 15, general form): For any prime p and a2 with D=a2^2-4p^2 < 0,");
print("  sf=squarefree(D) < 0, m^2=D/sf, p∤a2: the prime P above p in Q(sqrt(sf))");
print("  satisfies [P]^2 = 1 in Cl(Q(sqrt(sf))).");
print();
print("Part A: Numerical verification for 10 non-norm-form primes.");
print("  (5 sampled a2 values per prime, covering different CM fields sf.)");
print();

{
  my(primes_done, p, a2_list, a2, res, total_ok, total_tested);
  primes_done = 0;
  total_ok = 0;
  total_tested = 0;
  p = 5;

  while(primes_done < 10,
    p = nextprime(p + 1);
    if(is_norm_form(p), next);

    printf("p = %d (non-norm-form):\n", p);
    a2_list = sample_a2_values(p, 5);

    for(i = 1, #a2_list,
      a2 = a2_list[i];
      res = verify_order2_full(p, a2, 1);
      total_tested++;
      if(res == 1, total_ok++);
      if(res == -1,
        printf("    *** THEOREM VIOLATION at p=%d, a2=%d ***\n", p, a2)));

    primes_done++;
    printf("\n"));

  printf("Part A summary: %d/%d cases PASS.\n\n", total_ok, total_tested);
}

\\ ==============================================================
\\ Part B: Explicit algebraic verification for 3 larger primes
\\ Compute beta element, confirm minpoly = x^2 + a2*x + p^2,
\\ N(beta) = p^2, (beta) = P^2 in O_K.
\\ ==============================================================

print("Part B: Explicit algebraic verification for 3 larger non-norm-form primes.");
print("  For each: construct beta = (-a2 + m*sqrt(sf))/2 as Polmod in Q(sqrt(sf)),");
print("  confirm minpoly = x^2 + a2*x + p^2, N(beta) = p^2, (beta) = P^2 in O_K.");
print();

verify_algebraic(p, a2) = {
  my(D, sf, m2, m, beta, mp, nb, K, Pp, P, P2, I_beta, eq);
  D = a2^2 - 4*p^2;
  sf = sf_part(D);
  m2 = D / sf;
  m = sqrtint(m2);

  \\ beta = (-a2 + m*sqrt(sf)) / 2 as element of Q(sqrt(sf))
  beta = Mod((-a2 + m*x)/2, x^2 - sf);

  \\ minpoly of beta
  mp = minpoly(beta);

  \\ norm (constant term of monic minpoly)
  nb = polcoeff(mp, 0);

  \\ ideal (beta) in O_K
  K  = bnfinit(x^2 - sf, 1);
  Pp = idealprimedec(K, p);
  P  = Pp[1];
  P2 = idealpow(K, P, 2);
  I_beta = idealhnf(K, beta);

  \\ Is (beta) == P^2 or bar-P^2?
  eq = (I_beta == P2 || (#Pp > 1 && I_beta == idealpow(K, Pp[2], 2)));

  printf("  p=%d  a2=%d  sf=%d  m=%d\n", p, a2, sf, m);
  printf("    minpoly(beta) = %Ps\n", mp);
  printf("    expected      = x^2 + %d*x + %d\n", a2, p^2);
  printf("    match: %s\n", if(mp == x^2 + a2*x + p^2, "YES", "NO"));
  printf("    N(beta) = %d = p^2? %s\n", nb, if(nb == p^2, "YES", "NO"));
  printf("    (beta) = P^2? %s\n", if(eq, "YES", "NO"));
  printf("    [P]^2 = 1? %s\n\n",
    if(bnfisprincipal(K, P2)[1] == 0*bnfisprincipal(K, P2)[1], "YES", "NO"));
}

{
  \\ Pick 3 non-norm-form primes (larger than the Part A range)
  my(large_primes, p);
  large_primes = [];
  p = 200;
  while(#large_primes < 3,
    p = nextprime(p + 1);
    if(!is_norm_form(p), large_primes = concat(large_primes, [p])));

  for(i = 1, #large_primes,
    p = large_primes[i];
    \\ Use a mid-range a2 value (roughly a2 ~ p/2)
    my(a2_test);
    a2_test = -(p \ 2);
    \\ Adjust to ensure D < 0, sf < 0, issquare(D/sf)
    \\ Since |a2| < 2p, D = a2^2-4p^2 < 0 always. Just need issquare.
    \\ Try nearby values
    while(!issquare((a2_test^2 - 4*p^2) / sf_part(a2_test^2 - 4*p^2)),
      a2_test--);
    verify_algebraic(p, a2_test));
}

\\ ==============================================================
\\ Part C: Boundary check — does the theorem fail when p | a2?
\\ (This is the excluded case in condition (D).)
\\ ==============================================================

print("Part C: Confirm theorem does NOT apply (vacuously) when p | a2.");
print("  For p | a2, (beta) = (p) is principal trivially — but [P] itself");
print("  might not have order 2.  These are degenerate cases.");
print();

{
  my(p, a2, D, sf);
  p = 53;  \\ arbitrary non-norm-form prime
  a2 = p;  \\ a2 = p, so p | a2

  D = a2^2 - 4*p^2;
  sf = sf_part(D);
  printf("  p=%d, a2=%d (p|a2): D=%d, sf=%d\n", p, a2, D, sf);
  printf("  verify_order2 returns: %d (0=skipped as expected)\n",
    verify_order2_full(p, a2, 1));
  printf("  Reason: condition (D) p∤a2 fails; theorem hypothesis not satisfied.\n");
}

print();
print("Part D: Summary of CM fields encountered across all test cases.");
{
  my(sf_set, p, a2_list, a2, D, sf);
  sf_set = List([]);
  p = 5;
  while(p < 200,
    p = nextprime(p+1);
    if(!is_norm_form(p),
      a2_list = sample_a2_values(p, 5);
      for(i=1,#a2_list,
        a2 = a2_list[i];
        D = a2^2-4*p^2;
        if(D<0,
          sf = sf_part(D);
          if(sf < 0, listput(sf_set, sf))))));
  sf_set = Set(Vec(sf_set));
  printf("  Distinct CM discriminants (sf) encountered: %d unique values\n", #sf_set);
  printf("  Range: sf in [%d, %d]\n", vecmin(Vec(sf_set)), vecmax(Vec(sf_set)));
}

print();
print("DONE.");
