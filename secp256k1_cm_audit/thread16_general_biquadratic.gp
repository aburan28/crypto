\\ thread16_general_biquadratic.gp
\\
\\ Thread 16: Verify the universal [P]^2=1 theorem for non-norm-form primes.
\\
\\ THEOREM (proved algebraically in Thread 15):
\\   For any prime p and any biquadratic Weil polynomial T^4 + a2*T^2 + p^2
\\   with D = a2^2 - 4p^2 = sf*m^2 (sf squarefree, m > 0) and p ∤ a2,
\\   every prime P above p in K = Q(sqrt(sf)) satisfies [P]^2 = 1 in Cl(K).
\\
\\ PROOF SUMMARY (Thread 15, steps A-E):
\\   Let beta = (-a2 + m*sqrt(sf))/2.
\\   (A) beta is an algebraic integer satisfying x^2+a2*x+p^2 = 0.
\\   (B) beta in K = Q(sqrt(sf)).
\\   (C) N_{K/Q}(beta) = p^2, so (beta) is an O_K-ideal of norm p^2.
\\   (D) p does not divide a2, so (beta) != (p).
\\   (E) Therefore (beta) = P^2 or Pbar^2, so [P]^2 = 1. QED.
\\
\\ THREAD 16 GOAL: Confirm this theorem holds beyond the secp256k1 norm-form
\\   family (4p = 73 + 3k^2). We use the product construction:
\\     a2 = 2p - t^2  (from (T^2 - tT + p)(T^2 + tT + p) = T^4 + a2*T^2 + p^2)
\\   which gives a biquadratic Weil polynomial for any prime p and any trace t
\\   with |t| < 2*sqrt(p).
\\
\\ Run: gp --stacksize 128000000 -q thread16_general_biquadratic.gp

default(parisize, 128000000);
default(timer, 0);

\\ =====================================================================
\\ Utilities
\\ =====================================================================

sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

\\ Check if p is a secp256k1 norm-form prime: 4p = 73 + 3k^2 for some integer k.
is_norm_form(p) = {
  my(val, r, k);
  val = 4*p - 73;
  if(val < 0, return(0));
  r = val % 3;
  if(r != 0, return(0));
  k = sqrtint(val \ 3);
  return(k*k == val \ 3);
};

\\ =====================================================================
\\ Core verification function
\\ Returns: 1 = theorem holds, 0 = theorem fails,
\\          negative codes for degenerate/inapplicable cases.
\\ =====================================================================

verify_theorem(p, t, verbose) = {
  my(a2, D, sf, m2, m, K, Pp, P, P2, P2_prin, h, ord_P, is_prin);

  a2 = 2*p - t^2;
  if(a2 == 0,   return(-1));   \\ degenerate: |t|=sqrt(2p)
  if(a2 % p == 0, return(-2)); \\ p | a2: theorem doesn't directly apply

  D = a2^2 - 4*p^2;            \\ = t^2*(t^2 - 4p) < 0 for |t| < 2*sqrt(p)
  if(D >= 0, return(-3));       \\ need imaginary quadratic field
  sf = sf_part(D);

  m2 = D / sf;
  if(m2 <= 0 || denominator(m2) != 1, return(-4));
  m = sqrtint(m2);
  if(m*m != m2, return(-5));   \\ D/sf not a perfect square

  K = bnfinit(x^2 - sf, 1);
  h = K.clgp.no;
  Pp = idealprimedec(K, p);
  if(#Pp == 0, return(-6));    \\ p inert or ramifies differently

  P  = Pp[1];
  P2 = idealpow(K, P, 2);
  P2_prin = bnfisprincipal(K, P2);
  is_prin  = (P2_prin[1] == 0*P2_prin[1]);

  \\ Order of [P] in Cl(K): 1 if P principal, 2 if P^2 principal (but P not)
  ord_P = if(bnfisprincipal(K, P)[1] == 0*bnfisprincipal(K, P)[1], 1, 2);

  if(verbose,
    printf("p=%-7d t=%-4d a2=%-10d sf=%-10d m=%-5d h=%-4d ord([P])=%d  [P]^2=1:%s\n",
      p, t, a2, sf, m, h, ord_P, if(is_prin, "YES", "NO")));

  return(if(is_prin, 1, 0));
};

\\ =====================================================================
\\ Part 1: Small non-norm-form primes, t = 1
\\ =====================================================================

print("Thread 16: Universal [P]^2=1 — verification beyond norm-form primes");
print("=====================================================================");
print("Part 1: 15 small primes p not of form 4p=73+3k^2, t=1");
print("(Norm-form primes <=200: 19,37,79,109; all others non-norm-form)");
print();

{
  my(test_primes, p, r, total=0, pass=0);
  test_primes = [11, 17, 23, 29, 41, 43, 47, 53, 59, 61, 67, 71, 83, 89, 97];
  for(i = 1, #test_primes,
    p = test_primes[i];
    if(is_norm_form(p),
      printf("  SKIP: p=%d is norm-form\n", p); next);
    r = verify_theorem(p, 1, 1);
    if(r >= 0, total++);
    if(r == 1, pass++);
    if(r < 0, printf("  (inapplicable: p=%d t=1, code %d)\n", p, r));
  );
  printf("\nPart 1 summary: %d/%d passed\n", pass, total);
}

\\ =====================================================================
\\ Part 2: Larger primes with varied t values
\\ =====================================================================

print();
print("Part 2: Larger primes, varied traces");
print("--------------------------------------");

{
  my(cases, p, t, r, total=0, pass=0);
  \\ [p, t] pairs; all p non-norm-form
  cases = [[1009, 3], [2003, 5], [5003, 11], [10007, 13],
            [50021, 17], [100003, 19], [499979, 23], [999983, 29]];
  for(i = 1, #cases,
    p = cases[i][1];
    t = cases[i][2];
    if(is_norm_form(p),
      printf("  SKIP: p=%d is norm-form\n", p); next);
    r = verify_theorem(p, t, 1);
    if(r >= 0, total++);
    if(r == 1, pass++);
    if(r < 0, printf("  (inapplicable: p=%d t=%d, code %d)\n", p, t, r));
  );
  printf("\nPart 2 summary: %d/%d passed\n", pass, total);
}

\\ =====================================================================
\\ Part 3: Stress test — 30 consecutive non-norm-form primes around p=10^4
\\         all using t=7 (a "generic" non-trivial trace)
\\ =====================================================================

print();
print("Part 3: 30 consecutive primes near p=10000 with t=7");
print("-----------------------------------------------------");

{
  my(p, r, count=0, total=0, pass=0);
  p = nextprime(10000);
  while(count < 30,
    if(!is_norm_form(p),
      r = verify_theorem(p, 7, 1);
      if(r >= 0, total++; if(r==1, pass++));
      if(r < 0, printf("  (inapplicable: p=%d t=7, code %d)\n", p, r));
      count++;
    );
    p = nextprime(p + 1);
  );
  printf("\nPart 3 summary: %d/%d passed\n", pass, total);
}

\\ =====================================================================
\\ Part 4: Class number stress — find cases where h(K) >= 4
\\         (non-trivial class group; theorem is more interesting)
\\ =====================================================================

print();
print("Part 4: Cases with h(K) >= 4 (non-trivial class group)");
print("--------------------------------------------------------");

{
  my(p, t, r, count=0, total=0, pass=0, K, sf, a2, D, m2, m, h);
  p = nextprime(100);
  while(count < 15,
    if(!is_norm_form(p),
      \\ Try multiple t values to get higher h
      forstep(t = 1, floor(2*sqrt(p)) - 1, 1,
        a2 = 2*p - t^2;
        if(a2 == 0 || a2 % p == 0, next);
        D = a2^2 - 4*p^2;
        if(D >= 0, next);
        m2 = D / sf_part(D);
        if(m2 <= 0 || denominator(m2) != 1, next);
        m = sqrtint(m2);
        if(m*m != m2, next);
        sf = sf_part(D);
        K = bnfinit(x^2 - sf, 1);
        h = K.clgp.no;
        if(h >= 4,
          r = verify_theorem(p, t, 1);
          if(r >= 0, total++; if(r==1, pass++));
          count++;
          break);
      );
    );
    p = nextprime(p + 1);
    if(p > 50000, break);  \\ safety
  );
  printf("\nPart 4 summary: %d/%d passed (all with h>=4)\n", pass, total);
}

print();
print("=======================================================================");
print("CONCLUSION: If ALL parts show 100% pass rate,");
print("  the theorem is confirmed general (not just secp256k1 norm-form).");
print("  Proof is algebraic (A)-(E) in thread15_order2_algebraic.gp;");
print("  this script provides numerical evidence for the generality claim.");
print("=======================================================================");
print();
print("DONE.");
