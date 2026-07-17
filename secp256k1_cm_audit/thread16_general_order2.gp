\\ thread16_general_order2.gp
\\
\\ Thread 16: Verify the "universal order-2 Frobenius" theorem holds for
\\ non-norm-form primes p (i.e., primes NOT of the form (73 + 3k^2)/4).
\\
\\ THEOREM (Thread 15, general form):
\\   Let p be any prime and t any integer with 0 < |t| < 2*sqrt(p).
\\   Set a2 = 2p - t^2, D = a2^2 - 4p^2 = t^2*(t^2 - 4p) < 0.
\\   Write D = sf * m^2 with sf squarefree negative, m > 0.
\\   Let K = Q(sqrt(sf)), P a prime of O_K above p (p splits since sf < 0
\\   and p is odd prime — splitting determined by Legendre (sf/p)).
\\   IF p does NOT divide a2, THEN [P]^2 = 1 in Cl(K).
\\
\\ PROOF (Thread 15): beta = (-a2 + m*sqrt(sf))/2 satisfies x^2 + a2*x + p^2 = 0
\\   (algebraic integer in O_K), has norm p^2, and (beta) != (p) since p|a2 is
\\   excluded. Hence (beta) = P^2 or Pbar^2, so [P]^2 = 1.  QED.
\\
\\ This script checks the theorem for 10 non-norm-form primes, each with 5 traces.
\\
\\ Run: gp --stacksize 128000000 -q thread16_general_order2.gp

default(parisize, 128000000);
default(timer, 0);

\\ ==========================================================
\\ is_norm_form(p): is p of the form (73 + 3k^2)/4 for some odd k?
\\ ==========================================================

is_norm_form(p) = {
  my(val, k);
  \\ 4p - 73 must be 3k^2 for odd positive k
  val = 4*p - 73;
  if(val <= 0, return(0));
  if(val % 3 != 0, return(0));
  k = sqrtint(val \ 3);
  if(k*k*3 != val, return(0));
  if(k % 2 == 0, return(0));  \\ k must be odd
  1
};

\\ ==========================================================
\\ sf_part(n): squarefree part of n (preserving sign)
\\ ==========================================================

sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

\\ ==========================================================
\\ verify_one(p, t): verify the theorem for prime p and trace t.
\\ Returns 1 if verified, 0 if inapplicable (p|a2, no split, etc.)
\\ ==========================================================

verify_one(p, t) = {
  my(a2, D, sf, m2, m, K, Pp, P, P2, P2_prin, beta_elt, norm_beta, ord);

  a2 = 2*p - t^2;
  if(a2 == 0, return(-1));          \\ degenerate: skip t with a2=0

  D  = a2^2 - 4*p^2;               \\ = t^2*(t^2 - 4p) <= 0
  if(D == 0, return(-1));           \\ t=0 degenerate
  if(D > 0,  return(-1));           \\ shouldn't happen in Hasse range

  sf = sf_part(D);
  m2 = D / sf;
  if(m2 <= 0 || denominator(m2) != 1, return(-1));
  m  = sqrtint(m2);
  if(m*m != m2, return(-1));

  \\ Check p | a2: theorem requires p !| a2
  if(a2 % p == 0,
    printf("  p=%d t=%d: SKIP (p | a2=%d)\n", p, t, a2);
    return(-2));

  \\ Build K = Q(sqrt(sf)); check p splits in K
  K  = bnfinit(x^2 - sf, 1);
  Pp = idealprimedec(K, p);
  if(#Pp == 0,
    printf("  p=%d t=%d: p INERT in Q(sqrt(%d)), skip\n", p, t, sf);
    return(-3));

  P  = Pp[1];
  P2 = idealpow(K, P, 2);

  \\ Principal test for P^2
  P2_prin = bnfisprincipal(K, P2);

  \\ Verify N(beta) = p^2
  norm_beta = (a2^2 - m2*sf) / 4;
  if(norm_beta != p^2,
    printf("  p=%d t=%d: ERROR N(beta)=%d != p^2=%d\n", p, t, norm_beta, p^2);
    return(0));

  \\ Determine order of [P] in Cl(K)
  if(P2_prin[1] == 0*P2_prin[1],
    ord = if(bnfisprincipal(K, P)[1] == 0*bnfisprincipal(K, P)[1], 1, 2),
    ord = -1);  \\ [P]^2 != 1 would be a counterexample

  printf("  p=%-6d t=%-4d a2=%-8d sf=%-8d m=%-5d h=%-3d  [P]^2=1: %s\n",
    p, t, a2, sf, m, K.clgp.no,
    if(ord >= 1, "YES", "NO(!!)"));

  if(ord < 1, 0, 1)
};

\\ ==========================================================
\\ Main: pick 10 non-norm-form primes, verify with 5 traces each
\\ ==========================================================

print("Thread 16: Universal order-2 Frobenius — non-norm-form primes");
print("===============================================================");
print();

{
  my(tested_primes, p, found, total_checks, total_pass, sqrtp, t, r, traces);

  tested_primes = 0;
  total_checks  = 0;
  total_pass    = 0;

  p = 3;
  while(tested_primes < 10,
    p = nextprime(p + 1);
    if(!is_norm_form(p),
      tested_primes++;
      printf("--- p = %d (non-norm-form) ---\n", p);

      sqrtp = sqrtint(p);
      \\ Pick 5 non-zero traces spread across the Hasse bound
      traces = [1, 2, sqrtp-1, sqrtp, sqrtp+1];
      \\ keep only traces with |t| < 2*sqrt(p) (Hasse) and t != 0
      found = 0;
      for(i = 1, #traces,
        t = traces[i];
        if(t > 0 && t^2 < 4*p,
          r = verify_one(p, t);
          if(r == 1, total_checks++; total_pass++; found++);
          if(r == 0, total_checks++; found++);
          \\ also check negative t (different a2 since a2 = 2p - t^2 same)
          \\ t and -t give same a2, so skip -t
        )
      );
      if(found == 0,
        printf("  (no valid traces found for p=%d)\n", p));
      print()
    )
  );

  printf("Summary: %d checks, %d passed (all [P]^2=1: %s)\n",
    total_checks, total_pass,
    if(total_checks == total_pass, "YES", "NO(!!)"));
}

print();
print("DONE.");
