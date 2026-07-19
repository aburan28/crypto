\\ thread16_general_biquadratic.gp
\\
\\ Thread 16: Verify the universal order-2 Frobenius Theorem extends to
\\ NON-norm-form primes.
\\
\\ The Thread 15 proof (A)-(E) never used the norm-form condition 4p=73+3k^2.
\\ It only used:
\\   (i)  beta = (-a2 + m*sqrt(sf)) / 2 satisfies x^2 + a2*x + p^2 = 0
\\   (ii) p does not divide a2
\\   (iii) D = a2^2 - 4p^2 = sf * m^2 (sf squarefree)
\\
\\ NEW KEY INSIGHT (Lemma 1, proved here):
\\   sf = squarefree_part(a2^2 - 4p^2). Since a2^2-4p^2 ≡ a2^2 (mod p),
\\   we have sf * m^2 ≡ a2^2 (mod p), so sf ≡ (a2/m)^2 (mod p).
\\   Therefore sf is always a quadratic residue mod p (or 0 if p | D).
\\   So p NEVER inerts in K = Q(sqrt(sf)): it always splits or ramifies.
\\   This explains why the Theorem is self-consistent.
\\
\\ Run: gp -q thread16_general_biquadratic.gp

default(parisize, 128000000);
default(timer, 0);

\\ -------------------------------------------------------
\\ Utilities
\\ -------------------------------------------------------

sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

is_norm_form(p) = {
  my(v, k);
  v = 4*p - 73;
  if(v <= 0, return(0));
  if(v % 3 != 0, return(0));
  k = sqrtint(v \ 3);
  if(3*k*k != v, return(0));
  if(k % 2 == 0, return(0));
  1
};

\\ -------------------------------------------------------
\\ Verify one (p, a2) pair: returns 1 (pass), 0 (fail), -1 (skip)
\\ -------------------------------------------------------

verify_pair(p, a2) = {
  my(D, sf, m2, m, K, Pp, nP, P, P2, h,
     beta_elt, I_beta, P2_prin, P2_iso_beta, sf_qr, split_char, ord_str);

  if(a2 == 0, return(-1));
  if(a2 % p == 0, return(-1));
  if(abs(a2) >= 2*p, return(-1));

  D  = a2^2 - 4*p^2;
  if(D == 0, return(-1));

  sf = sf_part(D);
  m2 = D / sf;
  if(denominator(m2) != 1, return(-1));
  if(m2 <= 0, return(-1));
  m  = sqrtint(m2);
  if(m*m != m2, return(-1));

  \\ Lemma 1: sf must be QR or 0 mod p
  sf_qr = kronecker(sf, p);
  if(sf_qr == -1,
    printf("  FAIL(Lem1) p=%d a2=%d sf=%d: kron=-1\n", p, a2, sf);
    return(0));
  split_char = if(sf_qr == 0, "ram", "split");

  K  = bnfinit(x^2 - sf, 1);
  h  = K.clgp.no;
  Pp = idealprimedec(K, p);
  nP = #Pp;

  if(nP == 0, return(-1));

  P  = Pp[1];
  P2 = idealpow(K, P, 2);

  P2_prin = bnfisprincipal(K, P2);

  beta_elt = Mod((-a2 + m*x)/2, x^2 - sf);
  I_beta   = idealhnf(K, beta_elt);
  P2_iso_beta = (I_beta == P2);
  if(!P2_iso_beta && nP > 1,
    P2_iso_beta = (I_beta == idealpow(K, Pp[2], 2)));

  if(P2_prin[1] == 0*P2_prin[1],
    ord_str = if(bnfisprincipal(K, P)[1] == 0*bnfisprincipal(K, P)[1], "1", "2"),
    ord_str = ">2(FAIL)");

  printf(
    "  p=%-5d a2=%-7d sf=%-8d m=%-4d h=%-3d nP=%d %-6s  (b)=P2:%s  [P]^2=1:%s  ord:%s\n",
    p, a2, sf, m, h, nP, split_char,
    if(P2_iso_beta, "YES", "NO "),
    if(P2_prin[1] == 0*P2_prin[1], "YES", "NO!"),
    ord_str);

  if(P2_prin[1] != 0*P2_prin[1], return(0));
  1
};

\\ -------------------------------------------------------
\\ Collect 10 non-norm-form primes
\\ -------------------------------------------------------

print("Thread 16: Theorem generalisation to NON-norm-form primes");
print("===========================================================");
print();

{
  my(test_primes, pp);
  test_primes = [];
  pp = 7;
  while(#test_primes < 10,
    if(isprime(pp) && !is_norm_form(pp),
      test_primes = concat(test_primes, [pp]));
    pp = pp + 1);
  printf("Ten non-norm-form test primes: %Ps\n\n", test_primes);

  \\ For each prime, pick 4 a2 values: 2, -2, floor(p/3)+1, -(floor(p/3)+1)
  my(total, passed, failed, p, a2list, r);
  total = 0; passed = 0; failed = 0;

  for(i = 1, #test_primes,
    p = test_primes[i];
    printf("--- p = %d ---\n", p);
    a2list = [2, -2, p\3+1, -(p\3+1)];
    for(j = 1, #a2list,
      r = verify_pair(p, a2list[j]);
      if(r == 1, total++; passed++);
      if(r == 0, total++; failed++));
    print(""));

  printf("SUMMARY: %d pairs, %d passed, %d failed\n", total, passed, failed);
  if(failed == 0,
    print("ALL PASSED -- Theorem holds for non-norm-form primes."),
    print("FAILURES FOUND."));
}

\\ -------------------------------------------------------
\\ Wider sweep: all non-norm-form primes <= 200, a2 in {2,-2,p/3,-p/3}
\\ -------------------------------------------------------

print();
print("Wide sweep: all non-norm-form primes <= 200");
print("--------------------------------------------");

{
  my(total, passed, failed, a2list, r);
  total = 0; passed = 0; failed = 0;

  forprime(pp = 7, 200,
    if(!is_norm_form(pp),
      a2list = [2, -2, pp\3+1, -(pp\3+1)];
      for(j = 1, #a2list,
        r = verify_pair(pp, a2list[j]);
        if(r == 1, total++; passed++);
        if(r == 0, total++; failed++))));

  printf("Wide sweep: %d pairs, %d passed, %d failed\n", total, passed, failed);
  if(failed == 0,
    print("ALL PASSED -- Theorem is universal (norm-form and non-norm-form)."),
    print("FAILURES FOUND in wide sweep."));
}

\\ -------------------------------------------------------
\\ Lemma 1: explicit check over all non-norm-form p <= 200
\\ -------------------------------------------------------

print();
print("Lemma 1: sf(a2^2-4p^2) is always QR or ramified mod p (never inert)");
print("----------------------------------------------------------------------");

{
  my(all_pass, D, sf, kr);
  all_pass = 1;
  forprime(pp = 7, 200,
    if(!is_norm_form(pp),
      for(a2 = 1, pp-1,
        if(a2 % pp != 0 && 2*a2 < pp,
          D = a2^2 - 4*pp^2;
          sf = sf_part(D);
          kr = kronecker(sf, pp);
          if(kr == -1,
            printf("FAIL p=%d a2=%d sf=%d kron=-1\n", pp, a2, sf);
            all_pass = 0)))));
  if(all_pass,
    print("  All (p,a2) with p<=200, 1<=a2<p/2, p not norm-form: kron(sf,p)!=-1. LEMMA 1 OK."),
    print("  LEMMA 1 FAILED somewhere."));
}

print();
print("DONE.");
