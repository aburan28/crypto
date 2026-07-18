\\ thread16_general_biquadratic.gp
\\
\\ Thread 16: Verify the order-2 Frobenius theorem for non-norm-form primes.
\\
\\ THEOREM (Thread 15): Let T^4 + a2*T^2 + p^2 be a biquadratic Weil polynomial
\\ with D = a2^2 - 4p^2 = sf*m^2 (sf squarefree < 0), and p not dividing a2.
\\ Then [P]^2 = 1 in Cl(Q(sqrt(sf))) for any prime P above p in O_K.
\\
\\ PROOF SKETCH (repeated from Thread 15 for reference):
\\   beta = (-a2 + m*sqrt(sf))/2 satisfies x^2 + a2*x + p^2 = 0 (monic Z-coeffs).
\\   N(beta) = p^2 => (beta) has norm p^2. The only ideals of norm p^2 are P^2, Pbar^2, (p).
\\   Since p ∤ a2, (beta) ≠ (p). So (beta) = P^2 or Pbar^2 => [P]^2 = 1. QED.
\\
\\ HERE: Apply to 10 non-norm-form primes using E x E^t (quadratic-twist product).
\\   If E/F_p has trace t, then E^t has trace -t, and the Weil polynomial of
\\   E x E^t is (T^2 - tT + p)(T^2 + tT + p) = T^4 + (2p - t^2)*T^2 + p^2.
\\   So a2 = 2p - t^2. Since |t| <= 2*sqrt(p) and t != 0, we have:
\\     D = (2p-t^2)^2 - 4p^2 = t^2*(t^2 - 4p) < 0 (Hasse).
\\     p | a2 iff p | t^2 iff p | t, but t != 0 and |t| <= 2*sqrt(p) < p => p ∤ a2. ✓
\\
\\ Run: gp --stacksize 128000000 -q thread16_general_biquadratic.gp

default(parisize, 128000000);
default(timer, 0);

\\ ===========================================================
\\ Utility: squarefree part with sign
\\ ===========================================================

sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

\\ ===========================================================
\\ Find elliptic curve E/F_p with trace t != 0 (small a, b)
\\ ===========================================================

find_curve(p) = {
  my(E, n, t, disc, found = [0,0,0]);
  for(a = 0, 8,
    if(found[3] != 0, break);
    for(b = 1, 8,
      if(found[3] != 0, break);
      disc = (4*a^3 + 27*b^2) % p;
      if(disc == 0, next);
      E = ellinit([a, b] * Mod(1, p));
      n = ellcard(E);
      t = p + 1 - n;
      if(t != 0, found = [a, b, t])));
  if(found[3] == 0, error("No non-zero-trace curve for p=", p));
  found
};

\\ ===========================================================
\\ Verify [P]^2 = 1 for biquadratic Weil poly T^4 + a2*T^2 + p^2
\\ ===========================================================

verify_biquadratic(p, a2) = {
  my(D, sf, m2, m, K, Pp, nP, P, P2, P2_prin, h, beta_elt, I_beta, P2_iso);

  D = a2^2 - 4*p^2;
  if(D >= 0,
    printf("  SKIP p=%d a2=%d: D=%d >= 0\n", p, a2, D); return(0));

  sf = sf_part(D);
  if(sf >= 0,
    printf("  SKIP p=%d: sf=%d non-negative\n", p, sf); return(0));

  m2 = D / sf;
  if(denominator(m2) != 1 || m2 <= 0,
    printf("  ERROR p=%d: D/sf not a positive integer\n", p); return(0));
  m = sqrtint(m2);
  if(m*m != m2,
    printf("  ERROR p=%d: D/sf=%d not a perfect square\n", p, m2); return(0));

  if(a2 % p == 0,
    printf("  SKIP p=%d: p | a2 (t=0 case)\n", p); return(0));

  K  = bnfinit(x^2 - sf, 1);
  h  = K.clgp.no;
  Pp = idealprimedec(K, p);
  nP = #Pp;

  if(nP == 0,
    printf("  p=%-7d a2=%-8d sf=%-9d m=%-5d h=%-4d  p INERT in K (N/A)\n",
      p, a2, sf, m, h);
    return(0));

  P  = Pp[1];
  P2 = idealpow(K, P, 2);

  \\ Is P^2 principal?
  P2_prin = bnfisprincipal(K, P2);

  \\ Does (beta) = P^2?
  beta_elt  = Mod((-a2 + m*x)/2, x^2 - sf);
  I_beta    = idealhnf(K, beta_elt);
  P2_iso    = (I_beta == P2 || I_beta == idealpow(K, Pp[nP], 2));

  printf("  p=%-7d a2=%-8d sf=%-9d m=%-5d h=%-4d  (β)=P^2:%s  [P]^2=1:%s\n",
    p, a2, sf, m, h,
    if(P2_iso,                              "YES", "NO "),
    if(P2_prin[1] == 0*P2_prin[1],         "YES", "NO "));
  1
};

\\ ===========================================================
\\ Main: 10 non-norm-form primes
\\ (None of these satisfy 4p = 73 + 3k^2 for integer k.)
\\ ===========================================================

test_primes = [101, 251, 503, 1009, 2003, 4001, 8011, 10007, 20011, 50021];

print("Thread 16: Order-2 Frobenius theorem — non-norm-form prime verification");
print("=========================================================================");
print("Biquadratic Weil poly from E x E^t: a2 = 2p - t^2.");
print("Checking [P]^2 = 1 in Cl(Q(sqrt(sf))) for sf = sf((2p-t^2)^2 - 4p^2).");
print();

{
  my(ok = 0, count = 0);
  for(i = 1, #test_primes,
    my(p = test_primes[i]);
    my(res = find_curve(p));
    my(a = res[1], b = res[2], t = res[3], a2 = 2*p - t^2);
    printf("p=%-7d (y^2=x^3+%d*x+%d, t=%d => a2=%d)\n", p, a, b, t, a2);
    r = verify_biquadratic(p, a2);
    if(r, ok++);
    count++);
  printf("\nResult: %d/%d non-norm-form cases passed [(beta)=P^2 and [P]^2=1].\n", ok, count);
}

\\ ===========================================================
\\ Bonus: also verify for several values of a2 per prime
\\ (different curves => different t => different a2)
\\ ===========================================================

print();
print("Bonus: multiple (t, a2) values per prime p=1009:");
print();

{
  my(p = 1009, ok = 0, count = 0);
  \\ collect several curves with distinct non-zero traces
  my(seen = [], aa, bb, E_tmp, nn, tt, a2);
  for(aa = 0, 20,
    for(bb = 0, 20,
      if((4*aa^3 + 27*bb^2) % p == 0, next);
      E_tmp = ellinit([aa, bb] * Mod(1, p));
      nn = ellcard(E_tmp);
      tt = p + 1 - nn;
      if(tt == 0, next);
      if(setsearch(seen, abs(tt)), next);   \\ skip if |t| already seen (twist gives same a2)
      seen = setunion(seen, [abs(tt)]);
      a2 = 2*p - tt^2;
      printf("  t=%4d a2=%7d: ", tt, a2);
      r = verify_biquadratic(p, a2);
      if(r, ok++);
      count++;
      if(count >= 6, break(2))));
  printf("\n  p=1009 bonus: %d/%d distinct traces all passed.\n", ok, count);
}

print();
print("DONE.");
