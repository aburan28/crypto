\\ thread16_general_weil_order2.gp
\\
\\ Thread 16: Verify the biquadratic-Weil order-2 theorem for
\\ NON-norm-form primes.
\\
\\ BACKGROUND (Thread 15 algebraic proof):
\\   For any ordinary E/F_p with trace t (t != 0, p ∤ t):
\\     Weil poly of E x E^t  =  T^4 + a2*T^2 + p^2,   a2 = 2p - t^2
\\     D = a2^2 - 4*p^2     =  t^2*(t^2 - 4*p)   (<= 0 for |t| <= 2*sqrt(p))
\\     K = Q(sqrt(sf(D))),   P = prime above p in O_K
\\   Claim: [P]^2 = 1 in Cl(K).
\\   Proof (A)-(E):
\\   (A) beta = (-a2 + m*sqrt(sf)) / 2  in O_K  (root of x^2+a2*x+p^2)
\\   (B) N(beta) = p^2
\\   (C) (beta) != (p) since p ∤ a2 (equiv: p ∤ t for t != 0, |t| < 2*sqrt(p))
\\   (D) Therefore (beta) = P^2 or Pbar^2, so [P]^2 = 1.
\\
\\ EXPERIMENT (Thread 16):
\\   - Pick 10 random primes p NOT of norm-form (73+3k^2)/4.
\\   - For each, pick 3 random traces t (|t| < 2*sqrt(p), t!=0, p !| t).
\\   - Verify [P]^2 = 1 in all 30 cases.
\\
\\ Run: gp --stacksize 64000000 -q thread16_general_weil_order2.gp

default(parisize, 64000000);
default(timer, 0);

\\ ------------------------------------------------------------------
\\ Utility: squarefree part with sign
\\ ------------------------------------------------------------------

sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

\\ ------------------------------------------------------------------
\\ Is p a norm-form prime? (4p = 73 + 3k^2 for some odd k)
\\ ------------------------------------------------------------------

is_norm_form(p) = {
  my(v, k2);
  v = 4*p - 73;
  if(v <= 0 || v % 3 != 0, return(0));
  k2 = v \ 3;
  if(!issquare(k2), return(0));
  my(k = sqrtint(k2));
  (k % 2 == 1)
};

\\ ------------------------------------------------------------------
\\ Verify order-2 for a single (p, t) pair
\\   Returns 1 if [P]^2=1, 0 on any error/failure.
\\ ------------------------------------------------------------------

verify_one(p, t) = {
  my(a2, D, sf, m2, m, K, Pp, P, P2, P2_prin, beta_elt, norm_beta, label);

  if(t == 0 || (t % p) == 0,
    printf("  (p=%d, t=%d) SKIP: t=0 or p|t\n", p, t);
    return(-1));

  a2 = 2*p - t^2;
  D  = a2^2 - 4*p^2;       \\ = t^2*(t^2-4*p)

  if(D == 0,
    printf("  (p=%d, t=%d) SKIP: D=0 (supersingular-like)\n", p, t);
    return(-1));
  if(D > 0,
    printf("  (p=%d, t=%d) SKIP: D>0 (unexpected for |t|<2*sqrt(p))\n", p, t);
    return(-1));

  sf = sf_part(D);
  m2 = D / sf;
  if(denominator(m2) != 1 || m2 <= 0,
    printf("  (p=%d, t=%d) ERROR: D/sf not a positive integer\n", p, t);
    return(0));
  if(!issquare(m2),
    printf("  (p=%d, t=%d) ERROR: D/sf=%d not a perfect square\n", p, t, m2);
    return(0));
  m = sqrtint(m2);

  \\ Norm check
  norm_beta = (a2^2 - m2*sf) / 4;
  if(norm_beta != p^2,
    printf("  (p=%d, t=%d) ERROR: N(beta)=%d != p^2=%d\n", p, t, norm_beta, p^2);
    return(0));

  \\ p | a2 check: a2 = 2p - t^2; p|a2 iff p|t^2 iff p|t.
  if(a2 % p == 0,
    printf("  (p=%d, t=%d) SKIP: p|a2 (degenerate)\n", p, t);
    return(-1));

  \\ PARI number field and ideal computation
  K  = bnfinit(x^2 - sf, 1);
  Pp = idealprimedec(K, p);
  if(#Pp == 0,
    printf("  (p=%d, t=%d) SKIP: p inert in Q(sqrt(%d))\n", p, t, sf);
    return(-1));

  P  = Pp[1];
  P2 = idealpow(K, P, 2);
  P2_prin = bnfisprincipal(K, P2);

  my(is_principal = (P2_prin[1] == 0*P2_prin[1]));

  printf("  p=%-8d  t=%-6d  a2=%-10d  sf=%-8d  m=%-5d  h=%-4d  [P]^2=1:%s\n",
    p, t, a2, sf, m, K.clgp.no,
    if(is_principal, "YES", "NO(!)"));

  if(is_principal, 1, 0)
};

\\ ------------------------------------------------------------------
\\ Generate a non-norm-form prime >= lo with approximately nbits bits
\\ ------------------------------------------------------------------

next_non_norm_prime(lo) = {
  my(p = nextprime(lo));
  while(is_norm_form(p), p = nextprime(p+1));
  p
};

\\ ------------------------------------------------------------------
\\ Main experiment: 10 non-norm-form primes, 3 traces each
\\ ------------------------------------------------------------------

print("Thread 16: Biquadratic Weil order-2 — generalization to non-norm-form primes");
print("==============================================================================");
print();
print("Proof recap: For any ordinary E/F_p with trace t (t!=0, p!|t),");
print("  D = t^2*(t^2-4*p) < 0, sf = squarefree_part(D),");
print("  beta = (-a2 + m*sqrt(sf))/2 in O_{Q(sqrt(sf))},");
print("  N(beta)=p^2 and p!|a2 => (beta)=P^2 => [P]^2=1 in Cl(Q(sqrt(sf))).");
print();

{
  my(total = 0, ok = 0, skip = 0);
  my(primes_todo = 10, found = 0);
  my(p_start = 100);   \\ start search above 100

  \\ Use fixed small primes (deterministic: no Math.random)
  \\ 10 non-norm-form primes in various size ranges
  my(test_primes = [101, 127, 211, 307, 503, 751, 1009, 2003, 5003, 10007]);

  for(idx = 1, #test_primes,
    my(p_cand = test_primes[idx]);
    if(!isprime(p_cand), next());
    if(is_norm_form(p_cand),
      printf("p=%d IS norm-form, skipping\n", p_cand);
      next());

    found++;
    printf("\n[%d] p = %d (non-norm-form)  h_range(p)=[1..%d]\n",
      found, p_cand, ceil(4*sqrt(p_cand)*log(p_cand)/Pi));

    \\ Pick 3 traces: small, medium, near-maximal
    my(sq = sqrtint(p_cand));
    my(traces = [1, sq \ 3, sq]);

    for(ti = 1, #traces,
      my(t = traces[ti]);
      \\ Ensure |t| < 2*sqrt(p) strictly (Hasse) and t != 0 and p !| t
      if(t == 0 || t >= 2*sq + 2, t = 1);
      if(t % p_cand == 0, t++);

      my(r = verify_one(p_cand, t));
      if(r == 1, ok++; total++,
       r == 0, total++;  \\ failure
       \\ r == -1: skip (inert, D=0, etc.) => no count
      )
    )
  );

  printf("\n--------------------------------------------------------------\n");
  printf("Summary: %d cases verified, %d passed [P]^2=1, %d failed\n",
    total, ok, total - ok);
  printf("All passed: %s\n", if(total > 0 && ok == total, "YES", "NO or PARTIAL"));
}

print();
print("DONE.");
